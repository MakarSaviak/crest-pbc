#define _GNU_SOURCE

#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <linux/mempolicy.h>
#include <sched.h>
#include <signal.h>
#include <spawn.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/prctl.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

extern char **environ;

enum {
  /* Historical optimizer-affinity reference geometry.  These constants are
   * intentionally NOT process-MTD scheduler limits. */
  CREST_OPT_GROUPS = 48,
  CREST_OPT_GROUP_THREADS = 4,
  CREST_OPT_CPUS = CREST_OPT_GROUPS * CREST_OPT_GROUP_THREADS,
  CREST_MTD_MAX_CPUS = CPU_SETSIZE
};

static const char optimizer_affinity_optin_name[] =
    "CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY";

struct cpu_record {
  int cpu;
  int node;
  int package;
  int die;
  int core;
  int l3_first_cpu;
};

struct env_override {
  const char *name;
  const char *value;
};

static int env_name_matches(const char *entry, const char *name) {
  const size_t n = strlen(name);
  return strncmp(entry, name, n) == 0 && entry[n] == '=';
}

static void free_child_environment(char **child_env) {
  size_t i;
  if (child_env == NULL) return;
  for (i = 0; child_env[i] != NULL; ++i) free(child_env[i]);
  free(child_env);
}

static int build_child_environment(const char *places, int threads,
                                   char ***child_env_out) {
  char thread_text[32];
  const struct env_override overrides[] = {
    {"OMP_NUM_THREADS", thread_text},
    {"OMP_DYNAMIC", "FALSE"},
    {"OMP_NESTED", "FALSE"},
    {"OMP_MAX_ACTIVE_LEVELS", "1"},
    {"OMP_THREAD_LIMIT", thread_text},
    {"OMP_PROC_BIND", "close"},
    {"OMP_PLACES", places},
    {"OPENBLAS_NUM_THREADS", "1"},
    {"GOTO_NUM_THREADS", "1"},
    {"MKL_NUM_THREADS", "1"},
    {"BLIS_NUM_THREADS", "1"},
    {"VECLIB_MAXIMUM_THREADS", "1"},
    {"NUMEXPR_NUM_THREADS", "1"},
    {"CREST_MTD_INTERNAL_WORKER", "1"}
  };
  const size_t noverrides = sizeof(overrides) / sizeof(overrides[0]);
  const char *removed[] = {
    "GOMP_CPU_AFFINITY",
    "GOMP_SPINCOUNT",
    "OMP_SCHEDULE",
    "OMP_WAIT_POLICY",
    optimizer_affinity_optin_name
  };
  const size_t nremoved = sizeof(removed) / sizeof(removed[0]);
  size_t inherited = 0, kept = 0, i, j, out = 0;
  char **child_env;

  if (places == NULL || child_env_out == NULL || threads < 1 ||
      threads > CREST_MTD_MAX_CPUS) return EINVAL;
  if (snprintf(thread_text, sizeof(thread_text), "%d", threads) >=
      (int)sizeof(thread_text)) return EOVERFLOW;

  while (environ[inherited] != NULL) ++inherited;
  for (i = 0; i < inherited; ++i) {
    int replace = 0;
    for (j = 0; j < nremoved && !replace; ++j)
      replace = env_name_matches(environ[i], removed[j]);
    for (j = 0; j < noverrides && !replace; ++j)
      replace = env_name_matches(environ[i], overrides[j].name);
    if (!replace) ++kept;
  }
  child_env = (char **)calloc(kept + noverrides + 1, sizeof(char *));
  if (child_env == NULL) return ENOMEM;
  for (i = 0; i < inherited; ++i) {
    int replace = 0;
    for (j = 0; j < nremoved && !replace; ++j)
      replace = env_name_matches(environ[i], removed[j]);
    for (j = 0; j < noverrides && !replace; ++j)
      replace = env_name_matches(environ[i], overrides[j].name);
    if (!replace) {
      child_env[out] = strdup(environ[i]);
      if (child_env[out] == NULL) {
        free_child_environment(child_env);
        return ENOMEM;
      }
      ++out;
    }
  }
  for (j = 0; j < noverrides; ++j) {
    const size_t n = strlen(overrides[j].name) + strlen(overrides[j].value) + 2;
    child_env[out] = (char *)malloc(n);
    if (child_env[out] == NULL) {
      free_child_environment(child_env);
      return ENOMEM;
    }
    if (snprintf(child_env[out], n, "%s=%s", overrides[j].name,
                 overrides[j].value) >= (int)n) {
      free_child_environment(child_env);
      return EOVERFLOW;
    }
    ++out;
  }
  child_env[out] = NULL;
  *child_env_out = child_env;
  return 0;
}

#ifdef CREST_OPTIMIZER_AFFINITY_TESTING
int crest_optimizer_affinity_test_child_environment(void) {
  char **child_env = NULL;
  const char *parent_value = getenv(optimizer_affinity_optin_name);
  size_t i;
  int optin_entries = 0, worker_entries = 0, rc;

  if (parent_value == NULL || strcmp(parent_value, "1") != 0) return EINVAL;
  rc = build_child_environment("{0},{1},{2},{3}", 4, &child_env);
  if (rc != 0) return rc;
  for (i = 0; child_env[i] != NULL; ++i) {
    if (env_name_matches(child_env[i], optimizer_affinity_optin_name))
      ++optin_entries;
    if (strcmp(child_env[i], "CREST_MTD_INTERNAL_WORKER=1") == 0)
      ++worker_entries;
  }
  free_child_environment(child_env);
  if (optin_entries != 0) return EPERM;
  if (worker_entries != 1) return EDOM;
  return 0;
}
#endif

static int read_int_file_exact(const char *path, int *value) {
  FILE *fp = fopen(path, "r");
  int close_rc;
  if (value == NULL) return EINVAL;
  if (fp == NULL) return errno;
  if (fscanf(fp, "%d", value) != 1) {
    fclose(fp);
    return EIO;
  }
  close_rc = fclose(fp);
  return close_rc == 0 ? 0 : (errno ? errno : EIO);
}

static int cpu_numa_node(int cpu) {
  char path[PATH_MAX];
  DIR *dir;
  struct dirent *entry;
  int node = -1;

  if (snprintf(path, sizeof(path), "/sys/devices/system/cpu/cpu%d", cpu) >=
      (int)sizeof(path)) return -1;
  dir = opendir(path);
  if (dir == NULL) return -1;
  while ((entry = readdir(dir)) != NULL) {
    int candidate;
    char tail;
    if (sscanf(entry->d_name, "node%d%c", &candidate, &tail) == 1) {
      node = candidate;
      break;
    }
  }
  closedir(dir);
  return node;
}

static int cpu_l3_first_cpu(int cpu, int *first_cpu_out) {
  char base[PATH_MAX], path[PATH_MAX], buffer[256];
  DIR *dir;
  struct dirent *entry;
  int found = 0;
  if (first_cpu_out == NULL) return EINVAL;
  if (snprintf(base, sizeof(base), "/sys/devices/system/cpu/cpu%d/cache", cpu) >=
      (int)sizeof(base)) return EOVERFLOW;
  dir = opendir(base);
  if (dir == NULL) return errno;
  while ((entry = readdir(dir)) != NULL) {
    FILE *fp;
    int level = -1;
    char type[32] = "";
    char *end = NULL;
    long first;
    if (strncmp(entry->d_name, "index", 5) != 0) continue;
    if (snprintf(path, sizeof(path), "%s/%s/level", base, entry->d_name) >=
        (int)sizeof(path)) { closedir(dir); return EOVERFLOW; }
    fp = fopen(path, "r");
    if (fp == NULL) continue;
    if (fscanf(fp, "%d", &level) != 1) level = -1;
    fclose(fp);
    if (level != 3) continue;
    if (snprintf(path, sizeof(path), "%s/%s/type", base, entry->d_name) >=
        (int)sizeof(path)) { closedir(dir); return EOVERFLOW; }
    fp = fopen(path, "r");
    if (fp == NULL) continue;
    if (fscanf(fp, "%31s", type) != 1) type[0] = '\0';
    fclose(fp);
    if (strcasecmp(type, "Unified") != 0) continue;
    if (snprintf(path, sizeof(path), "%s/%s/shared_cpu_list", base, entry->d_name) >=
        (int)sizeof(path)) { closedir(dir); return EOVERFLOW; }
    fp = fopen(path, "r");
    if (fp == NULL) continue;
    if (fgets(buffer, sizeof(buffer), fp) == NULL) { fclose(fp); continue; }
    fclose(fp);
    errno = 0;
    first = strtol(buffer, &end, 10);
    if (errno != 0 || end == buffer || first < 0 || first >= CPU_SETSIZE) {
      closedir(dir);
      return EIO;
    }
    *first_cpu_out = (int)first;
    found = 1;
    break;
  }
  closedir(dir);
  return found ? 0 : ENODEV;
}

static int fill_cpu_record(struct cpu_record *rec, int cpu) {
  char path[PATH_MAX];
  int rc;
  rec->cpu = cpu;
  rec->node = cpu_numa_node(cpu);
  rec->l3_first_cpu = -1;
  if (rec->node < 0) return ENODEV;
  if (snprintf(path, sizeof(path),
               "/sys/devices/system/cpu/cpu%d/topology/physical_package_id", cpu) >=
      (int)sizeof(path)) return EOVERFLOW;
  rc = read_int_file_exact(path, &rec->package);
  if (rc != 0) return rc;
  if (snprintf(path, sizeof(path),
               "/sys/devices/system/cpu/cpu%d/topology/die_id", cpu) >=
      (int)sizeof(path)) return EOVERFLOW;
  rc = read_int_file_exact(path, &rec->die);
  if (rc != 0) return rc;
  if (snprintf(path, sizeof(path),
               "/sys/devices/system/cpu/cpu%d/topology/core_id", cpu) >=
      (int)sizeof(path)) return EOVERFLOW;
  return read_int_file_exact(path, &rec->core);
}

static int fill_mtd_cpu_record(struct cpu_record *rec, int cpu) {
  int rc = fill_cpu_record(rec, cpu);
  if (rc != 0) return rc;
  return cpu_l3_first_cpu(cpu, &rec->l3_first_cpu);
}

static int bind_memory_to_worker_node(const int cpus[CREST_OPT_GROUP_THREADS]) {
#if defined(SYS_set_mempolicy) && defined(SYS_get_mempolicy)
  enum { BITS_PER_WORD = (int)(8 * sizeof(unsigned long)) };
  unsigned long requested[(CPU_SETSIZE + BITS_PER_WORD - 1) / BITS_PER_WORD];
  unsigned long observed[(CPU_SETSIZE + BITS_PER_WORD - 1) / BITS_PER_WORD];
  const unsigned long maxnode = (unsigned long)(sizeof(requested) * 8);
  int node = cpu_numa_node(cpus[0]);
  int mode = -1;
  int i;
  if (node < 0 || (unsigned long)node >= maxnode) return ENODEV;
  for (i = 1; i < CREST_OPT_GROUP_THREADS; ++i)
    if (cpu_numa_node(cpus[i]) != node) return EXDEV;
  memset(requested, 0, sizeof(requested));
  requested[node / BITS_PER_WORD] |= 1UL << (node % BITS_PER_WORD);
  if (syscall(SYS_set_mempolicy, MPOL_BIND, requested, maxnode) != 0)
    return errno;
  memset(observed, 0, sizeof(observed));
  if (syscall(SYS_get_mempolicy, &mode, observed, maxnode, NULL, 0) != 0)
    return errno;
  if ((mode & ~MPOL_MODE_FLAGS) != MPOL_BIND) return EINVAL;
  if ((observed[node / BITS_PER_WORD] & (1UL << (node % BITS_PER_WORD))) == 0)
    return EINVAL;
  return 0;
#else
  (void)cpus;
  return ENOSYS;
#endif
}

static int compare_cpu_record(const void *left, const void *right) {
  const struct cpu_record *a = (const struct cpu_record *)left;
  const struct cpu_record *b = (const struct cpu_record *)right;
#define CMP_FIELD(field) do { if (a->field != b->field) \
  return (a->field > b->field) - (a->field < b->field); } while (0)
  CMP_FIELD(node);
  CMP_FIELD(package);
  CMP_FIELD(die);
  CMP_FIELD(core);
  CMP_FIELD(cpu);
#undef CMP_FIELD
  return 0;
}

static int slurm_allocation_is_full(void) {
  const char *names[] = {"SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK",
                         "SLURM_JOB_CPUS_PER_NODE"};
  size_t i;
  for (i = 0; i < sizeof(names) / sizeof(names[0]); ++i) {
    const char *value = getenv(names[i]);
    char *end = NULL;
    long parsed;
    if (value == NULL || *value == '\0') continue;
    errno = 0;
    parsed = strtol(value, &end, 10);
    if (errno == 0 && end != value && parsed == CREST_OPT_CPUS) return 1;
  }
  return 0;
}

/*
 * Build four-CPU groups within one NUMA node, then enumerate those groups in
 * NUMA round-robin order.  This is the process analogue of spread,close:
 * workers are spread while every worker's calculator threads stay local.
 * The prototype intentionally rejects anything other than a full 192-CPU
 * allocation.
 */
static int build_cpu_groups(int groups[CREST_OPT_GROUPS][CREST_OPT_GROUP_THREADS]) {
  struct cpu_record records[CPU_SETSIZE];
  cpu_set_t allowed;
  int nodes[CREST_OPT_CPUS];
  int node_group_start[CREST_OPT_CPUS];
  int node_group_count[CREST_OPT_CPUS];
  int node_cursor[CREST_OPT_CPUS];
  int ncpu = 0, nnode = 0, ngroup = 0;
  int cpu, i, j, round, progress;

  if (sysconf(_SC_NPROCESSORS_ONLN) != CREST_OPT_CPUS) return EDOM;
  if (!slurm_allocation_is_full()) return EPERM;
  CPU_ZERO(&allowed);
  if (sched_getaffinity(0, sizeof(allowed), &allowed) != 0) return errno;
  if (CPU_COUNT(&allowed) != CREST_OPT_CPUS) return EPERM;
  for (cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
    int rc;
    if (!CPU_ISSET(cpu, &allowed)) continue;
    if (ncpu >= CPU_SETSIZE) return EOVERFLOW;
    rc = fill_cpu_record(&records[ncpu], cpu);
    if (rc != 0) return rc;
    ++ncpu;
  }
  if (ncpu != CREST_OPT_CPUS) return EDOM;
  qsort(records, (size_t)ncpu, sizeof(records[0]), compare_cpu_record);
  for (i = 1; i < ncpu; ++i) {
    if (records[i].node == records[i - 1].node &&
        records[i].package == records[i - 1].package &&
        records[i].die == records[i - 1].die &&
        records[i].core == records[i - 1].core) {
      /* SMT siblings are not four independent physical cores. */
      return ENOTSUP;
    }
  }

  i = 0;
  while (i < ncpu) {
    const int node = records[i].node;
    int end = i + 1;
    while (end < ncpu && records[end].node == node) ++end;
    if (node < 0 || (end - i) % CREST_OPT_GROUP_THREADS != 0) return EXDEV;
    nodes[nnode] = node;
    node_group_start[nnode] = ngroup;
    node_group_count[nnode] = (end - i) / CREST_OPT_GROUP_THREADS;
    node_cursor[nnode] = 0;
    for (j = i; j < end; j += CREST_OPT_GROUP_THREADS) {
      int k;
      for (k = 0; k < CREST_OPT_GROUP_THREADS; ++k)
        groups[ngroup][k] = records[j + k].cpu;
      ++ngroup;
    }
    ++nnode;
    i = end;
  }
  if (ngroup != CREST_OPT_GROUPS) return EDOM;

  /* Reorder node-contiguous temporary groups into NUMA round-robin order. */
  {
    int ordered[CREST_OPT_GROUPS][CREST_OPT_GROUP_THREADS];
    int out = 0;
    (void)nodes;
    for (round = 0;; ++round) {
      progress = 0;
      for (i = 0; i < nnode; ++i) {
        if (round < node_group_count[i]) {
          int source = node_group_start[i] + node_cursor[i]++;
          for (j = 0; j < CREST_OPT_GROUP_THREADS; ++j)
            ordered[out][j] = groups[source][j];
          ++out;
          progress = 1;
        }
      }
      if (!progress) break;
    }
    if (out != CREST_OPT_GROUPS) return EDOM;
    memcpy(groups, ordered, sizeof(ordered));
  }
  return 0;
}

static int cached_cpu_groups(int groups[CREST_OPT_GROUPS][CREST_OPT_GROUP_THREADS]) {
  static int cache[CREST_OPT_GROUPS][CREST_OPT_GROUP_THREADS];
  static int cache_status = -1;
  static pid_t cache_pid = (pid_t)-1;
  const pid_t current_pid = getpid();
  if (cache_pid != current_pid) {
    cache_status = build_cpu_groups(cache);
    cache_pid = current_pid;
  }
  if (cache_status == 0) memcpy(groups, cache, sizeof(cache));
  return cache_status;
}


struct mtd_l3_domain {
  int start;
  int count;
  int package;
  int node;
  int key;
  int package_pos;
  int node_pos;
  int l3_pos;
};

static int compare_cpu_l3(const void *left, const void *right) {
  const struct cpu_record *a = (const struct cpu_record *)left;
  const struct cpu_record *b = (const struct cpu_record *)right;
#define CMP_FIELD(field) do { if (a->field != b->field) \
  return (a->field > b->field) - (a->field < b->field); } while (0)
  CMP_FIELD(package);
  CMP_FIELD(node);
  CMP_FIELD(l3_first_cpu);
  CMP_FIELD(core);
  CMP_FIELD(cpu);
#undef CMP_FIELD
  return 0;
}

static int compare_domain_order(const void *left, const void *right) {
  const struct mtd_l3_domain *a = (const struct mtd_l3_domain *)left;
  const struct mtd_l3_domain *b = (const struct mtd_l3_domain *)right;
#define CMP_FIELD(field) do { if (a->field != b->field) \
  return (a->field > b->field) - (a->field < b->field); } while (0)
  /* L3 round first, then NUMA position within socket, then socket.  On FIR this
   * yields one L3 from each NUMA/socket before the second L3 in any NUMA. */
  CMP_FIELD(l3_pos);
  CMP_FIELD(node_pos);
  CMP_FIELD(package_pos);
  CMP_FIELD(node);
  CMP_FIELD(key);
#undef CMP_FIELD
  return 0;
}

static int plan_mtd_cpu_records(struct cpu_record *records, int ncpu,
                                int workers, int threads, int *flat) {
  struct mtd_l3_domain domains[CPU_SETSIZE];
  unsigned char used[CPU_SETSIZE];
  int anchor_record[CPU_SETSIZE];
  int ndomain = 0;
  int i, j, w, lane;
  long long required;

  if (records == NULL || flat == NULL || ncpu < 1 || ncpu > CPU_SETSIZE ||
      workers < 1 || threads < 1) return EINVAL;
  required = (long long)workers * (long long)threads;
  if (required < 1 || required > CREST_MTD_MAX_CPUS) return E2BIG;
  if (ncpu < required) return ENOSPC;

  qsort(records, (size_t)ncpu, sizeof(records[0]), compare_cpu_l3);
  for (i = 1; i < ncpu; ++i) {
    if (records[i].package == records[i - 1].package &&
        records[i].die == records[i - 1].die &&
        records[i].core == records[i - 1].core) return ENOTSUP;
  }

  i = 0;
  while (i < ncpu) {
    int end = i + 1;
    while (end < ncpu && records[end].package == records[i].package &&
           records[end].l3_first_cpu == records[i].l3_first_cpu) ++end;
    domains[ndomain].start = i;
    domains[ndomain].count = end - i;
    domains[ndomain].package = records[i].package;
    domains[ndomain].node = records[i].node;
    domains[ndomain].key = records[i].l3_first_cpu;
    ++ndomain;
    i = end;
  }
  if (ndomain < 1) return ENODEV;

  /* Compute stable package/node/L3 positions used solely to balance seed
   * workers.  The later growth pass preserves L3, then NUMA, then socket
   * locality for every worker. */
  for (i = 0; i < ndomain; ++i) {
    int pkg_pos = 0, node_pos = 0, l3_pos = 0;
    int seen_pkg[CPU_SETSIZE], npkg = 0;
    int seen_node[CPU_SETSIZE], nnode = 0;
    for (j = 0; j < i; ++j) {
      int k, present = 0;
      for (k = 0; k < npkg; ++k) if (seen_pkg[k] == domains[j].package) present = 1;
      if (!present) seen_pkg[npkg++] = domains[j].package;
    }
    for (j = 0; j < npkg; ++j)
      if (seen_pkg[j] < domains[i].package) ++pkg_pos;
    for (j = 0; j < i; ++j) {
      int k, present = 0;
      if (domains[j].package != domains[i].package) continue;
      for (k = 0; k < nnode; ++k) if (seen_node[k] == domains[j].node) present = 1;
      if (!present) seen_node[nnode++] = domains[j].node;
    }
    for (j = 0; j < nnode; ++j)
      if (seen_node[j] < domains[i].node) ++node_pos;
    for (j = 0; j < i; ++j)
      if (domains[j].package == domains[i].package &&
          domains[j].node == domains[i].node) ++l3_pos;
    domains[i].package_pos = pkg_pos;
    domains[i].node_pos = node_pos;
    domains[i].l3_pos = l3_pos;
  }
  qsort(domains, (size_t)ndomain, sizeof(domains[0]), compare_domain_order);

  memset(used, 0, sizeof(used));
  for (w = 0; w < workers; ++w) {
    int attempt, chosen = -1;
    for (attempt = 0; attempt < ndomain && chosen < 0; ++attempt) {
      const struct mtd_l3_domain *d = &domains[(w + attempt) % ndomain];
      for (i = d->start; i < d->start + d->count; ++i) {
        if (!used[i]) { chosen = i; break; }
      }
    }
    if (chosen < 0) return ENOSPC;
    used[chosen] = 1;
    anchor_record[w] = chosen;
    flat[w * threads] = records[chosen].cpu;
  }

  /* Grow every worker one lane at a time.  This prevents an early worker from
   * monopolising a cache/NUMA domain and implements the desired hierarchy:
   * same L3 -> same NUMA -> same socket -> anywhere in the allocation. */
  for (lane = 1; lane < threads; ++lane) {
    for (w = 0; w < workers; ++w) {
      const struct cpu_record *anchor = &records[anchor_record[w]];
      int best = -1, best_relation = 99, best_cpu = INT_MAX;
      for (i = 0; i < ncpu; ++i) {
        int relation;
        if (used[i]) continue;
        if (records[i].package == anchor->package &&
            records[i].l3_first_cpu == anchor->l3_first_cpu) relation = 0;
        else if (records[i].node == anchor->node) relation = 1;
        else if (records[i].package == anchor->package) relation = 2;
        else relation = 3;
        if (relation < best_relation ||
            (relation == best_relation && records[i].cpu < best_cpu)) {
          best = i;
          best_relation = relation;
          best_cpu = records[i].cpu;
        }
      }
      if (best < 0) return ENOSPC;
      used[best] = 1;
      flat[w * threads + lane] = records[best].cpu;
    }
  }
  return 0;
}

static int build_mtd_cpu_plan(int workers, int threads, int *flat) {
  struct cpu_record records[CPU_SETSIZE];
  cpu_set_t allowed;
  int ncpu = 0;
  int cpu;
  long long required;

  if (flat == NULL || workers < 1 || threads < 1) return EINVAL;
  required = (long long)workers * (long long)threads;
  if (required < 1 || required > CREST_MTD_MAX_CPUS) return E2BIG;
  CPU_ZERO(&allowed);
  if (sched_getaffinity(0, sizeof(allowed), &allowed) != 0) return errno;
  if (CPU_COUNT(&allowed) < required) return ENOSPC;

  for (cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
    int rc;
    if (!CPU_ISSET(cpu, &allowed)) continue;
    rc = fill_mtd_cpu_record(&records[ncpu], cpu);
    if (rc != 0) return rc;
    ++ncpu;
  }
  return plan_mtd_cpu_records(records, ncpu, workers, threads, flat);
}

static int cached_mtd_cpu_plan(int workers, int threads, int *flat) {
  static int cache[CREST_MTD_MAX_CPUS];
  static int cache_workers = 0, cache_threads = 0, cache_status = -1;
  static pid_t cache_pid = (pid_t)-1;
  const pid_t current_pid = getpid();
  const long long required = (long long)workers * (long long)threads;
  if (flat == NULL || required < 1 || required > CREST_MTD_MAX_CPUS) return EINVAL;
  if (cache_pid != current_pid || cache_workers != workers || cache_threads != threads) {
    cache_status = build_mtd_cpu_plan(workers, threads, cache);
    cache_pid = current_pid;
    cache_workers = workers;
    cache_threads = threads;
  }
  if (cache_status == 0)
    memcpy(flat, cache, (size_t)required * sizeof(cache[0]));
  return cache_status;
}

int crest_mtd_validate_layout(int workers, int threads) {
  int flat[CREST_MTD_MAX_CPUS];
  return cached_mtd_cpu_plan(workers, threads, flat);
}

static int append_cpu_strings(const int *cpus, int threads,
                              char **cpu_list_out, char **places_out) {
  size_t cap, pos_cpu = 0, pos_places = 0;
  char *cpu_list, *places;
  int i;
  if (cpus == NULL || cpu_list_out == NULL || places_out == NULL || threads < 1)
    return EINVAL;
  cap = (size_t)threads * 24u + 1u;
  cpu_list = (char *)calloc(cap, 1);
  places = (char *)calloc(cap, 1);
  if (cpu_list == NULL || places == NULL) {
    free(cpu_list); free(places); return ENOMEM;
  }
  for (i = 0; i < threads; ++i) {
    int n;
    n = snprintf(cpu_list + pos_cpu, cap - pos_cpu, "%s%d", i ? "," : "", cpus[i]);
    if (n < 0 || (size_t)n >= cap - pos_cpu) { free(cpu_list); free(places); return EOVERFLOW; }
    pos_cpu += (size_t)n;
    n = snprintf(places + pos_places, cap - pos_places, "%s{%d}", i ? "," : "", cpus[i]);
    if (n < 0 || (size_t)n >= cap - pos_places) { free(cpu_list); free(places); return EOVERFLOW; }
    pos_places += (size_t)n;
  }
  *cpu_list_out = cpu_list;
  *places_out = places;
  return 0;
}

static int parse_cpu_list_exact(const char *text, int threads, int *cpus) {
  const char *p = text;
  int i, j;
  if (text == NULL || cpus == NULL || threads < 1) return EINVAL;
  for (i = 0; i < threads; ++i) {
    char *end = NULL;
    long value;
    errno = 0;
    value = strtol(p, &end, 10);
    if (errno != 0 || end == p || value < 0 || value >= CPU_SETSIZE) return EINVAL;
    cpus[i] = (int)value;
    for (j = 0; j < i; ++j) if (cpus[j] == cpus[i]) return EEXIST;
    p = end;
    if (i + 1 < threads) {
      if (*p != ',') return EINVAL;
      ++p;
    }
  }
  return *p == '\0' ? 0 : EINVAL;
}

static int bind_memory_to_worker_nodes(const int *cpus, int threads,
                                       int *node_count_out) {
#if defined(SYS_set_mempolicy) && defined(SYS_get_mempolicy)
  enum { BITS_PER_WORD = (int)(8 * sizeof(unsigned long)) };
  unsigned long requested[(CPU_SETSIZE + BITS_PER_WORD - 1) / BITS_PER_WORD];
  unsigned long observed[(CPU_SETSIZE + BITS_PER_WORD - 1) / BITS_PER_WORD];
  const unsigned long maxnode = (unsigned long)(sizeof(requested) * 8);
  int mode = -1, nodes = 0, i;
  if (cpus == NULL || threads < 1) return EINVAL;
  memset(requested, 0, sizeof(requested));
  for (i = 0; i < threads; ++i) {
    const int node = cpu_numa_node(cpus[i]);
    const unsigned long bit = 1UL << (node % BITS_PER_WORD);
    if (node < 0 || (unsigned long)node >= maxnode) return ENODEV;
    if ((requested[node / BITS_PER_WORD] & bit) == 0) {
      requested[node / BITS_PER_WORD] |= bit;
      ++nodes;
    }
  }
  if (syscall(SYS_set_mempolicy, MPOL_BIND, requested, maxnode) != 0) return errno;
  memset(observed, 0, sizeof(observed));
  if (syscall(SYS_get_mempolicy, &mode, observed, maxnode, NULL, 0) != 0) return errno;
  if ((mode & ~MPOL_MODE_FLAGS) != MPOL_BIND) return EINVAL;
  if (memcmp(requested, observed, sizeof(requested)) != 0) return EINVAL;
  if (node_count_out != NULL) *node_count_out = nodes;
  return 0;
#else
  (void)cpus; (void)threads; (void)node_count_out;
  return ENOSYS;
#endif
}

/*
 * Standalone --screen has no preceding process-MTD batch, so it needs a
 * separate, explicit experimental opt-in.  Integrated process-MTD arms the
 * same region only after the historical validated full-node 48x4 batch has
 * passed validation.  Dynamic process-MTD batches do not broaden that optimizer
 * gate.
 * Both calls are serial parent operations before the OpenMP team is created.
 */
static int optimizer_affinity_process_batch_ready;

int crest_optimizer_affinity_set_process_batch_ready(int ready) {
  if (ready != 0 && ready != 1) return EINVAL;
  if (ready != 0 && getenv("CREST_MTD_INTERNAL_WORKER") != NULL) return EPERM;
  optimizer_affinity_process_batch_ready = ready;
  return 0;
}

int crest_optimizer_affinity_resolve(int outer_threads, int inner_threads,
                                     int standalone_allowed,
                                     int *enabled_out) {
  const char *value;
  int explicit_requested = 0;
  int requested;

  if (enabled_out == NULL) return EINVAL;
  *enabled_out = 0;
  if (standalone_allowed != 0 && standalone_allowed != 1) return EINVAL;
  requested = optimizer_affinity_process_batch_ready;
  value = getenv(optimizer_affinity_optin_name);
  if (value != NULL) {
    if (strcmp(value, "1") != 0) return EINVAL;
    explicit_requested = 1;
    requested = 1;
  }
  if (!requested) return 0;
  if (getenv("CREST_MTD_INTERNAL_WORKER") != NULL) return EPERM;
  if (explicit_requested && !standalone_allowed) return EPERM;
  if (outer_threads != CREST_OPT_CPUS || inner_threads != 1) {
    /* Explicit opt-in is a strict request and therefore remains fail-closed.
     * Process-batch readiness is only an implicit opportunity to reuse the
     * historically validated 192x1 affinity plan.  Later optimizer stages may
     * legitimately auto-select smaller topologies as the ensemble shrinks; in
     * that case fall back to CREST's normal optimizer scheduling. */
    if (explicit_requested) return EDOM;
    return 0;
  }
  *enabled_out = 1;
  return 0;
}

/*
 * The process parent must remain unbound until the MTD scheduler has
 * discovered its complete allocation.  The later ensemble optimizer instead
 * uses one persistent OpenMP worker per physical CPU.  GNU OpenMP cannot
 * enable a place list after startup, so bind those already-created workers at
 * the OS level for this one region and restore every worker before it exits.
 *
 * This state is prepared and finished by the serial parent.  It is immutable
 * while the optimizer team is active.  CREST does not support concurrent or
 * nested ensemble-optimizer calls in this exact process-isolated prototype.
 */
struct optimizer_affinity_plan {
  int active;
  int threads;
  int cpu_by_rank[CREST_OPT_CPUS];
  cpu_set_t restore_set;
#ifdef CREST_OPTIMIZER_AFFINITY_TESTING
  int fail_bind_rank;
  int fail_restore_rank;
#endif
};

static struct optimizer_affinity_plan optimizer_plan;

static int store_optimizer_affinity_plan(int threads, const int *cpu_by_rank,
                                         const cpu_set_t *restore_set,
                                         int require_exact_union) {
  cpu_set_t rank_union;
  int i, j;

  if (optimizer_plan.active) return EALREADY;
  if (threads < 1 || threads > CREST_OPT_CPUS || cpu_by_rank == NULL ||
      restore_set == NULL)
    return EINVAL;
  CPU_ZERO(&rank_union);
  for (i = 0; i < threads; ++i) {
    if (cpu_by_rank[i] < 0 || cpu_by_rank[i] >= CPU_SETSIZE ||
        !CPU_ISSET(cpu_by_rank[i], restore_set))
      return EINVAL;
    for (j = 0; j < i; ++j)
      if (cpu_by_rank[i] == cpu_by_rank[j]) return EEXIST;
    CPU_SET(cpu_by_rank[i], &rank_union);
  }
  if (CPU_COUNT(&rank_union) != threads) return EDOM;
  if (require_exact_union && !CPU_EQUAL(&rank_union, restore_set)) return EPERM;

  optimizer_plan.threads = threads;
  memcpy(optimizer_plan.cpu_by_rank, cpu_by_rank,
         (size_t)threads * sizeof(cpu_by_rank[0]));
  optimizer_plan.restore_set = *restore_set;
#ifdef CREST_OPTIMIZER_AFFINITY_TESTING
  optimizer_plan.fail_bind_rank = -1;
  optimizer_plan.fail_restore_rank = -1;
#endif
  optimizer_plan.active = 1;
  return 0;
}

int crest_optimizer_affinity_prepare(int expected_threads) {
  int groups[CREST_OPT_GROUPS][CREST_OPT_GROUP_THREADS];
  int cpu_by_rank[CREST_OPT_CPUS];
  cpu_set_t current;
  int worker, lane, rc;

  if (getenv("CREST_MTD_INTERNAL_WORKER") != NULL) return EPERM;
  if (expected_threads != CREST_OPT_CPUS) return EDOM;
  if (optimizer_plan.active) return EALREADY;
  /* Re-run, rather than trust the scheduler cache: the parent must still own
   * the complete allocation immediately before the optimizer team starts. */
  rc = build_cpu_groups(groups);
  if (rc != 0) return rc;
  CPU_ZERO(&current);
  if (sched_getaffinity(0, sizeof(current), &current) != 0) return errno;
  for (worker = 0; worker < CREST_OPT_GROUPS; ++worker)
    for (lane = 0; lane < CREST_OPT_GROUP_THREADS; ++lane)
      cpu_by_rank[worker * CREST_OPT_GROUP_THREADS + lane] = groups[worker][lane];
  return store_optimizer_affinity_plan(expected_threads, cpu_by_rank, &current, 1);
}

int crest_optimizer_affinity_bind(int rank, int team_size) {
  cpu_set_t before, singleton, observed;
  int cpu;

  if (!optimizer_plan.active) return ENOENT;
  if (team_size != optimizer_plan.threads) return EDOM;
  if (rank < 0 || rank >= optimizer_plan.threads) return ERANGE;
  cpu = optimizer_plan.cpu_by_rank[rank];
  CPU_ZERO(&before);
  if (sched_getaffinity(0, sizeof(before), &before) != 0) return errno;
  if (!CPU_EQUAL(&before, &optimizer_plan.restore_set)) return EPERM;
#ifdef CREST_OPTIMIZER_AFFINITY_TESTING
  if (rank == optimizer_plan.fail_bind_rank) return EIO;
#endif
  CPU_ZERO(&singleton);
  CPU_SET(cpu, &singleton);
  if (sched_setaffinity(0, sizeof(singleton), &singleton) != 0) return errno;
  CPU_ZERO(&observed);
  if (sched_getaffinity(0, sizeof(observed), &observed) != 0) return errno;
  if (CPU_COUNT(&observed) != 1 || !CPU_ISSET(cpu, &observed)) return EIO;
  return 0;
}

int crest_optimizer_affinity_restore(int rank, int team_size) {
  cpu_set_t observed;

  if (!optimizer_plan.active) return ENOENT;
  if (team_size != optimizer_plan.threads) return EDOM;
  if (rank < 0 || rank >= optimizer_plan.threads) return ERANGE;
  if (sched_setaffinity(0, sizeof(optimizer_plan.restore_set),
                        &optimizer_plan.restore_set) != 0)
    return errno;
  CPU_ZERO(&observed);
  if (sched_getaffinity(0, sizeof(observed), &observed) != 0) return errno;
  if (!CPU_EQUAL(&observed, &optimizer_plan.restore_set)) return EIO;
#ifdef CREST_OPTIMIZER_AFFINITY_TESTING
  /* Report the injected failure only after restoring the real test thread. */
  if (rank == optimizer_plan.fail_restore_rank) return EIO;
#endif
  return 0;
}

int crest_optimizer_affinity_finish(void) {
  cpu_set_t observed;

  if (!optimizer_plan.active) return ENOENT;
  CPU_ZERO(&observed);
  if (sched_getaffinity(0, sizeof(observed), &observed) != 0) return errno;
  if (!CPU_EQUAL(&observed, &optimizer_plan.restore_set)) return EPERM;
  optimizer_plan.active = 0;
  optimizer_plan.threads = 0;
  CPU_ZERO(&optimizer_plan.restore_set);
  return 0;
}

#ifdef CREST_OPTIMIZER_AFFINITY_TESTING
/* Unit oracles use the caller's real allowed mask and a small prefix of its
 * CPUs.  This bypasses only the full-node/Slurm discovery precondition; bind,
 * verify, restore, re-entry, and failure coordination use the production
 * routines above.  This symbol is absent from normal builds. */
int crest_optimizer_affinity_test_prepare(int expected_threads,
                                          int fail_bind_rank,
                                          int fail_restore_rank) {
  int cpu_by_rank[CREST_OPT_CPUS];
  cpu_set_t current;
  int cpu, count = 0, rc;

  if (expected_threads < 1 || expected_threads > CREST_OPT_CPUS) return EINVAL;
  CPU_ZERO(&current);
  if (sched_getaffinity(0, sizeof(current), &current) != 0) return errno;
  for (cpu = 0; cpu < CPU_SETSIZE && count < expected_threads; ++cpu)
    if (CPU_ISSET(cpu, &current)) cpu_by_rank[count++] = cpu;
  if (count != expected_threads) return EDOM;
  rc = store_optimizer_affinity_plan(expected_threads, cpu_by_rank, &current, 0);
  if (rc != 0) return rc;
  optimizer_plan.fail_bind_rank = fail_bind_rank;
  optimizer_plan.fail_restore_rank = fail_restore_rank;
  return 0;
}
#endif

int crest_mtd_parent_pid(void) { return (int)getpid(); }

int crest_mtd_current_cpu(void) {
  const int cpu = sched_getcpu();
  if (cpu >= 0) return cpu;
  return errno != 0 ? -errno : -EIO;
}

int crest_mtd_get_thread_affinity(int cpu_capacity, int *cpu_count_out,
                                  int *cpus_out) {
  cpu_set_t observed;
  int cpu, count, stored = 0;

  if (cpu_capacity < 1 || cpu_count_out == NULL || cpus_out == NULL)
    return EINVAL;
  *cpu_count_out = -1;
  for (cpu = 0; cpu < cpu_capacity; ++cpu) cpus_out[cpu] = -1;
  CPU_ZERO(&observed);
  if (sched_getaffinity(0, sizeof(observed), &observed) != 0) return errno;
  count = CPU_COUNT(&observed);
  *cpu_count_out = count;
  for (cpu = 0; cpu < CPU_SETSIZE && stored < cpu_capacity; ++cpu)
    if (CPU_ISSET(cpu, &observed)) cpus_out[stored++] = cpu;
  return 0;
}

static int set_and_verify_thread_affinity(const cpu_set_t *requested) {
  cpu_set_t observed;

  if (requested == NULL || CPU_COUNT(requested) < 1) return EINVAL;
  if (sched_setaffinity(0, sizeof(*requested), requested) != 0) return errno;
  CPU_ZERO(&observed);
  if (sched_getaffinity(0, sizeof(observed), &observed) != 0) return errno;
  return CPU_EQUAL(requested, &observed) ? 0 : EIO;
}

int crest_mtd_secure_directory(const char *path) {
  struct stat info;
  if (path == NULL || *path == '\0') return EINVAL;
  if (lstat(path, &info) != 0) return errno;
  if (!S_ISDIR(info.st_mode) || info.st_uid != geteuid()) return EPERM;
  if (chmod(path, S_IRWXU) != 0) return errno;
  if (lstat(path, &info) != 0) return errno;
  if (!S_ISDIR(info.st_mode) || (info.st_mode & 0777) != 0700) return EPERM;
  return 0;
}

int crest_mtd_spawn_worker(const char *capsule, const char *workdir,
                           int worker_index, int worker_count, int thread_count,
                           int parent_pid, int *pid_out) {
  char executable[PATH_MAX];
  char index_arg[32], worker_count_arg[32], thread_arg[32], parent_arg[32];
  char stdout_path[PATH_MAX], stderr_path[PATH_MAX];
  char *argv[10];
  char **child_env = NULL;
  char *cpu_arg = NULL, *places = NULL;
  ssize_t executable_len;
  posix_spawn_file_actions_t actions;
  pid_t pid;
  int fdout = -1, fderr = -1, rc;
  int flat[CREST_MTD_MAX_CPUS];
  const int *worker_cpus;

  if (capsule == NULL || workdir == NULL || pid_out == NULL) return EINVAL;
  if (worker_count < 1 || thread_count < 1 || worker_index < 1 ||
      worker_index > worker_count) return EINVAL;
  if ((long long)worker_count * (long long)thread_count > CREST_MTD_MAX_CPUS)
    return E2BIG;
  rc = cached_mtd_cpu_plan(worker_count, thread_count, flat);
  if (rc != 0) return rc;
  worker_cpus = &flat[(worker_index - 1) * thread_count];
  rc = append_cpu_strings(worker_cpus, thread_count, &cpu_arg, &places);
  if (rc != 0) return rc;
  rc = build_child_environment(places, thread_count, &child_env);
  if (rc != 0) goto environment_done;

  executable_len = readlink("/proc/self/exe", executable, sizeof(executable) - 1);
  if (executable_len < 1 || executable_len >= (ssize_t)sizeof(executable) - 1) {
    rc = errno ? errno : ENAMETOOLONG;
    goto environment_done;
  }
  executable[executable_len] = '\0';

  if (snprintf(index_arg, sizeof(index_arg), "%d", worker_index) >= (int)sizeof(index_arg) ||
      snprintf(worker_count_arg, sizeof(worker_count_arg), "%d", worker_count) >= (int)sizeof(worker_count_arg) ||
      snprintf(thread_arg, sizeof(thread_arg), "%d", thread_count) >= (int)sizeof(thread_arg) ||
      snprintf(parent_arg, sizeof(parent_arg), "%d", parent_pid) >= (int)sizeof(parent_arg)) {
    rc = EOVERFLOW; goto environment_done;
  }
  if (snprintf(stdout_path, sizeof(stdout_path), "%s/worker.stdout", workdir) >=
      (int)sizeof(stdout_path)) { rc = ENAMETOOLONG; goto environment_done; }
  if (snprintf(stderr_path, sizeof(stderr_path), "%s/worker.stderr", workdir) >=
      (int)sizeof(stderr_path)) { rc = ENAMETOOLONG; goto environment_done; }

  fdout = open(stdout_path, O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0600);
  if (fdout < 0) { rc = errno; goto environment_done; }
  fderr = open(stderr_path, O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0600);
  if (fderr < 0) { rc = errno; close(fdout); fdout = -1; goto environment_done; }
  rc = posix_spawn_file_actions_init(&actions);
  if (rc != 0) goto done;
  rc = posix_spawn_file_actions_adddup2(&actions, fdout, STDOUT_FILENO);
  if (rc == 0) rc = posix_spawn_file_actions_adddup2(&actions, fderr, STDERR_FILENO);
  if (rc == 0) rc = posix_spawn_file_actions_addclose(&actions, fdout);
  if (rc == 0) rc = posix_spawn_file_actions_addclose(&actions, fderr);
  if (rc != 0) { posix_spawn_file_actions_destroy(&actions); goto done; }

  argv[0] = executable;
  argv[1] = (char *)"--crest-internal-mtd-worker-v1";
  argv[2] = (char *)capsule;
  argv[3] = (char *)workdir;
  argv[4] = index_arg;
  argv[5] = worker_count_arg;
  argv[6] = thread_arg;
  argv[7] = parent_arg;
  argv[8] = cpu_arg;
  argv[9] = NULL;
  rc = posix_spawn(&pid, executable, &actions, NULL, argv, child_env);
  posix_spawn_file_actions_destroy(&actions);
  if (rc == 0) *pid_out = (int)pid;

done:
  if (fdout >= 0) close(fdout);
  if (fderr >= 0) close(fderr);
environment_done:
  free_child_environment(child_env);
  free(cpu_arg);
  free(places);
  return rc;
}

int crest_mtd_prepare_worker(const char *workdir, int worker_index,
                             int worker_count, int thread_count,
                             int expected_parent_pid, const char *cpu_list) {
  cpu_set_t primary_set;
  int cpus[CREST_MTD_MAX_CPUS];
  char *expected_places = NULL, *unused_cpu_list = NULL;
  const char *actual_places;
  char thread_text[32];
  int node_count = 0, rc;

  if (workdir == NULL || cpu_list == NULL || worker_count < 1 || thread_count < 1 ||
      worker_index < 1 || worker_index > worker_count || thread_count > CREST_MTD_MAX_CPUS)
    return EINVAL;
  if ((int)getppid() != expected_parent_pid) return ESRCH;
  if (prctl(PR_SET_PDEATHSIG, SIGKILL) != 0) return errno;
  if ((int)getppid() != expected_parent_pid) return ESRCH;

  rc = parse_cpu_list_exact(cpu_list, thread_count, cpus);
  if (rc != 0) return rc;
  rc = append_cpu_strings(cpus, thread_count, &unused_cpu_list, &expected_places);
  free(unused_cpu_list);
  if (rc != 0) return rc;
  actual_places = getenv("OMP_PLACES");
  if (actual_places == NULL || strcmp(actual_places, expected_places) != 0) {
    free(expected_places);
    return EINVAL;
  }
  free(expected_places);
  if (getenv("CREST_MTD_INTERNAL_WORKER") == NULL) return EPERM;

  /* OMP_PLACES was installed before exec and describes singleton Linux CPU
   * IDs.  Never widen the libgomp primary thread away from place zero.  The
   * first application OpenMP team proves the remaining assigned CPUs. */
  CPU_ZERO(&primary_set);
  CPU_SET(cpus[0], &primary_set);
  rc = set_and_verify_thread_affinity(&primary_set);
  if (rc != 0) return rc;
  rc = bind_memory_to_worker_nodes(cpus, thread_count, &node_count);
  if (rc != 0) return rc;
  if (chdir(workdir) != 0) return errno;

  if (dprintf(STDOUT_FILENO,
              "Process MTD worker %d/%d CPU set=%s NUMA nodes=%d threads=%d\n",
              worker_index, worker_count, cpu_list, node_count, thread_count) < 0)
    return errno ? errno : EIO;

  if (snprintf(thread_text, sizeof(thread_text), "%d", thread_count) >=
      (int)sizeof(thread_text)) return EOVERFLOW;
  if (setenv("OMP_NUM_THREADS", thread_text, 1) != 0) return errno;
  if (setenv("OMP_THREAD_LIMIT", thread_text, 1) != 0) return errno;
  if (setenv("OMP_DYNAMIC", "FALSE", 1) != 0) return errno;
  if (setenv("OMP_NESTED", "FALSE", 1) != 0) return errno;
  if (setenv("OMP_MAX_ACTIVE_LEVELS", "1", 1) != 0) return errno;
  if (setenv("OMP_PROC_BIND", "close", 1) != 0) return errno;
  if (setenv("OPENBLAS_NUM_THREADS", "1", 1) != 0) return errno;
  if (setenv("GOTO_NUM_THREADS", "1", 1) != 0) return errno;
  if (setenv("MKL_NUM_THREADS", "1", 1) != 0) return errno;
  if (setenv("BLIS_NUM_THREADS", "1", 1) != 0) return errno;
  return 0;
}

int crest_mtd_wait_worker(int pid_value, int *exit_code, int *term_signal) {
  int status;
  pid_t rc;
  if (pid_value <= 0 || exit_code == NULL || term_signal == NULL) return EINVAL;
  *exit_code = -1;
  *term_signal = 0;
  do {
    rc = waitpid((pid_t)pid_value, &status, 0);
  } while (rc < 0 && errno == EINTR);
  if (rc < 0) return errno;
  *exit_code = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
  *term_signal = WIFSIGNALED(status) ? WTERMSIG(status) : 0;
  return 0;
}

int crest_mtd_kill_worker(int pid_value) {
  if (pid_value <= 0) return EINVAL;
  if (kill((pid_t)pid_value, SIGTERM) != 0 && errno != ESRCH) return errno;
  return 0;
}

int crest_mtd_reap_worker_bounded(int pid_value, int grace_milliseconds,
                                  int *exit_code, int *term_signal) {
  struct timespec pause_time = {0, 10000000L};
  int status, elapsed = 0, kill_elapsed = 0;
  pid_t rc;
  if (pid_value <= 0 || grace_milliseconds < 0 || exit_code == NULL ||
      term_signal == NULL)
    return EINVAL;
  *exit_code = -1;
  *term_signal = 0;
  while (elapsed <= grace_milliseconds) {
    rc = waitpid((pid_t)pid_value, &status, WNOHANG);
    if (rc == (pid_t)pid_value) {
      *exit_code = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
      *term_signal = WIFSIGNALED(status) ? WTERMSIG(status) : 0;
      return 0;
    }
    if (rc < 0) return errno;
    nanosleep(&pause_time, NULL);
    elapsed += 10;
  }
  if (kill((pid_t)pid_value, SIGKILL) != 0 && errno != ESRCH) return errno;
  while (kill_elapsed <= 500) {
    rc = waitpid((pid_t)pid_value, &status, WNOHANG);
    if (rc == (pid_t)pid_value) {
      *exit_code = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
      *term_signal = WIFSIGNALED(status) ? WTERMSIG(status) : 0;
      return 0;
    }
    if (rc < 0) return errno;
    nanosleep(&pause_time, NULL);
    kill_elapsed += 10;
  }
  return ETIMEDOUT;
}
