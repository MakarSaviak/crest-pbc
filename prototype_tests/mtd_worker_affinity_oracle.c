#define _GNU_SOURCE

#include <errno.h>
#include <omp.h>
#include <sched.h>
#include <spawn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

extern char **environ;

int crest_mtd_prepare_worker(const char *workdir, int worker_index,
                             int worker_count, int thread_count,
                             int expected_parent_pid, const char *cpu_list);
#ifdef CREST_MTD_AFFINITY_PROBE_AVAILABLE
int crest_mtd_get_thread_affinity(int cpu_capacity, int *cpu_count_out,
                                  int *cpus_out);
#endif

static void require(int condition, const char *message) {
  if (!condition) {
    fprintf(stderr, "MTD_WORKER_AFFINITY_ORACLE_FAIL: %s\n", message);
    exit(EXIT_FAILURE);
  }
}

static int read_thread_affinity(int *count_out, int *sole_cpu_out) {
  cpu_set_t observed;
  int cpu, count;

  if (count_out == NULL || sole_cpu_out == NULL) return EINVAL;
  *count_out = -1;
  *sole_cpu_out = -1;
  CPU_ZERO(&observed);
  if (sched_getaffinity(0, sizeof(observed), &observed) != 0) return errno;
  count = CPU_COUNT(&observed);
  *count_out = count;
  if (count == 1) {
    for (cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
      if (CPU_ISSET(cpu, &observed)) {
        *sole_cpu_out = cpu;
        break;
      }
    }
  }
  return 0;
}

static int child_main(int argc, char **argv) {
  int cpu[2], parent_pid, rc, count, sole_cpu;
#ifdef CREST_MTD_AFFINITY_PROBE_AVAILABLE
  int production_count = -1, production_cpus[2] = {-1, -1};
#endif
  int place_cpus[2] = {-1, -1};
  int team_size = 0;
  int thread_places[2] = {-1, -1};
  int thread_cpus[2] = {-1, -1};
  int affinity_status[2] = {-1, -1};
  int affinity_counts[2] = {-1, -1};
  int affinity_cpus[2] = {-1, -1};
  char *end = NULL;
  const char *workdir, *cpu_list;

  require(argc == 7, "invalid child argument count");
  errno = 0;
  cpu[0] = (int)strtol(argv[2], &end, 10);
  require(errno == 0 && end != argv[2] && *end == '\0', "invalid cpu0");
  errno = 0;
  cpu[1] = (int)strtol(argv[3], &end, 10);
  require(errno == 0 && end != argv[3] && *end == '\0', "invalid cpu1");
  errno = 0;
  parent_pid = (int)strtol(argv[4], &end, 10);
  require(errno == 0 && end != argv[4] && *end == '\0', "invalid parent PID");
  workdir = argv[5];
  cpu_list = argv[6];

  rc = crest_mtd_prepare_worker(workdir, 1, 1, 2, parent_pid, cpu_list);
  require(rc == 0, "production worker preparation failed");

  /* This is the deterministic regression point.  The old implementation
   * leaves {cpu0,cpu1} here; the fixed implementation must leave {cpu0}. */
  rc = read_thread_affinity(&count, &sole_cpu);
  require(rc == 0, "cannot read primary affinity after worker preparation");
  require(count == 1, "primary affinity is not singleton after worker preparation");
  require(sole_cpu == cpu[0], "primary affinity is not worker place zero");
#ifdef CREST_MTD_AFFINITY_PROBE_AVAILABLE
  rc = crest_mtd_get_thread_affinity(2, &production_count, production_cpus);
  require(rc == 0, "production primary-affinity probe failed");
  require(production_count == count && production_cpus[0] == sole_cpu,
          "production primary-affinity probe disagrees with direct syscall");
#endif

  omp_set_dynamic(0);
  omp_set_max_active_levels(1);
  omp_set_num_threads(2);
  require(omp_get_num_places() == 2, "libgomp did not expose two places");
  require(omp_get_place_num() == 0, "primary thread is not on OpenMP place zero");
  require(omp_get_place_num_procs(0) == 1 && omp_get_place_num_procs(1) == 1,
          "an OpenMP place is not singleton");
  omp_get_place_proc_ids(0, &place_cpus[0]);
  omp_get_place_proc_ids(1, &place_cpus[1]);
  require(place_cpus[0] == cpu[0] && place_cpus[1] == cpu[1],
          "OpenMP place CPU IDs differ from the absolute Linux CPU IDs");

#pragma omp parallel num_threads(2) default(none)                              \
    shared(team_size, thread_places, thread_cpus, affinity_status,             \
           affinity_counts, affinity_cpus)
  {
    const int tid = omp_get_thread_num();
#pragma omp single
    team_size = omp_get_num_threads();
    if (tid >= 0 && tid < 2) {
      thread_places[tid] = omp_get_place_num();
      thread_cpus[tid] = sched_getcpu();
      affinity_status[tid] =
          read_thread_affinity(&affinity_counts[tid], &affinity_cpus[tid]);
    }
  }

  require(team_size == 2, "OpenMP team size is not two");
  for (int i = 0; i < 2; ++i) {
    int declared_cpu;
    require(thread_places[i] >= 0 && thread_places[i] < 2,
            "thread reported an out-of-range OpenMP place");
    declared_cpu = place_cpus[thread_places[i]];
    require(affinity_status[i] == 0, "per-thread sched_getaffinity failed");
    require(affinity_counts[i] == 1, "an OpenMP thread mask is not singleton");
    require(affinity_cpus[i] == declared_cpu,
            "an OpenMP thread mask differs from its declared singleton place");
    require(thread_cpus[i] == declared_cpu,
            "sched_getcpu differs from the declared singleton place");
  }
  require(thread_places[0] != thread_places[1],
          "two OpenMP threads share one place");
  require(thread_cpus[0] != thread_cpus[1],
          "two OpenMP threads execute on one CPU");
  require((thread_cpus[0] == cpu[0] || thread_cpus[0] == cpu[1]) &&
              (thread_cpus[1] == cpu[0] || thread_cpus[1] == cpu[1]),
          "OpenMP team escaped the requested two-CPU set");

  puts("MTD_WORKER_AFFINITY_CHILD_ORACLE_PASS");
  return EXIT_SUCCESS;
}

static void set_child_environment(int cpu0, int cpu1) {
  char places[64];

  require(snprintf(places, sizeof(places), "{%d},{%d}", cpu0, cpu1) > 0,
          "cannot construct OMP_PLACES");
  require(setenv("OMP_NUM_THREADS", "2", 1) == 0, "cannot set OMP_NUM_THREADS");
  require(setenv("OMP_DYNAMIC", "FALSE", 1) == 0, "cannot set OMP_DYNAMIC");
  require(setenv("OMP_NESTED", "FALSE", 1) == 0, "cannot set OMP_NESTED");
  require(setenv("OMP_MAX_ACTIVE_LEVELS", "1", 1) == 0,
          "cannot set OMP_MAX_ACTIVE_LEVELS");
  require(setenv("OMP_THREAD_LIMIT", "2", 1) == 0, "cannot set OMP_THREAD_LIMIT");
  require(setenv("OMP_PROC_BIND", "close", 1) == 0, "cannot set OMP_PROC_BIND");
  require(setenv("OMP_PLACES", places, 1) == 0, "cannot set OMP_PLACES");
  require(setenv("CREST_MTD_INTERNAL_WORKER", "1", 1) == 0,
          "cannot mark internal worker");
  require(unsetenv("GOMP_CPU_AFFINITY") == 0, "cannot clear GOMP_CPU_AFFINITY");
  require(unsetenv("CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY") == 0,
          "cannot clear optimizer-affinity opt-in");
}

int main(int argc, char **argv) {
  cpu_set_t allowed;
  int cpus[2] = {-1, -1};
  int cpu, child_success, fallback_cpu = -1, status, rc;
  pid_t child, waited;
  char executable[4096], cpu0_text[32], cpu1_text[32], parent_text[32];
  char cpu_list[64], workdir_template[] = "/tmp/crest-mtd-affinity-oracle.XXXXXX";
  char *workdir;
  char *child_argv[8];
  ssize_t executable_len;

  if (argc > 1 && strcmp(argv[1], "--child") == 0)
    return child_main(argc, argv);
  require(argc == 1, "unexpected driver arguments");

  CPU_ZERO(&allowed);
  require(sched_getaffinity(0, sizeof(allowed), &allowed) == 0,
          "cannot read driver CPU allowance");
  for (cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
    if (!CPU_ISSET(cpu, &allowed)) continue;
    if (cpus[0] < 0) {
      cpus[0] = cpu;
    } else {
      if (fallback_cpu < 0) fallback_cpu = cpu;
      if (cpu > cpus[0] + 1) {
        cpus[1] = cpu;
        break;
      }
    }
  }
  if (cpus[1] < 0) cpus[1] = fallback_cpu;
  if (cpus[0] < 0 || cpus[1] < 0) {
    puts("MTD_WORKER_AFFINITY_ORACLE_SKIP requires_at_least_two_allowed_cpus");
    return 77;
  }

  executable_len = readlink("/proc/self/exe", executable, sizeof(executable) - 1);
  require(executable_len > 0 && executable_len < (ssize_t)sizeof(executable) - 1,
          "cannot resolve oracle executable");
  executable[executable_len] = '\0';
  workdir = mkdtemp(workdir_template);
  require(workdir != NULL, "cannot create private worker directory");

  require(snprintf(cpu0_text, sizeof(cpu0_text), "%d", cpus[0]) > 0,
          "cannot format cpu0");
  require(snprintf(cpu1_text, sizeof(cpu1_text), "%d", cpus[1]) > 0,
          "cannot format cpu1");
  require(snprintf(parent_text, sizeof(parent_text), "%d", (int)getpid()) > 0,
          "cannot format parent PID");
  require(snprintf(cpu_list, sizeof(cpu_list), "%d,%d", cpus[0], cpus[1]) > 0,
          "cannot format CPU list");
  set_child_environment(cpus[0], cpus[1]);

  child_argv[0] = executable;
  child_argv[1] = (char *)"--child";
  child_argv[2] = cpu0_text;
  child_argv[3] = cpu1_text;
  child_argv[4] = parent_text;
  child_argv[5] = workdir;
  child_argv[6] = cpu_list;
  child_argv[7] = NULL;
  rc = posix_spawn(&child, executable, NULL, NULL, child_argv, environ);
  require(rc == 0, "posix_spawn failed");
  do {
    waited = waitpid(child, &status, 0);
  } while (waited < 0 && errno == EINTR);
  require(waited == child, "waitpid failed");
  child_success = WIFEXITED(status) && WEXITSTATUS(status) == 0;
  require(rmdir(workdir) == 0, "cannot remove empty oracle worker directory");
  require(child_success, "affinity oracle child did not exit successfully");

  printf("MTD_WORKER_AFFINITY_ORACLE_PASS cpus=%d,%d\n", cpus[0], cpus[1]);
  return EXIT_SUCCESS;
}
