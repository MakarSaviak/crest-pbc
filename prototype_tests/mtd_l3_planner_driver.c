/* Local synthetic topology oracle for the dynamic process-MTD scheduler.
 * Compile this driver as one translation unit; it includes the production C
 * runtime so static planner helpers are exercised directly. */
#include "../src/mtd_process_runtime.c"

static void make_fir_records(struct cpu_record records[192]) {
  int cpu;
  for (cpu = 0; cpu < 192; ++cpu) {
    records[cpu].cpu = cpu;
    records[cpu].package = cpu / 96;
    records[cpu].node = cpu / 24;
    records[cpu].die = records[cpu].package;
    records[cpu].core = cpu % 96;
    records[cpu].l3_first_cpu = (cpu / 8) * 8;
  }
}

static int verify_plan(int workers, int threads, int expected_l3s,
                       int expected_nodes, int expected_active_per_used_l3,
                       int expected_l3_per_worker, int expected_nodes_per_worker) {
  struct cpu_record records[192];
  int flat[CREST_MTD_MAX_CPUS];
  int global_cpu[192] = {0}, l3_load[24] = {0};
  int node_used[8] = {0}, package_used[2] = {0};
  int w, lane, i, l3s = 0, nodes = 0, packages = 0;
  int rc;

  make_fir_records(records);
  rc = plan_mtd_cpu_records(records, 192, workers, threads, flat);
  if (rc != 0) {
    fprintf(stderr, "plan failed for N=%d K=%d: %d\n", workers, threads, rc);
    return 1;
  }
  for (w = 0; w < workers; ++w) {
    int worker_l3[24] = {0}, worker_node[8] = {0};
    int worker_l3s = 0, worker_nodes = 0;
    for (lane = 0; lane < threads; ++lane) {
      const int cpu = flat[w * threads + lane];
      const int l3 = cpu / 8;
      const int node = cpu / 24;
      const int package = cpu / 96;
      if (cpu < 0 || cpu >= 192 || global_cpu[cpu]) return 2;
      global_cpu[cpu] = 1;
      ++l3_load[l3];
      node_used[node] = 1;
      package_used[package] = 1;
      worker_l3[l3] = 1;
      worker_node[node] = 1;
    }
    for (i = 0; i < 24; ++i) worker_l3s += worker_l3[i];
    for (i = 0; i < 8; ++i) worker_nodes += worker_node[i];
    if (worker_l3s != expected_l3_per_worker ||
        worker_nodes != expected_nodes_per_worker) {
      fprintf(stderr, "worker locality failed for N=%d K=%d worker=%d: L3=%d NUMA=%d\n",
              workers, threads, w + 1, worker_l3s, worker_nodes);
      return 3;
    }
  }
  for (i = 0; i < 24; ++i) {
    if (l3_load[i] != 0) {
      ++l3s;
      if (l3_load[i] != expected_active_per_used_l3) {
        fprintf(stderr, "L3 load failed for N=%d K=%d L3=%d: %d\n",
                workers, threads, i, l3_load[i]);
        return 4;
      }
    }
  }
  for (i = 0; i < 8; ++i) nodes += node_used[i];
  for (i = 0; i < 2; ++i) packages += package_used[i];
  if (l3s != expected_l3s || nodes != expected_nodes || packages != 2) {
    fprintf(stderr, "spread failed for N=%d K=%d: L3=%d NUMA=%d sockets=%d\n",
            workers, threads, l3s, nodes, packages);
    return 5;
  }
  printf("FIR_PLAN_PASS N=%d K=%d L3=%d NUMA=%d active_per_L3=%d\n",
         workers, threads, l3s, nodes, expected_active_per_used_l3);
  return 0;
}

static const char *find_env(char **envp, const char *name) {
  size_t i, n = strlen(name);
  for (i = 0; envp != NULL && envp[i] != NULL; ++i)
    if (strncmp(envp[i], name, n) == 0 && envp[i][n] == '=') return envp[i] + n + 1;
  return NULL;
}

static int verify_dynamic_child_environment(void) {
  char **envp = NULL;
  int rc = build_child_environment("{0},{1},{2},{3},{4},{5},{6},{7}", 8, &envp);
  const char *omp_threads, *thread_limit, *blas_threads;
  if (rc != 0) return rc;
  omp_threads = find_env(envp, "OMP_NUM_THREADS");
  thread_limit = find_env(envp, "OMP_THREAD_LIMIT");
  blas_threads = find_env(envp, "OPENBLAS_NUM_THREADS");
  rc = (!omp_threads || strcmp(omp_threads, "8") != 0 ||
        !thread_limit || strcmp(thread_limit, "8") != 0 ||
        !blas_threads || strcmp(blas_threads, "1") != 0) ? EDOM : 0;
  free_child_environment(envp);
  if (rc == 0) puts("DYNAMIC_CHILD_ENV_PASS K=8 BLAS=1");
  return rc;
}

int main(void) {
  struct cpu_record records[192];
  int flat[CREST_MTD_MAX_CPUS];
  int rc;

  if (verify_plan(48, 1, 24, 8, 2, 1, 1)) return 11;
  if (verify_plan(48, 2, 24, 8, 4, 1, 1)) return 12;
  if (verify_plan(48, 3, 24, 8, 6, 1, 1)) return 13;
  if (verify_plan(48, 4, 24, 8, 8, 1, 1)) return 14;
  if (verify_plan(24, 8, 24, 8, 8, 1, 1)) return 15;
  if (verify_plan(6, 1, 6, 6, 1, 1, 1)) return 16;
  if (verify_plan(8, 16, 16, 8, 8, 2, 1)) return 17;
  if (verify_dynamic_child_environment()) return 18;

  make_fir_records(records);
  rc = plan_mtd_cpu_records(records, 192, 49, 4, flat);
  if (rc != ENOSPC) {
    fprintf(stderr, "invalid-capacity gate failed: rc=%d expected=%d\n", rc, ENOSPC);
    return 19;
  }
  puts("FIR_CAPACITY_GATE_PASS");
  return 0;
}
