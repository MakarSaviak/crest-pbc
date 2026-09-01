#define _GNU_SOURCE

#include <errno.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

int crest_optimizer_affinity_bind(int rank, int team_size);
int crest_optimizer_affinity_prepare(int expected_threads);
int crest_optimizer_affinity_restore(int rank, int team_size);
int crest_optimizer_affinity_finish(void);
int crest_optimizer_affinity_resolve(int outer_threads, int inner_threads,
                                     int standalone_allowed,
                                     int *enabled_out);
int crest_optimizer_affinity_set_process_batch_ready(int ready);
int crest_optimizer_affinity_test_child_environment(void);
int crest_optimizer_affinity_test_prepare(int expected_threads,
                                          int fail_bind_rank,
                                          int fail_restore_rank);

enum { TEST_THREADS = 4 };

struct worker_result {
  int rank;
  int bind_team_size;
  int bind_status;
  int restore_status;
  pthread_barrier_t *bound_barrier;
};

static void require(int condition, const char *message) {
  if (!condition) {
    fprintf(stderr, "OPTIMIZER_AFFINITY_C_FAIL: %s\n", message);
    exit(1);
  }
}

static void wait_at_barrier(pthread_barrier_t *barrier) {
  const int rc = pthread_barrier_wait(barrier);
  require(rc == 0 || rc == PTHREAD_BARRIER_SERIAL_THREAD,
          "pthread barrier failed");
}

static void run_gate_cases(void) {
  static const int sweep_workers[] = {1, 6, 24, 47, 48, 96, 191, 192, 193};
  static const int sweep_inner[] = {1, 2, 4, 8, 16};
  int enabled = -1;
  size_t i, j;

  require(unsetenv("CREST_MTD_INTERNAL_WORKER") == 0,
          "cannot clear inherited worker sentinel");
  require(unsetenv("CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY") == 0,
          "cannot clear inherited optimizer opt-in");
  require(crest_optimizer_affinity_set_process_batch_ready(0) == 0,
          "cannot clear process-batch readiness");
  require(crest_optimizer_affinity_resolve(4, 1, 0, &enabled) == 0 &&
              enabled == 0,
          "absent opt-in changed default optimizer behavior");
  puts("C_GATE_CASE_PASS label=absent_optin_unchanged enabled=0");

  require(setenv("CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY", "1", 1) == 0,
          "cannot establish standalone screen opt-in");
  require(crest_optimizer_affinity_test_child_environment() == 0,
          "constructed process-child environment inherited optimizer opt-in");
  puts("C_GATE_CASE_PASS label=child_environment_strips_parent_optin");

  enabled = -1;
  require(crest_optimizer_affinity_resolve(192, 1, 1, &enabled) == 0 &&
              enabled == 1,
          "standalone screen opt-in did not arm optimizer affinity");
  puts("C_GATE_CASE_PASS label=screen_optin_arms enabled=1");

  enabled = -1;
  require(crest_optimizer_affinity_resolve(191, 1, 1, &enabled) == EDOM &&
              enabled == 0,
          "explicit opt-in accepted a non-192x1 optimizer team");
  puts("C_GATE_CASE_PASS label=outer_191_rejected");

  enabled = -1;
  require(crest_optimizer_affinity_resolve(192, 2, 1, &enabled) == EDOM &&
              enabled == 0,
          "explicit opt-in accepted a 192x2 optimizer team");
  puts("C_GATE_CASE_PASS label=inner_2_rejected");

  enabled = -1;
  require(crest_optimizer_affinity_resolve(192, 1, 0, &enabled) == EPERM &&
              enabled == 0,
          "non-screen mode accepted standalone optimizer opt-in");
  puts("C_GATE_CASE_PASS label=non_screen_optin_rejected");

  require(setenv("CREST_MTD_INTERNAL_WORKER", "1", 1) == 0,
          "cannot establish worker sentinel");
  enabled = -1;
  require(crest_optimizer_affinity_resolve(192, 1, 1, &enabled) == EPERM &&
              enabled == 0,
          "process worker was allowed to resolve optimizer affinity");
  require(crest_optimizer_affinity_set_process_batch_ready(1) == EPERM,
          "process worker was allowed to arm batch readiness");
  require(crest_optimizer_affinity_prepare(192) == EPERM,
          "process worker was allowed to prepare optimizer affinity");
  require(unsetenv("CREST_MTD_INTERNAL_WORKER") == 0,
          "cannot clear worker sentinel");
  require(unsetenv("CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY") == 0,
          "cannot clear standalone screen opt-in");
  puts("C_GATE_CASE_PASS label=process_worker_rejected");

  require(crest_optimizer_affinity_set_process_batch_ready(1) == 0,
          "validated process-batch readiness could not be armed");
  enabled = -1;
  require(crest_optimizer_affinity_resolve(192, 1, 0, &enabled) == 0 &&
              enabled == 1,
          "successful process batch did not arm optimizer affinity");
  puts("C_GATE_CASE_PASS label=integrated_ready_non_screen_arms enabled=1");

  enabled = -1;
  require(crest_optimizer_affinity_resolve(47, 4, 0, &enabled) == 0 &&
              enabled == 0,
          "implicit process-batch readiness rejected a smaller auto optimizer topology");
  puts("C_GATE_CASE_PASS label=integrated_ready_47x4_falls_back enabled=0");

  enabled = -1;
  require(crest_optimizer_affinity_resolve(192, 1, 0, &enabled) == 0 &&
              enabled == 1,
          "implicit process-batch readiness was consumed by smaller-topology fallback");
  puts("C_GATE_CASE_PASS label=integrated_ready_persists_after_fallback enabled=1");

  require(setenv("CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY", "1", 1) == 0,
          "cannot establish explicit opt-in while process readiness is armed");
  enabled = -1;
  require(crest_optimizer_affinity_resolve(47, 4, 1, &enabled) == EDOM &&
              enabled == 0,
          "explicit optimizer-affinity opt-in lost strict topology rejection");
  require(unsetenv("CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY") == 0,
          "cannot clear explicit opt-in after precedence check");
  puts("C_GATE_CASE_PASS label=explicit_optin_precedence_47x4_rejected");

  for (i = 0; i < sizeof(sweep_workers) / sizeof(sweep_workers[0]); ++i) {
    for (j = 0; j < sizeof(sweep_inner) / sizeof(sweep_inner[0]); ++j) {
      const int exact = sweep_workers[i] == 192 && sweep_inner[j] == 1;
      enabled = -1;
      require(crest_optimizer_affinity_resolve(
                  sweep_workers[i], sweep_inner[j], 0, &enabled) == 0 &&
                  enabled == exact,
              "implicit process-batch topology sweep violated fallback semantics");
    }
  }
  puts("C_GATE_CASE_PASS label=integrated_ready_topology_sweep");

  require(setenv("CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY", "1", 1) == 0,
          "cannot establish explicit opt-in for topology sweep");
  for (i = 0; i < sizeof(sweep_workers) / sizeof(sweep_workers[0]); ++i) {
    for (j = 0; j < sizeof(sweep_inner) / sizeof(sweep_inner[0]); ++j) {
      const int exact = sweep_workers[i] == 192 && sweep_inner[j] == 1;
      const int expected_status = exact ? 0 : EDOM;
      enabled = -1;
      require(crest_optimizer_affinity_resolve(
                  sweep_workers[i], sweep_inner[j], 1, &enabled) == expected_status &&
                  enabled == exact,
              "explicit optimizer-affinity topology sweep lost strict semantics");
    }
  }
  require(unsetenv("CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY") == 0,
          "cannot clear explicit opt-in after topology sweep");
  puts("C_GATE_CASE_PASS label=explicit_topology_sweep_strict");

  require(crest_optimizer_affinity_set_process_batch_ready(0) == 0,
          "failed process-batch path could not clear readiness");
  enabled = -1;
  require(crest_optimizer_affinity_resolve(192, 1, 0, &enabled) == 0 &&
              enabled == 0,
          "failed process batch left optimizer affinity armed");
  puts("C_GATE_CASE_PASS label=failed_reset_batch_does_not_arm enabled=0");

  require(setenv("CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY", "true", 1) == 0,
          "cannot establish invalid optimizer opt-in");
  enabled = -1;
  require(crest_optimizer_affinity_resolve(192, 1, 1, &enabled) == EINVAL &&
              enabled == 0,
          "invalid optimizer opt-in did not fail closed");
  require(unsetenv("CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY") == 0,
          "cannot clear invalid optimizer opt-in");
  enabled = -1;
  require(crest_optimizer_affinity_resolve(192, 1, 2, &enabled) == EINVAL &&
              enabled == 0,
          "invalid standalone-allowed bit was accepted");
  require(crest_optimizer_affinity_set_process_batch_ready(2) == EINVAL,
          "invalid process-batch state was accepted");
  puts("C_GATE_CASE_PASS label=invalid_optin_rejected");
}

static void *run_worker(void *opaque) {
  struct worker_result *result = (struct worker_result *)opaque;
  result->bind_status = crest_optimizer_affinity_bind(
      result->rank, result->bind_team_size);
  wait_at_barrier(result->bound_barrier);
  result->restore_status = crest_optimizer_affinity_restore(
      result->rank, TEST_THREADS);
  return NULL;
}

static void run_region(const char *label, int fail_bind_rank,
                       int fail_restore_rank, int bind_team_size,
                       int expected_bind_errors, int expected_restore_errors,
                       int check_reentry) {
  pthread_t threads[TEST_THREADS];
  struct worker_result results[TEST_THREADS];
  pthread_barrier_t bound_barrier;
  int bind_errors = 0, restore_errors = 0;
  int i, rc;

  rc = crest_optimizer_affinity_test_prepare(
      TEST_THREADS, fail_bind_rank, fail_restore_rank);
  require(rc == 0, "test prepare failed");
  if (check_reentry) {
    rc = crest_optimizer_affinity_test_prepare(TEST_THREADS, -1, -1);
    require(rc == EALREADY, "active plan did not reject re-entry");
  }
  require(pthread_barrier_init(&bound_barrier, NULL, TEST_THREADS) == 0,
          "pthread barrier initialization failed");
  for (i = 0; i < TEST_THREADS; ++i) {
    results[i].rank = i;
    results[i].bind_team_size = bind_team_size;
    results[i].bind_status = -1;
    results[i].restore_status = -1;
    results[i].bound_barrier = &bound_barrier;
    require(pthread_create(&threads[i], NULL, run_worker, &results[i]) == 0,
            "pthread creation failed");
  }
  for (i = 0; i < TEST_THREADS; ++i)
    require(pthread_join(threads[i], NULL) == 0, "pthread join failed");
  require(pthread_barrier_destroy(&bound_barrier) == 0,
          "pthread barrier destruction failed");
  for (i = 0; i < TEST_THREADS; ++i) {
    if (results[i].bind_status != 0) ++bind_errors;
    if (results[i].restore_status != 0) ++restore_errors;
  }
  require(bind_errors == expected_bind_errors,
          "unexpected bind failure count");
  require(restore_errors == expected_restore_errors,
          "unexpected restore failure count");
  require(crest_optimizer_affinity_finish() == 0,
          "finish did not observe the restored parent mask");
  printf("C_CASE_PASS label=%s bind_errors=%d restore_errors=%d\n",
         label, bind_errors, restore_errors);
}

int main(void) {
  run_gate_cases();
  run_region("normal_first", -1, -1, TEST_THREADS, 0, 0, 1);
  run_region("bind_injected", 2, -1, TEST_THREADS, 1, 0, 0);
  run_region("team_mismatch", -1, -1, TEST_THREADS - 1,
             TEST_THREADS, 0, 0);
  run_region("restore_injected", -1, 1, TEST_THREADS, 0, 1, 0);
  run_region("normal_second", -1, -1, TEST_THREADS, 0, 0, 0);
  puts("OPTIMIZER_AFFINITY_C_ORACLE_PASS");
  return 0;
}
