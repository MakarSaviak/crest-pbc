#!/usr/bin/env bash
set -euo pipefail
umask 077
export LC_ALL=C

if (( $# != 1 )); then
  printf 'usage: %s EMPTY_OUTPUT_DIRECTORY\n' "$0" >&2
  exit 64
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source_root=$(cd -- "$script_dir/.." && pwd -P)
output_root=$1
if [[ -e "$output_root" ]]; then
  printf 'refusing existing output path: %s\n' "$output_root" >&2
  exit 65
fi
mkdir -m 700 -- "$output_root"
output_root=$(cd -- "$output_root" && pwd -P)

cc=${CC:-gcc}
fc=${FC:-gfortran}
# This standalone oracle compiles the complete runtime, whose process-MTD
# helpers are intentionally unreachable from the optimizer-only driver.
# Keep every other warning fatal.
cflags=(-std=c11 -O1 -g -Wall -Wextra -Wpedantic -Werror -Wno-unused-function)
sanflags=(-fsanitize=address,undefined -fno-omit-frame-pointer)
test_define=(-DCREST_OPTIMIZER_AFFINITY_TESTING)
runtime=$source_root/src/mtd_process_runtime.c
scheduler=$source_root/src/mtd_process_scheduler.F90
parallel_source=$source_root/src/algos/parallel.f90
c_driver=$script_dir/optimizer_affinity_c_driver.c
fortran_driver=$script_dir/optimizer_affinity_fortran_driver.F90
caller_checker=$script_dir/check_optimizer_affinity_callers.py

capture_command() {
  local stdout_path=$1
  local stderr_path=$2
  shift 2
  local status

  set +e
  "$@" >"$stdout_path" 2>"$stderr_path"
  status=$?
  set -e
  if [[ -s "$stdout_path" ]]; then
    sed -n '1,240p' "$stdout_path"
  fi
  if [[ -s "$stderr_path" ]]; then
    sed -n '1,240p' "$stderr_path" >&2
  fi
  if (( status != 0 )); then
    printf 'command failed with status %d: %s\n' "$status" "$*" >&2
    return "$status"
  fi
}

record_source_inputs() {
  local unsorted_paths=$output_root/.oracle-source-input-paths.unsorted
  local path_list=$output_root/oracle-source-input-paths.txt
  local hash_list=$output_root/oracle-source-inputs.sha256
  local -a source_inputs

  (
    cd -- "$source_root"
    find src -type f \( -iname '*.f' -o -iname '*.f90' \) -printf 'src/%P\n'
  ) >"$unsorted_paths"
  printf '%s\n' \
    src/mtd_process_runtime.c \
    prototype_tests/optimizer_affinity_c_driver.c \
    prototype_tests/optimizer_affinity_fortran_driver.F90 \
    prototype_tests/check_optimizer_affinity_callers.py \
    prototype_tests/run_optimizer_affinity_oracles.sh \
    >>"$unsorted_paths"
  sort -u -- "$unsorted_paths" >"$path_list"
  rm -f -- "$unsorted_paths"
  mapfile -t source_inputs <"$path_list"
  if (( ${#source_inputs[@]} == 0 )); then
    printf 'oracle source-input inventory is empty\n' >&2
    return 68
  fi
  (
    cd -- "$source_root"
    sha256sum -- "${source_inputs[@]}"
  ) >"$hash_list"
}

record_source_inputs
capture_command "$output_root/cc-version.stdout" "$output_root/cc-version.stderr" \
  "$cc" --version
capture_command "$output_root/fc-version.stdout" "$output_root/fc-version.stderr" \
  "$fc" --version

capture_command "$output_root/gate-wiring.stdout" "$output_root/gate-wiring.stderr" \
  awk '
  /^[[:space:]]*subroutine run_mtd_process_batch\(/ { inside = 1 }
  inside && /c_set_optimizer_affinity_process_batch_ready\(0_c_int\)/ {
    reset_count++; reset_line = NR
  }
  inside && /call validate_process_configuration/ { validate_line = NR }
  inside && /call print_process_finish/ { finish_line = NR }
  inside && /parent energy\/gradient call accounting overflow/ { overflow_line = NR }
  inside && /c_set_optimizer_affinity_process_batch_ready\(1_c_int\)/ {
    arm_count++; arm_line = NR
  }
  inside && /engrad_total = engrad_total\+batch_engrad_calls/ { accounting_line = NR }
  inside && /^[[:space:]]*end subroutine run_mtd_process_batch/ {
    end_line = NR; inside = 0
  }
  END {
    if (reset_count != 1 || arm_count != 1 ||
        !(reset_line < validate_line && finish_line < arm_line &&
          overflow_line < arm_line && arm_line < accounting_line &&
          accounting_line < end_line)) {
      print "OPTIMIZER_AFFINITY_GATE_WIRING_ORACLE_FAIL" > "/dev/stderr"
      exit 1
    }
    print "OPTIMIZER_AFFINITY_GATE_WIRING_ORACLE_PASS reset_line=" reset_line \
          " arm_line=" arm_line
  }
' "$scheduler"

capture_command "$output_root/screen-scope-wiring.stdout" \
  "$output_root/screen-scope-wiring.stderr" awk '
  /^subroutine crest_oloop\(/ { inside = 1 }
  inside && /affinity_standalone_allowed = 0_c_int/ {
    clear_count++; clear_line = NR
  }
  inside && /if \(env%crestver == crest_screen\) affinity_standalone_allowed = 1_c_int/ {
    screen_count++; screen_line = NR
  }
  inside && /affinity_cstat = c_optimizer_affinity_resolve/ {
    resolve_count++; resolve_line = NR
  }
  inside && /& affinity_standalone_allowed,affinity_gate_enabled\)/ {
    bit_count++; bit_line = NR
  }
  inside && /^end subroutine crest_oloop$/ { end_line = NR; inside = 0 }
  END {
    if (clear_count != 1 || screen_count != 1 || resolve_count != 1 ||
        bit_count != 1 ||
        !(clear_line < screen_line && screen_line < resolve_line &&
          resolve_line < bit_line && bit_line < end_line)) {
      print "OPTIMIZER_AFFINITY_SCREEN_SCOPE_WIRING_ORACLE_FAIL" > "/dev/stderr"
      exit 1
    }
    print "OPTIMIZER_AFFINITY_SCREEN_SCOPE_WIRING_ORACLE_PASS screen_line=" \
          screen_line " resolve_line=" resolve_line
  }
' "$parallel_source"

capture_command "$output_root/caller-hardening.stdout" \
  "$output_root/caller-hardening.stderr" python3 "$caller_checker" "$source_root"

capture_command "$output_root/compile-runtime-production.stdout" \
  "$output_root/compile-runtime-production.stderr" \
  "$cc" "${cflags[@]}" -c "$runtime" -o "$output_root/runtime-production.o"
capture_command "$output_root/runtime-production.nm" \
  "$output_root/runtime-production-nm.stderr" \
  nm "$output_root/runtime-production.o"

set +e
rg -q 'crest_optimizer_affinity_test_' "$output_root/runtime-production.nm" \
  >"$output_root/production-symbol-scan.stdout" \
  2>"$output_root/production-symbol-scan.stderr"
scan_status=$?
set -e
if (( scan_status == 0 )); then
  printf 'test-only optimizer-affinity symbol leaked into production object\n' >&2
  exit 66
elif (( scan_status != 1 )); then
  printf 'production symbol scan failed with status %d\n' "$scan_status" >&2
  exit 67
fi
printf 'OPTIMIZER_AFFINITY_PRODUCTION_SYMBOL_ORACLE_PASS test_symbols=absent\n' \
  >"$output_root/production-symbol-scan.result"
sed -n '1,20p' "$output_root/production-symbol-scan.result"

capture_command "$output_root/compile-runtime-test.stdout" \
  "$output_root/compile-runtime-test.stderr" \
  "$cc" "${cflags[@]}" "${test_define[@]}" -c "$runtime" \
  -o "$output_root/runtime-test.o"
capture_command "$output_root/compile-c-oracle.stdout" \
  "$output_root/compile-c-oracle.stderr" \
  "$cc" "${cflags[@]}" -pthread "$c_driver" "$output_root/runtime-test.o" \
  -o "$output_root/optimizer-affinity-c-oracle"
capture_command "$output_root/compile-fortran-oracle.stdout" \
  "$output_root/compile-fortran-oracle.stderr" \
  "$fc" -O1 -g -Wall -Wextra -Werror -fopenmp "$fortran_driver" \
  "$output_root/runtime-test.o" -o "$output_root/optimizer-affinity-fortran-oracle"

capture_command "$output_root/compile-runtime-sanitized.stdout" \
  "$output_root/compile-runtime-sanitized.stderr" \
  "$cc" "${cflags[@]}" "${sanflags[@]}" "${test_define[@]}" -c "$runtime" \
  -o "$output_root/runtime-test-sanitized.o"
capture_command "$output_root/compile-c-oracle-sanitized.stdout" \
  "$output_root/compile-c-oracle-sanitized.stderr" \
  "$cc" "${cflags[@]}" "${sanflags[@]}" -pthread "$c_driver" \
  "$output_root/runtime-test-sanitized.o" \
  -o "$output_root/optimizer-affinity-c-oracle-sanitized"
capture_command "$output_root/compile-fortran-oracle-sanitized.stdout" \
  "$output_root/compile-fortran-oracle-sanitized.stderr" \
  "$fc" -O1 -g -Wall -Wextra -Werror -fopenmp "${sanflags[@]}" \
  "$fortran_driver" "$output_root/runtime-test-sanitized.o" \
  -o "$output_root/optimizer-affinity-fortran-oracle-sanitized"

capture_command "$output_root/c-oracle.stdout" "$output_root/c-oracle.stderr" \
  env -u OMP_PROC_BIND -u OMP_PLACES -u GOMP_CPU_AFFINITY \
  -u CREST_MTD_INTERNAL_WORKER -u CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY \
  "$output_root/optimizer-affinity-c-oracle"
capture_command "$output_root/fortran-oracle.stdout" \
  "$output_root/fortran-oracle.stderr" \
  env -u OMP_PROC_BIND -u OMP_PLACES -u GOMP_CPU_AFFINITY \
  -u CREST_MTD_INTERNAL_WORKER -u CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY \
  OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OMP_THREAD_LIMIT=4 \
  "$output_root/optimizer-affinity-fortran-oracle"
capture_command "$output_root/c-oracle-sanitized.stdout" \
  "$output_root/c-oracle-sanitized.stderr" \
  env -u OMP_PROC_BIND -u OMP_PLACES -u GOMP_CPU_AFFINITY \
  -u CREST_MTD_INTERNAL_WORKER -u CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY \
  ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
  UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  "$output_root/optimizer-affinity-c-oracle-sanitized"
capture_command "$output_root/fortran-oracle-sanitized.stdout" \
  "$output_root/fortran-oracle-sanitized.stderr" \
  env -u OMP_PROC_BIND -u OMP_PLACES -u GOMP_CPU_AFFINITY \
  -u CREST_MTD_INTERNAL_WORKER -u CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY \
  ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
  UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OMP_THREAD_LIMIT=4 \
  "$output_root/optimizer-affinity-fortran-oracle-sanitized"

printf 'scope=all regular files in this evidence directory except oracle-output.sha256\n' \
  >"$output_root/oracle-output-manifest.scope.txt"
printf 'OPTIMIZER_AFFINITY_LOCAL_ORACLES_PASS output=%s\n' "$output_root" \
  >"$output_root/runner-summary.stdout"
: >"$output_root/runner-summary.stderr"

manifest_paths=$output_root/.oracle-output-paths
(
  cd -- "$output_root"
  find . -type f ! -name 'oracle-output.sha256' ! -name '.oracle-output-paths' \
    -printf '%P\n'
) >"$manifest_paths"
sort -u -o "$manifest_paths" -- "$manifest_paths"
mapfile -t evidence_files <"$manifest_paths"
rm -f -- "$manifest_paths"
(
  cd -- "$output_root"
  sha256sum -- "${evidence_files[@]}"
) >"$output_root/oracle-output.sha256"
sed -n '1,20p' "$output_root/runner-summary.stdout"
