#!/usr/bin/env bash
set -euo pipefail
umask 077
export LC_ALL=C

if (( $# != 2 )); then
  printf 'usage: %s EMPTY_OUTPUT_DIRECTORY PRE_FIX_SOURCE_TARBALL\n' "$0" >&2
  exit 64
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source_root=$(cd -- "$script_dir/.." && pwd -P)
output_root=$1
baseline_archive=$2
if [[ -e "$output_root" ]]; then
  printf 'refusing existing output path: %s\n' "$output_root" >&2
  exit 65
fi
if [[ ! -f "$baseline_archive" ]]; then
  printf 'missing pre-fix source tarball: %s\n' "$baseline_archive" >&2
  exit 66
fi
baseline_archive=$(cd -- "$(dirname -- "$baseline_archive")" && pwd -P)/$(basename -- "$baseline_archive")
mkdir -m 700 -- "$output_root"
output_root=$(cd -- "$output_root" && pwd -P)

cc=${CC:-gcc}
# The standalone oracle links the complete runtime directly, so static helpers
# unrelated to worker preparation are intentionally unreachable in this binary.
# Keep every other warning fatal.
cflags=(-std=gnu11 -O1 -g -Wall -Wextra -Wpedantic -Werror -Wno-unused-function -fopenmp)
runtime=$source_root/src/mtd_process_runtime.c
driver=$script_dir/mtd_worker_affinity_oracle.c
binary=$output_root/mtd-worker-affinity-oracle
baseline_binary=$output_root/mtd-worker-affinity-negative-control
baseline_runtime=$output_root/pre-fix-mtd_process_runtime.c

mapfile -t baseline_members < <(tar -tzf "$baseline_archive" | \
  awk '/\/src\/mtd_process_runtime[.]c$/')
if (( ${#baseline_members[@]} != 1 )); then
  printf 'expected one pre-fix mtd_process_runtime.c, found %d\n' \
    "${#baseline_members[@]}" >&2
  exit 67
fi
tar -xOzf "$baseline_archive" "${baseline_members[0]}" >"$baseline_runtime"

"$cc" --version >"$output_root/cc-version.stdout" \
  2>"$output_root/cc-version.stderr"
"$cc" "${cflags[@]}" -DCREST_MTD_AFFINITY_PROBE_AVAILABLE=1 \
  "$driver" "$runtime" -o "$binary" \
  >"$output_root/compile.stdout" 2>"$output_root/compile.stderr"
"$cc" "${cflags[@]}" "$driver" "$baseline_runtime" -o "$baseline_binary" \
  >"$output_root/negative-compile.stdout" \
  2>"$output_root/negative-compile.stderr"

set +e
env -u OMP_NUM_THREADS -u OMP_DYNAMIC -u OMP_NESTED \
  -u OMP_MAX_ACTIVE_LEVELS -u OMP_THREAD_LIMIT -u OMP_PROC_BIND \
  -u OMP_PLACES -u GOMP_CPU_AFFINITY -u CREST_MTD_INTERNAL_WORKER \
  -u CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY "$baseline_binary" \
  >"$output_root/negative-oracle.stdout" \
  2>"$output_root/negative-oracle.stderr"
negative_status=$?
set -e
if (( negative_status == 0 || negative_status == 77 )); then
  printf 'pre-fix negative control returned invalid status %d\n' \
    "$negative_status" >&2
  exit 68
fi
if ! grep -Fq \
  'MTD_WORKER_AFFINITY_ORACLE_FAIL: primary affinity is not singleton after worker preparation' \
  "$output_root/negative-oracle.stderr"; then
  sed -n '1,80p' "$output_root/negative-oracle.stdout"
  sed -n '1,80p' "$output_root/negative-oracle.stderr" >&2
  printf 'pre-fix negative control failed for an unrelated reason\n' >&2
  exit 69
fi
printf 'MTD_WORKER_AFFINITY_NEGATIVE_CONTROL_PASS expected_failure_status=%d\n' \
  "$negative_status" >"$output_root/negative-control.result"

set +e
env -u OMP_NUM_THREADS -u OMP_DYNAMIC -u OMP_NESTED \
  -u OMP_MAX_ACTIVE_LEVELS -u OMP_THREAD_LIMIT -u OMP_PROC_BIND \
  -u OMP_PLACES -u GOMP_CPU_AFFINITY -u CREST_MTD_INTERNAL_WORKER \
  -u CREST_EXPERIMENTAL_OPTIMIZER_AFFINITY "$binary" \
  >"$output_root/oracle.stdout" 2>"$output_root/oracle.stderr"
oracle_status=$?
set -e
if (( oracle_status != 0 )); then
  sed -n '1,80p' "$output_root/oracle.stdout"
  sed -n '1,80p' "$output_root/oracle.stderr" >&2
  exit "$oracle_status"
fi

(
  cd -- "$source_root"
  sha256sum -- src/mtd_process_runtime.c \
    src/mtd_process_scheduler.F90 \
    prototype_tests/mtd_worker_affinity_oracle.c \
    prototype_tests/run_mtd_worker_affinity_oracle.sh
) >"$output_root/oracle-source-inputs.sha256"
sha256sum -- "$baseline_archive" >"$output_root/pre-fix-source-tarball.sha256"
(
  cd -- "$output_root"
  sha256sum -- mtd-worker-affinity-oracle mtd-worker-affinity-negative-control \
    pre-fix-mtd_process_runtime.c \
    cc-version.stdout cc-version.stderr compile.stdout compile.stderr \
    negative-compile.stdout negative-compile.stderr negative-oracle.stdout \
    negative-oracle.stderr negative-control.result oracle.stdout oracle.stderr \
    oracle-source-inputs.sha256 pre-fix-source-tarball.sha256
) >"$output_root/oracle-output.sha256"

sed -n '1,20p' "$output_root/negative-control.result"
sed -n '1,40p' "$output_root/oracle.stdout"
