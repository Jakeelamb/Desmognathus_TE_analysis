#!/usr/bin/env bash
set -euo pipefail

comparison_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${comparison_dir}/../../.." && pwd)"
upstream_commit="b16ef6f0a164610327022c7d0b92de82b3596aad"
run_dir="${1:-${comparison_dir}/output/feature_bootstrap_200}"
replicates="${2:-200}"
parallel_jobs="${BONSAI_BOOTSTRAP_JOBS:-8}"
upstream_dir="${run_dir}/upstream"
upstream_venv="${upstream_dir}/.venv"

if [[ "${replicates}" != "200" ]]; then
  echo "The tracked campaign is declared as exactly 200 replicates; received ${replicates}." >&2
  exit 2
fi

uv run --frozen --no-sync python \
  "${comparison_dir}/validate_strategy_localization.py" prepare-bootstrap \
  --run-dir "${run_dir}" \
  --replicates "${replicates}" \
  --seed 20260826

if [[ ! -d "${upstream_dir}/.git" ]]; then
  git init -q "${upstream_dir}"
  git -C "${upstream_dir}" remote add origin \
    https://github.com/dhdegroot/Bonsai-data-representation.git
  git -C "${upstream_dir}" fetch -q --depth 1 origin "${upstream_commit}"
  git -C "${upstream_dir}" checkout -q --detach FETCH_HEAD
fi
test "$(git -C "${upstream_dir}" rev-parse HEAD)" = "${upstream_commit}"

if [[ ! -x "${upstream_venv}/bin/python" ]]; then
  uv venv --python 3.14 "${upstream_venv}"
fi
uv pip install --python "${upstream_venv}/bin/python" \
  --requirement "${upstream_dir}/requirements.txt"
"${upstream_venv}/bin/python" -c 'import numpy, pandas, scipy, sklearn'

run_one_bootstrap() {
  set -euo pipefail
  replicate_id="$1"
  replicate_dir="${run_dir}/bootstrap/replicate_${replicate_id}"
  final_tree="${replicate_dir}/results/final_bonsai/tree.nwk"
  log_path="${replicate_dir}/bonsai.log"
  if [[ -s "${final_tree}" ]]; then
    echo "SKIP replicate_${replicate_id}"
    return
  fi
  (
    cd "${upstream_dir}"
    env \
      OPENBLAS_NUM_THREADS=1 \
      OMP_NUM_THREADS=1 \
      MKL_NUM_THREADS=1 \
      NUMEXPR_NUM_THREADS=1 \
      "${upstream_venv}/bin/python" bonsai/bonsai_main.py \
      --config_filepath "${replicate_dir}/bonsai_config.yaml" \
      --step all
  ) >"${log_path}" 2>&1
  test -s "${final_tree}"
  echo "DONE replicate_${replicate_id}"
}

export run_dir upstream_dir upstream_venv
export -f run_one_bootstrap
seq -w 0 "${replicates}" | \
  xargs -P "${parallel_jobs}" -I '{}' bash -c 'run_one_bootstrap "$1"' _ '{}'

uv run --frozen --no-sync python \
  "${comparison_dir}/validate_strategy_localization.py" analyze-bootstrap \
  --run-dir "${run_dir}" \
  --output-dir "${run_dir}/validation" \
  --write

"${repo_root}/scripts/run_r.sh" Rscript \
  "${comparison_dir}/plot_bonsai_validation.R" \
  "${run_dir}/validation"

echo "PASS: completed 200-replicate Bonsai feature bootstrap under ${run_dir}"
