#!/usr/bin/env bash
set -euo pipefail

comparison_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${comparison_dir}/../../.." && pwd)"
upstream_commit="b16ef6f0a164610327022c7d0b92de82b3596aad"
run_dir="${1:-${comparison_dir}/output/bonsai_run}"
upstream_dir="${run_dir}/upstream"

if [[ -e "${run_dir}" ]]; then
  echo "Refusing to overwrite existing run directory: ${run_dir}" >&2
  exit 2
fi

mkdir -p "${run_dir}"
uv run --frozen --no-sync python \
  "${comparison_dir}/compare_bonsai_to_pca.py" prepare --run-dir "${run_dir}"

git init -q "${upstream_dir}"
git -C "${upstream_dir}" remote add origin \
  https://github.com/dhdegroot/Bonsai-data-representation.git
git -C "${upstream_dir}" fetch -q --depth 1 origin "${upstream_commit}"
git -C "${upstream_dir}" checkout -q --detach FETCH_HEAD
test "$(git -C "${upstream_dir}" rev-parse HEAD)" = "${upstream_commit}"

uv run --python 3.14 \
  --with-requirements "${upstream_dir}/requirements.txt" \
  python "${upstream_dir}/bonsai/bonsai_main.py" \
  --config_filepath "${run_dir}/bonsai_config.yaml" \
  --step all

uv run --frozen --no-sync python \
  "${comparison_dir}/compare_bonsai_to_pca.py" analyze \
  --tree "${run_dir}/results/final_bonsai/tree.nwk" \
  --metadata "${run_dir}/results/final_bonsai/metadata_bonsai.txt" \
  --output-dir "${run_dir}/comparison" \
  --write

"${repo_root}/scripts/run_r.sh" Rscript \
  "${comparison_dir}/plot_bonsai_comparison.R" \
  "${run_dir}/comparison"

echo "PASS: completed pinned Bonsai comparison under ${run_dir}"
