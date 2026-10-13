#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$project_root/scripts/r_environment.sh"
phylopath_version="1.3.1"
phylopath_sha256="fbbfe6ffe78fe4b6d59c98202e3c6e03332d379b4a18aed62008bc6944d5227f"

check_only=false
if [[ "${1:-}" == --check && "$#" == 1 ]]; then
  check_only=true
elif [[ "$#" != 0 ]]; then
  echo "Usage: scripts/setup_r.sh [--check]" >&2
  exit 2
fi

if resolve_r_environment; then
  cd "$project_root"
  isolate_r_environment
  if ! validate_r_environment; then
    echo "Existing environment was not changed. Create a separate environment with DESMOGNATHUS_R_ENV=<new-name-or-prefix> scripts/setup_r.sh." >&2
    exit 1
  fi
  exit 0
else
  status=$?
  [[ "$status" == 3 ]] || exit 1
fi
cd "$project_root"

if [[ "$check_only" == true ]]; then
  echo "Missing R environment '$environment_name'; --check does not create environments." >&2
  exit 1
fi
if ! command -v conda >/dev/null 2>&1; then
  echo "Conda is required to create R environment '$environment_name'." >&2
  exit 1
fi
creation_command=(conda env create)
environment_spec="$project_root/environment.yml"
if [[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]]; then
  creation_command=(conda create --yes)
  environment_spec="$project_root/environment-linux-64.explicit.txt"
  if [[ ! -f "$environment_spec" ]]; then
    echo "Missing Linux x86_64 package lock: $environment_spec" >&2
    exit 1
  fi
fi
if [[ "$environment_name" == */* ]]; then
  "${creation_command[@]}" --prefix "$environment_prefix" --file "$environment_spec"
else
  "${creation_command[@]}" --name "$environment_name" --file "$environment_spec"
fi
if ! resolve_r_environment; then
  echo "Conda did not create R environment '$environment_name'." >&2
  exit 1
fi
isolate_r_environment

setup_temp="$(mktemp -d)"
trap 'find "$setup_temp" -depth -delete' EXIT
source_tarball="$setup_temp/phylopath_${phylopath_version}.tar.gz"
source_urls=(
  "https://cran.r-project.org/src/contrib/phylopath_${phylopath_version}.tar.gz"
  "https://cran.r-project.org/src/contrib/Archive/phylopath/phylopath_${phylopath_version}.tar.gz"
)

downloaded=false
for source_url in "${source_urls[@]}"; do
  if curl -fsSL "$source_url" -o "$source_tarball"; then
    downloaded=true
    break
  fi
done
if [[ "$downloaded" != true ]]; then
  echo "Could not download phylopath $phylopath_version from CRAN." >&2
  exit 1
fi

observed_sha256="$(sha256sum "$source_tarball" | cut -d ' ' -f 1)"
if [[ "$observed_sha256" != "$phylopath_sha256" ]]; then
  echo "phylopath source checksum mismatch: $observed_sha256" >&2
  exit 1
fi

"$environment_prefix/bin/R" CMD INSTALL --no-multiarch \
  --library="$environment_prefix/lib/R/library" "$source_tarball"
echo "Installed checksum-verified CRAN phylopath $phylopath_version"
validate_r_environment
