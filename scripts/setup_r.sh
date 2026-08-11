#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
environment_name="desmognathus"
phylopath_version="1.3.1"
phylopath_sha256="fbbfe6ffe78fe4b6d59c98202e3c6e03332d379b4a18aed62008bc6944d5227f"

conda env update --file "$project_root/environment.yml" --prune

conda_base="$(conda info --base)"
environment_prefix="$conda_base/envs/$environment_name"
rscript="$environment_prefix/bin/Rscript"
if [[ ! -x "$rscript" ]]; then
  echo "Conda did not create the '$environment_name' environment." >&2
  exit 1
fi

export CONDA_PREFIX="$environment_prefix"
export PATH="$environment_prefix/bin:$PATH"
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null

installed_version="$($rscript --vanilla -e \
  "cat(if (requireNamespace('phylopath', quietly=TRUE)) as.character(packageVersion('phylopath')) else '')")"
if [[ "$installed_version" == "$phylopath_version" ]]; then
  echo "phylopath $phylopath_version is already installed"
  exit 0
fi

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
"$rscript" --vanilla -e \
  "stopifnot(as.character(packageVersion('phylopath')) == '$phylopath_version')"
echo "Installed checksum-verified CRAN phylopath $phylopath_version"
