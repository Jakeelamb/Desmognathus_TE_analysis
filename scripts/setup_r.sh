#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/activate_pixi.sh"

check_only=false
if [[ "${1:-}" == --check && "$#" == 1 ]]; then
  check_only=true
elif [[ "$#" != 0 ]]; then
  echo 'Usage: scripts/setup_r.sh [--check]' >&2
  exit 2
fi
resolve_pixi_environment || exit 1
if [[ "${PIXI_PROJECT_ROOT:-}" != "$project_root" || "${CONDA_PREFIX:-}" != "$environment_prefix" ]]; then
  echo 'Run setup-r inside this project Pixi environment: make setup-r (or make check-env).' >&2
  exit 1
fi
cd "$project_root"
isolate_pixi_environment
checker=("$environment_prefix/bin/python" -B "$project_root/scripts/check_environment.py" --r-only)
overlay="$project_root/.pixi/r-library"
if [[ "$check_only" == true || -e "$overlay" || -L "$overlay" ]]; then
  "${checker[@]}"
  exit "$?"
fi
# Validate the immutable runtime before downloading or creating an overlay.
"${checker[@]}" --without-phylopath
setup_lock="$project_root/.pixi/r-overlay.setup.lock"
mkdir "$setup_lock" 2>/dev/null || { echo "R overlay setup is already running (lock: $setup_lock)." >&2; exit 1; }
trap 'rmdir -- "$setup_lock"' EXIT
[[ ! -e "$overlay" && ! -L "$overlay" ]] || { echo 'R overlay appeared during setup; refusing to overwrite it.' >&2; exit 1; }
setup_temp="$(mktemp -d "$project_root/.pixi/phylopath-setup.XXXXXXXX")"
trap 'rm -rf -- "$setup_temp"; rmdir -- "$setup_lock"' EXIT
source_tarball="$setup_temp/phylopath_1.3.1.tar.gz"
source_urls=(
  'https://cran.r-project.org/src/contrib/phylopath_1.3.1.tar.gz'
  'https://cran.r-project.org/src/contrib/Archive/phylopath/phylopath_1.3.1.tar.gz'
)
downloaded=false
for source_url in "${source_urls[@]}"; do
  if curl -fsSL "$source_url" -o "$source_tarball"; then downloaded=true; break; fi
done
[[ "$downloaded" == true ]] || { echo 'Could not download phylopath 1.3.1 from CRAN.' >&2; exit 1; }
observed_sha256="$(sha256sum "$source_tarball" | cut -d ' ' -f 1)"
expected_sha256=fbbfe6ffe78fe4b6d59c98202e3c6e03332d379b4a18aed62008bc6944d5227f
[[ "$observed_sha256" == "$expected_sha256" ]] || { echo "phylopath source checksum mismatch: $observed_sha256" >&2; exit 1; }
mkdir "$setup_temp/library"
"$environment_prefix/bin/R" CMD INSTALL --no-multiarch --library="$setup_temp/library" "$source_tarball"
cat > "$setup_temp/library/.phylopath-source.json" <<'RECEIPT'
{"package":"phylopath","version":"1.3.1","r_version":"4.4.3","source_sha256":"fbbfe6ffe78fe4b6d59c98202e3c6e03332d379b4a18aed62008bc6944d5227f"}
RECEIPT
# Test the staged library before exposing it. Never overwrite an existing overlay.
R_LIBS_USER="$setup_temp/library" "${checker[@]}" --overlay "$setup_temp/library"
[[ ! -e "$overlay" && ! -L "$overlay" ]] || { echo 'R overlay appeared during setup; refusing to overwrite it.' >&2; exit 1; }
mv -T -- "$setup_temp/library" "$overlay"
echo 'Installed checksum-verified CRAN phylopath 1.3.1 in private Pixi R overlay.'
"${checker[@]}"
