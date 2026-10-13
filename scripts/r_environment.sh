#!/usr/bin/env bash
# Shared selection and isolation for the setup and R runner.
environment_name="${DESMOGNATHUS_R_ENV:-desmognathus}"

resolve_r_environment() {
  environment_prefix=""
  if [[ "$environment_name" == */* ]]; then
    environment_prefix="$(realpath -m -- "$environment_name")"
    [[ -d "$environment_prefix" ]] && return 0
    return 3
  fi
  if ! command -v conda >/dev/null 2>&1; then
    echo "Conda is required to resolve R environment '$environment_name'." >&2
    return 1
  fi
  local configuration line prefix directory
  local candidates=()
  if [[ "$environment_name" == base || "$environment_name" == root ]]; then
    if ! environment_prefix="$(conda info --base)"; then
      echo "Could not resolve the Conda base environment." >&2
      return 1
    fi
    [[ -d "$environment_prefix" ]] && return 0
    return 3
  fi
  # `conda env list` can rewrite the registry to remove stale entries. Read it
  # directly and query configuration so even --check leaves it untouched.
  if ! configuration="$(conda config --show envs_dirs)"; then
    echo "Could not read Conda environment directories." >&2
    return 1
  fi
  while IFS= read -r line; do
    if [[ "$line" =~ ^[[:space:]]*-[[:space:]]+(/.*)$ ]]; then
      directory="${BASH_REMATCH[1]}"
      candidates+=("$directory/$environment_name")
    fi
  done <<< "$configuration"
  local registry="$HOME/.conda/environments.txt"
  if [[ -f "$registry" ]]; then
    if [[ ! -r "$registry" ]]; then
      echo "Cannot read Conda registry '$registry'; use an explicit prefix." >&2
      return 1
    fi
    while IFS= read -r prefix || [[ -n "$prefix" ]]; do
      [[ "${prefix##*/}" == "$environment_name" ]] && candidates+=("$prefix")
    done < "$registry"
  fi
  for prefix in "${candidates[@]}"; do
    [[ -d "$prefix" ]] || continue
    if [[ -n "$environment_prefix" && "$environment_prefix" != "$prefix" ]]; then
      echo "Ambiguous R environment '$environment_name'; set DESMOGNATHUS_R_ENV to an explicit prefix." >&2
      return 1
    fi
    environment_prefix="$prefix"
  done
  [[ -n "$environment_prefix" ]] && return 0
  return 3
}

isolate_r_environment() {
  export CONDA_PREFIX="$environment_prefix"
  export PATH="$environment_prefix/bin:$PATH"
  export R_ENVIRON_USER=/dev/null
  export R_PROFILE_USER=/dev/null
  export R_LIBS_USER="$environment_prefix/lib/R/library"
  unset R_LIBS R_LIBS_SITE
}

validate_r_environment() {
  local rscript="$environment_prefix/bin/Rscript"
  if [[ ! -x "$rscript" || ! -x "$environment_prefix/bin/R" ]]; then
    echo "Missing R runtime in '$environment_prefix'." >&2
    return 1
  fi
  "$rscript" --vanilla - "$project_root/environment.yml" \
    "$project_root/scripts/publication/gbe_theme.R" <<'RSCRIPT' || return 1
args <- commandArgs(trailingOnly = TRUE)
spec <- readLines(args[[1]])
dependencies <- sub("^  - r-", "", grep("^  - r-", spec, value = TRUE))
for (dependency in dependencies) {
  fields <- strsplit(dependency, "=", fixed = TRUE)[[1]]
  package <- fields[[1]]
  if (package == "base") {
    expected <- fields[[2]]
    if (!startsWith(as.character(getRversion()), paste0(expected, "."))) {
      stop("R ", expected, ".x is required; found ", getRversion(), call. = FALSE)
    }
  } else {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Missing required R package: ", package, call. = FALSE)
    }
    if (length(fields) > 1L && as.character(packageVersion(package)) != fields[[2]]) {
      stop("Required ", package, " version: ", fields[[2]], "; found ",
           packageVersion(package), call. = FALSE)
    }
  }
}
if (!requireNamespace("phylopath", quietly = TRUE) ||
    as.character(packageVersion("phylopath")) != "1.3.1") {
  stop("phylopath 1.3.1 is required", call. = FALSE)
}
if (!capabilities("cairo")) stop("R Cairo support is required", call. = FALSE)
source(args[[2]])
stopifnot(identical(resolve_gbe_font(), "Arial"))
cat("R packages, runtime, Cairo and project Arial font check passed\n")
RSCRIPT
  local converter
  for converter in gs convert; do
    if [[ ! -x "$environment_prefix/bin/$converter" ]]; then
      echo "Missing figure converter: $environment_prefix/bin/$converter" >&2
      return 1
    fi
  done
  "$environment_prefix/bin/gs" --version >/dev/null || return 1
  "$environment_prefix/bin/convert" -version >/dev/null || return 1
  echo "R environment validated: $environment_prefix"
}
