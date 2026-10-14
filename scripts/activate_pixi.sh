#!/usr/bin/env bash
# Pixi activation hook and shared project-prefix guards.
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"

resolve_pixi_environment() {
  environment_name="${PIXI_ENVIRONMENT_NAME:-${PIXI_ENV_NAME:-default}}"
  case "$environment_name" in
    default|bonsai|image-qc|recovery) ;;
    *) echo "Unsupported Pixi environment: $environment_name" >&2; return 1 ;;
  esac
  environment_prefix="$project_root/.pixi/envs/$environment_name"
  if [[ -n "${PIXI_PROJECT_ROOT:-}" && "$(realpath -m -- "$PIXI_PROJECT_ROOT")" != "$project_root" ]]; then
    echo "Wrong Pixi project; activate $project_root/pixi.toml." >&2
    return 1
  fi
  if [[ -n "${PIXI_PROJECT_ROOT:-}" && -n "${CONDA_PREFIX:-}" && "$(realpath -m -- "$CONDA_PREFIX")" != "$environment_prefix" ]]; then
    echo "Wrong runtime prefix '$CONDA_PREFIX'; use this project's Pixi environment." >&2
    return 1
  fi
  if [[ ! -x "$environment_prefix/bin/python" || ! -x "$environment_prefix/bin/Rscript" || ! -x "$environment_prefix/bin/R" ]]; then
    echo "Missing Pixi runtime '$environment_prefix'. Run: make setup (checks never install)." >&2
    return 1
  fi
}

isolate_pixi_environment() {
  export CONDA_PREFIX="$environment_prefix"
  export PATH="$environment_prefix/bin:$PATH"
  local inherited_variable
  for inherited_variable in "${!PYTHON@}"; do unset "$inherited_variable"; done
  export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
  unset VIRTUAL_ENV UV_PROJECT_ENVIRONMENT R_HOME R_ENVIRON R_PROFILE R_DEFAULT_PACKAGES
  export R_ENVIRON=/dev/null R_PROFILE=/dev/null
  export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
  export R_LIBS_USER="$project_root/.pixi/r-library"
  # The frozen path runner resets .libPaths(); .Library.site retains this overlay.
  export R_LIBS_SITE="$project_root/.pixi/r-library"
  unset R_LIBS
}

# Pixi sources this hook during activation. Entrypoints source it for its functions.
if [[ "${PIXI_PROJECT_ROOT:-}" == "$project_root" && -n "${CONDA_PREFIX:-}" ]]; then
  resolve_pixi_environment || return 1
  isolate_pixi_environment
fi
