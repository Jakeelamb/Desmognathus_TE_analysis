# Local runtime and preservation

Use Pixi for this repository. `pixi.toml` is the task and dependency authority,
and `pixi.lock` is the active lockfile. This project policy replaces the older
UV and Conda setup instructions for local work.

- Run `make setup` to install all locked environments and prepare the private
  checksum-verified R library. Do not select a host Conda or UV environment.
- Run `make check-env`, `make test`, and `make lint` for the maintained checks.
  Make delegates to Pixi with `--locked --no-install`. New local workflows
  belong in the Pixi task registry.
- Use the default environment for the paper analyses and the `bonsai`
  environment for the pinned upstream Bonsai runner.
- Use `image-qc` for the existing-mask standard-slide audit. It adds image
  dependencies without changing the default or Bonsai dependency sets.
- Use `recovery` for microscopy recovery tooling. Run explicit `recovery-setup`
  before `microscopy-recovery`; never repair a failed tool cache in place.
  Keep original packets, authored ROI proposals and review decisions separate.
  Do not create real acceptance events in infrastructure tests or use diagnostic
  reports to promote a calibration. See `docs/MICROSCOPY_RECOVERY.md`.
- Keep `pyproject.toml`, `uv.lock`, `environment.yml`, and
  `environment-linux-64.explicit.txt` unchanged. They are retained dependency
  snapshots; scientific manifests hash several of them.
- Preserve data, Methods, figures, manifests, and existing analysis outputs.
  Exercise artifact-writing tasks in a complete scratch checkout when checking
  infrastructure changes. Do not regenerate the retained release as a side
  effect of a dependency migration.
- The report Rmd templates are maintained source code. Retained HTML records
  its historical template lineage; a fresh rendering must carry the current
  template and shared-contract fingerprints.
- Preserve hash-bound scientific scripts and historical runners. Check their
  manifests before editing them. The new Pixi Bonsai runner lives in
  `scripts/run_bonsai.py`; the historical shell runners remain evidence.
- Historical HPC job sources and command notes are indexed in
  `docs/HPC_SCRIPT_AUDIT.md`. Their exact copies are provenance evidence, not
  maintained Pixi tasks or proof of execution. The custom TE library, container
  images, raw microscopy inputs, and complete upstream run linkage remain
  reproducibility gaps. Do not invent replacements or execute HPC discovery
  workflows during local setup.

The verified platform is Linux x86-64 with glibc 2.39 or newer. Arial remains
an external font prerequisite checked by the actual project font resolver.
See `docs/REPRODUCIBILITY.md` for the verification scope and upstream boundary.
