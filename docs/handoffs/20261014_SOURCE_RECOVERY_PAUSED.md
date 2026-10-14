# Paused microscopy source recovery — 14 October 2026

**Paused at Jake's request before temporarily unplugging the computer.** Do not
resume jobs, source reads, tests or agent dispatch until Jake asks to resume.
This is a restart checkpoint, not a completed source-reader qualification or a
publication-readiness verdict.

## Resume from here

Paste this into a new Codex session:

> $pstack:Poteto Mode resume the paused Desmognathus microscopy source recovery
> from /home/jake/Projects/Desmognathus_TE/docs/handoffs/20261014_SOURCE_RECOVERY_PAUSED.md.
> Restore the saved roles, inspect retained evidence, and finish the pending
> independent verification. Preserve all data and outputs; reuse completed
> producer results rather than rerunning them.

The [team manifest](20261014_SOURCE_RECOVERY_TEAM.json) contains ownership,
statuses and individual resume briefs. It is a record for the coordinator;
it does not start agents or jobs on boot.

Start with local, read-only checks:

```bash
cd /home/jake/Projects/Desmognathus_TE
git status --short
git branch --show-current
git log -1 --format='%H %s' -- docs/handoffs/20261014_SOURCE_RECOVERY_PAUSED.md
sha256sum scripts/microscopy/source_reader.py tests/test_source_reader.py pixi.toml pixi.lock
findmnt --mountpoint /run/media/jake/easystore
```

Check the saved pins below before interpreting a previous receipt. If the drive
is absent, continue local code and receipt review; defer live source checks.
Do not recreate or repair caches, packets or output directories automatically.
The ignored evidence under `.pixi/` and `.git/` remains on this computer and is
not included in a portable Git clone.

## Repository checkpoint

| Item | Saved state |
| --- | --- |
| Checkout | `/home/jake/Projects/Desmognathus_TE` |
| Branch | `codex/microscopy-recovery-tooling` |
| Base HEAD before this checkpoint | `93282800655f6d3e8de71ebd27c9d371f98ce243` |
| Reader SHA-256 | `9cb60a82d9bf6e9eacf1a233e0ae49de1b17bfc08564078257f87305418f4251` |
| Reader tests SHA-256 | `c0601cb9ad479a763da2e37b18a41280d4881b9b3735756473b6088bb8ea2ef3` |
| `pixi.toml` SHA-256 | `c5e964b5dd56eeecbe2338d7154d0d28ec279d9c400a9002e0db93044a158d68` |
| `pixi.lock` SHA-256 | `792ffc8c6f71bc0b63c4fe183df816de5d4ec8c4c88e5af6d8f7269eb38c2cc4` |

Before saving the handoff, the only working-tree changes were
`scripts/microscopy/source_reader.py`, `tests/test_source_reader.py`,
`pixi.toml`, and `docs/MICROSCOPY_RECOVERY.md`. This checkpoint also adds the
handoff, team manifest and a short native-source recovery record. The checkpoint
commit can be identified with the `git log` command above; its title explicitly
marks final review pending.

## Completed work and reusable evidence

The source drive was mounted at `/run/media/jake/easystore` using NTFS3 with a
writable mount. Our source operations opened files for reading. This is not
OS-enforced read-only protection.

* Standards: original Process 413/414 VSI files and their ETS companions were
  found under `20250521/Slide glass_1321` and `Slide glass_1322`. Native metadata
  confirms 40× acquisition, 0.12 µm pixels and stored camera settings. Two small
  green-channel crops recovered context across former tile seams. Retained crop
  overlap and masks were unchanged. No replacement contour was authored here.
* Targets: four vendor RGB BigTIFF exports for Processes 314/315/320/321 were
  found under `Dried_blood_nopmount/Process_<ID>/Process_<ID>_raw.ome.btf`.
  Export reports link each input VSI basename to its output. Original target VSI
  files were not located; absolute input paths were not recovered.
* All **35 target crops** exactly match native green pixels from those exports:
  12 for Process 314, 8 for 315, 5 for 320 and 10 for 321. Shapes and uint8 dtypes
  also match. All four source stat fingerprints stayed unchanged. Large source
  payloads were **not fully hashed**. The job finished with exit 0 before pause;
  no partial or failed cases remain.
* A fresh `registered-reader-v1` recovery packet uses the current toolchain.
  All 82 case IDs, states, case metadata and raw/mask PNG bytes match the original
  pilot. Only toolchain binding and inventory path changed. No old decisions or
  proposals were transferred.
* The optional Bio-Formats 8.5.0 reader was exercised on genuine VSI crops and
  metadata. The latest implementation adds explicit TIFF-page and stat-only
  modes for the restored tiled RGB exports. It reads intersecting tiles rather
  than decoding an entire slide.

The [native-source record](../recovery_work_20261014/NATIVE_SOURCE_RECOVERY.md)
summarizes source and acquisition limits. Detailed receipts remain local:

| Local evidence | Purpose |
| --- | --- |
| `.git/drive_recovery_20261014/metadata/native-v1/` | Standard native metadata, commands and raw logs |
| `.git/drive_recovery_20261014/metadata/targets-v2/` | Export linkage, actual TIFF headers, metadata and decoder failures |
| `.pixi/runs/microscopy-recovery/source-context-v1/` | Two wider native standard crops and exact-overlap evidence |
| `.pixi/runs/microscopy-recovery/source-reader-live-v1/` | Maintained Bio-Formats reader's genuine VSI outputs |
| `.pixi/runs/microscopy-recovery/registered-reader-v1/equivalence.json` | Fresh packet versus original 82-case equivalence |
| `.pixi/runs/microscopy-recovery/target-source-pixels-v1/` | All 35 target comparisons, commands, byte ranges, source stats and logs |
| `.git/drive_recovery_20261014/tooling/verification-tiff-page.json` | Latest worker's 43-test, lint, format and diff receipt |
| `.git/drive_recovery_20261014/source-reader-review/tiff-page/PAUSED_HANDOFF.md` | Independent verifier's exact stopping point |
| `.git/drive_recovery_20261014/root/maintained_checks-before-btf.json` | Root checks for the earlier reader revision |
| `.git/drive_recovery_20261014/pause/` | Process audit and stopped gallery receipt |

Pinned local receipt SHA-256 values:

```text
89653b8c0b1de7c2b192e8d37e40a4d07c010767006c9e1a1db167c110898b18  target-source-pixels-v1/summary.json
5be2a9df6bb9f72c00d1839cc2bc929f6cb3a26c18ba499f4261224d379320d3  target-source-pixels-v1/run-pins.json
c6898da9e8d3b0d601d673345ef83b7788ae037dcce77b9c431da8c17fac573f  registered-reader-v1/equivalence.json
74061fa343ab08d4c6fc222416cb6e163e796fbe761310d1e698ed8cfa18e72c  tooling/verification-tiff-page.json
fcbfdc54d9a25f40ea2f590772120359a36fba80a1964fed10d3b14cbad0a76b  source-reader-review/tiff-page/PAUSED_HANDOFF.md
```

## Pending work, in order

1. Restore the independent verifier first. Finish review of the exact latest
   TIFF-page/stat-only revision. Inspect retained all-35 comparisons and packet
   equivalence before starting new commands. Independently check seam and image
   edge handling, missing/truncated tile rejection, selected byte ranges and
   complete crop initialization. The previous PASS covers an earlier reader.
2. Send accepted findings to the implementation worker with its original two-file
   ownership. Keep one coordinator for shared Pixi/docs and landing. Workers
   must accommodate other changes and never revert them.
3. Finish operator documentation for explicit TIFF-page/stat-only use. Preserve
   broken embedded OME metadata rather than rewriting the source files. Complete
   the small public recovery record from retained evidence.
4. Run appropriate final checks after code is frozen, through locked Pixi. Root's
   earlier `make check-env`, `make lint`, `make test` (429 passed, 38 skipped) and
   recovery suite (155 passed, no skips) preceded the TIFF-page revision. The
   latest worker reports 43 reader tests passed with no skips. Final independent
   qualification and root checks for the new revision remain pending.
5. Complete the preservation audit, review the final diff, and commit/push only
   small authorized code/docs. Keep raw imagery, native tools and private packets
   outside Git. Then return to bounded microscopy recovery and Methods gaps.

Do not rerun `target-source-pixels-v1/verify_targets.py` in its existing output
directory. Its producer finished. Repair reporting from saved artifacts if
needed. Any genuinely necessary new run must have a fresh output directory.

## Roles and restoration

The coordinator uses native Poteto routing. Current root was `gpt-6.1-sol` at
`xhigh`. Maximum concurrency is four including root. Existing agent IDs may
survive a short pause; a new session may require recreating them.

| Worker | Role | Pause state | Resume purpose |
| --- | --- | --- | --- |
| `ai_review_verification` | `pstack_verifier`, Sol/xhigh | Interrupted; no child process | Finish independent latest-reader review first |
| `source_reader_tooling` | `pstack_worker`, Sol/high | Completed; files frozen | Fix accepted findings only |
| `source_context_recovery` | `pstack_worker`, Sol/high | Completed; session 97278 exited 0 | Explain existing recovery receipts; no automatic producer rerun |
| `acquisition_evidence_audit` | `pstack_investigator`, Sol/medium | Earlier work completed; not live in final inventory | Resolve a specific remaining source/metadata question if needed |
| `ai_boundary_review_sol` | `pstack_deep`, Sol/xhigh | Earlier review completed; not live | Review genuinely new visual inputs if needed |
| `ai_boundary_review_luna` | `pstack_helper`, Luna/high | Earlier review completed; not live | Independent review of genuinely new visual inputs if needed |

Use `collaboration.followup_task` with a surviving worker and its saved resume
brief. Otherwise use `collaboration.spawn_agent` with the same `agent_type`,
`task_name`, `fork_turns: "none"` and the brief from the team manifest. The native
role selects the recorded model/effort; do not silently substitute models. Restore
only workers needed for the next step. Completed reviewers do not need to repeat
the unchanged 82 images.

Older cleanup/Bonsai worker names are preserved in the manifest as historical
participants, not running jobs: `bonsai_clean_env_fix`, `cleanup_code_plan`,
`cleanup_environment_plan`, `cleanup_python_entrypoints`,
`cleanup_r_entrypoints`. Their exact role settings are not retained here. Do not
restart cleanup or recreate those roles by guessing.

## Preservation and interpretation

Use Pixi exclusively for this repository. Read repository AGENTS instructions
and `docs/MICROSCOPY_RECOVERY.md` before resuming. Preserve dependency snapshots,
data, Methods, figures, manifests, hash-bound scripts and all existing outputs.
Exercise writing tasks in a complete scratch checkout. Do not run HPC discovery,
segmentation or release regeneration as part of this resume.

The original toolchain identity is
`5c35cb7aaba13706e1a11c0c819a65c841b964d735aa0ea9ab28a243e2c1aae5`;
the current registered toolchain is
`9ce2f46f132bf5f344d1bd39960f0d47dc9895b11dd999c68cba7a5cbd1de1e3`.
Original pilot packet SHA-256:
`fd9afdd0fc827d304553d31a003dc9cef6c46a0caaeed0680bbed106650baa94`.
Fresh packet SHA-256:
`46720c9388b7e85d1a4e455ffd244bfd899ee882209b0213511d1b4470073d33`.
Do not weaken runtime guards to load the old packet under the new identity.
Historical ROI trials and AI reviews remain attached to their original packet.

Failed Bio-Formats target OME parsing and earlier failed launcher/setup attempts
remain evidence. The actual TIFF has one tiled RGB page while copied OME metadata
references additional pyramid images/IFDs. TIFF-page mode makes this discrepancy
explicit. Full large-file hashes remain absent; stat equality is narrower.

Recovered pixels do not establish nuclear boundaries, camera linearity, a
physical zero offset or independent DNA content. Absolute genome-size calibration
remains unqualified. No masks, authored ROI proposals, acceptance decisions, IOD
calculations or scientific release outputs were changed in this source recovery
phase. The latest full preservation audit remains pending.

Jake's scope is Methods only. Preserve his confirmations: retained images and
standards use 40× dry brightfield green-channel images, one Feulgen staining run
with both standards, and the same stated fixation protocol. These confirmations
do not establish lifetime preservation history. Do not reopen animal-care
approval questions. Keep the possible specimen 3463/RAP3462 transcription issue
unresolved unless source evidence settles it; do not silently relabel data.

## Stop state

All current child workers were stopped or completed. The verifier had no active
process/session when interrupted. The crop producer had already exited 0. Root
stopped the task-owned reviewer HTTP server, PID 68509, after checking its exact
command. The retained process audit found no task reader, Java decoder or reader
test process. `fuser -m /run/media/jake/easystore` returned no users at that check.
The source drive was not unmounted. Resume requires an explicit user request.
