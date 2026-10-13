# Online record recovery, 2026-10-13

Jake's existing Chrome session supplied read-only access to the verified `jakelamb74@gmail.com` Drive account. This bypassed the account mismatch in the separate connector without changing that connection or attempting another sign-in.

| Recovered source | Retained evidence and finding |
|---|---|
| [Salamander slides Lab Work](https://docs.google.com/spreadsheets/d/1qmsaXn9pffPq68fphp9Wz5Tcvup43oas2AN8DKQh8QY/edit) | [Exact workbook review](ONLINE_LAB_WORK.md). `Sheet1!A84:K85` lists standards 32469/32470 as D. fuscus. The workbook has specimen collection records but no staining-event or physical-glass fields. |
| [Desmognathus dataset](https://docs.google.com/spreadsheets/d/1CaMYntCP1eDSa1kX297Y36b6hhb9Jcps_1SQqnE2M6E/edit) | [Link extraction](ONLINE_DATASET_LINKS.md). Genomic and phylogenetic resource links were recovered; this export supplied no microscopy batch mapping. |
| [Desmognathus GBE manuscript](https://docs.google.com/document/d/1IKIZmB_hmDO3XcjcJsQnZ_HNOGnx8SBy6rlRPuDXxVc/edit) | [Markdown snapshot](Desmognathus_GBE_manuscript_source_snapshot.md). This locates the loose full-paper draft; the Abstract, Results and Discussion remain incomplete. It contains older 60× oil-immersion and uniform ethanol-fixation wording that must be reconciled with the author's retained-40×-brightfield correction and unresolved specimen histories. |
| [Cell Stain and image methodology](https://docs.google.com/document/d/1uTvZyVQs7ZlLS_MY44zNUgSO0yTtLN1wliRnJw7lfI0/edit) | [Markdown snapshot](Cell_Stain_and_image_methodology_source_snapshot.md). General Feulgen densitometry and imaging guidance; no dated experimental event or target-to-standard mapping was established. |
| [Staining Procedure](https://docs.google.com/document/d/1JYC3_nhW3JNp66Th7hJNtBI5jL3G54xA4pPLNIQMd3Y/edit) | [Markdown snapshot](Staining_Procedure_source_snapshot.md). Recipes and intended fixation, hydrolysis, staining and imaging steps. This is protocol evidence, not proof of the actual retained-slide settings or staining date. |

All five source snapshots and URLs are pinned in [ONLINE_SOURCE_PINS.json](ONLINE_SOURCE_PINS.json). Browser exports do not independently establish complete native comments, notes, formatting or revision history. The bounded searches do not establish that no other records exist in Drive.

## One-run confirmation

Jake answered the question about all retained brightfield slides and reference specimens 32469/32470 with **“One staining run with both standards”**. [AUTHOR_STAINING_CONFIRMATION.json](AUTHOR_STAINING_CONFIRMATION.json) retains the exact exchange and its scope.

[STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv](STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv) records that assignment for 51 retained images, including the four current-denominator images, and two designated standards. The original identifiers and roles are preserved. `author_confirmed_single_run` is an internal label, not a recovered experimental run name. Staining dates remain blank. Existing local glass mappings match 22 rows; the remaining 31 glass identifiers are unknown. Scan dates are separately labeled and do not supply a staining date.

## Next scientific decision

The current release uses target D. fuscus Processes 314/315/320/321 as its equal-image IOD denominator. The author-designated standard images are Processes 413/414. Single-run membership alone does not resolve this mismatch, the previously recorded disagreement between standard IOD measurements, acquisition comparability, historical identity/preservation notes, or the conditional assembly-derived DNA-content anchor.

The next step is to inspect and, where necessary, remeasure the designated standards under the target QC rules before selecting a calibration denominator. No scientific producer was rerun and no figure, Methods source, online document or original mapping was changed during record recovery. [Independent verification](ONLINE_RECOVERY_VERIFICATION.md) checks the retained evidence and protected inputs.
