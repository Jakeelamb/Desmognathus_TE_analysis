# HPC source recovery audit — 13 October 2026

The missing HPC scripts were partly recoverable. Twelve historical shell jobs
survive in the April 2025 Git history shared by this repository, the official
backup and the archive. Five were already preserved during submission
preparation; this audit adds exact copies of the other seven. A separate project
note supplies commands for fastp, TE classifiers and library preparation, but
contains unfinished parameters and transcription errors. These sources improve
the record of the workflow; they do not establish which jobs produced every
retained output.

The [source registry](hpc_script_audit/manifest.json) records the original paths,
Git revisions and blobs, SHA256 values, and retained-copy locations. Copies are
historical evidence, not maintained Pixi tasks. No recovered job was executed
during this audit. Data, Methods and existing outputs were left unchanged.

The [independent verification receipt](hpc_script_audit/verification.json)
records exact source-copy checks and preservation of all 12,282 protected
entries. `make validate` passed for all 43 publication tables and the TE34,
Path24 and overlap21 panels.

## What was searched

The directory search covered `/home/jake/Projects`, `Documents`, `Downloads`,
`Desktop` and `.local/share/Trash/files`, to a maximum depth of seven. It excluded
environment, package, cache and build directories, unrelated mounts, browser
profiles and chat histories. It found 13 matching roots: this repository and
12 other directories. The nested Omarchy Git repository was examined separately.
This was a bounded local search, not an exhaustive search of every file on the
computer or of the original HPC filesystem.

A further filename-only search across Projects, Documents, Downloads and
Desktop checked the named upstream tools, library/container names and HPC job
extensions. It found no additional Desmognathus producer candidates outside the
13 roots. Eight other `.img` matches have filenames and locations suggesting
unrelated emulator/device firmware; their contents were not opened. Exact
patterns and exclusions are recorded in the shared scope inventory.

| Location | Finding |
|---|---|
| `Projects/Desmognathus_TE_official_bk_4_12_25_2pm` | Original jobs absent from the working tree, present in Git history; 17 reachable revisions and 131 distinct source blobs searched. |
| `Projects/Desmognathus_TE_archive` | Same early job blobs preserved in Git; its archive inventory also lacks the custom library and container images. |
| `Projects/Desmognathus_TE_Backup_OMarchy_safekeeping`, including nested `Desmognathus_TE` | Two January 2026 commits in each repository; no April 2025 history. All 31 compared source files are identical between the copies and concern downstream analysis. |
| `Projects/misc_desmog_scripts` | Exploratory notebooks and downstream processors. TEsorter annotations are read, rather than generated. |
| `Projects/Desmognathus_Results` and `Desmognathus_RPC_Results` | 89 result files; no recovered HPC jobs. Spreadsheet and RData metadata searches found no HPC command references. |
| Two `Projects/.cleanup_backups/Desmognathus*` directories | Local analysis, review and presentation backups; no additional original HPC jobs found. |
| `Documents/01-Projects/Desmognathus` | 15 project notes; `dnaPipeTE.md` contains the strongest additional command evidence. |
| Two `genometools-ltr-rs/artifacts/desmognathus_ectopic_ratios_*` directories | 26 and 215 files: tables, Parquet and logs from downstream LTR table conversion and ratio calculation; no HPC producer filename candidates. |
| Trashed `.gstack/projects/Jakeelamb-Desmognathus_TE_analysis` | One checkpoint note and three JSONL records; no recovered jobs. |

The shared [scope inventory](hpc_script_audit/scope.json) records searched
locations, exclusions and compact counts. Detailed search receipts remain locally under
`.git/hpc_script_audit_20261013/`: `local_directory_inventory.json`,
`original/recovered-inventory.json`, `backups/backup_audit.json`,
`backups/omarchy-source-parity.json`, `results/inventory.json`,
`results/extended-candidates.json`, and `root_extra_roots.json`. They retain
large inventories outside the shared source bundle.

## Historical jobs

All 12 source records use revision
`6a3d379fbbf32c58ddaed68035b4eee863b8a99a` on 5 April 2025 at
19:10:36 −06:00. The first five files are also available with identical bytes
at `cc64ffb134968250671ff0e3456255063a1aba97`, before their move under
`Projects/01_data_acquisition` or `Projects/02_Repeatlib_construction`. Commit
`d900a09592318790b4e9e5e4bb2730c1d57fc77f`, dated 12 April 2025,
removed the originals from the working tree. These are commit author timestamps,
not execution dates. Exact `commit:path` identities are in the registry.

| Source and relevant lines | Commands and dependencies | Recovery |
|---|---|---|
| [index_and_map.srun](hpc_script_audit/index_and_map.srun), 13–22 | Positional paired reads and mitochondrial reference → `bwa index`, `bwa mem -t 8` → SAMtools conversion/sort → nuclear BAM using `-f 12 -F 256`. Executes `bwa` despite the job label `bwamem2`. | New copy |
| [remove_mito.srun](hpc_script_audit/remove_mito.srun), 12–16 | Existing mapped SAM → sorted BAM → the same unmapped-pair filter. | New copy |
| [RepeatModeler_auto.sh](hpc_script_audit/RepeatModeler_auto.sh), 12–22 | Alpine scratch; Singularity 3.7.4; `tetools_1.87.sif`; assembly/species arguments → `BuildDatabase` → RepeatModeler with 32 threads, `-LTRStruct`, `-genomeSampleSizeMax 10000000000`. | New copy |
| [808080.sh](hpc_script_audit/808080.sh), 10–13 | `master_dfam_rmodeler_rmueller_TELib.fasta` → `_80.fasta`; CD-HIT-EST `-d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -T 64 -M 500000`. | New copy |
| [dnaPipeTE.sh](../outputs/submission_preparation_20261013/recovered_commands/dnaPipeTE.sh), 10–40 | SRX directory bound into `dnapipte.img`; mate 1 nuclear FASTQ; 16 CPUs, 15-Gb denominator, coverage 0.1, `RM_t 0.15`, two samples; library `dedupe_telib.fasta`; output `dnaPipeTE_15gb_10gc_RMt15_smpl_2`. | Previously recovered |
| [repeat_mask.sh](../outputs/submission_preparation_20261013/recovered_commands/repeat_mask.sh), 8–21 | Positional FASTA → RepeatMasker `-s -pa 4`, custom `dedupe_telib.fasta`, `${Fasta}_rmask_out`. | Previously recovered |
| [large_repeat_mask.sh](../outputs/submission_preparation_20261013/recovered_commands/large_repeat_mask.sh), 8–21 | Same configuration, `-pa 32` instead of 4. | Previously recovered |
| [LTR_pipeline.sh](../outputs/submission_preparation_20261013/recovered_commands/LTR_pipeline.sh), 24–69 | GCA-to-SRX lookup → dnaPipeTE `Annotation/LTR_annoted.fasta` → DANTE `Metazoa_v3.1`, 16 CPUs → genome suffix index → LTRharvest → sorted GFF3 → LTRdigest. | Previously recovered |
| [ltr_digest.sh](hpc_script_audit/ltr_digest.sh), 44–51 | Two-CPU digest-only variant; expects `${filename}.gff3_sorted.gff3`, differing from the pipeline's `${filename}_sorted.gff3`. | New copy |
| [Mass_map.sh](../outputs/submission_preparation_20261013/recovered_commands/Mass_map.sh), 27–44 | Paired nuclear FASTQs → existing genome `.fna` Bowtie2 index → species BAM. This is genome mapping, not mapping to individual LTR contigs. | Previously recovered |
| [combined_pipeline.sh](hpc_script_audit/combined_pipeline.sh), 10–18, 53–94 | Necturus project path; per-LTR index → Bowtie2 `-L 20 --very-sensitive-local` → SAMtools sort/index/depth `-aa`; parallel processing. Contains repeated `READ_FILE=$` assignments. Desmognathus use remains unverified. | New copy |
| [Calculate_depth.sh](hpc_script_audit/Calculate_depth.sh), 37–92 | Per-species coordinate CSVs and `.fa.depth.txt` files → `Data/coverage.txt`. Historical implementation requires caution below. | New copy |

All five previously recovered copies match the backup blobs byte for byte.
Current backup `scripts/ltr_pipeline.sh` and `scripts/process_ltr.sh` are later
rewrites; their altered paths and output handling are not evidence that those
versions ran on the HPC.

## Additional commands in the project note

The exact [dnaPipeTE.md note](hpc_script_audit/dnaPipeTE.md) was recovered from
`/home/jake/Documents/01-Projects/Desmognathus/dnaPipeTE.md`.
Its SHA256 is
`14c315608f67e2830155c1957ffb09f9ac8fd885132e5535d70c82e347d71b78`.
Its filesystem modification time is 13 February 2026; that does not date the
experiments described within it.

- Lines 17–41 record fastp with `--detect_adapter_for_pe`, BWA/SAMtools filtering,
  and BAM-to-FASTQ commands. The literal fastp option is `-t 32`; it must not be
  silently rewritten as a thread option. The BWA redirection is literally
  `%3E`, and the name-sort example writes `sorted.bam` while the next command
  reads `input.bam`. These are notes to reconcile, not a ready-to-run pipeline.
- Lines 89–96 repeat the historical RepeatModeler commands.
- Lines 104–116 record DeepTE with `-sp M` and a Metazoans model directory,
  `TEsorter $1 -db rexdb-metazoa -p 64`, and
  `TEClass/TEclass-2.1.3d/TEclassTest.pl $1`. The latter path does not settle
  which TEClass implementation/version should be described in the paper.
- Lines 122–142 describe two CD-HIT rounds, Dfam 3.8 extraction using
  `famdb.py ... families Amphibia -a -d -f fasta_name --include-class-in-name`,
  and concatenation. Names vary between `TELib_80.fasta`, `TELib_80`,
  `Master_TElib.fasta` and `Master_TELib.fasta`; no final rename to
  `dedupe_telib.fasta` is recorded.
- Lines 149–232 describe a dnaPipeTE parameter scan for two SRXs. Its example
  uses eight CPUs, variable coverage/threshold/sample count, and
  `master_dfam_rmodeler_rmueller_TELib.fasta`. It differs from the fixed
  production-style job above. Line 235 leaves final parameters blank.

The note supplies previously missing command descriptions. It is not a submitted
job record, and the scan configuration must not replace the retained analysis
configuration.

## What remains unresolved

1. **Library contents and lineage.** Neither `dedupe_telib.fasta` nor the
   referenced master-library FASTAs were recovered in the searched locations.
   The note and CD-HIT job narrow the search to the original
   `/nfs/home/jlamb/TE_libs` directory, but classifier output merging, manual
   curation and the final library rename remain unlinked.
2. **Runtime identity.** The jobs name Conda environments and container files,
   rather than providing version inventories or image checksums. Recover
   `tetools_1.87.sif`, `dnapipte.img`, environment records and actual tool
   versions before claiming a complete upstream environment.
3. **Execution-to-output correspondence.** The jobs and note preserve
   configurations. Submitted job IDs, argument lists, scheduler output and
   per-species provenance are still needed to identify which configuration
   produced each retained artifact. The standalone RepeatMasker jobs accept
   any FASTA; they do not establish whether the retained Trinity `.align`
   files came from those jobs or dnaPipeTE's internal annotation.
4. **Read handoff.** The original BWA jobs end at nuclear BAM. The note gives
   generic FASTQ conversion examples, but no verified command maps that BAM
   to the exact `${SRX}_nuclear_1.fastq` and mate 2 inputs.
5. **LTR mapping and depth.** `Mass_map.sh` uses genome indices. The per-LTR
   mapping candidate uses a Necturus path. Their relationship to the retained
   Desmognathus depth files requires job/input evidence.

`Calculate_depth.sh` is especially unsuitable as a replacement for the current
analysis: lines 58–70 use comma-delimited parsing and the first field as a
numeric coordinate, while lines 79–85 sum terminal means and omit inclusive
`+1` lengths. Its calculation differs from the retained length-weighted
terminal/internal definition. Preserve its source as evidence; do not rerun it
over the retained release or use it to revise scientific outputs.

The next recovery target is the original HPC home/scratch storage or an exact
backup of it: TE library files, container images, classifier outputs, submitted
scripts and job logs. No HPC discovery workflow is enabled by this audit, and
the active local runtime remains Pixi.
