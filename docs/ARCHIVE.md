# Pre-cleanup archive and recovery

The complete pre-cleanup checkout is stored at:

```text
/home/jake/Projects/Desmognathus_TE_archive
```

It was copied on 2026-08-10 before destructive cleanup and preserves the active
branch, full Git directory, dirty working tree, ignored data, HPC outputs,
intermediates, caches, generated galleries, retired analyses, and old
environments/specifications. The copy is on the same Btrfs filesystem and was
created with copy-on-write reflinks, so it is an independent namespace without
initially duplicating all physical blocks.

Archive controls:

- `PRE_CLEANUP_ARCHIVE.md` records the branch, commit, and preservation rule.
- `ARCHIVE_MANIFEST.csv` records every preserved working file's relative path,
  byte size, SHA-256, and retention class.
- `ARCHIVE_MANIFEST_SUMMARY.json` records aggregate counts and a manifest hash.

Do not develop the paper in the archive. To recover an individual file, copy it
from its exact relative path into a temporary location, compare its manifest
hash, and then decide whether it belongs in the active three-component model.
Do not bulk-restore an old top-level directory into this repository; doing so
would recreate the ambiguity this cleanup removes.

Nothing archive-only is implied to be reproducible. The archive is the evidence
and historical recovery layer for computations—especially HPC stages—that
cannot presently be recreated from complete inputs and software provenance.
