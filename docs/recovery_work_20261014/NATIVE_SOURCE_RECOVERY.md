# Native microscopy source recovery — paused checkpoint

This records work completed on 14 October 2026 before Jake requested a pause.
The [restart handoff](../handoffs/20261014_SOURCE_RECOVERY_PAUSED.md) is the task
authority. Latest TIFF-page reader verification and final operator documentation
remain incomplete.

## Recovered sources

Original Process 413 and 414 VSI files and their ETS companions were recovered
from the connected easystore archive. Native metadata records 40× acquisition,
0.12 µm pixels, RGB uint8 storage and exposure times of 45.296 and 43.869 ms,
respectively. Stored gamma was 1 in both. These settings do not establish a
physical intensity response, dark offset or an independent DNA-content standard.

Two native green-channel crops recover pixels across former tile seams:

| Case | Source | Recovered rectangle, half-open slide coordinates | New context |
| --- | --- | --- | --- |
| `N815f697d43b9efed92c4` | Process 414 | `[41975, 20416, 42212, 20595)` | 15,168 pixels north of the retained tile seam |
| `N3291d5e0a6bc3fe2af23` | Process 413 | `[45553, 12144, 45824, 12352)` | 17,344 pixels south of the retained tile seam |

The original packet overlap, full retained tile intersection and translated masks
match exactly. The Process 414 crop shows a round dark nucleus extending above
the original mask; the Process 413 crop still has touching/overlapping cell context
and an uncertain boundary. These are separate diagnostic crops. Original masks
and packet pixels remain unchanged; no replacement contour or decision was made.

Four vendor RGB BigTIFF exports were also recovered:

| Process | Actual RGB page dimensions | Released crops compared | Exposure recorded in copied metadata |
| --- | --- | --- | --- |
| 314 | 78,126 × 73,075 | 12 | 45.183 ms |
| 315 | 64,230 × 56,260 | 8 | 44.218 ms |
| 320 | 95,151 × 52,998 | 5 | 45.467 ms |
| 321 | 75,923 × 40,147 | 10 | 44.315 ms |

Files are under
`/run/media/jake/easystore/Dried_blood_nopmount/Process_<ID>/Process_<ID>_raw.ome.btf`.
Vendor export reports directly link each input VSI basename to its output. The
original target VSI files and their absolute input paths were not located.

All 35 retained target crops match the exported native green pixels exactly, with
equal uint8 dtype and shape. The read coordinates came only from the original
packet's tile origins and crop bounds. No alternative mapping, intensity
transformation, whole-slide decode, ROI proposal, decision or IOD calculation was
performed. Source size, mtime, ctime, device and inode remained unchanged. Full
large-payload SHA-256 values were not computed.

## Reader and metadata limits

Each target file has one actual uncompressed tiled RGB TIFF page. Embedded copied
OME metadata refers to additional pyramid images/IFDs, and Bio-Formats OME-TIFF
metadata reads failed for all four. Those failures are retained. The explicit
TIFF-page path uses actual page headers and selected tile byte ranges; it does
not repair source metadata or claim a valid OME-series interpretation.

The latest worker recorded 43 passing reader tests with no skips and successful
lint, format and whitespace checks. Independent verification of this exact
TIFF-page/stat-only revision was interrupted for the pause and has no final
verdict. Earlier root checks and independent PASS cover the previous reader.

Sources were on a writable NTFS3 mount. Reading them without source writes is
not equivalent to mounting them read-only. Native settings and exact crop
equivalence support source lineage; they do not validate biological boundaries,
camera linearity or absolute genome-size calibration.

Detailed crops, per-case byte-range receipts, commands, clock records and raw
logs remain in ignored local evidence directories listed in the handoff. They
are preserved on this computer and excluded from Git. No retained scientific
data, Methods, figures or release outputs were regenerated in this phase.
