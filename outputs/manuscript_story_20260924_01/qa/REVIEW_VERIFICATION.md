# Isolated review verification

Verified 24 September 2026 (US/Pacific).

- Five collages and twelve isolated panels rendered in PDF, PNG, and CMYK TIFF, with embedded Arial in all PDFs.
- `verify_run.py verify` checked 16,793 original files against the pre-render baseline: no changes.
- Private Chromium loaded all 17 figure images across the scrollable gallery and found 19 review cards, including both tables.
- All 54 distinct figure/table download URLs returned successful HTTP responses.
- A figure note and a table note were entered through the UI, saved by the separate server, and recovered after reload. The test text was then cleared through the UI and the empty state verified after reload. Cleared entries retain their normal revision history.
- Both workbook sheets were rendered and visually inspected. The workbook error scan returned zero matches.
- The new gallery was opened and its title/content verified in a separate visible Codex browser tab at http://127.0.0.1:8767/#figure-2. The original gallery remains at port 8766.

The persistent local review server uses user service `desmognathus-manuscript-story-review-20260924.service`. Its document root and note store are inside this run. `browser_figure_2.png` records the private-browser display check; it is distinct from opening the actual user-visible review tab.

Figure-specific numeric checks and export reports accompany the figures. These checks establish data preservation and review functionality, not final manuscript acceptance or resolution of the scientific limitations stated in the legends.
