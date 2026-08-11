# Genome Biology and Evolution requirements

**Verified:** 2026-07-29

**Target assumed here:** GBE **Article**; final classification remains an editorial decision.

**Evidence boundary:** Only current first-party GBE and Oxford University Press (OUP) sources are used. Journal-specific GBE instructions take precedence over publisher-wide OUP guidance.

## Submission and Methods

- An Article permits a maximum **10,000 words of text excluding references**, **250 words in the abstract**, **150 words in the significance statement**, **100 references**, and **8 combined tables plus figures**.
- Required Article order is: Abstract; Significance statement; Introduction; Results; Discussion; Materials and Methods; Data Availability.
- At initial submission, the main manuscript is one PDF with main-text tables and figures inline, page numbers, and continuous line numbers. Supplementary material is uploaded separately.
- The manuscript must be double-spaced (at least 6 mm between lines), with 25 mm margins; footnotes should be avoided and all PDF fonts embedded.
- At revision, text, tables, and figures are submitted as separate files.
- For field-collected samples, GBE expects authors to follow Nagoya Protocol principles and invites an appropriate statement in Materials and Methods or Acknowledgements.
- For animal experiments, authors must state whether institutional and national care-and-use guidelines were followed.
- GBE specifies the placement of Materials and Methods but does **not** prescribe internal Methods subsection names or an analysis-specific reporting checklist. The repository's concise section architecture can therefore follow the study design as long as the final manuscript retains the required Article order.

Source: [GBE Instructions to Authors](https://academic.oup.com/GBE/pages/General_Instructions) (Types of Manuscripts; Preparation of Manuscripts; Other Information; Online Submission).

### AI-use disclosure

- GBE says AI use is strongly discouraged.
- Any AI use—including content or image creation, complex code generation, data processing, proofreading, or language modification—must be disclosed both in the **cover letter** and in **Methods and/or Acknowledgements**.
- Because this repository has used AI-assisted research and coding, maintain an accurate human-reviewed record of what was used and for what purpose; the authors must decide the final disclosure language.

Source: [GBE Instructions to Authors](https://academic.oup.com/GBE/pages/General_Instructions) (Publication Ethics and the Use of AI).

## Main figures

| Item | Current GBE requirement |
|---|---|
| Count | Article maximum is **8 total display items: tables plus figures**, not 8 of each. |
| Initial submission | Embed main figures in the single manuscript PDF. |
| Revision | Submit each figure separately, on a separate page/file, at roughly final magnification. |
| Accepted formats | `.jpg`, `.gif`, `.tif`, `.pdf`, `.eps`. |
| Font | Use a sans-serif font such as Arial or Helvetica and keep font size uniform within each figure. |
| Minimum label size | Never below **6 pt at final magnification**. |
| Color figures and photographs | At least **300 pixels per inch at final printed size**. |
| Black-and-white line drawings | At least **600 pixels per inch at final printed size**. |
| Color space | Submit color art in **CMYK, not RGB**, for revision/publication-quality files. |
| Color charges | Color figures are published online in color at no author charge. |
| In-text form | Cite as “Figure 1” in text; abbreviate as “Fig. 1.” only in the legend. |
| Legends | Place at the end of the manuscript. Start with a sentence giving the key insight, then add details; define all symbols and abbreviations. |
| Alt text | Required for every image, figure, illustration, and photograph in the **main article**. Put it directly below its legend and prefix it with `Alt text:`. |

Source: [GBE Instructions to Authors](https://academic.oup.com/GBE/pages/General_Instructions) (Tables; Illustrations; Electronic Submission of Figures; Figure Legends; Color Figures; Figure Accessibility and Alt Text).

### What GBE does not prescribe

The current GBE instructions do not specify an exact ordinary-figure width or height, panel limit, panel-label style, palette, line weight, or plotting theme. The listed 89 mm single-column and 185 mm double-column widths apply explicitly to **revised tables**; using them as figure design targets is a practical project convention, not a stated GBE figure rule.

### OUP design guidance, subordinate to GBE

OUP labels its general artwork guide as tips rather than strict rules. Useful recommendations are:

- provide a multipanel figure as one file and keep its caption in the manuscript;
- prefer uncompressed TIF for raster art and PDF/EPS/SVG with embedded fonts for vector art; because GBE does not list SVG, PDF or EPS is the safer vector submission;
- do not upscale raster files to manufacture resolution;
- use text of at least 7 pt and lines 0.25–1 pt thick;
- avoid pale colors and red–green distinctions; never use color as the only information channel;
- keep font, color, and texture use consistent across figures.

The generic OUP guide allows RGB or CMYK, but the journal-specific GBE rule explicitly requires **CMYK**, so CMYK controls for GBE revision files.

Source: [OUP journals—guidance for preparing artwork (PDF)](https://static.primary.prod.gcms.the-infra.com/static/site/journals/document/images-author-guidance.pdf?node=1bf05d0b2fbd9c529a23&version=490455%3A30c2211aa70bba63a5ee).

OUP's accessibility guide further recommends, where feasible:

- 12–14 pt figure text;
- at least 3:1 contrast for graphical objects, 4.5:1 for normal text, and 3:1 for large text;
- redundant encodings such as labels, point shapes, patterns, or line types rather than color alone;
- concise alt text under 100 words, preferably about 25–30 words, without duplicating the legend.

These are accessibility recommendations, not replacements for GBE's hard 6 pt minimum.

Source: [OUP Making Figures Accessible, Journals Edition (PDF)](https://static.primary.prod.gcms.the-infra.com/static/umbrella/document/Figures_accessibility_journals_edition_v2.2_2024.pdf?node=aef6670590405d43ef10).

## Supplementary figures and files

- All supplementary information—including supplementary tables and supplementary figure legends—must be uploaded separately and not placed in the main manuscript.
- Initial submission asks for **one supplementary-information PDF** containing supplementary text, tables, and figures, formatted to fit a standard page.
- GBE states no numeric supplementary-figure cap on the cited instructions page.
- Publisher-wide OUP guidance additionally says supplementary material must be cited in the main text, self-explanatorily named, consistent with manuscript styling, browser-compatible, online-only, and no larger than 2 MB per file. It is not copyedited.
- Supplementary material should enhance the article without being necessary to understand it. Large files are better placed in a stable public repository; GBE asks the Data Availability statement to point to a public repository, ideally one assigning a DOI, for large supplementary files.
- GBE's explicit alt-text requirement is worded for images in the **main article**. Accessible design and alt text for supplementary figures remain good practice, but the page does not expressly mandate supplementary alt text.

Sources: [GBE Instructions to Authors](https://academic.oup.com/GBE/pages/General_Instructions) (Types of Manuscripts; Online Submission); [OUP Preparing and Submitting Your Manuscript](https://academic.oup.com/pages/for-authors/journals/preparing-and-submitting-your-manuscript) (Supplementary Material).

## Data, alignments, and code

- Public release of **all data underlying the paper**, where ethically possible, is a condition of GBE publication; GBE asks authors to follow FAIR principles.
- Underlying data must be in the manuscript/supporting files or deposited in a public repository whenever possible.
- New sequence data must be deposited in **GenBank/DDBJ/EMBL**.
- Every sequence alignment used must be supplied as supplementary data or deposited in a repository such as Figshare or Dryad.
- Code generated for the manuscript should be deposited on a persistent site such as **Zenodo**.
- A **Data availability** statement is mandatory in the Article end matter. It must explain access and provide links, accession numbers, and/or DOIs.
- Every public dataset must be fully cited in the reference list with a persistent accession or identifier. GBE requests Author(s), Year, Title, Repository/Publisher, and Identifier, with the temporary `[dataset]` tag at the beginning of the reference.

Source: [GBE Instructions to Authors](https://academic.oup.com/GBE/pages/General_Instructions) (Availability of Data and Materials; Data Availability Statement; Data Citation).

### Repository-specific consequences

These are implementation decisions derived from the requirements, not additional GBE rules:

- Under the current recommended four-figure main set, four display slots remain
  for main-text tables or a scientifically essential replacement. Recalculate
  this balance if figure placement changes; the hard cap is eight combined
  main-text figures plus tables.
- Retain editable/vector masters and produce CMYK publication derivatives at final size; verify embedded fonts, label size, resolution, and color space rather than relying on filename extensions.
- Keep one legend and one concise alt-text entry for every main figure.
- Package narrative supplementary methods, figures, legends, and compact tables into one submission PDF; archive larger datasets and scientific files in a versioned DOI-issuing repository.
- Archive the exact input/output artifacts underlying the published analyses: accession crosswalks, TE libraries and summaries, RepeatMasker landscapes, exact alignments and trees, microscopy measurements and review decisions, analysis-ready species tables, statistical outputs, scripts, environments, and a checksum/provenance manifest.
- Nothing in the cited instructions requires expensive analyses to be rerun for submission. For HPC-only dnaPipeTE, RepeatModeler, RepeatMasker, and related products that cannot be rerun, preserve and release the exact existing outputs, parameters, software/container versions, provenance, and downstream local transformations.

## Official sources

All sources verified 2026-07-29:

1. [Genome Biology and Evolution — Instructions to Authors](https://academic.oup.com/GBE/pages/General_Instructions).
2. [OUP — Preparing and submitting your manuscript](https://academic.oup.com/pages/for-authors/journals/preparing-and-submitting-your-manuscript).
3. [OUP journals—guidance for preparing artwork (PDF)](https://static.primary.prod.gcms.the-infra.com/static/site/journals/document/images-author-guidance.pdf?node=1bf05d0b2fbd9c529a23&version=490455%3A30c2211aa70bba63a5ee).
4. [OUP — Making Figures Accessible, Journals Edition (PDF)](https://static.primary.prod.gcms.the-infra.com/static/umbrella/document/Figures_accessibility_journals_edition_v2.2_2024.pdf?node=aef6670590405d43ef10).
