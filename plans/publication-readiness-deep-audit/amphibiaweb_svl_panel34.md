# AmphibiaWeb SVL audit for the TE resource panel (34 species)

Date: 2026-07-13

## Verdict

AmphibiaWeb does **not** provide a standardized adult-mean SVL for any of the 34 study species. Direct account text contains explicit SVL for five species: `organi`, `santeetlah`, `valentinei`, `welteri`, and `wrighti`. Four are ranges or sex-specific maxima; the `organi` account gives a 40–60 mm “average” range that conflicts with the curated supplemental maximum of 30.8 mm and should remain on QC hold. Eight additional accounts provide total or unspecified body length only, not SVL. Three narrative accounts contain no size measurement, and 18 records are taxon stubs with no narrative measurement account.

Do not relabel total length as SVL, do not convert a maximum or range midpoint into an adult mean without an explicit analysis decision, and do not treat a dataset merely linked from AmphibiaWeb's trait-database page as verified AmphibiaWeb account data.

## Scope and method

- Species universe: `path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv` (34 rows).
- Primary source: the AmphibiaWeb species-account query URLs listed below.
- Access limitation: the live site returned a Cloudflare human-verification page on 2026-07-13. Account text was therefore checked against Arquivo.pt captures from November 2023, with later Internet Archive CDX metadata used to confirm that the taxon URLs persisted through 2024–2025.
- `account_stub` means the AmphibiaWeb page supplies taxonomy/distribution scaffolding but no narrative account from which SVL can be extracted.
- `narrative_no_svl` means a substantive account exists but it contains no SVL statement.
- `length_only_not_svl` is retained as evidence but is intentionally blank for analysis-ready SVL.
- Supplemental primary-literature values are transcribed only from the repository's existing curated `path_analysis/data/templates/literature_trait_extraction.csv`; they are not attributed to AmphibiaWeb.

## 34-species extraction

| panel species | AmphibiaWeb taxon | account result | direct account evidence / interpretation | direct SVL (mm) | sex / statistic | analysis status | AmphibiaWeb source |
|---|---|---|---|---:|---|---|---|
| abditus | Desmognathus abditus | narrative_no_svl | Raffaelli narrative present; no SVL or body-length measurement |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=abditus) |
| adatsihi | Desmognathus adatsihi | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=adatsihi) |
| aeneus | Desmognathus aeneus | length_only_not_svl | Adults reported as 38–57 mm from snout to tail tip: explicitly total length |  | adult range | exclude as SVL | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=aeneus) |
| amphileucus | Desmognathus amphileucas | account_stub; taxonomy mismatch | AmphibiaWeb uses `amphileucas`; no narrative measurement account |  |  | missing; crosswalk required | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=amphileucas) |
| anicetus | Desmognathus anicetus | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=anicetus) |
| apalachicolae | Desmognathus apalachicolae | length_only_not_svl | French account gives 11 cm without naming SVL; treat as unspecified/likely total length |  | unspecified | exclude as SVL | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=apalachicolae) |
| aureatus | Desmognathus aureatus | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=aureatus) |
| auriculatus | Desmognathus auriculatus | length_only_not_svl | Account gives 80–164 mm length and says tail is half total body length; not SVL |  | range | exclude as SVL | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=auriculatus) |
| bairdi | Desmognathus bairdi | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=bairdi) |
| balsameus | Desmognathus balsameus | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=balsameus) |
| campi | Desmognathus campi | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=campi) |
| carolinensis | Desmognathus carolinensis | length_only_not_svl | French account gives 11 cm without naming SVL |  | unspecified | exclude as SVL | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=carolinensis) |
| catahoula | Desmognathus catahoula | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=catahoula) |
| cheaha | Desmognathus cheaha | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=cheaha) |
| conanti | Desmognathus conanti | length_only_not_svl | French account gives 12.7 cm without naming SVL |  | unspecified | exclude as SVL | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=conanti) |
| fuscus | Desmognathus fuscus | narrative_no_svl | Substantive account present; no SVL or body-length measurement |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=fuscus) |
| gvnigeusgwotli | Desmognathus gvnigeusgwotli | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=gvnigeusgwotli) |
| intermedius | Desmognathus intermedius | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=intermedius) |
| kanawha | Desmognathus kanawha | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=kanawha) |
| lycos | Desmognathus lycos | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=lycos) |
| marmoratus | Desmognathus marmoratus | length_only_not_svl | French account gives 15 cm and describes tail relative to total length; not SVL |  | unspecified | exclude as SVL | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=marmoratus) |
| mavrokoilius | Desmognathus mavrokoilius | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=mavrokoilius) |
| monticola | Desmognathus monticola | length_only_not_svl | French account gives 15 cm and tail as half total length; not SVL |  | unspecified | exclude as SVL | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=monticola) |
| ocoee | Desmognathus ocoee | narrative_no_svl | Raffaelli narrative present; no SVL or body-length measurement |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=ocoee) |
| orestes | Desmognathus orestes | length_only_not_svl | French account gives 11 cm without naming SVL |  | unspecified | exclude as SVL | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=orestes) |
| organi | Desmognathus organi | direct_svl_qc_hold | “snout to vent length on average measures 40 to 60 millimeters” | 40–60 | sex unspecified; ambiguous “average” range | **QC hold**: conflicts with curated 30.8 mm maximum | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=organi) |
| pascagoula | Desmognathus pascagoula | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=pascagoula) |
| perlapsus | Desmognathus perlapsus | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=perlapsus) |
| santeetlah | Desmognathus santeetlah | direct_svl | “up to 45 mm SVL in adult females and 55 mm SVL in adult males” | F max 45; M max 55 | sex-specific adult maximum | usable as maxima, not mean | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=santeetlah) |
| tilleyi | Desmognathus tilleyi | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=tilleyi) |
| valentinei | Desmognathus valentinei | direct_svl | Adults reported as 42–72 mm in males and 42–62 mm in females | M 42–72; F 42–62 | sex-specific adult range | usable as ranges, not mean | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=valentinei) |
| valtos | Desmognathus valtos | account_stub | No narrative measurement account |  |  | missing | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=valtos) |
| welteri | Desmognathus welteri | direct_svl | Adults reported with SVL between 50 and 95 mm | 50–95 | adult range; sex unspecified | usable as range, not mean | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=welteri) |
| wrighti | Desmognathus wrighti | direct_svl | Species reported as 17–20 mm SVL in both males and females | 17–20 | both sexes; range | usable as range, not mean | [account](https://amphibiaweb.org/cgi/amphib_query?where-genus=Desmognathus&where-species=wrighti) |

## Existing supplemental primary-literature SVL records

These are separate from AmphibiaWeb and are included only to expose already-curated alternatives for rows where the direct account is missing, total-length-only, or suspect.

| species | supplemental SVL | stage / statistic | source |
|---|---:|---|---|
| aeneus | 29.3 mm max | adult clutch-attending females, n=9 | Bruce 2014, doi:10.1643/CE-13-165 |
| apalachicolae | 45.43 mm mean | transformed specimens, n=23 | Graham et al. 2010, Table 4 |
| auriculatus | 48.46 mm mean | transformed specimens, n=33 | Graham et al. 2010, Table 4 |
| cheaha | 40–80 mm | adult range | Pyron et al. 2023, doi:10.1093/sysbio/syac065 |
| monticola | 44 mm mean; 19–72 mm observed | all specimens | Pyron et al. 2023, doi:10.1093/sysbio/syac065 |
| ocoee | 52 mm max | adult, n=1,210 | Bruce 2016, doi:10.1643/CE-14-204 |
| orestes | 38.9 ± 1.2 mm mean | adult, n=10 | Kozak et al. 2005, doi:10.1111/j.0014-3820.2005.tb01069.x |
| pascagoula | 47 mm mean; 32–56 mm range | adult | Pyron et al. 2022, doi:10.11646/zootaxa.5133.1.3 |
| perlapsus | 22–56 mm | adult, n=77 | Pyron and Beamer 2022, doi:10.11646/zootaxa.5190.2.3 |
| valentinei | 62 mm max | adult | Pyron et al. 2022, doi:10.11646/zootaxa.5133.1.3 |
| welteri | 95 mm max | adult | Felix 2001 thesis |
| wrighti | F 32.24 mm max; M 31.04 mm max | adults, n=132 | Kessler et al. 2024, doi:10.1655/Herpetologica-D-23-00049 |

## Taxonomy and QC flags

1. The panel uses `Desmognathus amphileucus`, whereas AmphibiaWeb's account is `Desmognathus amphileucas`. Preserve the study name and add an explicit source-name crosswalk rather than silently changing the panel taxonomy.
2. The AmphibiaWeb `organi` statement (40–60 mm SVL) is not credible as analysis-ready without returning to its cited source; it conflicts with the repository's curated 30.8 mm maximum from the Virginia Herpetological Society account. Keep it as reported evidence with a QC hold.
3. AmphibiaWeb's direct values are ranges/maxima, not a common estimand. A final body-size table must declare whether it uses adult mean, sex-specific maximum, pooled range, or another statistic; mixing these silently would be biologically and statistically invalid.
