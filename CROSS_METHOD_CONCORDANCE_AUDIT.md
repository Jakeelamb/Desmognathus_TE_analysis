# Cross-Method Concordance Audit

This note compares read-based dnaPipeTE summaries against raw RepeatMasker
classification summaries built directly from `merged_repeatmasker_data.csv`.

## Source Basis

- `results/data/dnaPipeTE_order_breakdown.csv`
- `results/data/dnaPipeTE_superfamily_breakdown.csv`
- `results/data/merged_repeatmasker_data.csv`

## Bottom Line

- The main order-level TE features show at least moderate cross-method agreement, which supports treating them as biological rather than purely method-specific summaries.

## Strongest Order-Level Concordance

- `DIRS`: Spearman rho `0.918`, Pearson r `0.963`, n `21`
- `LTR`: Spearman rho `0.857`, Pearson r `0.861`, n `21`
- `order_pielou`: Spearman rho `0.834`, Pearson r `0.856`, n `21`
- `PLE`: Spearman rho `0.773`, Pearson r `0.551`, n `21`
- `Helitron`: Spearman rho `0.726`, Pearson r `0.715`, n `21`
- `retro_dna_logratio`: Spearman rho `0.691`, Pearson r `0.701`, n `21`
- `TIR`: Spearman rho `0.677`, Pearson r `0.680`, n `21`
- `SINE`: Spearman rho `0.621`, Pearson r `0.531`, n `21`
- `ltr_line_logratio`: Spearman rho `0.595`, Pearson r `0.809`, n `21`
- `LINE`: Spearman rho `0.558`, Pearson r `0.821`, n `21`

## Highest Species-Level Disagreement

- `aureatus`: order Bray disagreement `0.211`, |LTR:LINE| diff `1.170`, |Pielou| diff `0.130`
- `intermedius`: order Bray disagreement `0.176`, |LTR:LINE| diff `0.963`, |Pielou| diff `0.109`
- `amphileucus`: order Bray disagreement `0.176`, |LTR:LINE| diff `1.009`, |Pielou| diff `0.105`
- `aeneus`: order Bray disagreement `0.174`, |LTR:LINE| diff `0.874`, |Pielou| diff `0.092`
- `fuscus`: order Bray disagreement `0.169`, |LTR:LINE| diff `0.893`, |Pielou| diff `0.118`

## Files

- `results/data/cross_method/cross_method_order_concordance.csv`
- `results/data/cross_method/cross_method_superfamily_concordance.csv`
- `results/data/cross_method/cross_method_species_disagreement.csv`
