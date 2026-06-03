# Reference Dataset Candidates

These are online datasets or documented example matrices that can extend the
cross-package validation suite. The most useful candidates are small enough to
vendor as CSV fixtures, have stable upstream package examples, and cover
indices already implemented in `DiversityAndDissimilarity.jl`.

## Highest Priority

| Dataset | Source | Size | Useful checks | Reference path |
|---|---|---:|---|---|
| scikit-bio diversity tutorial matrix | <https://scikit.bio/docs/latest/diversity.html> | 6 samples x 7 taxa | Richness, Bray-Curtis distance matrix, later phylogenetic examples if UniFrac/Faith PD are added | The documentation prints the count matrix, observed richness per sample, and Bray-Curtis distance matrix. This is the easiest immediate fixture because expected values are published directly in the docs. |
| `vegan::dune` | <https://rdrr.io/cran/vegan/man/dune.html> | 20 sites x 30 species | Richness, Shannon, Gini-Simpson, inverse Simpson, Bray-Curtis, Jaccard/Sorensen incidence, Hellinger after `decostand` | Use `vegan::diversity`, `vegan::specnumber`, `vegan::vegdist`, and `vegan::decostand`. The dataset is a canonical small community matrix used throughout vegan examples. |
| `vegan::varespec` | <https://rdrr.io/cran/vegan/man/varechem.html> | 24 sites x 44 species | Shannon, Simpson-family indices, Bray-Curtis, Hellinger, Euclidean/Manhattan after transformations | Use `vegan::diversity`, `vegan::vegdist`, and `vegan::decostand`. This is a good abundance/cover-valued companion to `dune`. |
| `vegan::BCI` | <https://rdrr.io/cran/vegan/man/BCI.html> | 50 plots x 225 species | Richness, Shannon, Fisher alpha, Chao1/coverage-style diagnostics, Bray-Curtis on larger count matrix | Use `vegan::diversity`, `vegan::fisher.alpha`, `vegan::specnumber`, and `vegan::vegdist`. Also overlaps with `biosampleR::BCI` and `biosampleR::calc_diversity_indices`. |

## Additional Candidates

| Dataset | Source | Size | Useful checks | Reference path |
|---|---|---:|---|---|
| `vegan::mite` | <https://rdrr.io/cran/vegan/man/mite.html> | 70 sites x 35 species | Bray-Curtis, Jaccard/Sorensen incidence, Hellinger, Shannon/Simpson-family row summaries | Use `vegan::vegdist`, `vegan::decostand`, and `vegan::diversity`. Good for sparse ecological count data. |
| `iNEXT::spider` | <https://rdrr.io/cran/iNEXT/man/spider.html> | 2 abundance vectors | Richness, sample coverage, Shannon effective diversity, Simpson/Hill-number checks | Use `iNEXT::DataInfo`, `iNEXT::estimateD`, or `iNEXT::iNEXT` for coverage-based and Hill-number reference outputs. |
| `biosampleR::BCI` | <https://rdrr.io/cran/biosampleR/man/calc_diversity_indices.html> | 50 plots x 225 species | Abundance, richness, Shannon, Simpson, Chao1 | Use `biosampleR::calc_diversity_indices(BCI)` as an independent R-package summary over the BCI count matrix. |

## Suggested Integration Plan

1. Vendor the scikit-bio tutorial matrix first because it has printed expected
   richness and Bray-Curtis values in the upstream documentation.
2. Export `dune`, `varespec`, `BCI`, and `mite` from vegan with a pinned vegan
   version and store both the raw CSVs and generated expected-value CSVs.
3. Generate expected values with small scripts in `validation/` rather than
   hard-coding package calls inside tests. This keeps CI independent of R and
   Python package availability.
4. Use `iNEXT::spider` for coverage/Hill-number validation once the package has
   matching coverage-based conventions documented.
