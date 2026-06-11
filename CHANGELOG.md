# Changelog

All notable changes to `DiversityAndDissimilarity.jl` will be documented in this file.

## Unreleased

- Corrected and stabilized executable documentation examples, and added
  ordinary tests that mirror source docstring and data-input examples.
- Added `Validated{T}` wrapper type and `validate()` function to separate
  input validation from computation. Calling `validate(community)` checks all
  preconditions once and returns a wrapper; subsequent calls on the wrapper
  skip per-call validation. The default raw-data pathway remains fully validated.
- Refactored all major computation methods (`richness`, `shannon_entropy`,
  `alpha_diversity`, `bray_curtis_distance`, `jaccard_distance`,
  `hellinger_distance`) to dispatch on `Validated{T}` for the computation
  kernel, with the raw-data dispatch validating and delegating.
- Updated Julia benchmark to report both safe (per-call validation) and
  pre-validated pathways as separate CSV rows, enabling fair cross-language
  comparison with Python/NumPy and R/vegan.
- Added Pielou evenness and Fisher alpha diversity metrics.
- Added labeled pairwise similarity/dissimilarity/distance helpers.
- Added vegan migration, benchmark, and examples documentation.
- Added benchmark report generation with hardware/runtime metadata and figures.
- Improved dense community-matrix performance for richness, Shannon entropy,
  alpha-diversity summaries, Bray-Curtis, Hellinger, and Jaccard.

## 0.1.0

- Initial public API for alpha diversity, entropy estimators, richness
  estimation, sample coverage, and pairwise similarity/dissimilarity indices.
