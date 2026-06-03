# Cross-Package Validation Data

This directory contains small, inspectable datasets that can be used to build a
larger cross-package validation suite for biodiversity software.

The aim is to separate three things that are often mixed together:

- the input data used for validation,
- the convention being tested,
- the package-specific call that implements that convention.

The files in `data/` are intentionally modest in size. They include uneven
abundances, all-zero absences within taxa, singletons, doubletons, shared
incidence structure, and samples with different total abundance. That makes
them useful for checking richness, Shannon-family estimators, Simpson-family
indices, coverage diagnostics, Bray-Curtis, Jaccard, Sorensen-Dice, Hellinger,
and related transformations across Julia, R, and Python.

## Files

- `data/community_counts.csv`: sample-by-taxon community matrix for alpha and
  beta diversity validation.
- `data/pairwise_vectors.csv`: hand-checkable paired vectors for pairwise
  similarity and dissimilarity formulas.
- `cross_package_manifest.csv`: package-call manifest mapping each validation
  task to the intended convention in DiversityIndices.jl, R vegan, SciPy, and
  scikit-bio-style APIs.

These are validation inputs, not ecological observations. They are simulated to
exercise convention-sensitive cases and should be treated as a reproducibility
corpus.
