# API Reference

This page is a compact index of the public API. Full docstrings are grouped by
topic on the diversity, dissimilarity, and framework pages.

## Core Types

```@docs
DiversityIndex
AlphaDiversityIndex
PairwiseIndex
```

## Data Preparation

```@docs
counts
community_matrix
proportions
alpha_diversity
```

## Diversity Indices

- [`ShannonEstimator`](@ref)
- [`Plugin`](@ref)
- [`MillerMadow`](@ref)
- [`HausserStrimmer`](@ref)
- [`Basharin`](@ref)
- [`AddGamma`](@ref)
- [`ChaoShen`](@ref)
- [`Richness`](@ref)
- [`Shannon`](@ref)
- [`Renyi`](@ref)
- [`Tsallis`](@ref)
- [`Simpson`](@ref)
- [`GiniSimpson`](@ref)
- [`GreenbergDiversityIndex`](@ref)
- [`LinguisticDiversityIndex`](@ref)
- [`InverseSimpson`](@ref)
- [`Hill`](@ref)
- [`Chao1`](@ref)
- [`ACE`](@ref)
- [`SampleCoverage`](@ref)
- [`PielouEvenness`](@ref)
- [`FisherAlpha`](@ref)

## Diversity Operations

- [`entropy`](@ref)
- [`entropy_variance`](@ref)
- [`entropy_confint`](@ref)
- [`diversity`](@ref)
- [`effective_diversity`](@ref)
- [`alpha_diversity`](@ref)
- [`richness`](@ref)
- [`shannon`](@ref)
- [`shannon_entropy`](@ref)
- [`shannon_variance`](@ref)
- [`shannon_confint`](@ref)
- [`shannon_diversity`](@ref)
- [`bootstrap`](@ref)
- [`jackknife`](@ref)
- [`renyi`](@ref)
- [`renyi_entropy`](@ref)
- [`renyi_diversity`](@ref)
- [`tsallis`](@ref)
- [`tsallis_entropy`](@ref)
- [`tsallis_diversity`](@ref)
- [`hill_number`](@ref)
- [`chao1`](@ref)
- [`ace`](@ref)
- [`sample_coverage`](@ref)
- [`pielou_evenness`](@ref)
- [`fisher_alpha`](@ref)
- [`simpson_index`](@ref)
- [`gini_simpson_index`](@ref)
- [`greenberg_diversity_index`](@ref)
- [`linguistic_diversity_index`](@ref)
- [`index_of_linguistic_diversity`](@ref)
- [`inverse_simpson_index`](@ref)

## Similarity And Dissimilarity Indices

- [`Jaccard`](@ref)
- [`SorensenDice`](@ref)
- [`Overlap`](@ref)
- [`BrayCurtis`](@ref)
- [`Ruzicka`](@ref)
- [`TotalVariation`](@ref)
- [`Manhattan`](@ref)
- [`Euclidean`](@ref)
- [`Canberra`](@ref)
- [`Hellinger`](@ref)
- [`Chord`](@ref)
- [`Bhattacharyya`](@ref)
- [`KullbackLeibler`](@ref)
- [`ShannonDifference`](@ref)
- [`JensenDifference`](@ref)
- [`JensenShannon`](@ref)
- [`MorisitaHorn`](@ref)

## Similarity And Dissimilarity Operations

- [`similarity`](@ref)
- [`dissimilarity`](@ref)
- [`distance`](@ref)
- [`labeled_similarity`](@ref)
- [`labeled_dissimilarity`](@ref)
- [`labeled_distance`](@ref)
- [`jaccard_similarity`](@ref)
- [`jaccard_index`](@ref)
- [`jaccard_distance`](@ref)
- [`sorensen_index`](@ref)
- [`sorensen_dice_index`](@ref)
- [`sorensen_distance`](@ref)
- [`sorensen_dice_distance`](@ref)
- [`sorensen_dice_dissimilarity`](@ref)
- [`bray_curtis_distance`](@ref)
- [`bray_curtis_dissimilarity`](@ref)
- [`overlap_similarity`](@ref)
- [`overlap_distance`](@ref)
- [`ruzicka_similarity`](@ref)
- [`quantitative_jaccard_similarity`](@ref)
- [`ruzicka_distance`](@ref)
- [`quantitative_jaccard_distance`](@ref)
- [`total_variation_distance`](@ref)
- [`manhattan_distance`](@ref)
- [`euclidean_distance`](@ref)
- [`canberra_distance`](@ref)
- [`hellinger_distance`](@ref)
- [`chord_distance`](@ref)
- [`bhattacharyya_coefficient`](@ref)
- [`bhattacharyya_distance`](@ref)
- [`kullback_leibler_divergence`](@ref)
- [`shannon_difference`](@ref)
- [`jensen_difference`](@ref)
- [`jensen_shannon_similarity`](@ref)
- [`jensen_shannon_divergence`](@ref)
- [`jensen_shannon_distance`](@ref)
- [`morisita_horn_similarity`](@ref)
- [`morisita_horn_distance`](@ref)

## Framework And Validation

- [`index_metadata`](@ref)
- [`index_family`](@ref)
- [`input_mode`](@ref)
- [`output_mode`](@ref)
- [`is_finite`](@ref)
- [`is_metric`](@ref)
- [`is_triangular`](@ref)
- [`is_nonnegative`](@ref)
- [`is_bounded`](@ref)
- [`is_pseudometric`](@ref)
- [`is_quasimetric`](@ref)
- [`is_metametric`](@ref)
- [`is_semimetric`](@ref)
- [`is_premetric`](@ref)
- [`is_supermetric`](@ref)
- [`is_similarity`](@ref)
- [`is_dissimilarity`](@ref)
- [`is_dissimiliarty`](@ref)
- [`index_range`](@ref)
- [`index_bounds`](@ref)
- [`requires_probabilities`](@ref)
- [`supports_matrix_kernel`](@ref)
- [`reference_cases`](@ref)
- [`validate_reference_cases`](@ref)
- [`estimator_report`](@ref)
- [`diversity_audit`](@ref)
- [`uncertainty_audit`](@ref)
