# Annotated Bibliography

This bibliography is organized around the methods implemented by
`DiversityAndDissimilarity.jl`. Local copies of downloaded materials are stored
under `notes/references/`. Papers or books that could not be downloaded
directly are listed in `notes/references/manual_downloads.md`; the BibTeX file
points to that manual list with the same keys used below.

## Diversity Indices

### Entropy, Effective Diversity, And Generalized Entropies

- `shannon1948mathematical`
  (`references/shannon1948_mathematical_theory_communication.pdf`) is the
  primary source for Shannon entropy. It supports `Shannon`, `Plugin`,
  `shannon_entropy`, `shannon`, and the entropy scale used by KL and
  Jensen-Shannon-style divergences.

- `hill1973diversity`
  (`references/manual_downloads.md#hill1973diversity`) introduces the Hill
  number notation that unifies richness, Shannon effective diversity, and
  inverse Simpson diversity. It supports `Hill`, `effective_diversity`,
  `hill_number`, and the package's order-0/1/2 identities.

- `jost2006entropy`
  (`references/manual_downloads.md#jost2006entropy`) clarifies the distinction
  between entropy and true diversity/effective numbers. It is the main
  interpretive reference for `shannon_diversity`, `renyi_diversity`,
  `tsallis_diversity`, and endpoint/bounds language for diversity values.

- `renyi1961measures`
  (`references/renyi1961_measures_entropy_information.pdf`) is the primary
  source for Renyi entropy. It supports `Renyi`, `renyi_entropy`, and
  `renyi_diversity`.

- `tsallis1988possible`
  (`references/tsallis1988_generalized_boltzmann_gibbs_statistics_springer.html`;
  PDF queued in `manual_downloads.md#tsallis1988possible`) is the primary
  source for Tsallis entropy. It supports `Tsallis`, `tsallis_entropy`, and
  `tsallis_diversity`; the package notes its deliberate base-scaling convention.

### Simpson Family, Evenness, Linguistic Diversity, And Fisher Alpha

- `simpson1949measurement`
  (`references/simpson1949_measurement_diversity_nature.html`; PDF queued in
  `manual_downloads.md#simpson1949measurement`) is the primary source for
  Simpson concentration. It supports `Simpson`, `GiniSimpson`,
  `InverseSimpson`, and the related convenience functions.

- `greenberg1956measurement`
  (`references/greenberg1956_measurement_linguistic_diversity.pdf`) defines
  Greenberg's linguistic diversity index as the probability that two randomly
  selected people have different mother tongues. It supports
  `GreenbergDiversityIndex`, `LinguisticDiversityIndex`,
  `greenberg_diversity_index`, and `index_of_linguistic_diversity`.

- `pielou1966measurement`
  (`references/manual_downloads.md#pielou1966measurement`) is the primary
  reference for Pielou evenness, implemented as Shannon entropy divided by the
  maximum entropy for the observed richness.

- `fisher1943relation`
  (`references/fisher_corbet_williams1943_logseries_species_abundance.pdf`)
  introduces Fisher's log-series alpha parameter. It supports `FisherAlpha` and
  `fisher_alpha`.

- `magurran2004measuring`
  (`references/manual_downloads.md#magurran2004measuring`) is a practical
  ecological reference for interpreting richness, Shannon, Simpson-family
  indices, evenness, and Fisher alpha. It also supports validation examples and
  user-facing language.

### Richness, Coverage, And Low-Sample Entropy Estimation

- `chao1984nonparametric`
  (`references/manual_downloads.md#chao1984nonparametric`) introduces the
  Chao1 lower-bound estimator for unseen classes. It supports `Chao1` and
  `chao1`.

- `chao1992coverage`
  (`references/manual_downloads.md#chao1992coverage`) develops sample-coverage
  richness estimation and supports the abundance-based coverage estimator
  `ACE`.

- `good1953population`
  (`references/good1953_population_frequencies_species.pdf`) is the classic
  Good-Turing unseen-frequency reference. It supports `SampleCoverage`,
  `sample_coverage`, and the unseen-mass component of `ChaoShen`.

- `chao2003nonparametric`
  (`references/chao_shen2003_unseen_species_shannon_springer.html`; PDF queued
  in `manual_downloads.md#chao2003nonparametric`) supports the `ChaoShen`
  estimator for Shannon entropy and its Good-Turing/Horvitz-Thompson
  correction for unseen species.

- `miller1955note`
  (`references/manual_downloads.md#miller1955note`) supports the
  `MillerMadow` first-order bias correction for Shannon entropy and the
  corresponding low-sample correction exposed for KL and Jensen-style
  divergences.

- `basharin1959statistical`
  (`references/basharin1959_entropy_estimation_download_failed.html`; PDF
  queued in `manual_downloads.md#basharin1959statistical`) supports the
  `Basharin` correction and entropy variance approximations.

- `laplace1812theorie`
  (`references/manual_downloads.md#laplace1812theorie`) and
  `jeffreys1946invariant`
  (`references/manual_downloads.md#jeffreys1946invariant`) support Bayesian
  pseudocount smoothing through `AddGamma`, especially `AddGamma(1)` for
  Laplace smoothing and `AddGamma(0.5)` for Jeffreys smoothing.

- `hausser2009entropy`
  (`references/hausser_strimmer2009_entropy_james_stein.pdf`) supports the
  `HausserStrimmer` shrinkage estimator.

- `paninski2003entropy`
  (`references/paninski2003_entropy_mutual_information.pdf`) surveys entropy
  estimation bias and supports the package's estimator documentation and
  warnings about plugin entropy under small samples.

- `nemenman2002entropy`
  (`references/nemenman2002_entropy_inference_revisited.pdf`) supports the
  Bayesian entropy-estimation discussion, including NSB-style motivation for
  Bayesian smoothing and information-functional estimation.

- `jiao2015minimax`, `wu2016minimax`, and `valiant2017estimating`
  (`references/jiao2015_minimax_information_functionals.pdf`,
  `references/wu_yang2016_minimax_entropy_large_alphabets.pdf`,
  `references/valiant_valiant2017_estimating_unseen.pdf`) support the minimax
  and unseen-species context for entropy, KL, JS, and related information
  functionals.

- `chao2014rarefaction`
  (`references/chao2014_rarefaction_extrapolation_hill_numbers.pdf`) supports
  Hill-number validation, coverage-aware interpretation, and reference
  calculations used around the iNEXT spider data.

## Dissimilarity Indices

### General Ecological Distances And Software Conventions

- `legendre2012numerical`
  (`references/manual_downloads.md#legendre2012numerical`) is the broad
  ecological distances reference. It supports matrix orientation,
  transformation language, and the ecological interpretation of
  Jaccard/Sorensen, Bray-Curtis, Ruzicka, Canberra, Hellinger, Chord, and
  related measures.

- `oksanen2022vegan`, `veganDiversityDocs`, and `veganVegdistDocs`
  (`references/vegan_diversity.html`, `references/vegan_vegdist.html`) support
  cross-package conventions with R `vegan`, especially Simpson naming,
  `diversity`, and `vegdist` distance choices.

- `the_scipy_community2020scipy`, `scipyPdistDocs`, and
  `scipyJensenShannonDocs`
  (`references/scipy_pdist.html`, `references/scipy_jensenshannon.html`) support
  Python/SciPy distance conventions, including unaveraged Canberra and
  square-root Jensen-Shannon distance.

### Incidence Similarities

- `jaccard1901distribution`
  (`references/manual_downloads.md#jaccard1901distribution`) supports
  `Jaccard`, `jaccard_similarity`, `jaccard_index`, and `jaccard_distance`.

- `sorensen1948method`
  (`references/manual_downloads.md#sorensen1948method`) supports
  `SorensenDice`, `sorensen_index`, `sorensen_dice_index`, and the associated
  distance/dissimilarity helpers.

- `simpson1960notes`
  (`references/manual_downloads.md#simpson1960notes`) supports `Overlap`,
  `overlap_similarity`, and `overlap_distance`.

### Abundance Similarities And Dissimilarities

- `bray1957ordination`
  (`references/bray_curtis1957_ordination_upland_forest.pdf`) is the primary
  source for Bray-Curtis dissimilarity and supports `BrayCurtis`,
  `bray_curtis_distance`, and `bray_curtis_dissimilarity`.

- `ruzicka1958application`
  (`references/manual_downloads.md#ruzicka1958application`) supports
  `Ruzicka`, `quantitative_jaccard_similarity`,
  `quantitative_jaccard_distance`, and `ruzicka_distance`.

- `lance1967general`
  (`references/manual_downloads.md#lance1967general`) and the SciPy distance
  documentation support `Canberra` and `canberra_distance`.

- `morisita1959measuring` and `horn1966measurement`
  (`references/manual_downloads.md#morisita1959measuring`,
  `references/manual_downloads.md#horn1966measurement`) support
  `MorisitaHorn`, `morisita_horn_similarity`, and `morisita_horn_distance`.

### Probability Metrics And Information Divergences

- `gibbs2002choosing`
  (`references/gibbs_su2002_probability_metrics.pdf`) supports total variation,
  Hellinger, KL, and relationships among probability metrics. It backs
  `TotalVariation`, `Hellinger`, and the package's metric-trait documentation.

- `bhattacharyya1943measure`
  (`references/manual_downloads.md#bhattacharyya1943measure`) is the primary
  Bhattacharyya source. `derpanis2008bhattacharyya`
  (`references/derpanis2008_bhattacharyya_measure.pdf`) is included as a short
  local technical note for the coefficient/distance relationship.

- `kullback1951information`
  (`references/kullback_leibler1951_information_sufficiency_download_failed.html`;
  PDF queued in `manual_downloads.md#kullback1951information`) is the primary
  source for `KullbackLeibler` and `kullback_leibler_divergence`.

- `lin1991divergence`
  (`references/manual_downloads.md#lin1991divergence`) introduces
  Shannon-entropy-based divergences, including Jensen-Shannon divergence. It
  supports `JensenDifference`, `JensenShannon`, `jensen_difference`, and
  `jensen_shannon_divergence`.

- `endres2003new`
  (`references/endres_schindelin2003_new_metric_probability_distributions.pdf`)
  proves the metric form of the square-root Jensen-Shannon distance. It supports
  `jensen_shannon_distance` and `is_metric(JensenShannon())`.

## Identities, Reference Datasets, And Validation Sources

- `whittaker1960vegetation`
  (`references/manual_downloads.md#whittaker1960vegetation`) supports the
  alpha/beta/gamma framing used in the architecture notes.

- `chao2014rarefaction`, `sackett2011disturbance`, `hsieh2016inext`, and
  `inextSpiderDocs`
  (`references/chao2014_rarefaction_extrapolation_hill_numbers.pdf`,
  `references/spider_dataset.md`, `references/inext_spider.html`) support the
  vendored spider reference values and the Hill-number/coverage validation
  notes.

- `scikitbio2024` and `scikitbioDiversityDocs`
  (`references/scikitbio_diversity.html`) support the scikit-bio validation
  candidate matrix and printed richness/Bray-Curtis examples.

- `veganDuneData`, `veganVarespecData`, `veganBciData`, and `veganMiteData`
  (`references/vegan_dune.html`, `references/vegan_varechem.html`,
  `references/vegan_bci.html`, `references/vegan_mite.html`) support the
  candidate R reference datasets described in `validation/`.

- `biosampleRCalcDiversityDocs`
  (`references/biosampler_calc_diversity_indices.html`) supports an independent
  BCI-based diversity summary candidate.

- `peng2011reproducible`, `wilkinson2016fair`, and `bezanson2017julia`
  (`references/manual_downloads.md#peng2011reproducible`,
  `references/manual_downloads.md#wilkinson2016fair`,
  `references/manual_downloads.md#bezanson2017julia`) support the notes
  manuscript framing around reproducible, reusable Julia-based computation.

## Method Coverage Checklist

| Package method or type | Bibliography keys |
|---|---|
| `Plugin`, `Shannon`, `shannon_entropy`, `shannon` | `shannon1948mathematical`, `paninski2003entropy` |
| `MillerMadow` | `miller1955note`, `paninski2003entropy` |
| `Basharin` | `basharin1959statistical` |
| `HausserStrimmer` | `hausser2009entropy` |
| `AddGamma`, Laplace, Jeffreys | `laplace1812theorie`, `jeffreys1946invariant` |
| `ChaoShen` | `chao2003nonparametric`, `good1953population` |
| `Richness` | `magurran2004measuring`, `chao1984nonparametric` |
| `Renyi` | `renyi1961measures` |
| `Tsallis` | `tsallis1988possible` |
| `Simpson`, `GiniSimpson`, `InverseSimpson` | `simpson1949measurement` |
| `GreenbergDiversityIndex`, `LinguisticDiversityIndex` | `greenberg1956measurement` |
| `Hill`, `effective_diversity`, `hill_number` | `hill1973diversity`, `jost2006entropy` |
| `Chao1` | `chao1984nonparametric` |
| `ACE` | `chao1992coverage` |
| `SampleCoverage` | `good1953population` |
| `PielouEvenness` | `pielou1966measurement` |
| `FisherAlpha` | `fisher1943relation` |
| `Jaccard` | `jaccard1901distribution`, `legendre2012numerical` |
| `SorensenDice` | `sorensen1948method`, `legendre2012numerical` |
| `Overlap` | `simpson1960notes`, `legendre2012numerical` |
| `BrayCurtis` | `bray1957ordination` |
| `Ruzicka` | `ruzicka1958application`, `legendre2012numerical` |
| `TotalVariation`, `Manhattan`, `Euclidean` | `gibbs2002choosing`, `scipyPdistDocs` |
| `Canberra` | `lance1967general`, `scipyPdistDocs`, `veganVegdistDocs` |
| `Hellinger`, `Chord` | `legendre2012numerical`, `gibbs2002choosing` |
| `Bhattacharyya` | `bhattacharyya1943measure`, `derpanis2008bhattacharyya` |
| `KullbackLeibler` | `kullback1951information`, `gibbs2002choosing` |
| `ShannonDifference` | `shannon1948mathematical`, `jost2006entropy` |
| `JensenDifference`, `JensenShannon` | `lin1991divergence`, `endres2003new` |
| `MorisitaHorn` | `morisita1959measuring`, `horn1966measurement` |
| Trait/bounds/metric descriptors | `legendre2012numerical`, `gibbs2002choosing`, `endres2003new` |
| Reference datasets | `chao2014rarefaction`, `sackett2011disturbance`, `hsieh2016inext`, `inextSpiderDocs`, `veganDuneData`, `scikitbioDiversityDocs` |
