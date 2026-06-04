# Reference Library

This directory stores local copies of the reference materials used by the
notes, documentation, validation fixtures, and method annotations for
`DiversityAndDissimilarity.jl`.

The curated guide is `../annotated_bibliography.md`. The BibTeX database is
`../references.bib`; each entry has a `file` field pointing to a local artifact
in this directory or to `manual_downloads.md` when a paper or book could not be
downloaded directly.

## Directory Conventions

- PDF files are downloaded papers or technical notes that were reachable from a
  stable public URL.
- HTML files are local snapshots of documentation pages, publisher metadata
  pages, or failed-download challenge pages.
- Markdown files contain hand-curated validation notes or download queues.
- `manual_downloads.md` lists library items that still need a PDF or print-copy
  note added manually.

When adding a new reference, use a descriptive lowercase filename:

```text
surname_year_short_topic.pdf
surname_year_short_topic.html
```

Then add or update the matching BibTeX entry in `../references.bib` with a
`file = {references/...}` field.

## Downloaded Diversity References

| File | BibTeX key | Use |
|---|---|---|
| `shannon1948_mathematical_theory_communication.pdf` | `shannon1948mathematical` | Shannon entropy, plugin estimator, information scale |
| `renyi1961_measures_entropy_information.pdf` | `renyi1961measures` | Renyi entropy and diversity |
| `greenberg1956_measurement_linguistic_diversity.pdf` | `greenberg1956measurement` | Greenberg / linguistic diversity index |
| `fisher_corbet_williams1943_logseries_species_abundance.pdf` | `fisher1943relation` | Fisher alpha and log-series diversity |
| `good1953_population_frequencies_species.pdf` | `good1953population` | Good-Turing sample coverage and unseen mass |
| `paninski2003_entropy_mutual_information.pdf` | `paninski2003entropy` | Entropy estimation bias |
| `hausser_strimmer2009_entropy_james_stein.pdf` | `hausser2009entropy` | James-Stein shrinkage entropy estimator |
| `nemenman2002_entropy_inference_revisited.pdf` | `nemenman2002entropy` | Bayesian/NSB entropy estimation |
| `jiao2015_minimax_information_functionals.pdf` | `jiao2015minimax` | Minimax information-functional estimation |
| `wu_yang2016_minimax_entropy_large_alphabets.pdf` | `wu2016minimax` | Large-alphabet minimax entropy estimation |
| `valiant_valiant2017_estimating_unseen.pdf` | `valiant2017estimating` | Estimating unseen species and entropy |
| `chao2014_rarefaction_extrapolation_hill_numbers.pdf` | `chao2014rarefaction` | Hill-number rarefaction/extrapolation and spider validation |

## Downloaded Dissimilarity References

| File | BibTeX key | Use |
|---|---|---|
| `bray_curtis1957_ordination_upland_forest.pdf` | `bray1957ordination` | Bray-Curtis dissimilarity |
| `gibbs_su2002_probability_metrics.pdf` | `gibbs2002choosing` | Total variation, Hellinger, KL, probability-metric bounds |
| `endres_schindelin2003_new_metric_probability_distributions.pdf` | `endres2003new` | Jensen-Shannon distance as a metric |
| `derpanis2008_bhattacharyya_measure.pdf` | `derpanis2008bhattacharyya` | Bhattacharyya coefficient/distance technical note |

## Local HTML Snapshots

These are local copies of online documentation or publisher pages. Some
publisher pages contain enough metadata for traceability but not the article
PDF; those are also listed in `manual_downloads.md`.

| File | BibTeX key | Use |
|---|---|---|
| `vegan_diversity.html` | `veganDiversityDocs` | R `vegan::diversity` conventions |
| `vegan_vegdist.html` | `veganVegdistDocs` | R `vegan::vegdist` conventions |
| `scipy_pdist.html` | `scipyPdistDocs` | SciPy pairwise distance conventions |
| `scipy_jensenshannon.html` | `scipyJensenShannonDocs` | SciPy Jensen-Shannon convention |
| `scikitbio_diversity.html` | `scikitbioDiversityDocs` | scikit-bio diversity examples |
| `inext_spider.html` | `inextSpiderDocs` | iNEXT spider dataset documentation |
| `vegan_dune.html` | `veganDuneData` | Candidate `vegan::dune` reference dataset |
| `vegan_varechem.html` | `veganVarespecData` | Candidate `vegan::varespec` reference dataset |
| `vegan_bci.html` | `veganBciData` | Candidate `vegan::BCI` reference dataset |
| `vegan_mite.html` | `veganMiteData` | Candidate `vegan::mite` reference dataset |
| `biosampler_calc_diversity_indices.html` | `biosampleRCalcDiversityDocs` | Candidate BCI diversity summaries |
| `chao_shen2003_unseen_species_shannon_springer.html` | `chao2003nonparametric` | Springer metadata/fulltext page; PDF still manual |
| `tsallis1988_generalized_boltzmann_gibbs_statistics_springer.html` | `tsallis1988possible` | Springer metadata/fulltext page; PDF still manual |
| `simpson1949_measurement_diversity_nature.html` | `simpson1949measurement` | Nature metadata/fulltext page; PDF still manual |

## Validation Notes

| File | Use |
|---|---|
| `spider_dataset.md` | Exact iNEXT spider abundance vectors and pinned diversity values |
| `scipy_conventions.md` | Convention comparison with `scipy.spatial.distance` |

## Failed Download Records

These files are intentionally retained as provenance for blocked direct
downloads. The corresponding citations are queued in `manual_downloads.md`.

| File | Related key |
|---|---|
| `basharin1959_entropy_estimation_download_failed.html` | `basharin1959statistical` |
| `kullback_leibler1951_information_sufficiency_download_failed.html` | `kullback1951information` |
| `good1953_good_turing_download_failed.html` | superseded by `good1953_population_frequencies_species.pdf` |

## Manual Library Downloads

See `manual_downloads.md` for the full queue. High-priority items are:

- `hill1973diversity`, `jost2006entropy`: Hill-number and effective-diversity
  identities.
- `chao1984nonparametric`, `chao1992coverage`, `chao2003nonparametric`:
  Chao1, ACE, and Chao-Shen estimators.
- `miller1955note`, `basharin1959statistical`, `laplace1812theorie`,
  `jeffreys1946invariant`: low-sample entropy and pseudocount corrections.
- `legendre2012numerical`: broad ecological distance conventions and metric
  properties.
- `jaccard1901distribution`, `sorensen1948method`, `ruzicka1958application`,
  `morisita1959measuring`, `horn1966measurement`: primary sources for
  incidence and abundance similarities.
- `kullback1951information`, `lin1991divergence`, `bhattacharyya1943measure`:
  primary information-divergence sources.
