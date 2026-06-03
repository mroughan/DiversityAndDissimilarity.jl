# Forest Gradient Workflow Case Study

This simulated case study uses `validation/data/community_counts.csv` as a
transparent stand-in for a community monitoring workflow. The point is to
exercise convention-sensitive biodiversity calculations, not to infer from
real ecological observations.

## Workflow

1. Read a sample-by-taxon matrix with site labels and habitat metadata.
2. Run `diversity_audit` for alpha summaries, estimator diagnostics, and a labeled Bray-Curtis matrix.
3. Run `uncertainty_audit` to bootstrap Shannon entropy and effective diversity for each site.
4. Inspect low-coverage warnings before using the pairwise matrix in downstream analyses.

## Dataset

- Samples: 8
- Taxa: 12
- Habitats: forest, wetland, grassland
- Total abundance range: 26.0 to 47.0

## Alpha And Uncertainty Summary

| Site | Habitat | Gradient | Richness | Coverage | Shannon H bootstrap | Effective diversity bootstrap |
|---|---|---|---:|---:|---:|---:|
| plot_A | forest | early | 7 | 0.979 | 2.31 [1.85, 2.509] | 4.96 [3.704, 5.6] |
| plot_B | forest | early | 8 | 0.912 | 2.443 [1.795, 2.657] | 5.438 [3.659, 6.147] |
| plot_C | forest | mid | 8 | 0.923 | 2.686 [2.054, 2.81] | 6.437 [4.183, 7.021] |
| plot_D | forest | mid | 7 | 0.923 | 2.439 [1.834, 2.599] | 5.423 [3.557, 6.201] |
| plot_E | wetland | late | 8 | 0.939 | 2.575 [2.017, 2.747] | 5.959 [3.952, 6.537] |
| plot_F | wetland | late | 8 | 0.951 | 2.504 [2.007, 2.712] | 5.671 [3.925, 6.469] |
| plot_G | grassland | stress | 8 | 0.909 | 2.457 [1.902, 2.634] | 5.491 [3.778, 6.454] |
| plot_H | grassland | stress | 7 | 0.941 | 2.352 [1.847, 2.557] | 5.107 [3.622, 5.714] |

## Audit Warnings

No audit warnings were produced.

## Interpretation

The workflow keeps alpha summaries, coverage diagnostics, bootstrap
intervals, and the Bray-Curtis distance matrix tied to the same labels and
input matrix. This is the core reproducibility benefit: the report exposes
which samples are convention-sensitive or coverage-sensitive before the
distance matrix is reused for ordination, clustering, or modelling.
