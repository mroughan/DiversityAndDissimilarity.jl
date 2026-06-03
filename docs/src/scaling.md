# Scaling And Performance

Most functions in `DiversityAndDissimilarity` are single-pass or small-number-of-pass
summaries over species abundance vectors. The main scaling variables are:

- ``n``: sample size, or total count in one assemblage.
- ``S``: observed support, the number of species/taxa/categories with positive
  abundance.
- ``K``: known finite support, including unobserved categories when supplied
  through `support`.
- ``M``: number of samples/sites, usually the number of rows in a community
  matrix.
- ``P``: number of species/taxa columns in a community matrix.
- ``B``: number of bootstrap replicates.

## Single Assemblages

For a dictionary or numeric abundance vector, alpha-diversity calculations are
linear in the number of categories that must be inspected. For a raw observation
vector, the package must first count observations, so cost is linear in sample
size and memory is proportional to observed support.

| Input or operation | Time | Extra memory | Scaling notes |
|---|---:|---:|---|
| Numeric abundance vector | ``O(P)`` | ``O(S)`` | Zeros are inspected but positive abundances drive most formulas. |
| Dictionary of abundances | ``O(S)`` | ``O(S)`` | Keys define observed support. |
| Raw observation vector | ``O(n)`` | ``O(S)`` | Observations are counted before proportions are computed. |
| `proportions` | ``O(S)`` or ``O(P)`` | ``O(S)`` | Matrix input returns a dense matrix of row proportions. |
| Richness, Shannon, Renyi, Tsallis, Simpson | ``O(S)`` | ``O(S)`` | After abundance extraction, formulas are linear in observed support. |
| `Chao1`, `ACE`, `SampleCoverage` | ``O(S)`` | ``O(S)`` | These count abundance frequencies such as singletons and doubletons. |

For large sparse supports, prefer dictionaries or table columns that contain
only the species you need. Dense community matrices are convenient, but they
scale with all ``M \\times P`` entries, including zeros.

## Known Support

Some Shannon estimators use `support` to represent a known finite category
space. Passing `support=K` has little memory cost; passing support labels also
checks that all observed categories are included.

| Estimator | Support requirement | Time in one assemblage | Notes |
|---|---|---:|---|
| [`Plugin`](@ref) | observed support | ``O(S)`` | Direct empirical entropy. |
| [`MillerMadow`](@ref) | observed support | ``O(S)`` | Adds an observed-support correction. |
| [`HausserStrimmer`](@ref) | known finite support | ``O(S + K)`` | Shrinks toward a uniform finite support. |
| [`Basharin`](@ref) | known finite support | ``O(S)`` | Correction depends on support size. |
| [`AddGamma`](@ref) | known finite support | ``O(S + K)`` | Adds pseudocount mass to observed and unobserved categories. |
| [`ChaoShen`](@ref) | unknown support allowed | ``O(S)`` | Uses sample coverage to account for unseen species. |

Use `support` deliberately for large ``K``. If ``K`` is much larger than the
observed support, finite-support estimators can move substantially and may do
more work than observed-support estimators.

## Matrix Fast Path

When computing `alpha_diversity` on a community matrix with the default
[`Plugin`](@ref) estimator and no `support` override, this package uses a
specialized cache-friendly kernel that accumulates all statistics in two
column-major passes. This is significantly faster than the fallback row-by-row
path, which is used for all other estimators.

```julia
alpha_diversity(community)                          # uses fast kernel (Plugin, no support)
alpha_diversity(community; estimator=MillerMadow()) # uses row-by-row fallback
alpha_diversity(community; support=10)              # uses row-by-row fallback
```

The fast kernel covers [`Richness`](@ref), [`Shannon`](@ref),
[`Simpson`](@ref)/[`GiniSimpson`](@ref)/[`InverseSimpson`](@ref),
[`Chao1`](@ref), [`ACE`](@ref), [`SampleCoverage`](@ref),
[`PielouEvenness`](@ref), and [`FisherAlpha`](@ref) simultaneously.
For non-Plugin or support-constrained workflows, benchmark with a
representative subset to check whether the row-by-row path is fast enough,
or pre-compute proportions and pass them directly.

## Community Matrices And Tables

Community matrices are interpreted as samples in rows and species/taxa in
columns. Alpha-diversity calculations are row-wise:

```julia
shannon_entropy(community)
chao1(community)
ace(community)
```

For a dense ``M \\times P`` community matrix, row-wise alpha diversity is
``O(MP)`` time. Result memory is usually ``O(M)``, except for helper operations
such as `proportions(community)`, which return a full ``M \\times P`` matrix.

Tables.jl-compatible inputs are converted to a dense community matrix. Pass
`species` explicitly when the table has metadata columns, and avoid selecting
columns that are not species/taxa abundances:

```julia
community_matrix(table; species=[:oak, :ash, :elm])
shannon_entropy(table; species=[:oak, :ash, :elm])
```

## Pairwise Similarity And Dissimilarity

Pairwise comparisons between two assemblages are linear in aligned support:
``O(P)`` for vectors and approximately ``O(S_1 + S_2)`` for dictionaries.
Low-sample corrections for KL and Jensen-Shannon-style divergences remain
linear in the aligned or supplied support. `AddGamma` and `HausserStrimmer`
scale with the supplied finite support ``K`` when `support` is larger than the
observed aligned support; `ChaoShen` adds one residual unseen-mass coordinate.

Passing a community matrix as a single argument computes all pairwise
comparisons across rows:

```julia
distance(BrayCurtis(), community)
distance(Hellinger(), community)
jensen_shannon_distance(community)
```

This returns a dense ``M \\times M`` matrix. Runtime is ``O(M^2P)`` for dense
matrices and memory is ``O(M^2)`` for the result. This is the scaling limit to
watch most carefully: doubling the number of sites roughly quadruples both
pairwise runtime and output size.

For large ``M``:

- compute only the pairs you need when a full matrix is unnecessary;
- filter or aggregate species columns before pairwise analysis when justified;
- consider storing downstream results in a sparse or long table format if many
  distances are discarded after thresholding;
- prefer bounded measures such as [`BrayCurtis`](@ref),
  [`TotalVariation`](@ref), or [`Hellinger`](@ref) when comparing many datasets
  because their ranges are easy to monitor for quality checks.
- use [`is_bounded`](@ref), [`is_finite`](@ref), [`is_metric`](@ref), and
  [`index_bounds`](@ref) to choose indices and attach bound-aware quality
  checks to large automated workflows.

## Resampling

Bootstrap and jackknife helpers repeat an estimator many times.

| Helper | Time | Memory | Notes |
|---|---:|---:|---|
| [`bootstrap`](@ref) | ``O(Bn)`` plus estimator cost | ``O(B)`` for replicates | Increase `nboot` for stable intervals; start small while developing. |
| [`jackknife`](@ref) | ``O(n)`` leave-one-out fits | ``O(n)`` for leave-one-out values | Best for moderate integer-count samples. |

On community matrices, resampling is applied row-wise, so the total cost scales
with the number of rows as well as the per-row sample size.

## Practical Defaults

For small to medium community datasets, direct calls such as
`shannon_entropy(community)` and `distance(BrayCurtis(), community)` are usually
the clearest choice.

For large datasets:

- use community matrices for dense, rectangular data and dictionaries for
  sparse single assemblages;
- keep raw observation vectors only when the observation sequence itself is
  needed, otherwise count them once;
- pass `species` for table inputs to avoid accidental metadata expansion;
- be cautious with full pairwise matrices when ``M`` is large;
- benchmark with a representative subset before increasing ``B`` for bootstrap
  intervals.
