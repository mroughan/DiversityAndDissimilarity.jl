# Data Input Formats

All functions in `DiversityAndDissimilarity.jl` accept several data shapes and
interpret them consistently. The word "species" follows ecology convention but
is a generic placeholder for any discrete category: taxa, words, languages,
demographic groups, or any other count-based label set.

## Overview

| Format | Typical use | Output shape |
|---|---|---|
| `Dict{K,V}` | Named categories with abundance counts | Scalar |
| Numeric vector | Unlabelled abundance counts | Scalar |
| Non-numeric vector | Raw observations (each element is one event) | Scalar |
| Numeric vector with `frequencies=false` | Numeric raw observations | Scalar |
| Matrix | Multiple samples — rows are samples, columns are taxa | Vector (one value per row) or pairwise matrix |
| Tables.jl-compatible table (DataFrame, etc.) | Labelled samples with named taxon columns | Same as matrix |

For community matrices and tables, pairwise index functions return a square
matrix of comparisons between every pair of rows; alpha-diversity functions
return one value per row.

---

## Dictionaries

A `Dict` maps each category label to its abundance count or relative frequency.
The key type can be anything (symbols, strings, integers).

```jldoctest datainput
julia> using DiversityAndDissimilarity

julia> assemblage = Dict(:oak => 12, :ash => 5, :elm => 3);

julia> richness(assemblage)
3

julia> shannon_entropy(assemblage)
1.3527241956246545

julia> dissimilarity(BrayCurtis(), Dict(:oak => 10, :ash => 2), Dict(:ash => 5, :elm => 3))
0.7
```

---

## Numeric Vectors

A numeric vector is treated as an ordered list of abundance values, one per
category. The identity of each category is not tracked — only the distribution
of counts matters.

```jldoctest datainput
julia> abundance = [12, 5, 3];

julia> richness(abundance)
3

julia> shannon_entropy(abundance)
1.3527241956246545
```

For pairwise comparisons, two numeric vectors must have the same length.
Position is used to align categories: element `i` in the left vector
corresponds to element `i` in the right vector.

```jldoctest datainput
julia> left = [12, 5, 0];

julia> right = [0, 5, 7];

julia> dissimilarity(BrayCurtis(), left, right)
0.6551724137931034
```

---

## Non-Numeric Observation Vectors

A non-numeric vector (e.g. `Vector{String}`) is treated as a sequence of raw
observations. Each element is one event; [`counts`](@ref) tallies how often
each value appears.

```jldoctest datainput
julia> observations = ["oak", "ash", "oak", "elm", "ash", "oak"];

julia> counts(observations)
Dict{String, Int64}("elm" => 1, "ash" => 2, "oak" => 3)

julia> richness(observations)
3

julia> shannon_entropy(observations)
1.4591479170272448
```

---

## Numeric Observation Vectors

When a numeric vector represents individual observations rather than abundance
counts — for example `[1, 2, 1, 3]` meaning four events whose labels happen to
be integers — pass `frequencies=false` to suppress the default abundance
interpretation.

```jldoctest datainput
julia> obs = [1, 2, 1, 3];

julia> richness(obs; frequencies=false)
3

julia> richness(obs)
3
```

Note that richness gives the same answer here because the number of distinct
values equals the number of nonzero positions. Shannon entropy differs:

```jldoctest datainput
julia> round(shannon_entropy(obs; frequencies=false); digits=4)
1.5

julia> round(shannon_entropy(obs); digits=4)
1.5
```

A clearer example where the distinction matters:

```jldoctest datainput
julia> raw = [1, 1, 2];          # three observations: two "1"s and one "2"

julia> richness(raw; frequencies=false)   # 2 distinct categories
2

julia> richness(raw)                       # 3 nonzero positions
3
```

---

## Community Matrices

A numeric matrix is interpreted as a **community matrix**: rows are samples
(sites, plots, time points, speakers, …) and columns are taxa or categories.
Every alpha-diversity function returns one value per row.

```jldoctest datainput
julia> community = [
           1 1 2 0 5
           3 0 1 1 0
       ];

julia> richness(community)
2-element Vector{Int64}:
 4
 3

julia> shannon_entropy(community)
2-element Vector{Float64}:
 1.6577427265048888
 1.3709505944546687

julia> alpha_diversity(community)
2-element Vector{@NamedTuple{…}}:
 (richness = 4, …)
 (richness = 3, …)
```

Pairwise index functions called with a single matrix argument return a square
pairwise matrix:

```jldoctest datainput
julia> bray_curtis_distance(community)
2×2 Matrix{Float64}:
 0.0       0.714286
 0.714286  0.0
```

### Orientation

The row = sample, column = taxon orientation is the standard in community
ecology (the "sites × species" matrix). If your data is the other way around,
transpose it before passing:

```julia
species_by_site_matrix = rand(1:10, 5, 20)   # 5 taxa, 20 sites
richness(species_by_site_matrix')             # transpose → 20 sites × 5 taxa
```

### Validation rules

Every community matrix passed to a function is checked for:
- no negative values;
- no non-finite values (`Inf`, `NaN`);
- positive total abundance in every row (no all-zero rows).

Errors include the row and column index so problems are easy to locate.

---

## Tables and DataFrames

Any [Tables.jl](https://github.com/JuliaData/Tables.jl)-compatible object —
including `DataFrames.DataFrame`, named tuples, and CSV rows — can be passed
directly. Numeric columns are treated as taxon columns by default.

```julia
using DataFrames

df = DataFrame(
    site   = ["plot-a", "plot-b"],
    oak    = [1, 3],
    ash    = [1, 0],
    elm    = [2, 1],
)

richness(df; species=[:oak, :ash, :elm])
shannon_entropy(df; species=[:oak, :ash, :elm])
alpha_diversity(df; species=[:oak, :ash, :elm])
```

Use [`community_matrix`](@ref) to inspect the numeric matrix that will be
extracted before analysis:

```julia
community_matrix(df; species=[:oak, :ash, :elm])
```

Use `label` (a column name) or `labels` (a vector) with pairwise helpers to
attach sample identifiers to the result:

```julia
labeled_distance(BrayCurtis(), df; label=:site, species=[:oak, :ash, :elm])
```

### Selecting taxon columns

By default every numeric column is treated as a taxon. If the table contains
numeric metadata (e.g. a plot ID column), pass `species` explicitly to select
only the taxon columns:

```julia
richness(df; species=[:oak, :ash, :elm])
```

`species` accepts a vector of column names as `Symbol` or `String`.

---

## Pre-validated Pipelines

When the same community matrix is passed to many functions in a pipeline,
each call repeats the full validation pass. Use [`validate`](@ref) to run
validation once and obtain a [`Validated`](@ref) wrapper:

```julia
v = validate(community)   # validates once; returns Validated{Matrix{Float64}}

richness(v)
shannon_entropy(v)
alpha_diversity(v)
bray_curtis_distance(v)
```

All major functions accept `Validated` arguments and skip per-call validation.
This matters most for inexpensive row-wise operations such as richness, where
validation cost is a significant fraction of total time.

!!! warning
    Constructing `Validated(data)` directly — without calling `validate` —
    bypasses all checks. Use `validate(data)` to obtain a wrapper from
    unverified input.

---

## Keyword Parameters

### `frequencies`

Controls how a **numeric vector** is interpreted.

| Value | Meaning |
|---|---|
| `true` (default) | Each element is an abundance count |
| `false` | Each element is one raw observation; the vector is first tallied with [`counts`](@ref) |

`frequencies` has no effect on dictionaries, non-numeric vectors, or matrices.

### `species`

Selects taxon columns when the input is a Tables.jl-compatible table. Accepts
a vector of column names (`Symbol` or `String`). When omitted, all numeric
columns are treated as taxon columns.

### `label` / `labels`

Used with pairwise matrix helpers ([`labeled_distance`](@ref),
[`labeled_dissimilarity`](@ref), [`labeled_similarity`](@ref)).

- `label`: a column name in a table; the values in that column become sample
  labels in the result.
- `labels`: a vector of label strings provided directly, one per row.

Pass either `label` or `labels`, not both.
