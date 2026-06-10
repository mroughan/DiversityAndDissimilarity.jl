const CountInput = Union{AbstractVector,AbstractDict}

_rows(data::AbstractMatrix) = (view(data, row, :) for row in axes(data, 1))
_check_matrix_frequencies(frequencies::Bool) =
    frequencies || throw(ArgumentError("community matrices are abundance matrices; use rows of observations with frequencies=false"))
_is_table(data) = Tables.istable(data) && !(data isa AbstractVector) && !(data isa AbstractDict) && !(data isa AbstractMatrix)

"""
    community_matrix(data; species=nothing)

Return a numeric community matrix with samples/sites in rows and species/taxa in
columns.

Matrices are returned as `Float64` copies. All abundances must be non-negative and
finite, and each row must have positive total abundance; otherwise an
`ArgumentError` is thrown.

Tables.jl-compatible inputs, including DataFrames, are converted column-wise. By
default, numeric columns are used as species columns. Pass `species` as a collection
of column names to choose species columns explicitly, which is recommended when a
table also contains numeric site identifiers or metadata.

```jldoctest
julia> using DiversityAndDissimilarity

julia> community_matrix([1 2 3; 4 5 6])
2×3 Matrix{Float64}:
 1.0  2.0  3.0
 4.0  5.0  6.0
```
"""
function community_matrix(data::AbstractMatrix{<:Real}; species=nothing)
    species === nothing || throw(ArgumentError("species selection is only used for Tables.jl-compatible inputs"))
    _validate_community_matrix(data)
    return float.(data)
end

function community_matrix(data; species=nothing)
    _is_table(data) || throw(ArgumentError("expected a numeric matrix or a Tables.jl-compatible table"))
    columns = Tables.columns(data)
    names = collect(Tables.columnnames(columns))
    selected = if species === nothing
        [name for name in names if eltype(Tables.getcolumn(columns, name)) <: Real]
    else
        collect(species)
    end
    isempty(selected) && throw(ArgumentError("at least one numeric species column is required"))

    nrows = _table_nrows(columns, names)
    result = Matrix{Float64}(undef, nrows, length(selected))
    for (column_index, name) in pairs(selected)
        name in names || throw(ArgumentError("species column $(repr(name)) was not found"))
        column = Tables.getcolumn(columns, name)
        eltype(column) <: Real || throw(ArgumentError("species column $(repr(name)) must be numeric"))
        length(column) == nrows || throw(ArgumentError("species columns must have the same length"))
        for row in 1:nrows
            result[row, column_index] = float(column[row])
        end
    end
    _validate_community_matrix(result)
    return result
end

function _table_nrows(columns, names)
    isempty(names) && throw(ArgumentError("table must have at least one column"))
    return length(Tables.getcolumn(columns, first(names)))
end

function _validate_community_matrix(data::AbstractMatrix{<:Real})
    for col in axes(data, 2)
        for row in axes(data, 1)
            value = data[row, col]
            value < 0 && throw(ArgumentError(
                "abundance at row $row, column $col is negative ($value); all abundances must be non-negative"))
            isfinite(value) || throw(ArgumentError(
                "abundance at row $row, column $col is not finite ($value); all abundances must be finite"))
        end
    end
    for row in axes(data, 1)
        sum(view(data, row, :)) > 0 || throw(ArgumentError(
            "row $row has zero total abundance; every community must contain at least one observation"))
    end
    return nothing
end

"""
    counts(x)

Return a `Dict` mapping each category to its count or validated abundance.

When `x` is a vector, each element is treated as an observation and the result
maps each unique value to the number of times it appears, with type
`Dict{eltype(x), Int64}`. When `x` is a `Dict`, it is treated as a
`category => abundance` mapping: abundances are validated to be non-negative and
finite, and returned as `Float64` values without changing the key type.

```jldoctest
julia> using DiversityAndDissimilarity

julia> sort(collect(counts(["oak", "ash", "oak"])), by=first)
2-element Vector{Pair{String, Int64}}:
 "ash" => 1
 "oak" => 2

julia> sort(collect(counts(Dict(:oak => 3, :ash => 1))), by=first)
2-element Vector{Pair{Symbol, Float64}}:
 :ash => 1.0
 :oak => 3.0
```
"""
function counts(xs::AbstractVector)
    result = Dict{eltype(xs), Int}()
    for x in xs
        result[x] = get(result, x, 0) + 1
    end
    return result
end

function counts(abundances::AbstractDict{K,<:Real}) where {K}
    result = Dict{K, Float64}()
    for (species, abundance) in pairs(abundances)
        abundance < 0 && throw(ArgumentError("abundances must be non-negative"))
        isfinite(abundance) || throw(ArgumentError("abundances must be finite"))
        result[species] = float(abundance)
    end
    return result
end

function _abundances(data::AbstractDict; frequencies::Bool=true)
    abundance = collect(float(v) for v in values(counts(data)) if v > 0)
    isempty(abundance) && throw(ArgumentError("at least one positive abundance is required"))
    return abundance
end

function _abundances(data::AbstractVector{<:Real}; frequencies::Bool=true)
    frequencies || return collect(float(v) for v in values(counts(data)))

    for value in data
        value < 0 && throw(ArgumentError("abundances must be non-negative"))
        isfinite(value) || throw(ArgumentError("abundances must be finite"))
    end

    abundance = collect(float(value) for value in data if value > 0)
    isempty(abundance) && throw(ArgumentError("at least one positive abundance is required"))
    return abundance
end

function _abundances(data::AbstractVector; frequencies::Bool=false)
    abundance = collect(float(v) for v in values(counts(data)) if v > 0)
    isempty(abundance) && throw(ArgumentError("at least one observation is required"))
    return abundance
end

function _abundances(data::AbstractMatrix{<:Real}; frequencies::Bool=true)
    _check_matrix_frequencies(frequencies)
    _validate_community_matrix(data)
    abundance = collect(float(value) for value in data if value > 0)
    isempty(abundance) && throw(ArgumentError("at least one positive abundance is required"))
    return abundance
end

"""
    proportions(x; frequencies=true)

Return relative abundances for `x` as a probability vector or matrix.

Numeric vectors are interpreted as abundance/frequency vectors by default; zero
entries are discarded and the returned vector contains only the positive-abundance
proportions. Pass `frequencies=false` to treat a numeric vector as raw observations
instead. Non-numeric vectors are always treated as observations.

Community matrices (samples in rows, taxa in columns) are normalized row-wise: the
returned matrix has the same shape as the input, each row sums to one, and zero
entries are preserved.

```jldoctest
julia> using DiversityAndDissimilarity

julia> proportions([3, 1])
2-element Vector{Float64}:
 0.75
 0.25

julia> proportions([3 1; 2 2])
2×2 Matrix{Float64}:
 0.75  0.25
 0.5   0.5
```
"""
function proportions(data::AbstractMatrix{<:Real}; frequencies::Bool=true)
    _check_matrix_frequencies(frequencies)
    _validate_community_matrix(data)
    result = Matrix{Float64}(undef, size(data))
    for row in axes(data, 1)
        row_total = 0.0
        for column in axes(data, 2)
            value = data[row, column]
            value < 0 && throw(ArgumentError("abundances must be non-negative"))
            isfinite(value) || throw(ArgumentError("abundances must be finite"))
            row_total += value
        end
        row_total > 0 || throw(ArgumentError("each community matrix row must have positive total abundance"))
        for column in axes(data, 2)
            result[row, column] = data[row, column] / row_total
        end
    end
    return result
end

function proportions(data; frequencies::Bool=true, species=nothing)
    if _is_table(data)
        _check_matrix_frequencies(frequencies)
        return proportions(community_matrix(data; species))
    end
    abundance = _abundances(data; frequencies)
    total = sum(abundance)
    total > 0 || throw(ArgumentError("total abundance must be positive"))
    return abundance ./ total
end

_community_input(data; species=nothing) = _is_table(data) ? community_matrix(data; species) : data

function _species_set(data::AbstractDict; frequencies::Bool=true)
    return Set(k for (k, v) in pairs(data) if v > 0)
end

function _species_set(data::AbstractVector{<:Real}; frequencies::Bool=true)
    frequencies || return Set(data)
    return Set(i for (i, v) in pairs(data) if v > 0)
end

function _species_set(data::AbstractVector; frequencies::Bool=false)
    return Set(data)
end

function _aligned_abundances(left::AbstractDict, right::AbstractDict; frequencies::Bool=true)
    left_counts = counts(left)
    right_counts = counts(right)
    species = union(keys(left_counts), keys(right_counts))
    return (
        [get(left_counts, key, 0.0) for key in species],
        [get(right_counts, key, 0.0) for key in species],
    )
end

function _aligned_abundances(left::AbstractVector{<:Real}, right::AbstractVector{<:Real}; frequencies::Bool=true)
    if frequencies
        length(left) == length(right) || throw(DimensionMismatch("abundance vectors must have the same length"))
        return (float.(left), float.(right))
    else
        return _aligned_abundances(counts(left), counts(right))
    end
end

function _aligned_abundances(left::AbstractVector, right::AbstractVector; frequencies::Bool=false)
    return _aligned_abundances(counts(left), counts(right))
end

"""
    Validated(data)

A thin wrapper asserting that `data` has already been validated as a community
matrix. Pass a `Validated` value to diversity and dissimilarity functions to
skip per-call input validation.

Obtain a `Validated` wrapper from raw data using [`validate`](@ref), which runs
all the usual checks and returns a type-stable `Float64` copy. Constructing
`Validated(data)` directly bypasses all checks entirely.

!!! warning "Use `validate`, not `Validated` directly"
    Constructing `Validated(data)` directly bypasses all input checks.
    Passing invalid data — negative abundances, non-finite values, or all-zero
    rows — may produce silently wrong results or unhelpful errors. Always use
    [`validate`](@ref) unless you are certain the data satisfies all
    preconditions and are deliberately trading safety for performance.
"""
struct Validated{T}
    data::T
end

"""
    validate(data; species=nothing)

Validate `data` and return a [`Validated`](@ref) wrapper.

Checks that all abundances are non-negative and finite and that each row has
positive total abundance. Matrices are returned as `Float64` copies. The
returned `Validated` object can be passed to diversity and dissimilarity
functions to skip per-call validation overhead, which is useful when the same
matrix is consumed by many functions in a pipeline.

The safe default pathway validates on every call:

```julia
richness(community)        # validates, then computes
shannon_entropy(community) # validates independently
```

With `validate`, validation runs once and the result is reused:

```julia
v = validate(community)    # validates once
richness(v)                # computation only
shannon_entropy(v)         # computation only
bray_curtis_distance(v)    # computation only
```

```jldoctest
julia> using DiversityAndDissimilarity

julia> v = validate([1 1 0; 0 1 1]);

julia> richness(v)
2-element Vector{Int64}:
 2
 2

julia> dissimilarity(BrayCurtis(), v)
2×2 Matrix{Float64}:
 0.0  0.5
 0.5  0.0
```
"""
function validate(data::AbstractMatrix{<:Real}; species=nothing)
    species === nothing || throw(ArgumentError("species selection is only used for Tables.jl-compatible inputs"))
    _validate_community_matrix(data)
    return Validated(float.(data))
end

function validate(data; species=nothing)
    _is_table(data) || throw(ArgumentError("expected a numeric matrix or a Tables.jl-compatible table"))
    return Validated(community_matrix(data; species))
end
