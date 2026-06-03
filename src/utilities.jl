const CountInput = Union{AbstractVector,AbstractDict}

_rows(data::AbstractMatrix) = (view(data, row, :) for row in axes(data, 1))
_check_matrix_frequencies(frequencies::Bool) =
    frequencies || throw(ArgumentError("community matrices are abundance matrices; use rows of observations with frequencies=false"))
_is_table(data) = Tables.istable(data) && !(data isa AbstractVector) && !(data isa AbstractDict) && !(data isa AbstractMatrix)

"""
    community_matrix(data; species=nothing)

Return a numeric community matrix with samples/sites in rows and species/taxa in
columns.

Matrices are returned as floating-point copies. Tables.jl-compatible inputs,
including DataFrames, are converted column-wise. By default, numeric columns are
used as species columns. Pass `species` as a collection of column names to choose
species columns explicitly, which is recommended when a table also contains
numeric site identifiers or metadata.
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
    for value in data
        value < 0 && throw(ArgumentError("abundances must be non-negative"))
        isfinite(value) || throw(ArgumentError("abundances must be finite"))
    end
    for row in axes(data, 1)
        sum(view(data, row, :)) > 0 || throw(ArgumentError("each community matrix row must have positive total abundance"))
    end
    return nothing
end

"""
    counts(x)

Return a `Dict` of category counts.

Vectors are treated as observations. Dictionaries are treated as category =>
abundance mappings and are validated without changing key types.
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

Return positive relative abundances for `x`.

Numeric vectors are interpreted as abundance/frequency vectors by default.
Use `frequencies=false` to treat a numeric vector as raw observations.
Community matrices are interpreted as samples in rows and taxa/categories in
columns; the returned matrix contains row-wise proportions.
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
