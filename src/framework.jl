"""
    index_family(index)

Return a symbolic family label for an index, such as `:entropy`,
`:richness`, `:evenness`, `:incidence`, `:abundance`, or `:probability`.
"""
index_family(::DiversityIndex) = :unknown
index_family(::Richness) = :richness
index_family(::Shannon) = :entropy
index_family(::Renyi) = :entropy
index_family(::Tsallis) = :entropy
index_family(::Simpson) = :dominance
index_family(::GiniSimpson) = :diversity
index_family(::Union{GreenbergDiversityIndex,LinguisticDiversityIndex}) = :linguistic_diversity
index_family(::InverseSimpson) = :effective_diversity
index_family(::Hill) = :effective_diversity
index_family(::Chao1) = :richness_estimator
index_family(::ACE) = :richness_estimator
index_family(::SampleCoverage) = :coverage
index_family(::PielouEvenness) = :evenness
index_family(::FisherAlpha) = :diversity
index_family(::Union{Jaccard,SorensenDice,Overlap}) = :incidence
index_family(::Union{BrayCurtis,Ruzicka,Canberra,MorisitaHorn}) = :abundance
index_family(::Union{TotalVariation,Manhattan,Euclidean,Hellinger,Chord,Bhattacharyya,KullbackLeibler,ShannonDifference,JensenDifference,JensenShannon}) = :probability

"""
    input_mode(index)

Return the expected input style for an index: `:single_assemblage`,
`:pairwise`, or `:either`.
"""
input_mode(::DiversityIndex) = :single_assemblage
input_mode(::PairwiseIndex) = :pairwise

"""
    output_mode(index)

Return the conventional output form, such as `:entropy`, `:diversity`,
`:similarity`, `:dissimilarity`, `:distance`, `:coefficient`, or `:estimate`.
"""
output_mode(::DiversityIndex) = :estimate
output_mode(::Union{Shannon,Renyi,Tsallis}) = :entropy
output_mode(::Union{Richness,Hill,InverseSimpson,FisherAlpha}) = :diversity
output_mode(::Simpson) = :dominance
output_mode(::Union{GiniSimpson,GreenbergDiversityIndex,LinguisticDiversityIndex}) = :diversity
output_mode(::SampleCoverage) = :coverage
output_mode(::PielouEvenness) = :evenness
output_mode(::Union{Jaccard,SorensenDice,Overlap,Ruzicka,MorisitaHorn}) = :similarity
output_mode(::Union{TotalVariation,Manhattan,Euclidean,Hellinger,Chord,JensenShannon}) = :distance
output_mode(::Union{BrayCurtis,Canberra,KullbackLeibler,ShannonDifference,JensenDifference}) = :dissimilarity
output_mode(::Bhattacharyya) = :coefficient

"""
    is_finite(index)

Return whether the index is expected to return finite values for valid finite
inputs. This is separate from [`is_bounded`](@ref): an index can be finite on
every finite data set while having no fixed finite upper bound.
"""
is_finite(::DiversityIndex) = true
is_finite(::KullbackLeibler) = false

"""
    is_metric(index)

Return whether the package's distance/dissimilarity form is a metric under the
usual assumptions for the index.
"""
is_metric(::DiversityIndex) = false
is_metric(::Union{Jaccard,TotalVariation,Manhattan,Euclidean,Canberra,Hellinger,Chord,JensenShannon}) = true

"""
    is_triangular(index)

Return whether the package's distance/dissimilarity form is known to obey the
triangle inequality. Returns `:unknown` when this package does not encode a
claim.
"""
is_triangular(index::DiversityIndex) = is_metric(index) ? true : :unknown
is_triangular(::Union{BrayCurtis,KullbackLeibler,JensenDifference}) = false
is_triangular(::ShannonDifference) = true

"""
    is_nonnegative(index)

Return whether the index output is known to be nonnegative.
"""
is_nonnegative(::DiversityIndex) = true

"""
    is_bounded(index)

Return whether the conventional output range has a finite upper bound.
"""
function is_bounded(index::DiversityIndex)
    bounds = index_range(index)
    return isfinite(bounds.lower) && isfinite(bounds.upper)
end

"""
    is_pseudometric(index)

Return whether the distance/dissimilarity form is known to be a pseudometric:
nonnegative, symmetric, zero for identical inputs, and triangular, while
allowing distinct inputs to have zero distance.
"""
is_pseudometric(index::DiversityIndex) = is_metric(index) ? true : :unknown
is_pseudometric(::ShannonDifference) = true
is_pseudometric(::Union{BrayCurtis,KullbackLeibler,JensenDifference}) = false

"""
    is_quasimetric(index)

Return whether the distance/dissimilarity form is known to be a quasimetric:
nonnegative, zero only for identical inputs, and triangular, without requiring
symmetry. Metrics are also quasimetrics under this convention.
"""
is_quasimetric(index::DiversityIndex) = is_metric(index) ? true : :unknown
is_quasimetric(::Union{BrayCurtis,KullbackLeibler,ShannonDifference,JensenDifference}) = false

"""
    is_metametric(index)

Return whether the distance/dissimilarity form is known to be a metametric
under this package's convention: nonnegative, symmetric, and zero for identical
inputs, without requiring identity of indiscernibles or the triangle
inequality. Returns `:unknown` where that classification is not encoded.
"""
is_metametric(index::DiversityIndex) = is_metric(index) ? true : :unknown
is_metametric(::Union{BrayCurtis,Overlap,Ruzicka,Canberra,ShannonDifference,JensenDifference}) = true
is_metametric(::KullbackLeibler) = false

"""
    is_semimetric(index)

Return whether the distance/dissimilarity form is known to be a semimetric
under this package's convention: nonnegative, symmetric, zero only for
identical inputs, but not necessarily triangular.
"""
is_semimetric(index::DiversityIndex) = is_metric(index) ? true : :unknown
is_semimetric(::Union{BrayCurtis,Ruzicka,Canberra}) = true
is_semimetric(::Union{Overlap,KullbackLeibler,ShannonDifference,JensenDifference}) = false

"""
    is_premetric(index)

Return whether the distance/dissimilarity form is known to be a premetric:
nonnegative and zero for identical inputs, without requiring symmetry or the
triangle inequality.
"""
is_premetric(index::DiversityIndex) = is_metric(index) ? true : :unknown
is_premetric(::Union{BrayCurtis,Overlap,Ruzicka,Canberra,KullbackLeibler,ShannonDifference,JensenDifference}) = true

"""
    is_supermetric(index)

Return whether the index is known to obey a supermetric or reverse-triangle
style condition. This property is uncommon for the indices implemented here, so
unknown cases return `:unknown`.
"""
is_supermetric(::DiversityIndex) = :unknown
is_supermetric(::Union{Jaccard,TotalVariation,Manhattan,Euclidean,Canberra,Hellinger,Chord,KullbackLeibler,ShannonDifference,JensenDifference,JensenShannon}) = false

"""
    is_similarity(index)

Return whether the primary output mode is a similarity or similarity
coefficient.
"""
is_similarity(index::DiversityIndex) = output_mode(index) in (:similarity, :coefficient)

"""
    is_dissimilarity(index)

Return whether the primary output mode is a dissimilarity or distance.
"""
is_dissimilarity(index::DiversityIndex) = output_mode(index) in (:dissimilarity, :distance)

"""
    is_dissimiliarty(index)

Deprecated misspelling of [`is_dissimilarity`](@ref), retained as a forgiving
alias.
"""
is_dissimiliarty(index::DiversityIndex) = is_dissimilarity(index)

"""
    is_symmetric(index)

Return whether the pairwise form satisfies `f(a, b) == f(b, a)` for all inputs.

Most indices are symmetric. [`KullbackLeibler`](@ref) is the notable exception:
`dissimilarity(KullbackLeibler(), a, b)` computes ``D_{KL}(a \\Vert b)`` which
generally differs from ``D_{KL}(b \\Vert a)``. Community distance matrices for
asymmetric indices are not symmetric matrices.
"""
is_symmetric(::DiversityIndex) = true
is_symmetric(::KullbackLeibler) = false

"""
    index_range(index)

Return a conventional numeric range for the index output when it is known.
"""
index_range(::DiversityIndex) = (lower=0.0, upper=Inf)
index_range(::Union{PielouEvenness,SampleCoverage,GiniSimpson,GreenbergDiversityIndex,LinguisticDiversityIndex,Jaccard,SorensenDice,Overlap,BrayCurtis,Ruzicka,TotalVariation,Canberra,Hellinger,Bhattacharyya,MorisitaHorn}) =
    (lower=0.0, upper=1.0)
index_range(::Manhattan) = (lower=0.0, upper=2.0)
index_range(::Union{Euclidean,Chord}) = (lower=0.0, upper=sqrt(2))
index_range(::Union{KullbackLeibler,ShannonDifference}) = (lower=0.0, upper=Inf)
index_range(index::JensenDifference) = (lower=0.0, upper=log(2) / log(index.base))
index_range(index::JensenShannon) = (
    lower=0.0,
    upper=index.distance ? sqrt(log(2) / log(index.base)) : log(2) / log(index.base),
)

"""
    index_bounds(index)

Return a named tuple describing the conventional numeric bounds and the usual
interpretation of those bounds. Unknown meanings are returned as `:unknown`.
"""
function index_bounds(index::DiversityIndex)
    range = index_range(index)
    return (
        lower=range.lower,
        upper=range.upper,
        lower_meaning=_lower_bound_meaning(index),
        upper_meaning=_upper_bound_meaning(index),
    )
end

function _lower_bound_meaning(index::DiversityIndex)
    mode = output_mode(index)
    if mode in (:similarity, :coefficient)
        return "minimal similarity; conventionally complete dissimilarity or no overlap"
    elseif mode in (:dissimilarity, :distance)
        return "minimal dissimilarity; identical or indistinguishable inputs"
    elseif mode == :entropy
        return "no uncertainty; all mass in one category"
    elseif mode == :diversity
        return "minimal diversity under the index convention"
    elseif mode == :dominance
        return "minimal dominance"
    elseif mode == :coverage
        return "no sampled probability mass covered"
    elseif mode == :evenness
        return "minimal evenness"
    else
        return :unknown
    end
end

function _upper_bound_meaning(index::DiversityIndex)
    mode = output_mode(index)
    if mode in (:similarity, :coefficient)
        return "maximal similarity; conventionally identical or complete overlap"
    elseif mode in (:dissimilarity, :distance)
        return isfinite(index_range(index).upper) ?
               "maximal dissimilarity under the index convention" :
               "unbounded dissimilarity; larger values mean greater separation"
    elseif mode == :entropy
        return "maximum uncertainty for the supplied or observed support"
    elseif mode == :diversity
        return isfinite(index_range(index).upper) ?
               "maximal diversity under the index convention" :
               "unbounded effective or richness-scale diversity"
    elseif mode == :dominance
        return "maximal dominance"
    elseif mode == :coverage
        return "complete sampled probability mass covered"
    elseif mode == :evenness
        return "maximal evenness"
    else
        return :unknown
    end
end

"""
    requires_probabilities(index)

Return whether an index is naturally defined on normalized probability vectors.
"""
requires_probabilities(::DiversityIndex) = false
requires_probabilities(::Union{Shannon,Renyi,Tsallis,Simpson,GiniSimpson,GreenbergDiversityIndex,LinguisticDiversityIndex,InverseSimpson,Hill,PielouEvenness,TotalVariation,Manhattan,Euclidean,Hellinger,Chord,Bhattacharyya,KullbackLeibler,ShannonDifference,JensenDifference,JensenShannon}) = true

"""
    supports_matrix_kernel(index)

Return whether this package has a specialized matrix implementation for the
index.
"""
supports_matrix_kernel(::DiversityIndex) = false
supports_matrix_kernel(::Union{Richness,Shannon,Jaccard,BrayCurtis,Hellinger}) = true

"""
    index_metadata(index)

Return convention-aware metadata for an index: family, input and output modes,
range, metric status, formula, aliases, and implementation notes.
"""
function index_metadata(index::DiversityIndex)
    return (
        type=typeof(index),
        family=index_family(index),
        input_mode=input_mode(index),
        output_mode=output_mode(index),
        range=index_range(index),
        bounds=index_bounds(index),
        is_finite=is_finite(index),
        is_metric=is_metric(index),
        is_triangular=is_triangular(index),
        is_nonnegative=is_nonnegative(index),
        is_bounded=is_bounded(index),
        is_pseudometric=is_pseudometric(index),
        is_quasimetric=is_quasimetric(index),
        is_metametric=is_metametric(index),
        is_semimetric=is_semimetric(index),
        is_premetric=is_premetric(index),
        is_supermetric=is_supermetric(index),
        is_similarity=is_similarity(index),
        is_dissimilarity=is_dissimilarity(index),
        is_symmetric=is_symmetric(index),
        requires_probabilities=requires_probabilities(index),
        supports_matrix_kernel=supports_matrix_kernel(index),
        formula=_index_formula(index),
        aliases=_index_aliases(index),
        notes=_index_notes(index),
    )
end

_index_formula(::DiversityIndex) = ""
_index_formula(::Richness) = "S_obs = number of positive-abundance taxa"
_index_formula(::Shannon) = "H_b = -sum_i p_i log_b(p_i)"
_index_formula(::Renyi) = "H_q = log_b(sum_i p_i^q)/(1-q)"
_index_formula(::Tsallis) = "T_q = (sum_i p_i^q - 1)/((1-q) log b)"
_index_formula(::Simpson) = "D = sum_i p_i^2"
_index_formula(::GiniSimpson) = "1 - sum_i p_i^2"
_index_formula(::GreenbergDiversityIndex) = "1 - sum_i p_i^2"
_index_formula(::LinguisticDiversityIndex) = "1 - sum_i p_i^2"
_index_formula(::InverseSimpson) = "1 / sum_i p_i^2"
_index_formula(::Hill) = "^qD = (sum_i p_i^q)^(1/(1-q))"
_index_formula(::Chao1) = "S_obs + f1(f1-1)/(2(f2+1))"
_index_formula(::ACE) = "abundance-based coverage estimator"
_index_formula(::SampleCoverage) = "C_hat = 1 - f1/n"
_index_formula(::PielouEvenness) = "J = H/log(S)"
_index_formula(::FisherAlpha) = "S = alpha log(1+n/alpha)"
_index_formula(::Jaccard) = "|A∩B| / |A∪B|"
_index_formula(::SorensenDice) = "2|A∩B| / (|A|+|B|)"
_index_formula(::Overlap) = "|A∩B| / min(|A|, |B|)"
_index_formula(::BrayCurtis) = "sum_i |x_i-y_i| / sum_i (x_i+y_i)"
_index_formula(::Ruzicka) = "sum_i min(x_i,y_i) / sum_i max(x_i,y_i)"
_index_formula(::TotalVariation) = "0.5 sum_i |p_i-q_i|"
_index_formula(::Hellinger) = "sqrt(sum_i (sqrt(p_i)-sqrt(q_i))^2 / 2)"
_index_formula(::KullbackLeibler) = "sum_i p_i log_b(p_i/q_i)"
_index_formula(::ShannonDifference) = "|H(p)-H(q)|"
_index_formula(::JensenDifference) = "H((p+q)/2) - (H(p)+H(q))/2"
_index_formula(::JensenShannon) = "sqrt((KL(p||m)+KL(q||m))/2) by default"

_index_aliases(::DiversityIndex) = String[]
_index_aliases(::Shannon) = ["vegan: diversity(index=\"shannon\")", "scikit-bio: shannon"]
_index_aliases(::GiniSimpson) = ["vegan: diversity(index=\"simpson\")", "scikit-bio: simpson"]
_index_aliases(::GreenbergDiversityIndex) = ["Linguistic Diversity Index", "LDI", "Greenberg's Diversity Index", "Gini-Simpson diversity"]
_index_aliases(::LinguisticDiversityIndex) = ["Greenberg's Diversity Index", "LDI", "Gini-Simpson diversity"]
_index_aliases(::InverseSimpson) = ["vegan: diversity(index=\"invsimpson\")", "scikit-bio: inv_simpson"]
_index_aliases(::FisherAlpha) = ["vegan: fisher.alpha", "scikit-bio: fisher_alpha"]
_index_aliases(::Jaccard) = ["vegan: vegdist(method=\"jaccard\", binary=TRUE)", "scikit-bio: jaccard"]
_index_aliases(::BrayCurtis) = ["vegan: vegdist(method=\"bray\")", "scikit-bio: braycurtis", "scipy: braycurtis"]
_index_aliases(::Hellinger) = ["vegan: vegdist(decostand(x, \"hellinger\"), method=\"euclidean\")"]

_index_notes(::DiversityIndex) = ""
_index_notes(::Simpson) = "This package's Simpson() is concentration. Use GiniSimpson() for vegan's index=\"simpson\" convention."
_index_notes(::GreenbergDiversityIndex) = "Equivalent to GiniSimpson(); interpreted as the probability that two randomly selected people have different mother tongues."
_index_notes(::LinguisticDiversityIndex) = "Equivalent to GreenbergDiversityIndex() and GiniSimpson(); the temporal Index of Linguistic Diversity is available as index_of_linguistic_diversity(current, baseline)."
_index_notes(::KullbackLeibler) = "Asymmetric (is_symmetric=false): dissimilarity(KullbackLeibler(), left, right) returns D_KL(left || right). Community matrices are not symmetric. Use estimator for Miller-Madow, pseudocount/shrinkage, or Good-Turing corrections."
_index_notes(::ShannonDifference) = "Measures the difference in entropy magnitudes |H(p)-H(q)|, not distributional divergence. Two assemblages with disjoint species but identical abundance profiles score zero. Use JensenShannon or KullbackLeibler for distributional divergence."
_index_notes(::JensenDifference) = "For Shannon entropy, the Jensen difference equals Jensen-Shannon divergence. Returns the raw divergence. Use JensenShannon (distance=true) for the metric square-root form. Use estimator for low-sample corrections."
_index_notes(::JensenShannon) = "JensenShannon(distance=true) returns the square root of the divergence (a metric). JensenShannon(distance=false) is identical to JensenDifference. Use estimator for low-sample corrections."

"""
    reference_cases()

Return curated cross-package and formula reference cases used to validate
conventions.
"""
function reference_cases()
    return [
        (name="vegan_shannon_natural_log", index=Shannon(; base=ℯ), data=[1, 1, 2], op=entropy, expected=1.0397207708399179, atol=1e-12),
        (name="vegan_gini_simpson", index=GiniSimpson(), data=[1, 1, 2], op=diversity, expected=0.625, atol=1e-12),
        (name="greenberg_linguistic_diversity", index=GreenbergDiversityIndex(), data=[1, 1, 2], op=diversity, expected=0.625, atol=1e-12),
        (name="vegan_inverse_simpson", index=InverseSimpson(), data=[1, 1, 2], op=diversity, expected=2.6666666666666665, atol=1e-12),
        (name="bray_curtis_formula", index=BrayCurtis(), left=[1, 2, 3], right=[2, 2, 0], op=dissimilarity, expected=0.4, atol=1e-12),
        (name="jaccard_incidence_formula", index=Jaccard(), left=[1, 1, 0, 1], right=[1, 0, 1, 1], op=similarity, expected=0.5, atol=1e-12),
        (name="sorensen_incidence_formula", index=SorensenDice(), left=[1, 1, 0, 1], right=[1, 0, 1, 1], op=similarity, expected=2 / 3, atol=1e-12),
        (name="pielou_formula", index=PielouEvenness(), data=[1, 1, 2], op=diversity, expected=entropy(Shannon(; base=ℯ), [1, 1, 2]) / log(3), atol=1e-12),
    ]
end

"""
    validate_reference_cases(; cases=reference_cases())

Evaluate reference cases and return one result named tuple per case.
"""
function validate_reference_cases(; cases=reference_cases())
    return map(cases) do case
        observed = if haskey(case, :left)
            case.op(case.index, case.left, case.right)
        else
            case.op(case.index, case.data)
        end
        passed = isapprox(observed, case.expected; atol=case.atol, rtol=0)
        (name=case.name, observed=observed, expected=case.expected, atol=case.atol, passed=passed)
    end
end

"""
    estimator_report(data; support=nothing, base=2, frequencies=true, species=nothing)

Return a compact report comparing Shannon entropy estimators and basic
coverage diagnostics for one assemblage.
"""
function estimator_report(data; support=nothing, base::Real=2, frequencies::Bool=true, species=nothing)
    if _is_table(data)
        return estimator_report(community_matrix(data; species); support, base, frequencies)
    end
    if data isa AbstractMatrix{<:Real}
        _check_matrix_frequencies(frequencies)
        return [estimator_report(row; support=_matrix_support(support), base, frequencies) for row in _rows(data)]
    end
    abundance = _abundances(data; frequencies)
    total = sum(abundance)
    f1 = count(==(1), abundance)
    f2 = count(==(2), abundance)
    estimator_specs = Any[
        (name=:plugin, estimator=Plugin(), support=nothing),
        (name=:miller_madow, estimator=MillerMadow(), support=nothing),
        (name=:chao_shen, estimator=ChaoShen(), support=nothing),
    ]
    if support !== nothing
        append!(estimator_specs, [
            (name=:basharin, estimator=Basharin(), support=support),
            (name=:hausser_strimmer, estimator=HausserStrimmer(), support=support),
            (name=:add_one, estimator=AddGamma(1), support=support),
        ])
    end
    estimates = map(estimator_specs) do spec
        value = spec.support === nothing ?
            entropy(Shannon(; base, estimator=spec.estimator), data; frequencies) :
            entropy(Shannon(; base, estimator=spec.estimator), data; frequencies, support=spec.support)
        (name=spec.name, estimator=typeof(spec.estimator), entropy=value, effective_diversity=base^value)
    end
    warnings = String[]
    f1 > 0 && push!(warnings, "singletons present; unseen taxa may affect richness and entropy estimates")
    support === nothing && push!(warnings, "no known finite support supplied")
    return (
        total_abundance=total,
        observed_richness=length(abundance),
        singletons=f1,
        doubletons=f2,
        sample_coverage=1 - f1 / total,
        estimates=estimates,
        warnings=warnings,
    )
end

"""
    diversity_audit(data; species=nothing, labels=nothing, label=nothing,
                    pairwise_index=BrayCurtis())

Return validation diagnostics, alpha summaries, estimator diagnostics, and an
optional labeled pairwise matrix for a diversity workflow.
"""
function diversity_audit(data; species=nothing, labels=nothing, label=nothing,
        pairwise_index::DiversityIndex=BrayCurtis())
    matrix = _is_table(data) ? community_matrix(data; species) : community_matrix(data)
    row_totals = vec(sum(matrix; dims=2))
    alpha = alpha_diversity(matrix)
    warnings = String[]
    any(<=(0), row_totals) && push!(warnings, "one or more samples have non-positive total abundance")
    any(summary -> summary.sample_coverage < 0.8, alpha) &&
        push!(warnings, "one or more samples have sample coverage below 0.8")
    pairwise = labeled_distance(pairwise_index, data; labels, label, species)
    return (
        n_samples=size(matrix, 1),
        n_taxa=size(matrix, 2),
        row_totals=row_totals,
        alpha=alpha,
        estimator_report=estimator_report(matrix),
        pairwise=pairwise,
        warnings=warnings,
    )
end

"""
    uncertainty_audit(data; species=nothing, labels=nothing, label=nothing,
                      index=Shannon(), nboot=1000, level=0.95,
                      quantities=(:entropy, :diversity), rng=nothing)

Return bootstrap uncertainty summaries for Shannon entropy and effective
diversity alongside sample labels and coverage diagnostics. For matrices and
Tables.jl-compatible inputs, one report is returned per sample.
"""
function uncertainty_audit(data; species=nothing, labels=nothing, label=nothing,
        index::Shannon=Shannon(), nboot::Integer=1000, level::Real=0.95,
        quantities=(:entropy, :diversity), rng=nothing)
    matrix = _is_table(data) ? community_matrix(data; species) : community_matrix(data)
    sample_labels = _pairwise_labels(data; labels, label)
    rows = collect(_rows(matrix))
    reports = map(eachindex(rows)) do i
        row = rows[i]
        estimates = map(quantities) do quantity
            boot = bootstrap(index, row; nboot, level, quantity, rng)
            (
                quantity=quantity,
                estimate=boot.estimate,
                stderr=boot.stderr,
                lower=boot.lower,
                upper=boot.upper,
                level=boot.level,
            )
        end
        (
            label=sample_labels[i],
            total_abundance=sum(row),
            richness=richness(row),
            sample_coverage=sample_coverage(row),
            estimates=estimates,
        )
    end
    warnings = String[]
    any(report -> report.sample_coverage < 0.8, reports) &&
        push!(warnings, "one or more samples have sample coverage below 0.8")
    nboot < 200 && push!(warnings, "low nboot; intervals are suitable for workflow checks, not final inference")
    return (
        index=typeof(index),
        nboot=nboot,
        level=float(level),
        labels=sample_labels,
        reports=reports,
        warnings=warnings,
    )
end
