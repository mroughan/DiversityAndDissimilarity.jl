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
input_mode(::Union{Jaccard,SorensenDice,Overlap,BrayCurtis,Ruzicka,TotalVariation,Manhattan,Euclidean,Canberra,Hellinger,Chord,Bhattacharyya,KullbackLeibler,ShannonDifference,JensenDifference,JensenShannon,MorisitaHorn}) = :pairwise

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
output_mode(::Union{BrayCurtis,TotalVariation,Manhattan,Euclidean,Canberra,Hellinger,Chord,JensenShannon}) = :distance
output_mode(::Union{KullbackLeibler,ShannonDifference,JensenDifference}) = :dissimilarity
output_mode(::Bhattacharyya) = :coefficient

"""
    is_metric(index)

Return whether the package's distance/dissimilarity form is a metric under the
usual assumptions for the index.
"""
is_metric(::DiversityIndex) = false
is_metric(::Union{Jaccard,TotalVariation,Manhattan,Euclidean,Hellinger,Chord,JensenShannon}) = true

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
        is_metric=is_metric(index),
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
_index_notes(::KullbackLeibler) = "Asymmetric: dissimilarity(KullbackLeibler(), left, right) returns D_KL(left || right) and may be Inf when right has zero probability where left is positive."
_index_notes(::JensenDifference) = "For Shannon entropy, the Jensen difference equals Jensen-Shannon divergence."
_index_notes(::JensenShannon) = "JensenShannon(distance=true) returns the square root of the divergence."

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
