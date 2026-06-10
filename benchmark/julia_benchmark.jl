using DiversityAndDissimilarity

function simulated_community(nsites::Integer=400, ntaxa::Integer=200; total::Integer=10_000)
    matrix = Matrix{Float64}(undef, nsites, ntaxa)
    for site in 1:nsites
        weights = [1 / (1 + mod(taxon + 3site, ntaxa))^1.15 for taxon in 1:ntaxa]
        scale = total / sum(weights)
        for taxon in 1:ntaxa
            value = round(Int, scale * weights[taxon])
            matrix[site, taxon] = max(value, taxon == site % ntaxa + 1 ? 1 : 0)
        end
    end
    return matrix
end

function best_time(f; repeats::Integer=5, inner::Integer=1)
    inner >= 1 || throw(ArgumentError("inner must be positive"))
    best = Inf
    for _ in 1:repeats
        elapsed = @elapsed begin
            for _ in 1:inner
                f()
            end
        end
        best = min(best, elapsed / inner)
    end
    return best
end

function warmup(community, validated)
    richness(community)
    shannon_entropy(community)
    alpha_diversity(community)
    bray_curtis_distance(community)
    jaccard_distance(community)
    hellinger_distance(community)
    richness(validated)
    shannon_entropy(validated)
    alpha_diversity(validated)
    bray_curtis_distance(validated)
    jaccard_distance(validated)
    hellinger_distance(validated)
    return nothing
end

function main()
    nsites = parse(Int, get(ENV, "DIVERSITY_BENCH_NSITES", "400"))
    ntaxa = parse(Int, get(ENV, "DIVERSITY_BENCH_NTAXA", "200"))
    total = parse(Int, get(ENV, "DIVERSITY_BENCH_TOTAL", "10000"))
    repeats = parse(Int, get(ENV, "DIVERSITY_BENCH_REPEATS", "5"))
    inner = parse(Int, get(ENV, "DIVERSITY_BENCH_INNER", "1"))
    community = simulated_community(nsites, ntaxa; total)
    # validate once; the Validated pathway skips per-call validation
    validated = validate(community)
    warmup(community, validated)

    println("language,package,task,nsites,ntaxa,repeats,inner,best_seconds")

    # Safe pathway: validates on every call (default user experience)
    for (task, f) in [
        ("richness",                  () -> richness(community)),
        ("shannon_entropy",           () -> shannon_entropy(community)),
        ("alpha_diversity",           () -> alpha_diversity(community)),
        ("bray_curtis_distance_matrix", () -> bray_curtis_distance(community)),
        ("jaccard_distance_matrix",   () -> jaccard_distance(community)),
        ("hellinger_distance_matrix", () -> hellinger_distance(community)),
    ]
        println(join(("Julia", "DiversityAndDissimilarity.jl", task, nsites, ntaxa, repeats,
            inner, best_time(f; repeats, inner)), ","))
    end

    # Pre-validated pathway: validate once, then computation only (fair comparison with
    # Python/NumPy and R/vegan, which do not re-validate on each call)
    for (task, f) in [
        ("richness_prevalidated",                  () -> richness(validated)),
        ("shannon_entropy_prevalidated",           () -> shannon_entropy(validated)),
        ("alpha_diversity_prevalidated",           () -> alpha_diversity(validated)),
        ("bray_curtis_distance_matrix_prevalidated", () -> bray_curtis_distance(validated)),
        ("jaccard_distance_matrix_prevalidated",   () -> jaccard_distance(validated)),
        ("hellinger_distance_matrix_prevalidated", () -> hellinger_distance(validated)),
    ]
        println(join(("Julia", "DiversityAndDissimilarity.jl", task, nsites, ntaxa, repeats,
            inner, best_time(f; repeats, inner)), ","))
    end
end

main()
