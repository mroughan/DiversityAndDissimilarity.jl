using DiversityIndices

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

function warmup(community)
    richness(community)
    shannon_entropy(community)
    alpha_diversity(community)
    bray_curtis_distance(community)
    jaccard_distance(community)
    hellinger_distance(community)
    return nothing
end

function main()
    nsites = parse(Int, get(ENV, "DIVERSITY_BENCH_NSITES", "400"))
    ntaxa = parse(Int, get(ENV, "DIVERSITY_BENCH_NTAXA", "200"))
    total = parse(Int, get(ENV, "DIVERSITY_BENCH_TOTAL", "10000"))
    repeats = parse(Int, get(ENV, "DIVERSITY_BENCH_REPEATS", "5"))
    inner = parse(Int, get(ENV, "DIVERSITY_BENCH_INNER", "1"))
    community = simulated_community(nsites, ntaxa; total)
    warmup(community)

    println("language,package,task,nsites,ntaxa,repeats,inner,best_seconds")
    println(join(("Julia", "DiversityIndices.jl", "richness", nsites, ntaxa, repeats,
        inner, best_time(() -> richness(community); repeats, inner)), ","))
    println(join(("Julia", "DiversityIndices.jl", "shannon_entropy", nsites, ntaxa, repeats,
        inner, best_time(() -> shannon_entropy(community); repeats, inner)), ","))
    println(join(("Julia", "DiversityIndices.jl", "alpha_diversity", nsites, ntaxa, repeats,
        inner, best_time(() -> alpha_diversity(community); repeats, inner)), ","))
    println(join(("Julia", "DiversityIndices.jl", "bray_curtis_distance_matrix", nsites, ntaxa, repeats,
        inner, best_time(() -> bray_curtis_distance(community); repeats, inner)), ","))
    println(join(("Julia", "DiversityIndices.jl", "jaccard_distance_matrix", nsites, ntaxa, repeats,
        inner, best_time(() -> jaccard_distance(community); repeats, inner)), ","))
    println(join(("Julia", "DiversityIndices.jl", "hellinger_distance_matrix", nsites, ntaxa, repeats,
        inner, best_time(() -> hellinger_distance(community); repeats, inner)), ","))
end

main()
