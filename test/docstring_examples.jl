using DiversityAndDissimilarity
using Test

@testset "source docstring examples" begin
    @testset "input utilities" begin
        @test community_matrix([1 2 3; 4 5 6]) == [1.0 2.0 3.0; 4.0 5.0 6.0]
        @test sort(collect(counts(["oak", "ash", "oak"])), by=first) ==
            ["ash" => 1, "oak" => 2]
        @test sort(collect(counts(Dict(:oak => 3, :ash => 1))), by=first) ==
            [:ash => 1.0, :oak => 3.0]
        @test proportions([3, 1]) == [0.75, 0.25]
        @test proportions([3 1; 2 2]) == [0.75 0.25; 0.5 0.5]

        validated = validate([1 1 0; 0 1 1])
        @test richness(validated) == [2, 2]
        @test dissimilarity(BrayCurtis(), validated) == [0.0 0.5; 0.5 0.0]
    end

    @testset "metadata traits" begin
        @test index_family(Shannon()) == :entropy
        @test index_family(BrayCurtis()) == :abundance
        @test index_family(Jaccard()) == :incidence
        @test input_mode(Shannon()) == :single_assemblage
        @test input_mode(BrayCurtis()) == :pairwise
        @test output_mode(Shannon()) == :entropy
        @test output_mode(Jaccard()) == :similarity
        @test output_mode(BrayCurtis()) == :dissimilarity
        @test output_mode(Hellinger()) == :distance
        @test is_finite(JensenShannon())
        @test !is_finite(KullbackLeibler())
        @test is_metric(Jaccard())
        @test !is_metric(BrayCurtis())
        @test is_metric(JensenShannon())
        @test is_triangular(Jaccard())
        @test !is_triangular(BrayCurtis())
        @test is_triangular(Overlap()) == :unknown
        @test is_nonnegative(Shannon())
        @test is_nonnegative(KullbackLeibler())
        @test is_bounded(Jaccard())
        @test !is_bounded(KullbackLeibler())
        @test is_bounded(JensenShannon())
        @test is_pseudometric(Jaccard())
        @test is_pseudometric(ShannonDifference())
        @test !is_pseudometric(BrayCurtis())
        @test is_pseudometric(Overlap()) == :unknown
        @test is_quasimetric(Jaccard())
        @test !is_quasimetric(BrayCurtis())
        @test is_quasimetric(Overlap()) == :unknown
        @test is_metametric(BrayCurtis())
        @test !is_metametric(KullbackLeibler())
        @test is_metametric(Jaccard())
        @test is_semimetric(BrayCurtis())
        @test !is_semimetric(Overlap())
        @test is_semimetric(Jaccard())
        @test is_premetric(BrayCurtis())
        @test is_premetric(KullbackLeibler())
        @test is_premetric(Shannon()) == :unknown
        @test !is_supermetric(Jaccard())
        @test is_supermetric(BrayCurtis()) == :unknown
        @test is_similarity(Jaccard())
        @test is_similarity(Bhattacharyya())
        @test !is_similarity(BrayCurtis())
        @test is_dissimilarity(BrayCurtis())
        @test is_dissimilarity(Hellinger())
        @test !is_dissimilarity(Jaccard())
        @test is_symmetric(BrayCurtis())
        @test !is_symmetric(KullbackLeibler())
        @test index_range(Jaccard()) == (lower=0.0, upper=1.0)
        @test index_range(Shannon()) == (lower=0.0, upper=Inf)
        @test index_range(JensenDifference()) == (lower=0.0, upper=1.0)
        @test requires_probabilities(Hellinger())
        @test !requires_probabilities(BrayCurtis())
        @test supports_matrix_kernel(BrayCurtis())
        @test !supports_matrix_kernel(Canberra())
    end

    @testset "metadata reports and audits" begin
        jaccard_bounds = index_bounds(Jaccard())
        @test (jaccard_bounds.lower, jaccard_bounds.upper) == (0.0, 1.0)
        @test jaccard_bounds.lower_meaning ==
            "minimal similarity; conventionally complete dissimilarity or no overlap"
        @test jaccard_bounds.upper_meaning ==
            "maximal similarity; conventionally identical or complete overlap"

        bray_bounds = index_bounds(BrayCurtis())
        @test bray_bounds.lower_meaning ==
            "minimal dissimilarity; identical or indistinguishable inputs"
        @test bray_bounds.upper_meaning ==
            "maximal dissimilarity under the index convention"

        kl_bounds = index_bounds(KullbackLeibler())
        @test kl_bounds.upper == Inf
        @test kl_bounds.upper_meaning ==
            "unbounded dissimilarity; larger values mean greater separation"

        metadata = index_metadata(BrayCurtis())
        @test metadata.family == :abundance
        @test metadata.output_mode == :dissimilarity
        @test !metadata.is_metric
        @test metadata.is_semimetric
        @test metadata.is_symmetric

        results = validate_reference_cases()
        @test all(result -> result.passed, results)
        @test results[1].name == "vegan_shannon_natural_log"

        report = estimator_report([1, 1, 2, 0, 5])
        @test report.observed_richness == 4
        @test report.singletons == 2
        @test report.estimates[1].name == :plugin

        community = [1 1 2 0 5; 3 0 1 1 0]
        audit = diversity_audit(community; labels=["a", "b"])
        @test audit.n_samples == 2
        @test audit.n_taxa == 5
        @test audit.pairwise.labels == ["a", "b"]

        uncertainty = uncertainty_audit(community; labels=["a", "b"], nboot=50, rng=nothing)
        @test length(uncertainty.reports) == 2
        @test uncertainty.labels == ["a", "b"]
        @test uncertainty.reports[1].label == "a"
    end
end

@testset "data-input documentation examples" begin
    assemblage = Dict(:oak => 12, :ash => 5, :elm => 3)
    @test richness(assemblage) == 3
    @test round(shannon_entropy(assemblage); digits=4) == 1.3527
    @test dissimilarity(
        BrayCurtis(),
        Dict(:oak => 10, :ash => 2),
        Dict(:ash => 5, :elm => 3),
    ) == 0.8

    observations = ["oak", "ash", "oak", "elm", "ash", "oak"]
    @test sort(collect(counts(observations)), by=first) ==
        ["ash" => 2, "elm" => 1, "oak" => 3]
    @test richness(observations) == 3
    @test round(shannon_entropy(observations); digits=4) == 1.4591

    obs = [1, 2, 1, 3]
    @test richness(obs; frequencies=false) == 3
    @test richness(obs) == 4
    @test round(shannon_entropy(obs; frequencies=false); digits=4) == 1.5
    @test round(shannon_entropy(obs); digits=4) == 1.8424

    community = [1 1 2 0 5; 3 0 1 1 0]
    @test richness(community) == [4, 3]
    @test round.(shannon_entropy(community); digits=4) == [1.6577, 1.371]
    @test getproperty.(alpha_diversity(community), :richness) == [4, 3]
end
