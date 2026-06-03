using DataFrames
using DiversityIndices
using Test

function test_valid_probabilities(data; frequencies=true)
    probabilities = proportions(data; frequencies)
    @test all(>=(0), probabilities)
    @test all(isfinite, probabilities)
    @test sum(probabilities) ≈ 1
    return probabilities
end

@testset "abundance handling" begin
    @test counts(["a", "b", "a"]) == Dict("a" => 2, "b" => 1)
    @test proportions(Dict(:a => 2, :b => 2)) == [0.5, 0.5]
    @test proportions([1 1 2 0 5; 3 0 1 1 0]) ≈ [
        1 / 9 1 / 9 2 / 9 0 5 / 9
        3 / 5 0 1 / 5 1 / 5 0
    ]
    @test test_valid_probabilities(Dict(:a => 2, :b => 2)) == [0.5, 0.5]
    @test test_valid_probabilities([1, 2, 3]) ≈ [1 / 6, 2 / 6, 3 / 6]
    @test test_valid_probabilities([0.1, 0.2, 0.3]) ≈ [1 / 6, 2 / 6, 3 / 6]
    @test sort(test_valid_probabilities(["a", "b", "a"])) == [1 / 3, 2 / 3]
    @test sort(test_valid_probabilities([1, 2, 1, 3]; frequencies=false)) == [0.25, 0.25, 0.5]
    @test_throws ArgumentError proportions([0, 0, 0])
    @test_throws ArgumentError proportions([-1, 2])
    @test_throws ArgumentError proportions([Inf, 1])
    @test_throws ArgumentError proportions([NaN, 1])
    @test_throws ArgumentError proportions([1 2; 0 0])
    @test_throws ArgumentError proportions([1 -1; 1 1])
end

@testset "alpha diversity" begin
    data = Dict(:a => 2, :b => 2)

    @test diversity(Richness(), data) == 2
    @test Shannon().base == 2
    @test entropy(Shannon(), data) ≈ 1
    @test diversity(Shannon(), data) ≈ 2
    @test effective_diversity(Shannon(), data) ≈ 2
    @test entropy(Shannon(; base=2), data) ≈ 1
    @test diversity(Shannon(; base=2), data) ≈ 2
    @test effective_diversity(Shannon(; base=2), data) ≈ 2
    @test Shannon().estimator isa Plugin
    @test entropy(Shannon(; estimator=MillerMadow()), data) ≈ 1 + 1 / (2 * 4 * log(2))
    @test entropy(Shannon(; estimator=Basharin()), data; support=4) ≈ 1 + 3 / (4 * log(2))
    @test entropy(Shannon(; estimator=AddGamma(1)), data; support=4) ≈
        -2 * (3 / 8) * log2(3 / 8) - 2 * (1 / 8) * log2(1 / 8)
    @test entropy(Shannon(; estimator=HausserStrimmer()), data; support=4) ≈
        -2 * (1 / 3) * log2(1 / 3) - 2 * (1 / 6) * log2(1 / 6)
    @test entropy(Shannon(; estimator=ChaoShen()), data) ≈ 1 / (1 - 0.5^4)
    @test shannon_entropy(data; estimator=MillerMadow()) ≈ entropy(Shannon(; estimator=MillerMadow()), data)
    @test shannon(data; estimator=MillerMadow()) ≈ shannon_entropy(data; estimator=MillerMadow())
    @test shannon_diversity(data; estimator=MillerMadow()) ≈
        effective_diversity(Shannon(; estimator=MillerMadow()), data)
    @test_throws ArgumentError entropy(Shannon(; estimator=AddGamma(1)), data; support=1)

    shannon_estimators = (
        (Plugin(), nothing),
        (MillerMadow(), nothing),
        (Basharin(), 4),
        (AddGamma(1), 4),
        (HausserStrimmer(), 4),
        (ChaoShen(), nothing),
    )
    for (estimator, support) in shannon_estimators
        entropy_value = support === nothing ?
            entropy(Shannon(; estimator), data) :
            entropy(Shannon(; estimator), data; support)
        diversity_value = support === nothing ?
            diversity(Shannon(; estimator), data) :
            diversity(Shannon(; estimator), data; support)
        effective = support === nothing ?
            effective_diversity(Shannon(; estimator), data) :
            effective_diversity(Shannon(; estimator), data; support)
        wrapper = support === nothing ?
            shannon_diversity(data; estimator) :
            shannon_diversity(data; estimator, support)

        @test isfinite(entropy_value)
        @test diversity_value ≈ effective
        @test effective ≈ 2^entropy_value
        @test wrapper ≈ effective
    end

    @test Renyi(2).base == 2
    @test entropy(Renyi(2), data) ≈ 1
    @test diversity(Renyi(2), data) ≈ 2
    @test entropy(Renyi(1), data) ≈ entropy(Shannon(), data)
    @test effective_diversity(Renyi(2), data) ≈ diversity(Hill(2), data)
    @test renyi_entropy(data, 2) ≈ entropy(Renyi(2), data)
    @test renyi(data, 2) ≈ renyi_entropy(data, 2)
    @test renyi_diversity(data, 2) ≈ effective_diversity(Renyi(2), data)
    @test Tsallis(2).base == 2
    @test entropy(Tsallis(2), data) ≈ 1 / (2log(2))
    @test diversity(Tsallis(2), data) ≈ 2
    @test entropy(Tsallis(1), data) ≈ entropy(Shannon(), data)
    @test effective_diversity(Tsallis(2), data) ≈ diversity(Hill(2), data)
    @test tsallis_entropy(data, 2) ≈ entropy(Tsallis(2), data)
    @test tsallis(data, 2) ≈ tsallis_entropy(data, 2)
    @test tsallis_diversity(data, 2) ≈ effective_diversity(Tsallis(2), data)
    @test diversity(Simpson(), data) ≈ 0.5
    @test diversity(GiniSimpson(), data) ≈ 0.5
    @test diversity(InverseSimpson(), data) ≈ 2
    @test diversity(Hill(0), data) == diversity(Richness(), data)
    @test diversity(Hill(1), data) ≈ effective_diversity(Shannon(), data)
    @test diversity(Hill(2), data) ≈ diversity(InverseSimpson(), data)
    @test hill_number(data, 2) ≈ diversity(Hill(2), data)
    @test diversity(Chao1(), [1, 1, 2, 0, 5]) ≈ 4.5
    @test chao1([1, 1, 2, 0, 5]) ≈ diversity(Chao1(), [1, 1, 2, 0, 5])
    @test diversity(ACE(), [1, 1, 2, 0, 5]) ≈ 6.612244897959183
    @test ace([1, 1, 2, 0, 5]) ≈ diversity(ACE(), [1, 1, 2, 0, 5])
    @test sample_coverage([1, 1, 2, 0, 5]) ≈ 7 / 9
    @test diversity(SampleCoverage(), [1, 1, 2, 0, 5]) ≈ sample_coverage([1, 1, 2, 0, 5])
    @test pielou_evenness(data) ≈ 1
    @test diversity(PielouEvenness(), data) ≈ pielou_evenness(data)
    fisher_value = fisher_alpha([1, 1, 2, 0, 5])
    @test fisher_value * log1p(9 / fisher_value) ≈ 4
    @test diversity(FisherAlpha(), [1, 1, 2, 0, 5]) ≈ fisher_value

    summary = alpha_diversity(data)
    @test summary.richness == richness(data)
    @test summary.shannon_entropy ≈ shannon_entropy(data)
    @test summary.shannon_diversity ≈ shannon_diversity(data)
    @test summary.inverse_simpson ≈ inverse_simpson_index(data)
    @test summary.pielou_evenness ≈ pielou_evenness(data)
    @test summary.fisher_alpha ≈ fisher_alpha(data)
end

@testset "community matrices" begin
    matrix = [
        1 1 2 0 5
        3 0 1 1 0
    ]

    @test richness(matrix) == [4, 3]
    @test entropy(Shannon(), matrix) ≈ [entropy(Shannon(), matrix[1, :]), entropy(Shannon(), matrix[2, :])]
    @test diversity(Shannon(), matrix) ≈ [diversity(Shannon(), matrix[1, :]), diversity(Shannon(), matrix[2, :])]
    @test diversity(Simpson(), matrix) ≈ [0.3827160493827161, 0.44000000000000006]
    @test hill_number(matrix, 2) ≈ inverse_simpson_index(matrix)
    @test chao1(matrix) ≈ [4.5, 4.0]
    @test ace(matrix) ≈ [6.612244897959183, 6.666666666666667]
    @test sample_coverage(matrix) ≈ [7 / 9, 3 / 5]
    @test pielou_evenness(matrix) ≈ [pielou_evenness(matrix[1, :]), pielou_evenness(matrix[2, :])]
    @test fisher_alpha(matrix) ≈ [fisher_alpha(matrix[1, :]), fisher_alpha(matrix[2, :])]
    @test_throws ArgumentError richness(matrix; frequencies=false)
    @test entropy(Shannon(; estimator=AddGamma(1)), matrix; support=5) ≈
        [entropy(Shannon(; estimator=AddGamma(1)), matrix[1, :]; support=5),
            entropy(Shannon(; estimator=AddGamma(1)), matrix[2, :]; support=5)]
    @test entropy(Shannon(; estimator=AddGamma(1)), matrix; support=[:a, :b, :c, :d, :e]) ≈
        entropy(Shannon(; estimator=AddGamma(1)), matrix; support=5)

    @test distance(BrayCurtis(), matrix) ≈ [
        0.0 bray_curtis_distance(matrix[1, :], matrix[2, :])
        bray_curtis_distance(matrix[1, :], matrix[2, :]) 0.0
    ]
    @test jaccard_similarity(matrix) ≈ [
        1.0 jaccard_similarity(matrix[1, :], matrix[2, :])
        jaccard_similarity(matrix[1, :], matrix[2, :]) 1.0
    ]
    @test length(alpha_diversity(matrix)) == 2
    @test alpha_diversity(matrix)[1].richness == richness(matrix)[1]
    @test alpha_diversity(matrix)[1] == alpha_diversity(matrix[1, :])
end

@testset "Shannon uncertainty" begin
    data = [4, 3, 2, 1]
    plugin_variance = entropy_variance(Shannon(), data)
    basharin_variance = entropy_variance(Shannon(; estimator=Basharin()), data; support=4)

    @test plugin_variance >= 0
    @test basharin_variance > plugin_variance
    @test shannon_variance(data; estimator=Basharin(), support=4) ≈ basharin_variance
    @test entropy_variance(Shannon(; estimator=ChaoShen()), data) >= 0
    @test_throws ArgumentError entropy_variance(Shannon(; estimator=AddGamma(1)), data; support=4)

    interval = entropy_confint(Shannon(; estimator=Basharin()), data; support=4)
    @test interval.lower <= interval.estimate <= interval.upper
    @test interval.stderr ≈ sqrt(interval.variance)
    @test shannon_confint(data; estimator=Basharin(), support=4).estimate ≈
        entropy(Shannon(; estimator=Basharin()), data; support=4)

    boot = bootstrap(Shannon(), data; nboot=25)
    @test boot.lower <= boot.upper
    @test length(boot.replicates) == 25
    @test isfinite(boot.stderr)

    boot_diversity = bootstrap(Shannon(), data; nboot=25, quantity=:diversity)
    @test boot_diversity.estimate ≈ shannon_diversity(data)

    jack = jackknife(Shannon(), data)
    @test jack.lower <= jack.upper
    @test length(jack.leave_one_out) == 4
    @test isfinite(jack.bias_corrected)

    matrix = [
        4 3 2 1
        3 3 0 0
    ]
    @test length(entropy_confint(Shannon(), matrix)) == 2
    @test length(bootstrap(Shannon(), matrix; nboot=10)) == 2
    @test length(jackknife(Shannon(), matrix)) == 2
end

@testset "Tables and DataFrames integration" begin
    table = DataFrame(
        site=["a", "b"],
        block=[101, 102],
        oak=[1, 3],
        ash=[1, 0],
        elm=[2, 1],
        pine=[0, 1],
        birch=[5, 0],
    )
    matrix = [
        1 1 2 0 5
        3 0 1 1 0
    ]

    @test community_matrix(table; species=[:oak, :ash, :elm, :pine, :birch]) == float.(matrix)
    @test richness(table; species=[:oak, :ash, :elm, :pine, :birch]) == richness(matrix)
    @test shannon_entropy(table; species=[:oak, :ash, :elm, :pine, :birch]) ≈ shannon_entropy(matrix)
    @test inverse_simpson_index(table; species=[:oak, :ash, :elm, :pine, :birch]) ≈ inverse_simpson_index(matrix)
    @test chao1(table; species=[:oak, :ash, :elm, :pine, :birch]) ≈ chao1(matrix)
    @test pielou_evenness(table; species=[:oak, :ash, :elm, :pine, :birch]) ≈ pielou_evenness(matrix)
    @test bray_curtis_distance(table; species=[:oak, :ash, :elm, :pine, :birch]) ≈ bray_curtis_distance(matrix)
    @test distance(Jaccard(), table; species=[:oak, :ash, :elm, :pine, :birch]) ≈ distance(Jaccard(), matrix)
    labeled = labeled_distance(BrayCurtis(), table; label=:site, species=[:oak, :ash, :elm, :pine, :birch])
    @test labeled.labels == ["a", "b"]
    @test labeled.matrix ≈ bray_curtis_distance(matrix)
    @test labeled_similarity(Jaccard(), matrix; labels=["a", "b"]).labels == ["a", "b"]
    @test_throws ArgumentError community_matrix(table; species=[:site, :oak])
end

@testset "observation vectors" begin
    observations = ["a", "b", "a", "c"]

    @test richness(observations) == 3
    @test richness([1, 2, 1, 3]; frequencies=false) == 3
    @test richness([1, 2, 1, 3]) == 4
end

@testset "metadata and audit framework" begin
    metadata = index_metadata(BrayCurtis())
    @test metadata.family == :abundance
    @test metadata.input_mode == :pairwise
    @test metadata.supports_matrix_kernel
    @test !is_metric(BrayCurtis())
    @test is_metric(Jaccard())
    @test output_mode(Shannon()) == :entropy
    @test index_range(PielouEvenness()) == (lower=0.0, upper=1.0)
    @test all(result -> result.passed, validate_reference_cases())

    report = estimator_report([1, 1, 2, 0, 5]; support=6)
    @test report.observed_richness == 4
    @test report.singletons == 2
    @test any(estimate -> estimate.name == :add_one, report.estimates)

    linguistic = diversity(LinguisticDiversityIndex(), [1, 1, 2])
    @test linguistic ≈ diversity(GiniSimpson(), [1, 1, 2])
    @test greenberg_diversity_index([1, 1, 2]) ≈ 0.625
    @test linguistic_diversity_index([1, 1, 2]) ≈ 0.625
    @test index_of_linguistic_diversity([1, 1, 2], [1, 3]) ≈ 0.625 / 0.375
    @test index_metadata(LinguisticDiversityIndex()).family == :linguistic_diversity

    audit = diversity_audit([1 1 2 0 5; 3 0 1 1 0]; labels=["a", "b"])
    @test audit.n_samples == 2
    @test audit.n_taxa == 5
    @test audit.pairwise.labels == ["a", "b"]
    @test size(audit.pairwise.matrix) == (2, 2)

    uncertainty = uncertainty_audit([1 1 2 0 5; 3 0 1 1 0]; labels=["a", "b"], nboot=10)
    @test uncertainty.labels == ["a", "b"]
    @test length(uncertainty.reports) == 2
    @test uncertainty.reports[1].label == "a"
    @test any(result -> result.quantity == :entropy, uncertainty.reports[1].estimates)
    @test any(contains("low nboot"), uncertainty.warnings)
end

@testset "pairwise indexes" begin
    left = Dict(:a => 2, :b => 1)
    right = Dict(:b => 3, :c => 1)
    left_vector = [1, 2, 0]
    right_vector = [0, 1, 3]
    p = [1 / 3, 2 / 3, 0]
    q = [0, 1 / 4, 3 / 4]
    hellinger_expected = sqrt(sum(abs2, sqrt.(p) .- sqrt.(q)) / 2)
    jensen_shannon_expected =
        (
            sum(pi > 0 ? pi * log2(pi / ((pi + qi) / 2)) : 0 for (pi, qi) in zip(p, q)) +
            sum(qi > 0 ? qi * log2(qi / ((pi + qi) / 2)) : 0 for (pi, qi) in zip(p, q))
        ) / 2

    @test similarity(Jaccard(), left, right) ≈ 1 / 3
    @test dissimilarity(Jaccard(), left, right) ≈ 2 / 3
    @test similarity(SorensenDice(), left, right) ≈ 1 / 2
    @test dissimilarity(SorensenDice(), left, right) ≈ 1 / 2
    @test dissimilarity(BrayCurtis(), left, right) ≈ 5 / 7
    @test similarity(Overlap(), left, right) ≈ 1 / 2
    @test dissimilarity(Overlap(), left, right) ≈ 1 / 2
    @test similarity(Ruzicka(), left_vector, right_vector) ≈ 1 / 6
    @test dissimilarity(Ruzicka(), left_vector, right_vector) ≈ 5 / 6
    @test total_variation_distance(left_vector, right_vector) ≈ 3 / 4
    @test manhattan_distance(left_vector, right_vector) ≈ 2 * total_variation_distance(left_vector, right_vector)
    @test euclidean_distance(left_vector, right_vector) ≈ sqrt(61 / 72)
    @test canberra_distance(left_vector, right_vector) ≈ 7 / 9
    @test hellinger_distance(left_vector, right_vector) ≈ hellinger_expected
    @test chord_distance(left_vector, right_vector) ≈ sqrt(2) * hellinger_expected
    @test bhattacharyya_coefficient(left_vector, right_vector) ≈ sqrt(1 / 6)
    @test bhattacharyya_distance(left_vector, right_vector) ≈ -log(sqrt(1 / 6))
    @test jensen_shannon_divergence(left_vector, right_vector) ≈ jensen_shannon_expected
    @test jensen_shannon_distance(left_vector, right_vector) ≈ sqrt(jensen_shannon_expected)
    @test dissimilarity(JensenShannon(; distance=false), left_vector, right_vector) ≈ jensen_shannon_expected
    @test morisita_horn_similarity(left_vector, right_vector) ≈ 24 / 85
    @test morisita_horn_distance(left_vector, right_vector) ≈ 1 - 24 / 85

    @test jaccard_index(left, right) ≈ similarity(Jaccard(), left, right)
    @test jaccard_similarity(left, right) ≈ jaccard_index(left, right)
    @test distance(Jaccard(), left, right) ≈ dissimilarity(Jaccard(), left, right)
    @test sorensen_index(left, right) ≈ sorensen_dice_index(left, right)
    @test sorensen_distance(left, right) ≈ sorensen_dice_dissimilarity(left, right)
    @test sorensen_dice_distance(left, right) ≈ sorensen_dice_dissimilarity(left, right)
    @test bray_curtis_distance(left, right) ≈ bray_curtis_dissimilarity(left, right)
    @test bray_curtis_dissimilarity(left, right) ≈ dissimilarity(BrayCurtis(), left, right)
    @test overlap_distance(left, right) ≈ dissimilarity(Overlap(), left, right)
    @test ruzicka_similarity(left_vector, right_vector) ≈ quantitative_jaccard_similarity(left_vector, right_vector)
    @test ruzicka_distance(left_vector, right_vector) ≈ quantitative_jaccard_distance(left_vector, right_vector)
    @test_throws ArgumentError total_variation_distance([0, 0], [1, 0])
    @test_throws ArgumentError euclidean_distance([-1, 1], [1, 0])
end

@testset "reference values from vegan and scikit-bio" begin
    # scikit-bio community diversity tutorial example:
    # https://scikit.bio/docs/latest/diversity.html
    data = [
        [23, 64, 14, 0, 0, 3, 1],
        [0, 3, 35, 42, 0, 12, 1],
        [0, 5, 5, 0, 40, 40, 0],
        [44, 35, 9, 0, 1, 0, 0],
        [0, 2, 8, 0, 35, 45, 1],
        [0, 0, 25, 35, 0, 19, 0],
    ]

    @test richness.(data) == [5, 5, 4, 4, 5, 3]
    @test bray_curtis_distance(data[1], data[2]) ≈ 0.78787879 atol = 1e-8
    @test bray_curtis_distance(data[1], data[4]) ≈ 0.30927835 atol = 1e-8
    @test bray_curtis_distance(data[3], data[5]) ≈ 0.09392265 atol = 1e-8

    @test jaccard_distance(data[1], data[2]) ≈ 1 / 3
    @test jaccard_similarity(data[1], data[2]) ≈ 2 / 3

    # vegan documents Shannon with configurable logarithm base, Simpson as
    # 1 - sum(p_i^2), inverse Simpson as 1 / sum(p_i^2), and Bray-Curtis as
    # sum(abs(x_i - y_i)) / sum(x_i + y_i):
    # https://vegandevs.github.io/vegan/reference/diversity.html
    # https://vegandevs.github.io/vegan/reference/vegdist.html
    x = [1, 1, 2]
    p = [1 / 4, 1 / 4, 1 / 2]

    @test entropy(Shannon(; base=ℯ), x) ≈ -sum(pi -> pi * log(pi), p)
    @test gini_simpson_index(x) ≈ 1 - sum(abs2, p)
    @test inverse_simpson_index(x) ≈ inv(sum(abs2, p))
    @test pielou_evenness(x) ≈ entropy(Shannon(; base=ℯ), x) / log(3)
    @test fisher_alpha(x) * log1p(sum(x) / fisher_alpha(x)) ≈ 3
    @test bray_curtis_distance([1, 2, 3], [2, 2, 0]) ≈ 4 / 10
end
