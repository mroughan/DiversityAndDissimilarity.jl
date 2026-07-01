using DataFrames
using DiversityAndDissimilarity
using LinearAlgebra
using Test

include("docstring_examples.jl")

function test_valid_probabilities(data; frequencies=true)
    probabilities = proportions(data; frequencies)
    @test all(>=(0), probabilities)
    @test all(isfinite, probabilities)
    @test sum(probabilities) ≈ 1
    return probabilities
end

struct UnknownFrameworkIndex <: DiversityIndex end

@testset "abundance handling" begin
    @test counts(["a", "b", "a"]) == Dict("a" => 2, "b" => 1)
    @test proportions(Dict(:a => 2, :b => 2)) == [0.5, 0.5]
    @test proportions([1 1 2 0 5; 3 0 1 1 0]) ≈ [
                                                  1 / 9     1 / 9     2 / 9     0        5 / 9
                                                  3 / 5         0     1 / 5     1 / 5    0
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
    @test AddGamma().gamma == 1.0
    @test AddGamma(0).gamma == 0.0
    @test_throws ArgumentError AddGamma(-0.5)
    @test_throws ArgumentError entropy(Shannon(; estimator=AddGamma(1)), data; support=1)
    @test entropy(Shannon(; estimator=AddGamma(1)), Dict(:a => 2, :b => 1); support=[:a, :b, :c]) ≈
        -sum(p -> p * log2(p), [1 / 2, 1 / 3, 1 / 6])
    @test_throws ArgumentError entropy(Shannon(; estimator=AddGamma(1)), Dict(:a => 2); support=[:b])
    @test_throws ArgumentError entropy(Shannon(; estimator=ChaoShen()), [1, 1, 1])

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
    @test ACE().threshold == 10
    @test ACE(; threshold=3).threshold == 3
    @test_throws ArgumentError ACE(; threshold=0)
    @test_throws ArgumentError ACE(0)
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
    @test entropy(Shannon(; estimator=MillerMadow()), matrix) ≈
        [entropy(Shannon(; estimator=MillerMadow()), matrix[1, :]),
            entropy(Shannon(; estimator=MillerMadow()), matrix[2, :])]
    @test shannon(matrix) ≈ shannon_entropy(matrix)
    @test shannon_diversity(matrix) ≈ effective_diversity(Shannon(), matrix)
    @test renyi_entropy(matrix, 2) ≈ entropy(Renyi(2), matrix)
    @test renyi(matrix, 2) ≈ renyi_entropy(matrix, 2)
    @test renyi_diversity(matrix, 2) ≈ effective_diversity(Renyi(2), matrix)
    @test tsallis_entropy(matrix, 2) ≈ entropy(Tsallis(2), matrix)
    @test tsallis(matrix, 2) ≈ tsallis_entropy(matrix, 2)
    @test tsallis_diversity(matrix, 2) ≈ effective_diversity(Tsallis(2), matrix)
    @test effective_diversity(Simpson(), matrix) ≈ inverse_simpson_index(matrix)
    @test simpson_index(matrix) ≈ diversity(Simpson(), matrix)
    @test gini_simpson_index(matrix) ≈ diversity(GiniSimpson(), matrix)
    @test greenberg_diversity_index(matrix) ≈ diversity(GreenbergDiversityIndex(), matrix)
    @test linguistic_diversity_index(matrix) ≈ diversity(LinguisticDiversityIndex(), matrix)
    @test inverse_simpson_index(matrix) ≈ diversity(InverseSimpson(), matrix)
    @test length(shannon_confint(matrix)) == 2

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
    @test alpha_diversity(matrix; estimator=MillerMadow())[1].shannon_entropy ≈
        shannon_entropy(matrix[1, :]; estimator=MillerMadow())
    @test_throws ArgumentError alpha_diversity(matrix; threshold=0)

    validated = validate(matrix)
    @test validated.data == float.(matrix)
    @test richness(validated) == richness(matrix)
    @test alpha_diversity(validated) == alpha_diversity(matrix)
    @test alpha_diversity(validated; estimator=MillerMadow())[2].shannon_entropy ≈
        shannon_entropy(matrix[2, :]; estimator=MillerMadow())
    @test distance(BrayCurtis(), validated) ≈ distance(BrayCurtis(), matrix)
    @test similarity(Jaccard(), validated) ≈ similarity(Jaccard(), matrix)
    @test dissimilarity(Hellinger(), validated) ≈ dissimilarity(Hellinger(), matrix)
    @test distance(JensenShannon(), validated) ≈ distance(JensenShannon(), matrix)
    @test_throws ArgumentError validate(matrix; species=[:oak])
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
    wide_interval = entropy_confint(Shannon(), data; level=0.999999)
    @test wide_interval.lower < wide_interval.estimate < wide_interval.upper
    @test DiversityAndDissimilarity._normal_quantile(1e-6) < 0
    @test DiversityAndDissimilarity._normal_quantile(0.5) ≈ 0 atol = 1e-12
    @test_throws ArgumentError entropy_confint(Shannon(), data; level=1)
    @test_throws ArgumentError DiversityAndDissimilarity._normal_quantile(0)

    boot = bootstrap(Shannon(), data; nboot=25)
    @test boot.lower <= boot.upper
    @test length(boot.replicates) == 25
    @test isfinite(boot.stderr)

    boot_diversity = bootstrap(Shannon(), data; nboot=25, quantity=:diversity)
    @test boot_diversity.estimate ≈ shannon_diversity(data)
    @test_throws ArgumentError bootstrap(Shannon(), data; nboot=1)
    @test_throws ArgumentError bootstrap(Shannon(), data; nboot=3, quantity=:bad)
    @test_throws ArgumentError bootstrap(Shannon(), [1.5, 2.5]; nboot=3)

    jack = jackknife(Shannon(), data)
    @test jack.lower <= jack.upper
    @test length(jack.leave_one_out) == 4
    @test isfinite(jack.bias_corrected)
    jack_diversity = jackknife(Shannon(), data; quantity=:diversity)
    @test jack_diversity.estimate ≈ shannon_diversity(data)
    @test_throws ArgumentError jackknife(Shannon(), [1])
    @test_throws ArgumentError jackknife(Shannon(), data; quantity=:bad)
    @test_throws ArgumentError entropy_variance(Shannon(; estimator=ChaoShen()), [1, 1, 1])

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
    @test size(community_matrix(table)) == (2, 6)
    @test community_matrix(table)[:, 1] == [101.0, 102.0]
    @test proportions(table; species=[:oak, :ash, :elm, :pine, :birch]) ≈ proportions(matrix)
    @test validate(table; species=[:oak, :ash, :elm, :pine, :birch]).data == float.(matrix)
    @test_throws ArgumentError community_matrix([1 2; 3 4]; species=[:oak])
    @test richness(table; species=[:oak, :ash, :elm, :pine, :birch]) == richness(matrix)
    @test shannon_entropy(table; species=[:oak, :ash, :elm, :pine, :birch]) ≈ shannon_entropy(matrix)
    @test shannon_variance(table; species=[:oak, :ash, :elm, :pine, :birch]) ≈ shannon_variance(matrix)
    @test shannon_confint(table; species=[:oak, :ash, :elm, :pine, :birch])[1].estimate ≈
        shannon_confint(matrix)[1].estimate
    @test length(bootstrap(Shannon(), table; species=[:oak, :ash, :elm, :pine, :birch], nboot=3)) == 2
    @test length(jackknife(Shannon(), table; species=[:oak, :ash, :elm, :pine, :birch])) == 2
    @test alpha_diversity(table; species=[:oak, :ash, :elm, :pine, :birch]) == alpha_diversity(matrix)
    @test inverse_simpson_index(table; species=[:oak, :ash, :elm, :pine, :birch]) ≈ inverse_simpson_index(matrix)
    @test chao1(table; species=[:oak, :ash, :elm, :pine, :birch]) ≈ chao1(matrix)
    @test pielou_evenness(table; species=[:oak, :ash, :elm, :pine, :birch]) ≈ pielou_evenness(matrix)
    @test bray_curtis_distance(table; species=[:oak, :ash, :elm, :pine, :birch]) ≈ bray_curtis_distance(matrix)
    @test distance(Jaccard(), table; species=[:oak, :ash, :elm, :pine, :birch]) ≈ distance(Jaccard(), matrix)
    labeled = labeled_distance(BrayCurtis(), table; label=:site, species=[:oak, :ash, :elm, :pine, :birch])
    @test labeled.labels == ["a", "b"]
    @test labeled.matrix ≈ bray_curtis_distance(matrix)
    @test labeled_dissimilarity(BrayCurtis(), table; label=:site, species=[:oak, :ash, :elm, :pine, :birch]).matrix ≈
        bray_curtis_distance(matrix)
    @test labeled_similarity(Jaccard(), table; label=:site, species=[:oak, :ash, :elm, :pine, :birch]).matrix ≈
        jaccard_similarity(matrix)
    @test labeled_similarity(Jaccard(), matrix; labels=["a", "b"]).labels == ["a", "b"]
    @test labeled_distance(BrayCurtis(), matrix).labels == [1, 2]
    @test_throws ArgumentError labeled_distance(BrayCurtis(), table; labels=["a", "b"], label=:site,
        species=[:oak, :ash, :elm, :pine, :birch])
    @test_throws ArgumentError labeled_distance(BrayCurtis(), matrix; label=:site)
    @test_throws ArgumentError labeled_distance(BrayCurtis(), matrix; labels=["a"])
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
    @test is_symmetric(Shannon())
    @test is_symmetric(BrayCurtis())
    @test !is_symmetric(KullbackLeibler())
    @test index_metadata(KullbackLeibler()).is_symmetric == false
    @test output_mode(Shannon()) == :entropy
    @test output_mode(BrayCurtis()) == :dissimilarity
    @test output_mode(Canberra()) == :dissimilarity
    @test output_mode(Hellinger()) == :distance
    @test index_range(PielouEvenness()) == (lower=0.0, upper=1.0)
    @test index_bounds(Jaccard()).lower_meaning ==
        "minimal similarity; conventionally complete dissimilarity or no overlap"
    @test index_bounds(BrayCurtis()).upper == 1.0
    @test !is_finite(KullbackLeibler())
    @test is_finite(JensenShannon())
    @test is_nonnegative(KullbackLeibler())
    @test is_bounded(JensenDifference())
    @test !is_bounded(KullbackLeibler())
    @test is_triangular(Jaccard())
    @test !is_triangular(KullbackLeibler())
    @test is_pseudometric(ShannonDifference())
    @test is_quasimetric(Jaccard())
    @test is_metametric(Overlap())
    @test is_semimetric(BrayCurtis())
    @test is_premetric(KullbackLeibler())
    @test !is_supermetric(JensenShannon())
    @test is_similarity(Jaccard())
    @test !is_similarity(BrayCurtis())
    @test is_dissimilarity(BrayCurtis())
    @test is_dissimiliarty(BrayCurtis())
    @test index_metadata(BrayCurtis()).is_semimetric
    @test index_metadata(Jaccard()).bounds.upper_meaning ==
        "maximal similarity; conventionally identical or complete overlap"
    @test index_bounds(Shannon()).lower_meaning ==
        "no uncertainty; all mass in one category"
    @test index_bounds(Shannon()).upper_meaning ==
        "maximum uncertainty for the supplied or observed support"
    @test index_bounds(Richness()).lower_meaning ==
        "minimal diversity under the index convention"
    @test index_bounds(Richness()).upper_meaning ==
        "unbounded effective or richness-scale diversity"
    @test index_bounds(GiniSimpson()).upper_meaning ==
        "maximal diversity under the index convention"
    @test index_bounds(Simpson()).lower_meaning == "minimal dominance"
    @test index_bounds(Simpson()).upper_meaning == "maximal dominance"
    @test index_bounds(SampleCoverage()).lower_meaning ==
        "no sampled probability mass covered"
    @test index_bounds(SampleCoverage()).upper_meaning ==
        "complete sampled probability mass covered"
    @test index_bounds(PielouEvenness()).lower_meaning == "minimal evenness"
    @test index_bounds(PielouEvenness()).upper_meaning == "maximal evenness"
    @test index_bounds(Chao1()).lower_meaning == :unknown
    @test index_bounds(Chao1()).upper_meaning == :unknown
    @test is_triangular(BrayCurtis()) == false
    @test is_supermetric(Richness()) == :unknown
    @test all(result -> result.passed, validate_reference_cases())

    catalogue = DiversityIndex[
        Richness(), Shannon(), Renyi(2), Tsallis(2), Simpson(), GiniSimpson(),
        GreenbergDiversityIndex(), LinguisticDiversityIndex(), InverseSimpson(),
        Hill(2), Chao1(), ACE(), SampleCoverage(), PielouEvenness(), FisherAlpha(),
        Jaccard(), SorensenDice(), Overlap(), BrayCurtis(), Ruzicka(), TotalVariation(),
        Manhattan(), Euclidean(), Canberra(), Hellinger(), Chord(), Bhattacharyya(),
        KullbackLeibler(), ShannonDifference(), JensenDifference(), JensenShannon(),
        MorisitaHorn(),
    ]
    for index in catalogue
        metadata = index_metadata(index)
        @test metadata.type == typeof(index)
        @test metadata.family == index_family(index)
        @test metadata.input_mode == input_mode(index)
        @test metadata.output_mode == output_mode(index)
        @test metadata.range == index_range(index)
        @test metadata.bounds == index_bounds(index)
        @test metadata.is_finite == is_finite(index)
        @test metadata.is_metric == is_metric(index)
        @test metadata.is_triangular == is_triangular(index)
        @test metadata.is_nonnegative == is_nonnegative(index)
        @test metadata.is_bounded == is_bounded(index)
        @test metadata.is_pseudometric == is_pseudometric(index)
        @test metadata.is_quasimetric == is_quasimetric(index)
        @test metadata.is_metametric == is_metametric(index)
        @test metadata.is_semimetric == is_semimetric(index)
        @test metadata.is_premetric == is_premetric(index)
        @test metadata.is_supermetric == is_supermetric(index)
        @test metadata.is_similarity == is_similarity(index)
        @test metadata.is_dissimilarity == is_dissimilarity(index)
        @test metadata.is_symmetric == is_symmetric(index)
        @test metadata.requires_probabilities == requires_probabilities(index)
        @test metadata.supports_matrix_kernel == supports_matrix_kernel(index)
        @test metadata.formula isa AbstractString
        @test metadata.aliases isa Vector{String}
        @test metadata.notes isa AbstractString
    end

    expected_families = (
        Richness() => :richness,
        Shannon() => :entropy,
        Simpson() => :dominance,
        GreenbergDiversityIndex() => :linguistic_diversity,
        InverseSimpson() => :effective_diversity,
        Chao1() => :richness_estimator,
        SampleCoverage() => :coverage,
        PielouEvenness() => :evenness,
        FisherAlpha() => :diversity,
        Jaccard() => :incidence,
        BrayCurtis() => :abundance,
        KullbackLeibler() => :probability,
    )
    for (index, family) in expected_families
        @test index_family(index) == family
    end
    @test input_mode(Shannon()) == :single_assemblage
    @test input_mode(Jaccard()) == :pairwise
    @test output_mode(Bhattacharyya()) == :coefficient
    @test index_range(Manhattan()) == (lower=0.0, upper=2.0)
    @test index_range(Euclidean()) == (lower=0.0, upper=sqrt(2))
    @test index_range(JensenDifference(; base=4)).upper ≈ 0.5
    @test index_range(JensenShannon(; base=4)).upper ≈ sqrt(0.5)
    @test index_range(JensenShannon(; base=4, distance=false)).upper ≈ 0.5
    @test requires_probabilities(Shannon())
    @test !requires_probabilities(Richness())
    @test supports_matrix_kernel(Shannon())
    @test !supports_matrix_kernel(Canberra())
    @test contains(index_metadata(KullbackLeibler()).notes, "Asymmetric")
    @test "vegan: vegdist(method=\"bray\")" in index_metadata(BrayCurtis()).aliases
    @test !isempty(index_metadata(SorensenDice()).formula)

    unknown = UnknownFrameworkIndex()
    unknown_metadata = index_metadata(unknown)
    @test unknown_metadata.family == :unknown
    @test unknown_metadata.input_mode == :single_assemblage
    @test unknown_metadata.output_mode == :estimate
    @test unknown_metadata.range == (lower=0.0, upper=Inf)
    @test unknown_metadata.bounds.lower_meaning == :unknown
    @test unknown_metadata.bounds.upper_meaning == :unknown
    @test unknown_metadata.is_finite
    @test !unknown_metadata.is_metric
    @test unknown_metadata.is_triangular == :unknown
    @test unknown_metadata.is_pseudometric == :unknown
    @test unknown_metadata.is_quasimetric == :unknown
    @test unknown_metadata.is_metametric == :unknown
    @test unknown_metadata.is_semimetric == :unknown
    @test unknown_metadata.is_premetric == :unknown
    @test unknown_metadata.is_supermetric == :unknown
    @test !unknown_metadata.is_similarity
    @test !unknown_metadata.is_dissimilarity
    @test unknown_metadata.is_symmetric
    @test !unknown_metadata.requires_probabilities
    @test !unknown_metadata.supports_matrix_kernel
    @test isempty(unknown_metadata.formula)
    @test isempty(unknown_metadata.aliases)
    @test isempty(unknown_metadata.notes)

    report = estimator_report([1, 1, 2, 0, 5]; support=6)
    @test report.observed_richness == 4
    @test report.singletons == 2
    @test any(estimate -> estimate.name == :add_one, report.estimates)

    table_report = estimator_report(
        DataFrame(site=["a", "b"], oak=[1, 3], ash=[1, 0], elm=[2, 1]);
        species=[:oak, :ash, :elm],
    )
    @test length(table_report) == 2
    @test table_report[1].observed_richness == 3
    @test table_report[2].observed_richness == 2

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
    table_audit = diversity_audit(
        DataFrame(site=["a", "b"], oak=[1, 3], ash=[1, 0], elm=[2, 1]);
        label=:site,
        species=[:oak, :ash, :elm],
        pairwise_index=Jaccard(),
    )
    @test table_audit.pairwise.labels == ["a", "b"]
    @test table_audit.pairwise.matrix ≈ jaccard_distance([1 1 2; 3 0 1])

    uncertainty = uncertainty_audit([1 1 2 0 5; 3 0 1 1 0]; labels=["a", "b"], nboot=10)
    @test uncertainty.labels == ["a", "b"]
    @test length(uncertainty.reports) == 2
    @test uncertainty.reports[1].label == "a"
    @test any(result -> result.quantity == :entropy, uncertainty.reports[1].estimates)
    @test any(contains("low nboot"), uncertainty.warnings)
    table_uncertainty = uncertainty_audit(
        DataFrame(site=["a", "b"], oak=[1, 3], ash=[1, 0], elm=[2, 1]);
        label=:site,
        species=[:oak, :ash, :elm],
        quantities=(:entropy,),
        nboot=3,
    )
    @test table_uncertainty.labels == ["a", "b"]
    @test length(table_uncertainty.reports[1].estimates) == 1
end

@testset "pairwise indexes" begin
    left = Dict(:a => 2, :b => 1)
    right = Dict(:b => 3, :c => 1)
    left_vector = [1, 2, 0]
    right_vector = [0, 1, 3]
    p = [1 / 3, 2 / 3, 0]
    q = [0, 1 / 4, 3 / 4]
    hellinger_expected = sqrt(sum(abs2, sqrt.(p) .- sqrt.(q)) / 2)
    shannon_p = -sum(pi > 0 ? pi * log2(pi) : 0 for pi in p)
    shannon_q = -sum(qi > 0 ? qi * log2(qi) : 0 for qi in q)
    shannon_difference_expected = abs(shannon_p - shannon_q)
    kl_left_right_expected = sum(pi > 0 ? pi * log2(pi / qi) : 0 for (pi, qi) in zip(p, q))
    kl_right_left_expected = Inf
    finite_kl_left = [1 / 2, 1 / 2]
    finite_kl_right = [3 / 4, 1 / 4]
    finite_kl_expected = sum(pi * log2(pi / qi) for (pi, qi) in zip(finite_kl_left, finite_kl_right))
    jensen_shannon_expected =
        (
            sum(pi > 0 ? pi * log2(pi / ((pi + qi) / 2)) : 0 for (pi, qi) in zip(p, q)) +
            sum(qi > 0 ? qi * log2(qi / ((pi + qi) / 2)) : 0 for (pi, qi) in zip(p, q))
        ) / 2
    pairwise_matrix = [
        1 2 0
        0 1 3
        2 0 1
    ]

    @test KullbackLeibler(; base=10).base == 10.0
    @test KullbackLeibler(; estimator=AddGamma(1), support=4).support == 4
    @test ShannonDifference(; base=10).base == 10.0
    @test JensenDifference(; base=10).base == 10.0
    @test JensenShannon(; distance=false).distance == false
    for constructor in (KullbackLeibler, ShannonDifference, JensenDifference, JensenShannon)
        @test_throws ArgumentError constructor(; base=1)
        @test_throws ArgumentError constructor(; base=0)
    end

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
    @test kullback_leibler_divergence(left_vector, right_vector) ≈ kl_left_right_expected
    @test dissimilarity(KullbackLeibler(), left_vector, right_vector) ≈ kl_left_right_expected
    @test kullback_leibler_divergence(right_vector, left_vector) == kl_right_left_expected
    @test kullback_leibler_divergence(finite_kl_left, finite_kl_right) ≈ finite_kl_expected
    @test kullback_leibler_divergence(finite_kl_left, finite_kl_right) !=
        kullback_leibler_divergence(finite_kl_right, finite_kl_left)
    @test kullback_leibler_divergence([1, 1], [1, 1]) ≈ 0
    kl_matrix = kullback_leibler_divergence([1 1; 3 1])
    @test kl_matrix[1, 2] ≈ finite_kl_expected
    @test kl_matrix[2, 1] ≈ kullback_leibler_divergence(finite_kl_right, finite_kl_left)
    @test kl_matrix[1, 2] != kl_matrix[2, 1]
    @test shannon_difference(left_vector, right_vector) ≈ shannon_difference_expected
    @test dissimilarity(ShannonDifference(), left_vector, right_vector) ≈ shannon_difference_expected
    @test similarity(ShannonDifference(), left_vector, right_vector) ≈
        1 - shannon_difference_expected / log2(3)
    @test jensen_difference(left_vector, right_vector) ≈ jensen_shannon_expected
    @test dissimilarity(JensenDifference(), left_vector, right_vector) ≈ jensen_shannon_expected
    @test similarity(JensenDifference(), left_vector, right_vector) ≈
        1 - jensen_shannon_expected
    @test jensen_shannon_divergence(left_vector, right_vector) ≈ jensen_shannon_expected
    @test jensen_shannon_distance(left_vector, right_vector) ≈ sqrt(jensen_shannon_expected)
    @test jensen_shannon_similarity(left_vector, right_vector) ≈
        1 - sqrt(jensen_shannon_expected)
    @test dissimilarity(JensenShannon(; distance=false), left_vector, right_vector) ≈ jensen_shannon_expected
    @test dissimilarity(JensenDifference(), left_vector, right_vector) ≈
        dissimilarity(JensenShannon(; distance=false), left_vector, right_vector)
    @test dissimilarity(KullbackLeibler(), [1, 0], [0, 1]) == Inf
    @test dissimilarity(KullbackLeibler(), [0, 1], [1, 0]) == Inf
    @test dissimilarity(Bhattacharyya(), [1, 0], [0, 1]) == Inf
    @test similarity(ShannonDifference(), [5], [10]) ≈ 1
    @test kullback_leibler_divergence(left_vector, right_vector; estimator=MillerMadow()) <=
        kullback_leibler_divergence(left_vector, right_vector)
    @test isfinite(kullback_leibler_divergence(left_vector, right_vector; estimator=AddGamma(1)))
    @test isfinite(kullback_leibler_divergence(left_vector, right_vector; estimator=AddGamma(0.5)))
    @test isfinite(kullback_leibler_divergence(left_vector, right_vector; estimator=HausserStrimmer()))
    @test isfinite(kullback_leibler_divergence(left_vector, right_vector; estimator=ChaoShen()))
    @test isfinite(jensen_shannon_divergence(left_vector, right_vector; estimator=MillerMadow()))
    @test isfinite(jensen_shannon_divergence(left_vector, right_vector; estimator=AddGamma(1), support=4))
    @test isfinite(jensen_shannon_divergence(left_vector, right_vector; estimator=AddGamma(0.5), support=4))
    @test isfinite(jensen_shannon_divergence(left_vector, right_vector; estimator=HausserStrimmer(), support=4))
    @test isfinite(jensen_shannon_divergence(left_vector, right_vector; estimator=ChaoShen()))
    support = [:a, :b, :c, :d]
    @test isfinite(kullback_leibler_divergence(left, right; estimator=AddGamma(1), support))
    @test isfinite(kullback_leibler_divergence(left, right; estimator=ChaoShen(), support))
    @test isfinite(jensen_difference(left, right; estimator=AddGamma(0.5), support))
    @test isfinite(jensen_shannon_divergence(left, right; estimator=ChaoShen(), support))
    @test_throws ArgumentError kullback_leibler_divergence(left, right; estimator=AddGamma(1), support=[:a, :b])
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
    @test distance(Euclidean(), left_vector, right_vector) ≈ euclidean_distance(left_vector, right_vector)
    @test jaccard_index(pairwise_matrix) ≈ similarity(Jaccard(), pairwise_matrix)
    @test jaccard_similarity(pairwise_matrix) ≈ jaccard_index(pairwise_matrix)
    @test jaccard_distance(pairwise_matrix) ≈ dissimilarity(Jaccard(), pairwise_matrix)
    @test sorensen_dice_index(pairwise_matrix) ≈ similarity(SorensenDice(), pairwise_matrix)
    @test sorensen_index(pairwise_matrix) ≈ sorensen_dice_index(pairwise_matrix)
    @test sorensen_dice_dissimilarity(pairwise_matrix) ≈ dissimilarity(SorensenDice(), pairwise_matrix)
    @test sorensen_dice_distance(pairwise_matrix) ≈ sorensen_dice_dissimilarity(pairwise_matrix)
    @test sorensen_distance(pairwise_matrix) ≈ sorensen_dice_distance(pairwise_matrix)
    @test bray_curtis_dissimilarity(pairwise_matrix) ≈ dissimilarity(BrayCurtis(), pairwise_matrix)
    @test bray_curtis_distance(pairwise_matrix) ≈ bray_curtis_dissimilarity(pairwise_matrix)
    @test overlap_similarity(pairwise_matrix) ≈ similarity(Overlap(), pairwise_matrix)
    @test overlap_distance(pairwise_matrix) ≈ dissimilarity(Overlap(), pairwise_matrix)
    @test ruzicka_similarity(pairwise_matrix) ≈ similarity(Ruzicka(), pairwise_matrix)
    @test quantitative_jaccard_similarity(pairwise_matrix) ≈ ruzicka_similarity(pairwise_matrix)
    @test ruzicka_distance(pairwise_matrix) ≈ dissimilarity(Ruzicka(), pairwise_matrix)
    @test quantitative_jaccard_distance(pairwise_matrix) ≈ ruzicka_distance(pairwise_matrix)
    @test total_variation_distance(pairwise_matrix) ≈ dissimilarity(TotalVariation(), pairwise_matrix)
    @test manhattan_distance(pairwise_matrix) ≈ dissimilarity(Manhattan(), pairwise_matrix)
    @test euclidean_distance(pairwise_matrix) ≈ dissimilarity(Euclidean(), pairwise_matrix)
    @test canberra_distance(pairwise_matrix) ≈ dissimilarity(Canberra(), pairwise_matrix)
    @test hellinger_distance(pairwise_matrix) ≈ dissimilarity(Hellinger(), pairwise_matrix)
    @test chord_distance(pairwise_matrix) ≈ dissimilarity(Chord(), pairwise_matrix)
    @test bhattacharyya_coefficient(pairwise_matrix) ≈ similarity(Bhattacharyya(), pairwise_matrix)
    @test bhattacharyya_distance(pairwise_matrix) ≈ dissimilarity(Bhattacharyya(), pairwise_matrix)
    @test shannon_difference(pairwise_matrix) ≈ dissimilarity(ShannonDifference(), pairwise_matrix)
    @test jensen_difference(pairwise_matrix) ≈ dissimilarity(JensenDifference(), pairwise_matrix)
    @test jensen_shannon_similarity(pairwise_matrix) ≈ similarity(JensenShannon(), pairwise_matrix)
    @test jensen_shannon_divergence(pairwise_matrix) ≈ dissimilarity(JensenShannon(; distance=false), pairwise_matrix)
    @test jensen_shannon_distance(pairwise_matrix) ≈ dissimilarity(JensenShannon(), pairwise_matrix)
    @test morisita_horn_similarity(pairwise_matrix) ≈ similarity(MorisitaHorn(), pairwise_matrix)
    @test morisita_horn_distance(pairwise_matrix) ≈ dissimilarity(MorisitaHorn(), pairwise_matrix)
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

@testset "reference values from Magurran (2004)" begin
    # Five-species community used as a recurring worked example in:
    # Magurran, A.E. (2004) "Measuring Biological Diversity", Blackwell, Table 4.1.
    # Independent verification via R vegan:
    #   x <- c(135, 76, 45, 22, 13)
    #   vegan::diversity(x, "shannon")    # 1.329752 (nats)
    #   vegan::diversity(x, "simpson")    # 0.684947 (= GiniSimpson)
    #   vegan::diversity(x, "invsimpson") # 3.174054
    x = [135, 76, 45, 22, 13]

    @test richness(x) == 5

    # Simpson concentration D = Σp_i² = (135²+76²+45²+22²+13²)/291²; exact fraction.
    D_num = 135^2 + 76^2 + 45^2 + 22^2 + 13^2  # 26679
    D_den = sum(x)^2                             # 84681
    @test simpson_index(x) ≈ D_num / D_den
    @test gini_simpson_index(x) ≈ 1 - D_num / D_den
    @test inverse_simpson_index(x) ≈ D_den / D_num

    # Shannon entropy in nats; verifiable via vegan::diversity(c(135,76,45,22,13),"shannon").
    @test entropy(Shannon(; base=ℯ), x) ≈ 1.3296981241766592 atol = 1e-10

    # Pielou evenness J = H_nats / ln(S).
    @test pielou_evenness(x) ≈ entropy(Shannon(; base=ℯ), x) / log(5)

    # Fisher's alpha satisfies the implicit equation S = α log(1 + n/α).
    α = fisher_alpha(x)
    @test α * log1p(sum(x) / α) ≈ 5
end

@testset "reference values from iNEXT spider dataset" begin
    # Spider abundance data from two hemlock-forest canopy treatments at Harvard Forest.
    # Original study: Sackett et al. (2011) Can. J. Forest Res. 41(2):394-409.
    #   doi:10.1139/X10-207
    # Distributed with the iNEXT R package:
    #   Chao et al. (2014) Ecol. Monogr. 84(1):45-67. doi:10.1890/13-0133.1
    #   Hsieh, Ma & Chao (2016) Methods Ecol. Evol. 7(12):1451-1456.
    #
    # Vectors extracted from iNEXT package source:
    #   https://github.com/JohnsonHsieh/iNEXT
    #
    # To verify in R:
    #   library(iNEXT); data(spider)
    #   library(vegan)
    #   vegan::diversity(spider$Girdled, "shannon")  # nats
    #   estimateR(spider$Girdled)                    # Chao1, ACE
    #
    # See also: notes/references/spider_dataset.md
    girdled = [46, 22, 17, 15, 15, 9, 8, 6, 6, 4, 2, 2, 2, 2,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    logged  = [88, 22, 16, 15, 13, 10, 8, 8, 7, 7, 7, 5, 4, 4, 4,
               3, 3, 3, 3, 2, 2, 2, 2,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    # --- Girdled site: 26 species, 168 individuals ---
    @test sum(girdled) == 168
    @test richness(girdled) == 26

    # Exact Chao1 (bias-corrected): 26 + f₁(f₁-1)/(2(f₂+1)) = 26 + 12·11/(2·5)
    # f₁=12 singletons, f₂=4 doubletons; result 39.2 is exact.
    @test chao1(girdled) ≈ 39.2

    # Exact sample coverage: 1 − f₁/n = 1 − 12/168 = 13/14
    @test sample_coverage(girdled) ≈ 13 / 14

    # Shannon entropy (nats); verified via vegan::diversity(spider$Girdled, "shannon").
    @test entropy(Shannon(; base=ℯ), girdled) ≈ 2.4898652 atol = 1e-5

    # Hill(1) = exp(H_nats) is equivalent to Shannon effective diversity.
    @test hill_number(girdled, 1) ≈ exp(entropy(Shannon(; base=ℯ), girdled))

    # Hill(2) = InverseSimpson; exact from integer counts.
    # D = Σnᵢ²/n² = (46²+22²+…)/168²; InverseSimpson = 1/D = 168²/Σnᵢ²
    D_num_g = sum(v^2 for v in girdled)
    @test hill_number(girdled, 2) ≈ 168^2 / D_num_g

    # ACE from package (Chao & Lee 1992); verifiable via estimateR() in vegan.
    @test ace(girdled) ≈ 48.41522903033908 atol = 1e-6

    # Pielou evenness J = H_nats / ln(S).
    @test pielou_evenness(girdled) ≈ entropy(Shannon(; base=ℯ), girdled) / log(26)

    # --- Logged site: 37 species, 252 individuals ---
    @test sum(logged) == 252
    @test richness(logged) == 37

    # Exact Chao1: 37 + 14·13/(2·5) = 37 + 182/10 = 55.2
    @test chao1(logged) ≈ 55.2

    # Exact sample coverage: 1 − 14/252 = 17/18
    @test sample_coverage(logged) ≈ 17 / 18

    # Shannon entropy (nats); verified via vegan::diversity(spider$Logged, "shannon").
    @test entropy(Shannon(; base=ℯ), logged) ≈ 2.6686855819621367 atol = 1e-5

    # Hill(2) = InverseSimpson; exact from integer counts.
    D_num_l = sum(v^2 for v in logged)
    @test hill_number(logged, 2) ≈ 252^2 / D_num_l

    # ACE; verifiable via estimateR() in vegan.
    @test ace(logged) ≈ 52.68499427262313 atol = 1e-6
end

@testset "SciPy convention comparison" begin
    # Documents where this package agrees with scipy.spatial.distance and where
    # conventions differ. See notes/references/scipy_conventions.md for details.
    # SciPy reference: https://docs.scipy.org/doc/scipy/reference/spatial.distance.html
    # Virtanen et al. (2020) SciPy 1.0. Nature Methods 17:261-272.

    # --- Bray-Curtis: identical formula ---
    # Both packages use sum(|u-v|) / sum(u+v).
    @test bray_curtis_dissimilarity([1, 0, 0], [0, 1, 0]) ≈ 1.0
    @test bray_curtis_dissimilarity([1, 1, 0], [0, 1, 0]) ≈ 1 / 3

    # --- Canberra: this package averages by nonzero terms; SciPy does NOT ---
    # For [1,2,3] vs [2,1,0]:
    #   Each term |u-v|/(|u|+|v|): 1/3 + 1/3 + 3/3 = 5/3
    #   SciPy (unaveraged): 5/3 ≈ 1.6667
    #   This package (averaged over m=3 nonzero terms): (5/3)/3 = 5/9 ≈ 0.5556
    # Convention follows Legendre & Legendre (2012) "Numerical Ecology" §7.4.
    @test canberra_distance([1, 2, 3], [2, 1, 0]) ≈ 5 / 9
    # Relationship to unaveraged SciPy form: ours * m = scipy (where m = nonzero terms)
    m = count(>( 0), [1, 2, 3] .+ [2, 1, 0])  # 3 nonzero pairs
    scipy_canberra = 5 / 3
    @test canberra_distance([1, 2, 3], [2, 1, 0]) * m ≈ scipy_canberra

    # --- Jensen-Shannon distance: base-2 (bits) vs natural log (nats) ---
    # SciPy uses natural log and returns sqrt(JSD_nats).
    # This package defaults to base=2 and returns sqrt(JSD_bits).
    # Conversion: ours = scipy × sqrt(log₂ e) ≈ scipy × 1.2011.
    # Numerically: sqrt(JSD_bits) = sqrt(JSD_nats / ln 2) = sqrt(JSD_nats) × 1/sqrt(ln 2).
    #
    # SciPy reference values (from scipy.spatial.distance.jensenshannon docs):
    #   [1.0, 0.0, 0.0] vs [0.0, 1.0, 0.0] → 0.83255 (scipy), 1.0 (ours, base=2)
    #   [1.0, 0.0]      vs [0.5, 0.5]       → 0.46450 (scipy), 0.55792 (ours, base=2)
    #
    scipy_js_disjoint = sqrt(log(2))           # ≈ 0.8326
    scipy_js_half     = sqrt(0.31128 / log(2)) * sqrt(log(2))  # same as sqrt(0.2158)

    # Our base=2 result for disjoint distributions (max JSD = 1 bit → distance = 1).
    @test jensen_shannon_distance([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]) ≈ 1.0

    # Using base=ℯ matches SciPy exactly.
    @test jensen_shannon_distance([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]; base=ℯ) ≈
        sqrt(log(2)) atol = 1e-10
    @test jensen_shannon_distance([1.0, 0.0], [0.5, 0.5]; base=ℯ) ≈ 0.4645 atol = 1e-4

    # The ratio between our (base=2) and scipy (nats) values equals sqrt(log₂ e).
    our_val   = jensen_shannon_distance([1.0, 0.0], [0.5, 0.5])
    scipy_val = jensen_shannon_distance([1.0, 0.0], [0.5, 0.5]; base=ℯ)
    @test our_val / scipy_val ≈ sqrt(log2(ℯ)) atol = 1e-10
end

@testset "mathematical identities — entropy and diversity" begin
    # --- Shannon entropy (Cover & Thomas 2006, "Elements of Information Theory") ---

    # Uniform distribution achieves maximum entropy log₂(n).
    # H(p₁,...,pₙ) ≤ log₂(n) with equality iff pᵢ = 1/n (Cover & Thomas §2.6).
    for n in [2, 3, 5, 10, 50]
        @test shannon_entropy(ones(n)) ≈ log2(n)
    end

    # Single non-zero mass: H = 0 (no uncertainty).
    @test shannon_entropy([1, 0, 0, 0]) ≈ 0.0
    @test shannon_entropy([7]) ≈ 0.0

    # Adding zero-probability species is a no-op (continuity axiom, Shannon 1948).
    x = [1, 3, 2]
    @test shannon_entropy(x) ≈ shannon_entropy([x; 0])
    @test shannon_entropy(x) ≈ shannon_entropy([0; x; 0])
    @test richness(x) == richness([x; 0]) == richness([0; x; 0])

    # Base change: H_b = H_e / ln(b).
    @test shannon_entropy(x; base=2) ≈ shannon_entropy(x; base=ℯ) / log(2)
    @test shannon_entropy(x; base=10) ≈ shannon_entropy(x; base=ℯ) / log(10)

    # MillerMadow corrects upward: H_MM ≥ H_Plugin (correction is always +).
    # Miller (1955) "Note on the bias of information estimates".
    for x in [[1, 1], [1, 2, 3], [4, 3, 2, 1]]
        @test shannon_entropy(x; estimator=MillerMadow()) ≥ shannon_entropy(x)
    end

    # --- Simpson family identities ---

    # Complementarity: D + (1-D) = 1 and D × (1/D) = 1 for any distribution.
    for x in [[1, 2, 3], [10, 10], [5, 3, 1, 1], [1, 1, 1, 1, 1]]
        @test simpson_index(x) + gini_simpson_index(x) ≈ 1.0
        @test simpson_index(x) * inverse_simpson_index(x) ≈ 1.0
    end

    # Uniform distribution: GiniSimpson = 1 - 1/n (maximum for given richness).
    for n in [2, 4, 10]
        @test gini_simpson_index(ones(n)) ≈ 1 - 1 / n
    end

    # --- Renyi entropy at special orders (Rényi 1961, Hill 1973) ---
    x = [1, 2, 3, 4]

    # q = 0: H₀ = log₂(S) since Σpᵢ⁰ = S (number of non-zero categories).
    @test renyi_entropy(x, 0) ≈ log2(richness(x))

    # q = 1 (limit): reduces to Shannon entropy.
    @test renyi_entropy(x, 1) ≈ shannon_entropy(x)

    # q = 2: H₂ = −log₂(Σpᵢ²) = −log₂(D_Simpson).
    @test renyi_entropy(x, 2) ≈ -log2(simpson_index(x))

    # Uniform distribution: Rényi(q) = log₂(n) for all q.
    for q in [0.0, 0.5, 1.0, 2.0, 5.0]
        @test renyi_entropy(ones(6), q) ≈ log2(6)
    end

    # Rényi is non-increasing in q for non-uniform distributions (Hill 1973).
    @test renyi_entropy(x, 0) ≥ renyi_entropy(x, 1) - 1e-10
    @test renyi_entropy(x, 1) ≥ renyi_entropy(x, 2) - 1e-10
    @test renyi_entropy(x, 2) ≥ renyi_entropy(x, 5) - 1e-10

    # --- Tsallis entropy at special orders ---
    # Tsallis (1988) "Possible generalization of Boltzmann-Gibbs statistics".

    # q = 1: reduces to Shannon (limit definition, scaled to same base).
    @test tsallis_entropy(x, 1) ≈ shannon_entropy(x)

    # q = 2 with base b: T₂ = GiniSimpson / log(b).
    # From T_q = (Σpᵢ^q − 1)/((1−q)log b): at q=2, T₂ = (D−1)/(−log b) = (1−D)/log b.
    for b in [2.0, ℯ, 10.0]
        @test tsallis_entropy(x, 2; base=b) ≈ gini_simpson_index(x) / log(b)
    end

    # --- Hill numbers (Hill 1973, "Diversity and evenness") ---

    # Hill(0) = Richness (counts all species equally).
    @test hill_number(x, 0) == richness(x)

    # Hill(1) = exp(Shannon in nats) = Shannon effective diversity.
    @test hill_number(x, 1) ≈ exp(shannon_entropy(x; base=ℯ))

    # Hill(2) = Inverse Simpson.
    @test hill_number(x, 2) ≈ inverse_simpson_index(x)

    # Hill–Rényi correspondence: Hill(q) = base^Renyi_base(q) for any base.
    for q in [0.5, 1.5, 2.0, 3.0]
        @test hill_number(x, q) ≈ 2.0^renyi_entropy(x, q)
        @test hill_number(x, q) ≈ ℯ^renyi_entropy(x, q; base=ℯ)
    end

    # Uniform distribution: Hill(q) = n for all q.
    for n in [3, 5, 10]
        for q in [0.0, 0.5, 1.0, 2.0, 5.0]
            @test hill_number(ones(n), q) ≈ float(n)
        end
    end

    # Hill numbers are non-increasing in q for non-uniform distributions.
    @test hill_number(x, 0) ≥ hill_number(x, 1) - 1e-10
    @test hill_number(x, 1) ≥ hill_number(x, 2) - 1e-10

    # --- Richness estimator bounds ---

    # Chao1 = Richness when no singletons (f₁ = 0); see Chao (1984) eq. 2.
    @test chao1([2, 4, 3, 6]) == float(richness([2, 4, 3, 6]))

    # Chao1 ≥ Richness and ACE ≥ Richness (both are upward-biased estimators).
    for x in [[1, 1, 2, 3], [4, 3, 2, 1], [2, 2, 2], [1, 5, 3]]
        @test chao1(x) ≥ richness(x) - 1e-10
        @test ace(x) ≥ richness(x) - 1e-10
    end

    # SampleCoverage = 1 when no singletons (all species observed at least twice).
    @test sample_coverage([2, 3, 4]) ≈ 1.0
    @test sample_coverage([2, 2, 2, 2]) ≈ 1.0

    # PielouEvenness = 1 for uniform distributions.
    for n in [2, 5, 10]
        @test pielou_evenness(ones(n)) ≈ 1.0
    end
end

@testset "mathematical identities — pairwise indices" begin
    # --- Incidence index relationships ---

    # Ruzicka = Jaccard on binary (0/1) data.
    # When data are binary, min(x,y)=x∧y and max(x,y)=x∨y, so Ruzicka reduces
    # to |A∩B|/|A∪B| = Jaccard (Ruzicka 1958 is the abundance generalisation).
    for (a, b) in [([1, 1, 0, 1], [1, 0, 1, 1]), ([1, 0, 1], [0, 1, 1])]
        @test similarity(Ruzicka(), a, b) ≈ similarity(Jaccard(), a, b)
    end

    # Jaccard–Sørensen conversion: S = 2J/(1+J) and J = S/(2−S).
    # Follows from their definitions with |A|+|B| = 2|A∩B| + |A△B|.
    for (a, b) in [([1, 1, 0, 1], [1, 0, 1, 1]), ([1, 0, 0, 1], [0, 1, 0, 1])]
        J = similarity(Jaccard(), a, b)
        S = similarity(SorensenDice(), a, b)
        @test S ≈ 2J / (1 + J)
        @test J ≈ S / (2 - S)
    end

    # --- Probability distance relationships ---

    # Manhattan = 2 × TotalVariation.
    # TV = ½Σ|pᵢ−qᵢ|, Manhattan = Σ|pᵢ−qᵢ|, so Manhattan = 2 TV.
    for (a, b) in [([1, 2, 3], [3, 2, 1]), ([4, 1, 0, 5], [1, 3, 2, 4])]
        @test manhattan_distance(a, b) ≈ 2 * total_variation_distance(a, b)
    end

    # Chord = √2 × Hellinger.
    # Hellinger = √(Σ(√pᵢ−√qᵢ)²/2), Chord = √(Σ(√pᵢ−√qᵢ)²), so Chord = √2 H.
    # See Legendre & Gallagher (2001) "Ecologically meaningful transformations".
    for (a, b) in [([1, 2, 3], [3, 1, 2]), ([4, 1, 0, 5], [1, 3, 2, 4])]
        @test chord_distance(a, b) ≈ sqrt(2) * hellinger_distance(a, b)
    end

    # Hellinger² = 1 − Bhattacharyya coefficient.
    # H²(p,q) = Σ(√pᵢ−√qᵢ)²/2 = (Σpᵢ+Σqᵢ−2Σ√(pᵢqᵢ))/2 = 1−BC(p,q).
    # See Deza & Deza (2016) "Encyclopedia of Distances", ch. 14.
    for (a, b) in [([1, 3, 2], [2, 1, 3]), ([5, 1, 4], [2, 4, 4])]
        @test hellinger_distance(a, b)^2 ≈ 1 - bhattacharyya_coefficient(a, b)
    end

    # JensenShannon (metric, distance=true) = √(JensenDifference divergence).
    # JensenShannon is the square-root metric form; JensenDifference is the raw divergence.
    for (a, b) in [([1, 2, 3], [3, 1, 2]), ([1, 0, 4], [0, 3, 1])]
        @test jensen_shannon_distance(a, b)^2 ≈ jensen_difference(a, b)
    end

    # Jensen–Shannon divergence = ½KL(p‖m) + ½KL(q‖m), m = (p+q)/2.
    # Lin (1991) "Divergence measures based on the Shannon entropy", IEEE Trans. Inf. Theory.
    a, b = [1, 2, 3], [3, 1, 2]
    pa = proportions(a); pb = proportions(b)
    m = (pa .+ pb) ./ 2
    js_via_kl = (kullback_leibler_divergence(pa, m) + kullback_leibler_divergence(pb, m)) / 2
    @test jensen_difference(a, b) ≈ js_via_kl

    # KL(p‖Uniform) = log₂(n) − H(p).
    # KL(p‖U) = Σpᵢ log(pᵢ/(1/n)) = −H(p) + log₂(n) (Cover & Thomas §2.3).
    for x in [[1, 2, 3], [4, 1, 5, 2], [3, 3, 3]]
        n = richness(x)
        uniform = ones(n)
        @test kullback_leibler_divergence(x, uniform) ≈ log2(n) - shannon_entropy(x)
    end

    # --- Self-identity: d(x,x)=0 for distances, s(x,x)=1 for similarities ---
    # Follows from definitions; failure indicates a normalization or sign bug.
    for x in [[1, 2, 3], [5, 0, 1, 3], [2, 2, 2]]
        for index in [Jaccard(), SorensenDice(), Overlap(), BrayCurtis(), Ruzicka(),
                      TotalVariation(), Manhattan(), Euclidean(), Hellinger(), Chord(),
                      JensenShannon(), MorisitaHorn()]
            if is_similarity(index)
                @test similarity(index, x, x) ≈ 1.0
                @test dissimilarity(index, x, x) ≈ 0.0
            else
                @test dissimilarity(index, x, x) ≈ 0.0
            end
        end
        @test similarity(Bhattacharyya(), x, x) ≈ 1.0
        @test dissimilarity(Bhattacharyya(), x, x) ≈ 0.0
        @test dissimilarity(KullbackLeibler(), x, x) ≈ 0.0
        @test dissimilarity(ShannonDifference(), x, x) ≈ 0.0
        @test dissimilarity(JensenDifference(), x, x) ≈ 0.0
    end

    # --- Symmetry for all symmetric indices ---
    # is_symmetric(index) == true iff f(a,b) = f(b,a); KullbackLeibler is the exception.
    for (a, b) in [([1, 2, 3], [3, 1, 2]), ([4, 1, 0, 2], [1, 3, 2, 0])]
        for index in [Jaccard(), SorensenDice(), BrayCurtis(), Ruzicka(),
                      TotalVariation(), Hellinger(), JensenShannon(), MorisitaHorn()]
            @test dissimilarity(index, a, b) ≈ dissimilarity(index, b, a)
        end
    end

    # --- Triangle inequality for metric indices ---
    # For any triplet (a,b,c): d(a,c) ≤ d(a,b) + d(b,c).
    # Indices with is_metric == true must satisfy this exactly (up to float tolerance).
    a, b, c = [4, 2, 0, 1], [1, 3, 2, 0], [0, 1, 4, 3]
    for index in [Jaccard(), TotalVariation(), Manhattan(), Euclidean(), Hellinger(),
                  Chord(), JensenShannon()]
        dab = dissimilarity(index, a, b)
        dbc = dissimilarity(index, b, c)
        dac = dissimilarity(index, a, c)
        @test dac ≤ dab + dbc + 1e-10
    end

    # --- Pairwise matrix properties ---
    # Diagonal must be 0 (distances) or 1 (similarities); off-diagonal symmetric.
    mat = [1 1 2 0 5; 3 0 1 1 0; 2 3 0 2 1]
    for index in [BrayCurtis(), Jaccard(), Hellinger(), JensenShannon()]
        D = dissimilarity(index, mat)
        @test all(≈(0.0), diag(D))
        @test D ≈ D'
    end
    S = similarity(Jaccard(), mat)
    @test all(≈(1.0), diag(S))
    @test S ≈ S'
end

@testset "analytic bounds" begin
    # Properties that hold for any valid abundance vector.
    # These serve as lightweight property tests over representative inputs.
    function check_alpha_bounds(x)
        n = richness(x)
        H = shannon_entropy(x)
        # Shannon entropy ≤ log₂(n), equality iff uniform (Cover & Thomas §2.6).
        @test H ≤ log2(n) + 1e-10
        # Pielou evenness ∈ [0, 1], achieves 1 iff uniform (Pielou 1966).
        J = pielou_evenness(x)
        @test 0 ≤ J ≤ 1 + 1e-10
        # GiniSimpson ∈ [0, 1−1/n] for n≥2 (Simpson 1949).
        G = gini_simpson_index(x)
        @test 0 ≤ G ≤ 1 - 1 / n + 1e-10
        # SampleCoverage ∈ [0, 1] (Good 1953, "The population frequencies of species").
        @test 0 ≤ sample_coverage(x) ≤ 1
        # Richness estimators ≥ observed richness (Chao 1984, §2).
        @test chao1(x) ≥ n - 1e-10
        # ACE requires sample coverage > 0; skip vectors where all rare species are singletons.
        if count(==(1), x) < sum(v for v in x if v ≤ 10)
            @test ace(x) ≥ n - 1e-10
        end
        # Hill numbers are non-increasing in q for non-uniform distributions (Hill 1973).
        if n > 1
            @test hill_number(x, 0) ≥ hill_number(x, 1) - 1e-10
            @test hill_number(x, 1) ≥ hill_number(x, 2) - 1e-10
        end
    end

    function check_pairwise_bounds(a, b)
        # Bounded [0,1] incidence and abundance similarity/dissimilarity indices.
        # For these, similarity(index, a, b) + dissimilarity(index, a, b) = 1 by construction.
        for index in [Jaccard(), SorensenDice(), Overlap(), BrayCurtis(), Ruzicka(),
                      TotalVariation(), Hellinger()]
            s = similarity(index, a, b)
            d = dissimilarity(index, a, b)
            @test 0 ≤ s ≤ 1 + 1e-10
            @test 0 ≤ d ≤ 1 + 1e-10
            @test s + d ≈ 1.0
        end
        # Jensen-Shannon divergence ∈ [0, 1] in bits (Lin 1991 Theorem 3).
        @test 0 ≤ jensen_difference(a, b) ≤ 1 + 1e-10
        # Jensen-Shannon distance (√divergence) ∈ [0, 1] in bits.
        @test 0 ≤ jensen_shannon_distance(a, b) ≤ 1 + 1e-10
        # Bhattacharyya coefficient ∈ [0, 1].
        # Note: Bhattacharyya *distance* = -log(coefficient) ∈ [0, ∞], not bounded to 1.
        # dissimilarity + similarity ≠ 1 for this index.
        @test 0 ≤ bhattacharyya_coefficient(a, b) ≤ 1 + 1e-10
    end

    # Use assemblages where ACE is well-defined: at least one non-singleton rare species.
    assemblages = [[1, 2, 3], [10, 10, 10], [1, 2, 3, 4, 100], [5, 3, 1], [2, 2]]
    for x in assemblages
        check_alpha_bounds(x)
    end

    pairs = [
        ([1, 2, 3], [3, 2, 1]),
        ([4, 1, 0, 5], [1, 3, 2, 4]),
        ([10, 1, 0], [0, 10, 1]),   # overlapping support (avoids infinite Bhattacharyya)
        ([1, 1, 1], [1, 1, 1]),
    ]
    for (a, b) in pairs
        check_pairwise_bounds(a, b)
    end
end
