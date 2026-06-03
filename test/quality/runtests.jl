using Aqua
using DiversityAndDissimilarity
using JET
using Test

@testset "package quality" begin
    Aqua.test_all(DiversityAndDissimilarity)
    JET.test_package(DiversityAndDissimilarity; target_modules=(DiversityAndDissimilarity,))
end
