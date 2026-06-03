using Aqua
using DiversityIndices
using JET
using Test

@testset "package quality" begin
    Aqua.test_all(DiversityIndices)
    JET.test_package(DiversityIndices; target_modules=(DiversityIndices,))
end
