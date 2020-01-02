
# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019

using Test

printstyled("Running tests:\n", color=:blue)

tests = ["univariate", "multivariate"]

for t in tests
    @testset "Test $t multitaper" begin
      include("$(t)_test.jl")
    end
end

