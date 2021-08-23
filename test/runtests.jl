
# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019


using Multitaper, Test 
using DelimitedFiles, LinearAlgebra, Statistics

printstyled("Running tests:\n", color=:blue)

tests = ["utils", "univariate", "crossings", "mdmwps", "multivariate", "twodim", "bronez"]

for t in tests
    @testset "$t multitaper" begin
      include("$(t)_test.jl")
    end
end

