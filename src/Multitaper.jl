
# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019

__precompile__()

module Multitaper

  using StatsFuns, FastGaussQuadrature, Statistics, FFTW 
  using SparseArrays, LinearAlgebra, SpecialFunctions, RecipesBase, Arpack, Distributions

  # Utiltity

  include("mtbase.jl")
    export conv

  include("StructsTypes.jl")
    export ecoef, mtparams, mtspec, mtdemod, mtacf, mtacvf, mtceps
    export mtccvf, mtccf, mttransf, mtdemod, mtcoh, mtcspecqspec 
    export Demodulate

  include("Utils.jl")
    export tanhtrans, atanhtrans, EJN

  # Univariate 

  include("Univariate.jl")
    export multispec, mt_acvf, welsh, blockerr

  include("PlotsRecipes/univariaterecipes.jl")

  # Univariate multitaper on data with gaps 

  include("MDmwps.jl")
    export MDslepian, MDmwps

  include("Demodulation.jl")
    export demodulate

  include("crossings.jl")

  # Multivariate

  include("Multivariate.jl")
    export multispec, mt_ccvf
  
  include("PlotsRecipes/multivariaterecipes.jl")

end

