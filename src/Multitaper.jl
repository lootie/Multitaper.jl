
# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019

__precompile__()

module Multitaper

  using StatsFuns, FastGaussQuadrature, Statistics, FFTW 
  using LinearAlgebra, SpecialFunctions, RecipesBase, Arpack
  using Distributions

  # Utiltity

  include("StructsTypes.jl")
    export Ecoef, MtParams, MtSpec, Demodulate, MtAcf, MtAcvf, MtCeps
    export MtCcvf, MtCcf, MtTransf, MtCoh 
    export Demodulate

  include("Utils.jl")
    export tanhtrans, atanhtrans, ejn

  # Univariate 

  include("Univariate.jl")
    export multispec, mt_acvf, welsh, blockerr

  include("PlotsRecipes/univariaterecipes.jl")

  # Univariate multitaper on data with gaps 

  include("MDmwps.jl")
    export mdslepian, mdmwps

  include("Demodulation.jl")
    export demodulate

  include("crossings.jl")

  # Multivariate

  include("Multivariate.jl")
    export multispec, mt_ccvf
  
  include("PlotsRecipes/multivariaterecipes.jl")

  # Two-dimensional 
  include("2DRectangle.jl")
    export multispec2_Rectangle

end

