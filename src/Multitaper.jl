
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
    export tanhtrans, atanhtrans, ejn, unwrapphase

  # Univariate 

  include("Univariate.jl")
    export dpss_tapers, multispec, mt_acvf, welch, blockerr

  include("PlotsRecipes/univariaterecipes.jl")

  # Univariate multitaper on data with gaps 

  include("MDmwps.jl")
    export mdslepian, mdmultispec

  include("Demodulation.jl")
    export demodulate

  include("crossings.jl")
    export  Pgram_upcrossings, MT_Upcrossings, uctable

  # Multivariate

  include("Multivariate.jl")
    export multispec, mt_ccvf
  
  include("PlotsRecipes/multivariaterecipes.jl")

  # Two-dimensional 
  include("2DRectangle.jl")
    export rectsleps, multispec2_Rectangle

end

