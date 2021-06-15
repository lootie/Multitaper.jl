
# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019

__precompile__()

module Multitaper

  using StatsFuns, FastGaussQuadrature, Statistics, FFTW 
  using LinearAlgebra, SpecialFunctions, RecipesBase, Arpack
  using Distributions

  # Utiltity

  include("StructsTypes.jl")
    export EigenCoefficient, MTParameters, MTSpectrum, Demodulate
    export MTAutocorrelationFunction, MTAutocovarianceFunction, MTCepstrum
    export MtCrossCovarianceFunction, MTCrossCorrelationFunction, MTTransferFunction
    export MTCoherence 

  include("Utils.jl")
    export tanhtrans, atanhtrans, ejn, unwrapphase

  # Univariate 

  include("Univariate.jl")
    export dpss_tapers, multispec, mt_acvf, mt_acf, mt_cepstrum, welch, blockerr

  include("PlotsRecipes/univariaterecipes.jl")

  # Multitaper on data with gaps 

  include("MDmwps.jl")
    export mdslepian, mdmultispec, gpss

  include("Demodulation.jl")
    export demodulate

  include("crossings.jl")
    export  Pgram_upcrossings, MT_Upcrossings, uctable

  # Unequally sampled multitaper
  include("bronez.jl")
    export bspec

  # Multivariate

  include("Multivariate.jl")
    export multispec, mt_ccf, mt_ccvf, mtm_svd, LFV_spectrum
  
  include("PlotsRecipes/multivariaterecipes.jl")

  # Two-dimensional 
  include("2DRectangle.jl")
    export rectsleps, multispec2_Rectangle

end

