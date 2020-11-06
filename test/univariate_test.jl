

# Example data for computing spectra

fn  = @__DIR__()*"/../Examples/data/soirecruitment.dat"
dat = readdlm(fn)

# Spectrum tests
println()
NW = 4.0
K  = 6
dt = 1.0/12
univ = multispec(dat[:, 3], dt = dt, NW = NW, K = K, jk = true, Ftest = true, guts =
        true, pad = 2.0, a_weight = true)
 
@testset "Spectrum:" begin 

  # univariate frequency
  @test univ.f[1:5]         == 12*range(0, 1, length=454 + 453)[1:5]

  # univariate MT spectrum
  @test univ.S[1:5]         ≈ [0.05776450047987208, 0.0510836499696889, 
                               0.053947076576308435, 
                               0.03960810660934419, 0.03238905443102084]

  # univariate MT phase
  @test univ.phase          == nothing 

  # univariate parameters
  @test univ.params.M       == 906

  # univariate eigencoefficients
  @test univ.coef.coef[1:2] ≈ [0.0763450617941093 + 0.0im, 
                              -0.018461837449426248 - 0.07253045824096403im] 

  @test univ.coef.wts       == nothing

  # univariate F-test
  @test univ.Fpval[1:5]     ≈ [0.9999999498486979, 0.3106832963307704, 
                               0.46150154664925935, 
                               0.7590806611936394, 0.023822662121927296]

  # univariate jackknife
  @test univ.jkvar[1:5]     ≈ [0.48594970945868243, 0.48941505597497825, 
                               0.23173477815035817,
                               0.2766449068431087, 0.16729659627580204]

end

uu = Multitaper.testTsq(ones(K),EigenCoefficient(univ.coef.coef[[1,3,5],:],
      univ.coef.wts))

@testset "T^2 test" begin
  
  # T squared test
  @test uu ≈ 0.5645801176110739

end

univ =  welch(dat[:, 3], 3, 0.5, dt = dt, NW = NW, K = K, guts = false, pad = 2.0)

@testset "Welch" begin
  
  # Univariate frequency
  @test univ[1].f[1:3]  ≈ 0.0:0.02654867256637168:0.05309734513274336

  # Univariate Welch Spectrum
  @test univ[1].S[1:5] ≈[0.031851476832981986, 0.02690415047577252, 
                  0.028692054543592537, 0.032023262389080605, 0.039447614292936535] 

  # Welch Parameters
  @test univ[1].params == MTParameters(4.0, 6, 453, 0.08333333333333333, 452, 3, 0.5)

  # Effective bandwidth
  @test univ[2] ≈ 24.15894039735099

end

# Time-domain type tests
acvf = mt_acvf(dat[:,3], NW = NW, K = K, dt = dt)

acf = mt_acf(dat[:,3], NW = NW, K = K, dt = dt)

ceps = mt_cepstrum(dat[:,3], NW = NW, K = K, dt = dt)

cet_path = @__DIR__()*"/../Examples/data/CETmonthly.dat"
cet = readdlm(cet_path)
cdm = demodulate(cet[:,3], 1.0, 2.0, 15*12, true, 1.0/12, 0.0)

@testset "Autocovariance" begin

  # univariate acvf
  @test collect(acvf.lags[1:3]) ≈ [0.0, 0.0833333333333, 0.166666666667] 
  @test acvf.acvf[1:5]      ≈  [0.012255897713702007, 0.007260759850546666, 
                                0.004502534551242927, 
                                0.002412704328674198, 0.0002545617261240325]
  @test acvf.params.K       == K

end

@testset "Autocorrelation" begin

  # univariate acvf
  @test collect(acf.lags[1:3]) ≈ [0.0, 0.0833333333333, 0.166666666667] 
  @test acf.acf[1:5]     ≈ [1.0, 0.5924298668411039, 0.3673769687396398, 
                            0.196860677612935, 0.020770549173189858] 
  @test acf.params       == MTParameters(4.0, 6, 453, 0.08333333333333333, 453, 1, 
                            nothing)

end

@testset "Cepstrum" begin

  @test collect(ceps.lags[1:3]) ≈ [0.0, 0.0833333333333, 0.166666666667] 
  @test ceps.ceps[1:5]      ≈ [-5.143067306451898, 0.3666296616675006, 
                               0.16488294222178018, 
                               0.1880899501961062, 0.1398326259962087]
  @test ceps.params.NW      == NW

end

@testset "Demodulation" begin

  @test cdm.mag[1:5] ≈ [6.786697185396128, 6.79296024005296, 6.799489409493339, 
                        6.805499133878922, 6.809510248764921]
  @test cdm.phase[1:5] ≈ [140.9216771230031, 140.91171594928556, 140.90139636523648,
  140.8842740916691, 140.8859505605725]

end
