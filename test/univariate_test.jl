
module univariate_test

# Example data for computing spectra

using Multitaper, Test 
using DelimitedFiles, LinearAlgebra, Statistics

fn  = "../Examples/data/soirecruitment.dat"
dat = readdlm(fn)

# Spectrum tests
println("Spectrum:")
NW = 4.0
K  = 6
dt = 1.0/12
univ = multispec(dat[:, 3], dt = dt, NW = NW, K = K, jk = true, Ftest = true, guts = true, pad = 2.0)
  
println("1. Testing univariate frequency")
@test univ.f[1:5]         == 12*LinRange(0, 1, 454 + 453)[1:5]

println("2. Testing univariate MT spectrum")
@test univ.S[1:5]         ≈ [0.05776450047987208, 0.0510836499696889, 0.053947076576308435,
                             0.03960810660934419, 0.03238905443102084]

println("3. Testing univariate MT phase")
@test univ.phase          == nothing 

println("4. Testing univariate parameters")
@test univ.params.M       == 906

println("5. Testing univariate eigencoefficients")
@test univ.coef.coef[1:2] ≈ [0.0763450617941093 + 0.0im, -0.018461837449426248 - 0.07253045824096403im]

println("6. Testing univariate F-test")
@test univ.Fpval[1:5]     ≈ [0.9999999498486979, 0.3106832963307704, 0.46150154664925935, 
                             0.7590806611936394, 0.023822662121927296]

println("7. Testing univariate jackknife")
@test univ.jkvar[1:5]     ≈ [0.48594970945868243, 0.48941505597497825, 0.23173477815035817,
                             0.2766449068431087, 0.16729659627580204]

# Currently not testing this functionality
println("8. Testing univariate T-squared test")
@test univ.Tsq_pval       == nothing

# Time-domain type tests
println("Autocorrelation:")
acvf = mt_acvf(dat[:,3], NW = NW, K = K, dt = dt)
acvfspec = mt_acvf(univ)

println("9-11. Testing univariate acvf")
@test collect(acvf.lags[1:3]) ≈ [0.0, 0.0833333333333, 0.166666666667] 
@test acvf.acvf[1:5]      ≈  [0.012255897713702007, 0.007260759850546666, 0.004502534551242927, 
                              0.002412704328674198, 0.0002545617261240325]
@test acvf.params.K       == K

println("Cepstrum:")
ceps = mt_acvf(dat[:,3], NW = NW, K = K, dt = dt, typ = :ceps)
cepsspec = mt_acvf(univ, typ = :ceps)

println("12-14. Testing univariate cepstrum")
@test collect(ceps.lags[1:3]) ≈ [0.0, 0.0833333333333, 0.166666666667] 
@test ceps.ceps[1:5]      ≈ [-5.143067306451898, 0.3666296616675006, 0.16488294222178018, 
                             0.1880899501961062, 0.1398326259962087]
@test ceps.params.NW      == NW

println("15-16. Testing complex demodulation")
cet = readdlm("../Examples/data/CETmonthly.dat")
cdm = demodulate(cet[:,3], 1.0/12, 1.0, 2.0, 15*12, true)
@test cdm.mag[1:5] ≈ [6.786697185396128, 6.79296024005296, 6.799489409493339, 
                      6.805499133878922, 6.809510248764921]
@test cdm.phase[1:5] ≈ [-30.0, -219.08828405071444, -219.09860363476352, 
                      -219.1157259083309, -219.1140494394275]
end
