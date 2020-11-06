
fn  = @__DIR__()*"/../Examples/data/soirecruitment.dat"
dat = readdlm(fn)

# Spectrum tests
NW = 4.0
K  = 6
dt = 1.0/12
multiv = multispec(dat[:, 3:4], dt = dt, NW = NW, K = K, jk = true, Ftest = true,
                    guts = true, pad = 2.0, Tsq = [0.1, 0.2])

@testset "Coherency" begin

  #  multivariate frequency
  @test multiv[1][1].f[1:5]         == 12*LinRange(0, 1, 454 + 453)[1:5]

  #  univariate MT spectra
  @test multiv[1][1].S[1:5]         ≈ [0.05776450047987208, 0.0510836499696889,
                               0.053947076576308435,
                               0.03960810660934419, 0.03238905443102084]

  #  univariate MT phase
  @test multiv[1][2].phase          == nothing 

  #  multivariate parameters
  @test multiv[1][2].params.M       == 906

  #  multivariate eigencoefficients
  @test multiv[1][1].coef.coef[1:2] ≈ [0.0763450617941093 + 0.0im,
                              -0.018461837449426248 - 0.07253045824096403im]

  #  multivariate F-test
  @test multiv[1][1].Fpval[1:5]     ≈ [0.9999999498486979, 0.3106832963307704,
                              0.46150154664925935, 
                              0.7590806611936394, 0.023822662121927296]

  #  univariate jackknife
  @test multiv[1][1].jkvar[1:5]     ≈ [0.48594970945868243, 0.48941505597497825, 
                                       0.23173477815035817, 
                                       0.2766449068431087, 0.16729659627580204] 

  # multivariate jackknife: coherence variance
  @test multiv[2][1,2].jkvar[1][1:5] ≈ [11.663257249685838, 3.652684030823364, 0.8240763280267569, 1.4806959512973812, 1.6987467111535954]

  # multivariate jackknife: phase variance
  @test multiv[2][1,2].jkvar[2][1:5] ≈ [0.0, 2.2808214820275725, 1.252986500086001, 1.6043511134415998, 2.433432774031582]

  #  multivariate T-squared test
  @test multiv[1][1].Tsq_pval       ≈ [0.07046082220259418]

  #  the coherency
  @test multiv[2][1,2].coh[1:5]     ≈ [7.046455236530036, 6.697940886784345, 
                                       6.226973678450996, 
                                       5.18871118000256, 5.082716554324644]

  #  the phase
  @test multiv[2][1,2].phase[1:5]   ≈ [180.0, 178.71279387445034, 178.30193928288293,
  177.60013658787375, 179.26526345068578]

end

multiv2 = multispec(dat[:, 3], dat[:,4], outp = :transf, dt = dt, NW = NW, K = K, jk =
                    true, guts = true, pad = 2.0)

@testset "Transfer Function" begin
  #  multivariate frequency
  @test multiv2.f[1:5] ≈ 0.0:0.013245033112582781:0.052980132450331126

  #  transfer function
  @test multiv2.transf[1:5] ≈ [8977.335209581075, 9414.676069658892, 9808.02657929821, 8655.730134294896, 9433.008704433021] 

  #  phase  
  @test multiv2.phase[1:5] ≈ [180.0, 182.39849609134032, 181.98452791869974, 182.58379890869585, 181.47460681717862]

  #  parameters
  @test multiv2.params == MTParameters(4.0, 6, 453, 0.08333333333333333, 906, 1, nothing)

  #  coefficients, jackknife variance already tested via univariate tests

end

S = multispec(dat[:,3], dat[:,4], outp = :spec, offset = 0.25,
              NW = NW, K = K, dt = dt, ctr = true, pad = 2.0,
              guts = false, jk = false, Tsq = nothing)
 
@testset "Offset Cross-spectrum" begin

  #  offset frequency
  @test S.f[1:5] == 2.8609271523178808:0.013245033112582781:2.9139072847682117

  # offset cross-spectrum
  @test S.S[1:5] ≈ [1.5074259822190432, 1.128276855356535, 1.1648381894848983, 
                    1.2577989895592048, 1.7636931138403107]

  # offset phase
  @test S.phase[1:5] ≈ [-163.89510356195487, -130.37341672575664, 
                -159.68242136376077, -192.09704297403184, -127.52800412648241]
  
end

# Time-domain type tests
ccvf = mt_ccvf(dat[:,3], dat[:,4], NW = NW, K = K, dt = dt, pad = 2.0)

@testset "Cross-covariance" begin

  #  multivariate ccvf
  @test ccvf.lags[1:3]              ≈ -37.75:0.16666666666666666:-37.416666666666664
  @test ccvf.ccvf[1:5]              ≈ [-0.0034710006247021317, 0.009173222535744282, 0.017410059264281845, 0.0454050228426607, -0.039622614768931544]
  @test ccvf.params.K               == K

end

ccf = mt_ccvf(dat[:,3], dat[:,4], typ = :ccf, NW = NW, K = K, dt = dt, pad = 2.0)

@testset "Cross-correlation" begin

  #  multivariate ccvf
  @test         ccf.lags[1:3]      ≈ -37.75:0.16666666666666666:-37.416666666666664
  
  @test ccf.ccf[1:5]               ≈ [-0.0002030558867807063, 0.0005366397295858057, 0.0010185002554174937, 0.002656224580312745, -0.0023179497926959173]
  @test ccf.params                 == MTParameters(4.0, 6, 453, 0.08333333333333333, 
                                    906, 1, nothing)

end
