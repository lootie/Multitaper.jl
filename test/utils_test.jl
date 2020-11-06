
N = 1024
sig = [0.05, 0.5, 0.9, 0.99, 0.999, 0.9999]
K = 7
NW = 4.0

@testset "Utility functions test" begin 
  
  # Expected jackknife variance
  @test ejn(K) ≈ 0.2204585537918871 

  # Bandwidth
  @test Multitaper.bw(NW, N, 1.0) ≈ 0.00390625
 
  # Logarithmic range
  @test Multitaper.logRange(1, 5, 5) == [0.9, 0.99, 0.999, 0.9999, 0.99999]  
 
  # Magnitude squared coherence significance
  @test Multitaper.mscsig(0.9, K) ≈ 0.999999

  # Hyperbolic tangent transformation
  @test tanhtrans(7.0, 7) ≈  0.9654629992362074

end


