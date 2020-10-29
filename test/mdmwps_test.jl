
# Example data for computing spectra with gaps

fn  = @__DIR__()*"/../Examples/data/temp.txt"
dat = readdlm(fn,'\t')

# Spectrum with gaps:
tt = dat[:,1]
xx = dat[:,2:3]

bw = 5/length(tt);
k = Int64(round(2*bw*length(tt))) - 1;
kk = Int64(2*bw*length(xx) - 1)
nz = 0;
alpha = 0.95;

lams,ei = Multitaper.mdslepian(bw,k,tt)
out1 = mdmultispec(tt,xx[:,1],jk=true,dof=true)
out = mdmultispec(tt,xx)

@testset "Slepians with gaps" begin
  # Slepians with gaps: eigenvalues
  @test lams ≈ [0.9999999999814994, 0.9999999979919965, 0.9999999112387221,
              0.9999972080829853, 0.9999348274751064, 0.9990806667306628, 0.9901498911725759,
              0.9185956480485494, 0.8611080144951507]
end

@testset "Spectrum analysis with gaps" begin

  # spectrum with gaps
  @test out[1][1].S[end-5:end] ≈ [0.14565610885883232, 0.13580424778417852,
  0.1293182864247158, 0.1395702616312291, 0.1664040092175493, 0.09202992618124456]

  # dof
  @test out1[2][end-5:end] ≈ [17.444104148892166, 17.408259466237617,
  17.382026447855413, 17.422494161582264, 17.506998428080035, 17.170675010806065]

  # F-test with gaps
  @test out1[1].Fpval[end-5:end] ≈ [0.16807947246930965, 0.09822111646921738,
  0.0341711343274661, 0.4933223598055142, 0.9312676978366197, 0.7307073177440282]

  # the number of lines found
  @test findall(out1[1].Fpval .<0.05)[1:10] == [16, 26, 45, 50, 51, 68, 85, 94, 113, 126] 

end

@testset "Coherence with gaps" begin

  @test out[2][1,2].coh[end-5:end] ≈ [0.21020101113566358, 0.43985588013870025,
  0.39732804488526385, 0.10864516699443383, 0.12246458417213957, 0.08184136299900113]

end
