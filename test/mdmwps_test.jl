
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
out2 = mdmultispec(tt,xx[:,1],xx[:,2])

@testset "Slepians with gaps" begin
  # Slepians with gaps: eigenvalues
  @test lams ≈ [0.999999999981499, 0.9999999979919967, 0.9999999112387216, 0.9999972080829845, 0.9999348274751061, 0.9990806667306631, 0.9901498911725768, 0.9185956480485502, 0.8611080144951503]
end

@testset "Spectrum analysis with gaps" begin

  # spectrum with gaps
  @test out[1][1].S[end-5:end] ≈ [0.06543035539617903, 0.06120799684389478, 0.06496856845578666, 0.07705250113383467, 0.08705314879916222, 0.04244661937652378]

  # coherence with gaps
  @test out2.coh[1:5] ≈[0.3914653089257123, 0.5814720910611886, 0.18221545009607953, 0.2847615929878291, 0.23658188530825575]  

  # dof
  @test out1[2][end-5:end] ≈ [16.90312736502505, 16.843935582791087, 16.896950430892392, 17.038629637400643, 17.13101303801798, 16.479401696185825]

  # F-test with gaps
  @test out1[1].Fpval[end-5:end] ≈ [0.1059516846656845, 0.038428506229240655, 0.505329149922955, 0.9332725255153663, 0.7313065910084713, 0.8575739191868774]

  # the number of lines found
  @test findall(out1[1].Fpval .<0.05)[1:10] == [15, 44, 49, 50, 67, 84, 93, 112, 125, 129]
  end

@testset "Coherence with gaps" begin

  @test out[2][1,2].coh[end-5:end] ≈ [0.4398558660852374, 0.39732802025166336, 0.10864515185504586, 0.12246456200932032, 0.08184138313232747, 0.03477254755728354]

end


