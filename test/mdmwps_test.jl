
# Example data for computing spectra with gaps

fn  = "../Examples/data/temp.txt"
dat = readdlm(fn,',')

# Spectrum tests
println("Spectrum with gaps:")
tt = dat[:,1]
xx = dat[:,2]

bw = 5/length(tt);
k = Int64(round(2*bw*length(tt))) - 1;
kk = Int64(2*bw*length(xx) - 1)
nz = 0;
alpha = 0.95;

lams,ei = Multitaper.mdslepian(bw,k,tt)
out = mdmwps(tt,xx,jk=true,Tsq=nothing,dof=true,alpha=0.95)

@testset "Slepians with gaps, spectrum analysis with gaps." begin
 
  println("Testing Slepians with gaps: eigenvalues")
  @test lams ≈ [0.9999999999814994, 0.9999999979919965, 0.9999999112387221,
              0.9999972080829853, 0.9999348274751064, 0.9990806667306628, 0.9901498911725759,
              0.9185956480485494, 0.8611080144951507]

  println("Testing spectrum with gaps")
  @test out[1].S[end-5:end] ≈ [ 0.1358023997817664, 0.12939199720491235, 0.13953326445088568, 
              0.16563989397260973, 0.18402382557990213, 0.08885119971059109]

  println("Testing dof")
  @test out[2][end-5:end] ≈ [17.40825230897549, 17.382337634406642, 17.422357651541745, 
              17.50493555937301, 17.550216141919826, 17.145772544877815]

  println("Testing F-test with gaps")
  @test out[1].Fpval[end-5:end] ≈ [ 0.09822119696194964, 0.03416880618582274, 
              0.49332546035901725, 0.9312765181262197, 0.7249977043925413, 0.8515331773112675]

  println("Testing the number of lines found")
  @test findall(out[1].Fpval .<0.05)[1:10] == [15, 25, 44, 49, 50, 67, 84, 93, 112, 125]

end
