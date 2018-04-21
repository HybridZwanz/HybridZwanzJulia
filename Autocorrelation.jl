@everywhere include("LJChain.jl")
using LJChain,CGStats,PyPlot

# Set up simulation parameters.
C=Chain(NBeads=10,SEps=100,LEps=1,SM=1,LM=100,Beta=1,ConstrainedBeads=2)

S=Simulator(FullSteps=Int(1e4),ApproxSteps=Int(1e4),dt=1e-3,OrthoSamples=10000,MacroSamples=1000,LagStep=10,NumLags=1000)

# LagStep = spacing of lag times for autocorrelation
# NumLags = number of lag times to compute

FOut	 = open("./F.csv","w+");
POut     = open("./P.csv","w+");
FFOut	 = open("./FFull.csv","w+");
PFOut     = open("./PFull.csv","w+");

# Generate IC
(P,Q,F,H,PCG,QCG) = InitConds(C);

println("Equilibrate ...")
# Run in to equilibrate
for i = 1:S.FullSteps
 P,Q,F = LeapfrogStep(P,Q,F,C,S.dt)
end
println("... done")

Pfull=P;
Qfull=Q;

# Compute autocorrelation stats
(StressStats,StressMean,MomStats,MomMean) = ConstrainedAutocorr(P,Q,C,S);

# Use PyPlot to plot autocorrelation of momenta
A=MomStats.-kron((MomMean.^2)',ones(S.NumLags));
plot(S.dt*S.LagStep*[0:S.NumLags-1;],A./A[1])
writecsv(POut, [S.dt*S.LagStep*[0:S.NumLags-1;] A./A[1]]);
savefig("OrthocorrM.svg")

# Use PyPlot to plot autocorrelation of forces
A=StressStats.-kron((StressMean.^2)',ones(S.NumLags));
plot(S.dt*S.LagStep*[0:S.NumLags-1;],A./A[1])
writecsv(FOut,[S.dt*S.LagStep*[0:S.NumLags-1;] A./A[1]]);
savefig("OrthocorrF.svg")



# Compute autocorrelation stats
(StressStats,StressMean,MomStats,MomMean) = FullAutocorr(Pfull,Qfull,C,S);

A=MomStats.-kron((MomMean.^2)',ones(S.NumLags));
plot(S.dt*S.LagStep*[0:S.NumLags-1;],A./A[1])
writecsv(PFOut, [S.dt*S.LagStep*[0:S.NumLags-1;] A./A[1]]);
savefig("FullcorrM.svg")

# Use PyPlot to plot autocorrelation of forces
A=StressStats.-kron((StressMean.^2)',ones(S.NumLags));
plot(S.dt*S.LagStep*[0:S.NumLags-1;],A./A[1])
writecsv(FFOut,[S.dt*S.LagStep*[0:S.NumLags-1;] A./A[1]]);
savefig("FullcorrF.svg")
