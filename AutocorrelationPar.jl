@everywhere include("LJChain.jl")
using LJChain,CGStats,PyPlot

# Set up simulation parameters.
C=Chain(NBeads=100,SEps=1,LEps=100,SM=10,LM=1,Beta=1,ConstrainedBeads=2)
#C=Chain(NBeads=10,SEps=100,LEps=1,SM=1,LM=100,Beta=1,ConstrainedBeads=2)

#S=Simulator(FullSteps=Int(1e3),ApproxSteps=Int(1e3),dt=1e-3,OrthoSamples=1e3,MacroSamples=20,LagStep=10,NumLags=1000)
S=Simulator(FullSteps=Int(1e4),ApproxSteps=Int(1e4),dt=1e-3,OrthoSamples=1e4,MacroSamples=1e3,LagStep=20,NumLags=1000)

# # LagStep = spacing of lag times for autocorrelation
# # NumLags = number of lag times to compute
#
FOut	   = open("./FParOrtho.csv","w+");
FFOut	  = open("./FParFull.csv","w+");
PFOut   = open("./PParFull.csv","w+");
AppFOut = open("./FParApp.csv","w+");
AppMOut = open("./PParApp.csv","w+");


# Generate data
(RR,PP,Veff,Fmean,Fstd) = VEffConstrained(C,S);

C=Chain(NBeads=100,SEps=1,LEps=100,SM=10,LM=1,Beta=1,ConstrainedBeads=100)

# Compute autocorrelation stats for full dynamics
(StressAutocorr,MomAutocorr) = FullAutocorrPar(C,S,RR,Fmean);

# Use PyPlot to plot autocorrelation of momenta
plot(S.dt*S.LagStep*[0:S.NumLags-1;],MomAutocorr)
writecsv(PFOut, [S.dt*S.LagStep*[0:S.NumLags-1;] MomAutocorr]);
close(PFOut)
savefig("FullcorrM.svg")

# Use PyPlot to plot autocorrelation of forces
plot(S.dt*S.LagStep*[0:S.NumLags-1;],StressAutocorr)
writecsv(FFOut,[S.dt*S.LagStep*[0:S.NumLags-1;] StressAutocorr]);
close(FFOut)
savefig("FullcorrF.svg")



# Compute autocorrelation stats for constrained dynamics
(StressAutocorr,MomAutocorr) = ConstrainedAutocorrPar(C,S);

# Don't plot autocorrelation of momenta (as it's not well-defined)

# Use PyPlot to plot autocorrelation of forces
plot(S.dt*S.LagStep*[0:S.NumLags-1;],StressAutocorr)
writecsv(FOut,[S.dt*S.LagStep*[0:S.NumLags-1;] StressAutocorr]);
close(FOut)
savefig("OrthocorrF.svg")



# Compute autocorrelation stats for approximate dynamics
(StressAutocorr,MomAutocorr) = ApproximateAutocorrPar(C,S,RR,Fmean,Fstd*0.5);

# Use PyPlot to plot autocorrelation of momenta
plot(S.Stochdt*S.LagStep*[0:S.NumLags-1;],MomAutocorr)
writecsv(AppMOut, [S.Stochdt*S.LagStep*[0:S.NumLags-1;] MomAutocorr]);
close(AppMOut)
savefig("ApproximatecorrM.svg")

# Use PyPlot to plot autocorrelation of forces
plot(S.Stochdt*S.LagStep*[0:S.NumLags-1;],StressAutocorr)
writecsv(AppFOut,[S.Stochdt*S.LagStep*[0:S.NumLags-1;] StressAutocorr]);
close(AppFOut)
savefig("ApproximatecorrF.svg")
