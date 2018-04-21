@everywhere include("LJChain.jl")
using LJChain,CGStats

# Set up simulation parameters.
C=Chain(NBeads=10,SEps=1,LEps=100,SM=1,LM=10,Beta=1,ConstrainedBeads=2)
# C=Chain(NBeads=10,SEps=100,LEps=.01,SM=.1,LM=10)

S=Simulator(FullSteps=Int(1e4),ApproxSteps=Int(5e4),dt=1e-3,OrthoSamples=1e5,MacroSamples=1000)
# S=Simulator(FullSteps=Int(1e4),dt=1e-6,OrthoSamples=1e4,MacroSamples=2)

# Generate data
(RR,PP,Veff,Fmean,Fstd) = VEffConstrained(C,S);

# Write data to file
ROut	 = open("./RPar.csv","w+");
POut     = open("./PPar.csv","w+");
VeffOut  = open("./VeffPar.csv","w+");
FmeanOut = open("./FmeanPar.csv","w+");
FstdOut  = open("./FstdPar.csv","w+");
writecsv(ROut,RR);
writecsv(POut,PP);
writecsv(VeffOut,Veff);
writecsv(FmeanOut,Fmean);
writecsv(FstdOut,Fstd);

## Read in old data
#RR = readcsv("./RPar.csv"); RR = reshape(RR,length(RR));
#Fmean = readcsv("./FmeanPar.csv"); Fmean = reshape(Fmean,length(Fmean));
#Fstd = readcsv("./FstdPar.csv"); Fstd = reshape(Fstd,length(Fmean));

# Generate a random trajectory, using data-generated potentials
P0 = zeros(C.NBeads);
Q0 = [2.0:3.0:30.0;];
PApprox,QApprox,FApprox = ApproximateTraj(P0,Q0,RR,Fmean,Fstd,C,S)

## Write data to file
PApproxOut = open("./PApprox.csv","w+");
QApproxOut = open("./QApprox.csv","w+");
FApproxOut = open("./FApprox.csv","w+");
writecsv(PApproxOut,PApprox);
writecsv(QApproxOut,QApprox);
writecsv(FApproxOut,FApprox);
