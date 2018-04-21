@everywhere include("LJChain.jl")
using LJChain,CGStats

# Set up simulation parameters.
C=Chain(NBeads=10,SEps=1,LEps=100,SM=10,LM=1,Beta=1,ConstrainedBeads=2)
# C=Chain(NBeads=10,SEps=100,LEps=.01,SM=.1,LM=10)

S=Simulator(FullSteps=Int(1e6),ApproxSteps=Int(1e6),dt=1e-3,OrthoSamples=1e5,MacroSamples=2000)
# S=Simulator(FullSteps=Int(1e4),dt=1e-6,OrthoSamples=1e4,MacroSamples=2)

# Generate data
(P0,Q0,F,H,PCG,QCG) = InitConds(C);
P0CG = C.Phi*P0
Q0CG = C.Psi*Q0

# Generate data
(RR,PP,Veff,Fmean,Fstd) = VEffConstrained(C,S);

# Generate a random trajectory, using data-generated potentials
PFull,QFull,FFull = FullTraj(P0,Q0,C,S)
PApprox,QApprox,FApprox = ApproximateTraj(P0CG,Q0CG,RR,Fmean,Fstd,C,S)
PApprox2,QApprox2,FApprox2 = ApproximateTraj(P0CG,Q0CG,RR,Fmean,Fstd*0,C,S)


## Write data to file
PFullOut = open("./PFull.csv","w+");
QFullOut = open("./QFull.csv","w+");
FFullOut = open("./FFull.csv","w+");
writecsv(PFullOut,(C.Phi*PFull)');
writecsv(QFullOut,(C.Psi*QFull)');
writecsv(FFullOut,(C.PhiF*FFull)');
PApproxOut = open("./PApprox.csv","w+");
QApproxOut = open("./QApprox.csv","w+");
FApproxOut = open("./FApprox.csv","w+");
writecsv(PApproxOut,PApprox);
writecsv(QApproxOut,QApprox);
writecsv(FApproxOut,FApprox);
PApproxOut2 = open("./PApprox2.csv","w+");
QApproxOut2 = open("./QApprox2.csv","w+");
FApproxOut2 = open("./FApprox2.csv","w+");
writecsv(PApproxOut2,PApprox2);
writecsv(QApproxOut2,QApprox2);
writecsv(FApproxOut2,FApprox2);
