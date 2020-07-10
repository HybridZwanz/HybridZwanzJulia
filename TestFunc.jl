# NB: to run this script from the terminal is not necessary to give the following command "julia -p n scritpt.jl"
#     because I explicitly added in the script the command-line option "-p n".

# using Distributed
#
# # Set up parallel computing
# addprocs(Sys.CPU_THREADS);
#
# @everywhere include("LJChain.jl")
# @everywhere include("CGStats.jl")

include("LJChain.jl")
include("CGStats.jl")
using Main.LJChain,Main.CGStats,SharedArrays,Parameters,Statistics,LinearAlgebra,StatsBase,DelimitedFiles # add Pyplot if necessary


C = Chain();
S = Simulator();

# S = Simulator(dt=1e-4,Stochdt=1e-5,FullSteps=Int(1e5),OrthoSamples=Int(1e5),ApproxSteps=Int(1e5),LagStep=10,NumLags=100,MacroSamples=50);

(P,Q,F,H,PCG,QCG) = InitConds(C);

####################################################################################################

# (RR,PP,Veff,Fmean,Fstd) = VEffConstrained(C,S);

# P = randn(C.NBeads);
# Q = [0.0:C.BeadsDist:C.BeadsDist*(C.NBeads-1);];
# Q = [2.0:3.0:30.0;];

# CG = readdlm("ABCD_Systems/Outcomes/sC_100.csv",Float64);
# RR = CG[:,1];
# Fmean = CG[:,4];
# Fstd = CG[:,5];

# Pstore,Qstore,Fstore = ApproximateTraj(P0CG,Q0CG,RR,Fmean,Fstd,C,S);
# (StressStats,StressMean,StressSqMean,MomStats,MomMean,MomSqMean) = ApproximateAutocorr(P,Q,C,S,RR,Fmean,Fstd);

# (P,Q,F) = ApproxDynStep(P,Q,F,MeanForce,FlucStep,C,S.Stochdt);
# PNext = ApproxFlucStep(P,Q,RR,Fstd,C,S.Stochdt)

#####################################################################################################

# (Pstore,Qstore,Fstore) = FullTraj(P0,Q0,C,S);
#
# Pi = Pstore[:,end];
# Qi = Qstore[:,end];
#
# PCG = C.PhiConstr*Pi;
# QCG = C.PsiConstr*Qi;
#
# (StressMoments,Pf,Qf) = ConstrainedSample(Pi,Qi,C,S);

#####################################################################################################

# P0 = zeros(C.NPart);
# Q0 = Vector(1.0:1.0:30.0);
#
# @time begin
#     (Pstore,Qstore,Fstore,Kmean,Hmean) = LangevinFullTraj(P0,Q0,C,S);
# end
#
# Pi = Pstore[:,end];
# Qi = Qstore[:,end];
#
# PCG = C.PhiConstr*Pi;
# QCG = C.PsiConstr*Qi;
#
# println("Computing orthogonal dynamics through constrained Langevin dynamics ...");
#
# @time begin
#     S = Simulator();
#     (StressMoments,Pf,Qf) = ConstrSampleLangDyn(Pi,Qi,C,S);
#     println("  dt, time step size = ",S.dt," seconds.");
#     println("  Total simulation time = ",S.dt*S.OrthoSamples," seconds.");
# end;
# #
# @time begin
#     S = Simulator();
#     (StressMoments,Pf,Qf) = ConstrSampleLangDyn(Pi,Qi,C,S);
#     println("  dt, time step size = ",S.dt," seconds.");
#     println("  Total simulation time = ",S.dt*S.OrthoSteps*S.OrthoSamples," seconds.");
# end;

#####################################################################################################
#
# PhiConstr = kron(Diagonal(ones(C.NBeads))[1:C.ConstrainedBeads,:],ones(Float64,length(C.Patt))');
# PsiConstr = kron(Diagonal(ones(C.NBeads))[1:C.ConstrainedBeads,:],C.MassPatt'./C.BeadMass);
# ZetaConstr = Diagonal(ones(C.NPart))-PsiConstr'PhiConstr;
#
# #########################################################################################################
#
# P0,Q0,F0 = InitConds(C);
#
# (StressMoments,Portho,Qortho) = ConstrSampleLangDyn(P0,Q0,C,S);
#
# MomentSamples = SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads-1,S.Moments);
# ii=1;
# p=1;
# MomentSamples[ii,:,:] = remotecall_fetch(ConstrSampleLangDyn,p,P0,Q0,C,S,0)[1];
# # NB:ConstrSampleLangDyn returns a Tuple and remotecall_fetch()[1] returns the first value of this Tuple
#
# ############################################################################################################
#
# EpsPatt = (sqrt.(C.LEps)*C.LTotPatt+sqrt.(C.SEps)*C.STotPatt).*circshift(sqrt.(C.LEps)*C.LTotPatt+sqrt.(C.SEps)*C.STotPatt,-1)
