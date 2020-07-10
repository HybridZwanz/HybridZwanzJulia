# NB: to run this script from the terminal is not necessary to give the following command "julia -p n scritpt.jl"
#     because I explicitly added in the script the command-line option "-p n".

using Distributed

@time begin

# Set up parallel computing
addprocs(Sys.CPU_THREADS);

@everywhere include("LJChain.jl")
@everywhere include("CGStats.jl")
using Main.LJChain,Main.CGStats,SharedArrays,Parameters,Statistics,LinearAlgebra,StatsBase,DelimitedFiles

#####################################################################################################################
# Set up simulation parameters.

C = Chain(); # system C ≣ [Patt=[0,1,0],NBeads=10,ConstrainedBeads=2,SEps=1,LEps=100,SM=10,LM=1,Dist=1,Beta=1,Gamma=1]

# S = Simulator(); # [dt=1e-4,Stochdt=1e-4,FullSteps=Int(1e5),OrthoSteps=1,ApproxSteps=Int(1e5),OrthoSamples=Int(1e5),MacroSamples=50,LagStep=400,NumLags=30]
S = Simulator(dt=1e-3, Stochdt=1e-4, FullSteps=Int(1e6), ApproxSteps=Int(1e5));

#################################################################################
# Generate initial data

(P0, Q0, F0, H0, P0CG, Q0CG) = InitConds(C);

#################################################################################
# Generate effective potential data

# (RR,PP,Veff,Fmean,Fstd) = VEffConstrainedLangDyn(C,S);
CG = readdlm("ABCD_Systems/Outcomes/sC_2000.csv", Float64);
RR = CG[:, 1];
Fmean = CG[:, 4];
Fstd = CG[:, 5];

########################################################################################
# Generate a random trajectory, using data-generated potentials

## FGD
# (PFull, QFull, FFull, Kmean, Hmean) = LangevinFullTraj(P0, Q0, C, S);

## MMZD
(PApprox, QApprox, FApprox) = ApproximateTraj(P0CG, Q0CG, RR, Fmean, Fstd, C, S);

## DCGD
# (PApprox2, QApprox2, FApprox2) = DetCGCanEnsTraj(P0CG, Q0CG, RR, Fmean, C, S);

#########################################################################################
# Write data to file

## FGD
# open("Store_Trajectory/Outcomes/PFull.csv", "w+") do x
#     writedlm(x, (C.Phi*PFull)');
# end
# open("Store_Trajectory/Outcomes/QFull.csv", "w+") do x
#     writedlm(x, (C.Psi*QFull)');
# end
# open("Store_Trajectory/Outcomes/FFull.csv", "w+") do x
#     writedlm(x, (C.PhiF*FFull)');
# end

## MMZD
open("Store_Trajectory/Outcomes/PApprox.csv", "w+") do x
    writedlm(x, PApprox);
end
open("Store_Trajectory/Outcomes/QApprox.csv", "w+") do x
    writedlm(x, QApprox);
end
open("Store_Trajectory/Outcomes/FApprox.csv", "w+") do x
    writedlm(x, FApprox);
end

## DCGD
# open("Store_Trajectory/Outcomes/PApprox_DCGD.csv", "w+") do x
#     writedlm(x, PApprox2);
# end
# open("Store_Trajectory/Outcomes/QApprox_DCGD.csv", "w+") do x
#     writedlm(x, QApprox2);
# end
# open("Store_Trajectory/Outcomes/FApprox_DCGD.csv", "w+") do x
#     writedlm(x, FApprox2);
# end

#################################################################################
end # @time
