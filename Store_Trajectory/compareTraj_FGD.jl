
@time begin

include("../LJChain.jl")
include("../CGStats.jl")
using Main.LJChain,Main.CGStats,DelimitedFiles

#####################################################################################################################
## Set up simulation parameters.

C = Chain(SM=1, LM=1, Beta=1.096); # system C ≣ [Patt=[0,1,0],NBeads=10,ConstrainedBeads=2,SEps=1,LEps=100,SM=10,LM=1,Dist=1,Beta=1,Gamma=1]

# S = Simulator(); # [dt=1e-4,Stochdt=1e-4,FullSteps=Int(1e5),OrthoSteps=1,ApproxSteps=Int(1e5),OrthoSamples=Int(1e5),MacroSamples=50,LagStep=400,NumLags=30]
S = Simulator(dt=1e-3, Stochdt=1e-3, FullSteps=Int(1e7), ApproxSteps=Int(1e7));
# S = Simulator(dt=1e-3, Stochdt=1e-3, FullSteps=Int(1e5), ApproxSteps=Int(1e5));

######################################################################################################################
## Generate initial data

(P0, Q0, F0, H0, P0CG, Q0CG) = InitConds(C);

########################################################################################
## Generate a random trajectory, using data-generated potentials

(PFull, QFull, FFull, Kmean, Hmean) = LangFullTraj_NVT(P0, Q0, C, S);

#########################################################################################
## Write data to file

open("Store_Trajectory/Outcomes/P_FGD_M1m1_1e7(B1.096).csv", "w+") do x
    writedlm(x, (C.Phi*PFull)');
end
open("Store_Trajectory/Outcomes/Q_FGD_M1m1_1e7(B1.096).csv", "w+") do x
    writedlm(x, (C.Psi*QFull)');
end
open("Store_Trajectory/Outcomes/F_FGD_M1m1_1e7(B1.096).csv", "w+") do x
    writedlm(x, (C.PhiF*FFull)');
end

#################################################################################
end # @time
