
@time begin

include("../LJChain.jl")
include("../CGStats.jl")
using Main.LJChain,Main.CGStats,DelimitedFiles

#####################################################################################################################
## Set up simulation parameters.

C = Chain(SM=1, LM=1, Beta=1.096); # system C ≣ [Patt=[0,1,0],NBeads=10,ConstrainedBeads=2,SEps=1,LEps=100,SM=10,LM=1,Dist=1,Beta=1,Gamma=1]

# S = Simulator(); # [dt=1e-4,Stochdt=1e-4,FullSteps=Int(1e5),OrthoSteps=1,ApproxSteps=Int(1e5),OrthoSamples=Int(1e5),MacroSamples=50,LagStep=400,NumLags=30]
S = Simulator(dt=1e-3, Stochdt=1e-3, FullSteps=Int(1e7), ApproxSteps=Int(1e7), m_MF=421.52, q_MF=-1237.5);
# S = Simulator(dt=1e-3, Stochdt=1e-3, FullSteps=Int(1e5), ApproxSteps=Int(1e5));

######################################################################################################################
## Generate initial data

(P0, Q0, F0, H0, P0CG, Q0CG) = InitConds(C);

#########################################################################################
## Generate effective potential data

# ## NB: thesis data
# data extracted from the Kriging interpolation of the original Mean Force data
# force_data = readdlm("Store_Trajectory/force_to_interpolate.dat", Float64);
#
# # to generate the coefficients for the cubic spline interpolation
# # NB: the grid points *must* be equally spaced
# RR = force_data[:, 1];
# dx = RR[2] - RR[1];
#
# Fmean = force_data[:, 2];
# Fmean_d = force_data[:, 3];
#
# Fmean_coeff1 = (2. * (Fmean - circshift(Fmean, -1)) + (Fmean_d + circshift(Fmean_d, -1)) * dx) / dx;
# Fmean_coeff2 = circshift(Fmean, -1) - Fmean + (-Fmean_d - Fmean_coeff1) * dx;
#
# Fdata = zeros(length(Fmean) - 1, 5);
# Fdata[:, 1] = RR[1:end - 1];
# Fdata[:, 2] = Fmean[1:end - 1];
# Fdata[:, 3] = Fmean_d[1:end - 1];
# Fdata[:, 4] = Fmean_coeff1[1:end - 1];
# Fdata[:, 5] = Fmean_coeff2[1:end - 1];

## NB: paper data
Fdata = readdlm("ABCD_Systems/Outcomes/Coeff_MF_SplineI_sC_M1m1.csv", Float64);


########################################################################################
## Generate a random trajectory, using data-generated potentials

(PApprox2, QApprox2, FApprox2) = DetCGTraj_NVT(P0CG, Q0CG, Fdata, C, S);

#########################################################################################
## Write data to file

open("Store_Trajectory/Outcomes/PApprox_DCGD_M1m1_1e7(B1.096)_CORRETTA.csv", "w+") do x
    writedlm(x, PApprox2);
end
open("Store_Trajectory/Outcomes/QApprox_DCGD_M1m1_1e7(B1.096)_CORRETTA.csv", "w+") do x
    writedlm(x, QApprox2);
end
open("Store_Trajectory/Outcomes/FApprox_DCGD_M1m1_1e7(B1.096)_CORRETTA.csv", "w+") do x
    writedlm(x, FApprox2);
end

#################################################################################
end # @time
