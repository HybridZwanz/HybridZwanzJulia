
@time begin

include("../LJChain.jl")
include("../CGStats.jl")
using Main.LJChain, Main.CGStats, DelimitedFiles, Random

# Random.seed!(1234)

####################################################################################################################
## Set up simulation parameters.

C = Chain(SM=1, LM=1, Beta=1.096); # system C ≣ [Patt=[0,1,0],NBeads=10,ConstrainedBeads=2,SEps=1,LEps=100,SM=10,LM=1,Dist=1,Beta=1,Gamma=1]

# S = Simulator(); # [dt=1e-4,Stochdt=1e-4,FullSteps=Int(1e5),OrthoSteps=1,ApproxSteps=Int(1e5),OrthoSamples=Int(1e5),MacroSamples=50,LagStep=400,NumLags=30]
S = Simulator(dt=1e-3, Stochdt=1e-3, FullSteps=Int(1e7), ApproxSteps=Int(1e7), m_FF=-105.73, q_FF=320.07, coeffEXP_A=6.827e11, coeffEXP_B=8.5891);
# S = Simulator(dt=1e-3, Stochdt=1e-3, FullSteps=Int(1e5), ApproxSteps=Int(1e5));

#####################################################################################################################
## Generate initial data

(P0, Q0, F0, H0, P0CG, Q0CG) = InitConds(C);

#####################################################################################################################
## Generate effective potential data

# ## NB: thesis data
# # (RR,PP,Veff,Fmean,Fstd) = VEffConstrainedLangDyn(C,S);
#
# # CG = readdlm("ABCD_Systems/Outcomes/sC_2000.csv", Float64);
# # RR = CG[:, 1];
# # Fmean = CG[:, 4];
# # Fstd = CG[:, 5];
#
# ## to generate the coefficients for the cubic spline interpolation for mean force and fluctuating force
# # NB: the grid points *must* be equally spaced
#
# ## Mean Force
#
# force_data = readdlm("Store_Trajectory/force_to_interpolate.dat", Float64);
#
# # i punti della griglia *devono* essere equispaziati
# RR = force_data[:, 1];
# dx = RR[2] - RR[1];
#
# Fmean = force_data[:, 2];
# Fmean_d = force_data[:, 3];
#
# Fmean_coeff1 = (2. * (Fmean - circshift(Fmean, -1)) + (Fmean_d + circshift(Fmean_d, -1)) * dx) / dx;
# Fmean_coeff2 = circshift(Fmean, -1) - Fmean + (-Fmean_d - Fmean_coeff1) * dx;
#
# FdataMF = zeros(length(Fmean) - 1, 5);
# FdataMF[:, 1] = RR[1:end - 1];
# FdataMF[:, 2] = Fmean[1:end - 1];
# FdataMF[:, 3] = Fmean_d[1:end - 1];
# FdataMF[:, 4] = Fmean_coeff1[1:end - 1];
# FdataMF[:, 5] = Fmean_coeff2[1:end - 1];
#
# ## Fluctuating Force
#
# flucts_data = readdlm("Store_Trajectory/flucts_to_interpolate.dat", Float64);
#
# # i punti della griglia *devono* essere equispaziati
# RR = flucts_data[:, 1];
# dx = RR[2] - RR[1];
#
# FF = flucts_data[:, 2];
# FF_d = flucts_data[:, 3];
#
# FF_coeff1 = (2. * (FF - circshift(FF, -1)) + (FF_d + circshift(FF_d, -1)) * dx) / dx;
# FF_coeff2 = circshift(FF, -1) - FF + (-FF_d - FF_coeff1) * dx;
#
# FdataFF = zeros(length(FF) - 1, 5);
# FdataFF[:, 1] = RR[1:end - 1];
# FdataFF[:, 2] = FF[1:end - 1];
# FdataFF[:, 3] = FF_d[1:end - 1];
# FdataFF[:, 4] = FF_coeff1[1:end - 1];
# FdataFF[:, 5] = FF_coeff2[1:end - 1];

## NB: paper data
FdataMF = readdlm("ABCD_Systems/Outcomes/Coeff_MF_SplineI_sC_M1m1.csv", Float64);
FdataFF = readdlm("ABCD_Systems/Outcomes/Coeff_FF_SplineI_sC_M1m1.csv", Float64);

#########################################################################################
## Generate a random trajectory, using data-generated potentials

(PApprox, QApprox, FApprox) = ApproximateTraj_MMZD(P0CG, Q0CG, FdataMF, FdataFF, C, S);

#########################################################################################
## Write data to file

open("Store_Trajectory/Outcomes/PApprox_MMZD_M1m1_1e7(B1.096)_CORRETTA.csv", "w+") do x
    writedlm(x, PApprox);
end
open("Store_Trajectory/Outcomes/QApprox_MMZD_M1m1_1e7(B1.096)_CORRETTA.csv", "w+") do x
    writedlm(x, QApprox);
end
open("Store_Trajectory/Outcomes/FApprox_MMZD_M1m1_1e7(B1.096)_CORRETTA.csv", "w+") do x
    writedlm(x, FApprox);
end

#################################################################################
end # @time
