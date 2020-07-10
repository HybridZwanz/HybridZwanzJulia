include("../../LJChain.jl")
include("../../CGStats.jl")
# include("LJChain.jl")
# include("CGStats.jl")

using Main.LJChain,Main.CGStats,StatsBase,DelimitedFiles,Statistics

C = Chain();
S = Simulator();

################################################ pdf(P) ################################################################

# P = readdlm("/home/pmams3/Desktop/PFull.csv", Float64);

# P_FGD = readdlm("/home/mspino/MarcoMasterProject/Julia_Code/Store_Trajectory/Outcomes/P_FGD_1e7(B1.096)_CORRETTA.csv", Float64);
P_FGD = readdlm("Store_Trajectory/Outcomes/P_FGD_1e7.csv", Float64);
P_FGD = P_FGD[:];
h1 = fit(Histogram, P_FGD);

x01 = Vector(h1.edges[1]);
H1 = zeros(length(h1.weights), 2);
H1[:, 1] = (x01[2:end] + x01[1:end-1]) / 2;
bin_size = x01[2] - x01[1];
H1[:, 2] = h1.weights ./ (length(P_FGD[:, 1]) * bin_size)
# open("Data_Analysis/PDF/Outcomes/HistP_FGD_1e7(B1.096)_CORRETTA.csv", "w+") do x
open("Data_Analysis/PDF/Outcomes/HistP_FGD_1e7.csv", "w+") do x
   writedlm(x, H1);
end



# P_DCGD = readdlm("/home/mspino/MarcoMasterProject/Julia_Code/Store_Trajectory/Outcomes/PApprox_DCGD_1e7(B1.096)_CORRETTA.csv", Float64);
# P_DCGD = P_DCGD[:];
# h2 = fit(Histogram, P_DCGD);
#
# x02 = Vector(h2.edges[1]);
# H2 = zeros(length(h2.weights), 2);
# H2[:, 1] = (x02[2:end] + x02[1:end-1]) / 2;
# bin_size = x02[2] - x02[1];
# H2[:, 2] = h2.weights ./ (length(P_DCGD[:, 1]) * bin_size);
# open("Data_Analysis/PDF/Outcomes/HistP_DCGD_1e7(B1.096)_CORRETTA.csv", "w+") do x
#     writedlm(x, H2);
# end
#
#
# # P_MMZD = readdlm("Store_Trajectory/Outcomes/PApprox_MMZD_1.csv", Float64);
# P_MMZD = readdlm("/home/mspino/MarcoMasterProject/Julia_Code/Store_Trajectory/Outcomes/PApprox_MMZD_1e7(B1.096).csv", Float64);
# for i in 1:size(P_MMZD, 1)
#     # println(mean(P_MMZD[i, :]));
#     P_MMZD[i, :] .-= mean(P_MMZD[i, :]);
# end
#
# P_MMZD = P_MMZD[:];
# h3 = fit(Histogram, P_MMZD);
#
# x03 = Vector(h3.edges[1]);
# H3 = zeros(length(h3.weights), 2);
# H3[:, 1] = (x03[2:end] + x03[1:end-1]) / 2;
# bin_size = x03[2] - x03[1];
# H3[:, 2] = h3.weights ./ (length(P_MMZD[:, 1]) * bin_size);
# open("Data_Analysis/PDF/Outcomes/HistP_MMZD_1e7(B1.096).csv", "w+") do x
#     writedlm(x, H3);
# end


############################################# pdf(R) ###############################################################

# Q_FGD = readdlm("/home/mspino/MarcoMasterProject/Julia_Code/Store_Trajectory/Outcomes/Q_FGD_1e7(B1.096)_CORRETTA.csv", Float64);
# lf = size(Q_FGD,1);
# dQ4 = zeros(lf, C.NBeads);
# for i in 1:lf
#     dQ4[i, :] = circshift(Q_FGD[i, :], -1) - Q_FGD[i, :];
#     dQ4[i, C.NBeads] += C.Period;
# end
#
# dQ4 = dQ4[:];
# h4 = fit(Histogram, dQ4);
#
# x04 = Vector(h4.edges[1]);
# H4 = zeros(length(h4.weights), 2);
# bin_size = x04[2] - x04[1];
# H4[:, 1] = (x04[2:end] +  x04[1:end-1]) / 2;
# H4[:, 2] = h4.weights ./ (length(dQ4[:, 1]) * bin_size);
# open("Data_Analysis/PDF/Outcomes/HistQ_FGD_1e7(B1.096)_CORRETTA.csv", "w+") do x
#     writedlm(x, H4);
# end
#
#
# Q_DCGD = readdlm("/home/mspino/MarcoMasterProject/Julia_Code/Store_Trajectory/Outcomes/QApprox_DCGD_1e7(B1.096)_CORRETTA.csv", Float64);
# ld = size(Q_DCGD,1);
# dQ5 = zeros(ld, C.NBeads);
# for i in 1:ld
#     dQ5[i, :] = circshift(Q_DCGD[i, :], -1) - Q_DCGD[i, :];
#     dQ5[i, C.NBeads] += C.Period;
# end
#
# dQ5 = dQ5[:];
# h5 = fit(Histogram, dQ5);
#
# x05 = Vector(h5.edges[1]);
# H5 = zeros(length(h5.weights), 2);
# bin_size = x05[2] - x05[1];
# H5[:, 1] = (x05[2:end] +  x05[1:end-1]) / 2;
# H5[:, 2] = h5.weights ./ (length(dQ5[:, 1]) * bin_size);
# open("Data_Analysis/PDF/Outcomes/HistQ_DCGD_1e7(B1.096)_CORRETTA.csv", "w+") do x
#     writedlm(x, H5);
# end
#
#
# # Q_MMZD = readdlm("Store_Trajectory/Outcomes/QApprox_MMZD_1.csv", Float64);
# Q_MMZD = readdlm("/home/mspino/MarcoMasterProject/Julia_Code/Store_Trajectory/Outcomes/QApprox_MMZD_1e7(B1.096).csv", Float64);
# lm = size(Q_MMZD,1);
# dQ6 = zeros(lm, C.NBeads);
# for i in 1:lm
#     dQ6[i, :] = circshift(Q_MMZD[i, :], -1) - Q_MMZD[i, :];
#     dQ6[i, C.NBeads] += C.Period;
# end
#
# dQ6 = dQ6[:];
# h6 = fit(Histogram, dQ6);
#
# x06 = Vector(h6.edges[1]);
# H6 = zeros(length(h6.weights), 2);
# bin_size = x06[2] - x06[1];
# H6[:, 1] = (x06[2:end] +  x06[1:end-1]) / 2;
# H6[:, 2] = h6.weights ./ (length(dQ6[:, 1]) * bin_size);
# open("Data_Analysis/PDF/Outcomes/HistQ_MMZD_1e7(B1.096).csv", "w+") do x
#     writedlm(x, H6);
# end
