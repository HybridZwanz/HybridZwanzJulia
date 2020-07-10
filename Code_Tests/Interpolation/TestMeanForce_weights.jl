# NB: to run this script from the terminal is not necessary to give the following command "julia -p n scritpt.jl"
#     because I explicitly added in the script the command-line option "-p n".

using Distributed

@time begin

# Set up parallel computing
# addprocs(Sys.CPU_THREADS);

@everywhere include("../../LJChain.jl")
@everywhere include("../../CGStats.jl")
using Main.LJChain,
      Main.CGStats,
      SharedArrays,
      Parameters,
      Statistics,
      LinearAlgebra,
      StatsBase,
      SparseArrays,
      DelimitedFiles # add Pyplot if necessary


C = Chain();
S = Simulator(
    dt = 1e-3,
    Stochdt = 1e-4,
    FullSteps = Int(1e6),
    OrthoSamples = Int(1e6),
    ApproxSteps = Int(1e5),
    LagStep = 10,
    NumLags = 100,
    MacroSamples = 50,
);

# CG = readdlm("../../ABCD_Systems/Outcomes/sC_20_300p.csv", Float64);
# RR = CG[:, 1];
# Fmean = CG[:, 4];
# Fstd = CG[:, 5];

RR = readdlm("../../ABCD_Systems/Outcomes/RR_sB.csv", Float64);
Fmean = readdlm("../../ABCD_Systems/Outcomes/Fmean_sB.csv", Float64);


N = 500;
σ = 0.1;
Rmin = 1.5;
Rmax = 4.5;
range = Rmax - Rmin;
delta = range / N

dQ = [Rmin:delta:Rmax;]

### Mean Force
Sigma = zeros(size(dQ));
for i in 1:length(dQ)
    # if dQ[i] < 0.9
    #     m = (- 9620.78 + 11462.61) / (1.09-0.95);
    #     q = -9620.78 - m * 1.09;
    #     Sigma[i] = m * dQ[i] + q;
    # elseif dQ[i] > 3.01
    #     W = spdiagm(0 => exp.(-(RR .- 3.01) .^2 ./ σ^2));
    #     J = ones(3, length(RR));
    #     for j = 1:2
    #         J[j+1, :] = (3.01 .- RR).^j;
    #     end
    #     K = J * W;
    #     c = (K * J') \ (K * Fmean);
    #     Sigma[i] = c[1];
    # else
        W = spdiagm(0 => exp.(-(RR .- dQ[i]) .^2 ./ σ^2));

        J = ones(3, length(RR));
        for j = 1:2
            J[j+1, :] = (dQ[i] .- RR).^j;
        end
        K = J * W;
        c = (K * J') \ (K * Fmean);
        Sigma[i] = c[1];
    # end
end

### Fluctuating Force
# Gvec = zeros(size(dQ,1));
# for i=1:length(dQ)
#     if dQ[i] < 0.9
#         m = (6464.26 - 8859.34) / (1.2-0.90);
#         q = 8859.34 - m * 0.90;
#         Gvec[i] = m * dQ[i] + q;
#     elseif dQ[i] > 3.01
#         W = spdiagm(0 => exp.(-(RR .- 3.01) .^2 ./ σ^2));
#         J = ones(3, length(RR));
#         for j = 1:2
#             J[j+1, :] = (3.01 .- RR).^j;
#         end
#         K = J * W;
#         coeffs = (K * J') \ (K * Fstd);
#         Gvec[i] = coeffs[1];
#     else
#         W = spdiagm(0 => exp.(-(RR .- dQ[i]) .^2 ./ σ^2));
#
#         J = ones(3, length(RR));
#         for j=1:2
#             J[j+1,:] = (dQ[i] .- RR).^j;
#         end
#         K = J * W;
#         coeffs = (K * J') \ (K * Fstd);
#         Gvec[i] = coeffs[1];
#     end
# end

###
#
# l = length(dQ[:, 1]);
# RB = reshape(dQ[:, 1], l);
# i = sortperm(RB);
# RB = RB[i];
# # println("R =", RB);
#
# Sigma1 = reshape(Sigma, l);
# Sigma1 = Sigma1[i];
# # println("Sigma = ", Sigma1);
#
# Gvec1 = reshape(Gvec, l);
# Gvec1 = Gvec1[i];
# # println("FF = ", Gvec1);
#
#open("./Code_Tests/Interpolation/Outcomes/RFmean_algInt_test.csv", "w+") do x
# open("force_300p.dat", "w+") do x
#     writedlm(x, [RB Sigma1])
# end
#
# open("flucts_300p.dat", "w+") do x
#     writedlm(x, [RB Gvec1])
# end

open("force_test.dat", "w+") do x
    writedlm(x, [RB Sigma1])
end


######################################################################

end # @time
