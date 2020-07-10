include("../../LJChain.jl")
include("../../CGStats.jl")

using Main.LJChain,Main.CGStats,SharedArrays,Parameters,Statistics,LinearAlgebra,StatsBase, DelimitedFiles



###########################################################################################################

C = Chain();      # System C ([0,1,0], SM (mass of 0 particles)=10, LM (mass of 1 particles)=1, γ=1, β=1)
S = Simulator();  # dt=1e-4, FullSteps=Int(1e5),  OrthoSteps=1, OrthoSamples=Int(1e5)


P=readdlm("Store_FullandConstrainedTrajectory/P_compMBdBETA(1e-2).csv",Float64);

## dt=1e-4
h1 = fit(Histogram,P[:,1]);

x01 = Vector(h1.edges[1]);

H1 = zeros(length(h1.weights),2);
H1[:,1] = (x01[2:end]+x01[1:end-1])/2;
H1[:,2] = h1.weights./length(P[:,1]);

open("Data_Analysis/Maxwell-Boltzmann_distribution/HistMaxBolt_dt1e-4_B1e-2.csv","w+") do x
    writedlm(x,H1);
end

## dt=1e-3
h2 = fit(Histogram,P[:,2]);

x02 = Vector(h2.edges[1]);

H2 = zeros(length(h2.weights),2);
H2[:,1] = (x02[2:end]+x02[1:end-1])/2;
H2[:,2] = h2.weights./length(P[:,2]);

open("Data_Analysis/Maxwell-Boltzmann_distribution/HistMaxBolt_dt1e-3_B1e-2.csv","w+") do x
    writedlm(x,H2);
end

## dt=1e-2
h3 = fit(Histogram,P[:,3]);

x03 = Vector(h3.edges[1]);

H3 = zeros(length(h3.weights),2);
H3[:,1] = (x03[2:end]+x03[1:end-1])/2;
H3[:,2] = h3.weights./length(P[:,3]);

open("Data_Analysis/Maxwell-Boltzmann_distribution/HistMaxBolt_dt1e-2_B1e-2.csv","w+") do x
    writedlm(x,H3);
end


# x=(x0(2:end)+x0(1:end-1))/2;
# plot(x,px0,'r')
# [px0,x0]=histcounts(Q0,100,'Normalization','probability');
# N,edges
#
# """
#     MBd = MaxwellBoltzmannDist(P::Array{Float64,1},C::Chain)
#
# Computes the Maxwell-Boltzmann Distribution of bead with momenta P.
# """
# function MaxwellBoltzmannDist(P::Array{Float64,1},C::Chain)
# 	a = C.Beta/(2*C.BeadMass);
# 	b = sqrt(16*a^3/pi);
#   return b*(P.^2).*exp.(-a*P.^2)
# end
