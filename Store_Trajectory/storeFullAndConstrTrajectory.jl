# NB: to run this script from the terminal is not necessary to give the following command "julia -p n scritpt.jl"
#     because I explicitly added in the script the command-line option "-p n".

using Distributed

# Set up parallel computing
addprocs(Sys.CPU_THREADS);

@everywhere include("LJChain.jl")
@everywhere include("CGStats.jl")
using Main.LJChain,Main.CGStats,SharedArrays,Parameters,Statistics,LinearAlgebra,StatsBase,DelimitedFiles

################################################################################################################

C = Chain(Beta=1e-2);          # System C ([0,1,0], SM (mass of 0 particles)=10, LM (mass of 1 particles)=1, γ=1, β=1)
S = Simulator(FullSteps=Int(1e7));               # dt=1e-4, FullSteps=Int(1e5),  OrthoSteps=1, OrthoSamples=Int(1e5)

P0 = zeros(C.NPart);
Q0 = Vector(1.0:1.0:30.0);

#################################################################################################################
## Store the data to compute Maxwell-Boltzmann distribution ##

# to change the values of the time step
dtSamples = [1e-4,1e-3,1e-2];
nsamples = length(dtSamples);


j=0;
np = nprocs();
PFull = SharedArray{Float64}(C.NPart,S.FullSteps,nsamples);
P     = SharedArray{Float64}(S.FullSteps,nsamples);

@sync begin
    while j < nsamples
        for p in 1:np
            if p != myid() || np == 1
                if j >= nsamples
                    break;
                end
                global j += 1;
                println("Start dt sample ",j);
                pp = p;
                jj = j;
                @async begin
                    println(" CPU_THREADS ",pp-1,", Worker ",pp,", dt=",dtSamples[jj]);
                    S = Simulator(dt=dtSamples[jj],FullSteps=Int(1e7));
                    PFull[:,:,jj] = remotecall_fetch(LangevinFullTraj,pp,P0,Q0,C,S)[1];
                    for i =1:S.FullSteps
                        P[i,jj] = (C.PhiConstr*PFull[:,i,jj])[1];
                    end
                end
            end
        end
    end
end

println(" ...Momenta of single bead computed!");

# Write data to file
Pb = zeros(S.FullSteps,nsamples);
Pb[:,1:3] = P;
open("Store_FullandConstrainedTrajectory/Outcomes/P_compMBdBETA(1e-2).csv","w+") do x
    writedlm(x,Pb);
end




#############################################################################################################
#
# nsamples = 1;
#
# Pstore_DET = SharedArray{Float64}(C.NPart,S.FullSteps,nsamples);
# Qstore_DET = SharedArray{Float64}(C.NPart,S.FullSteps,nsamples);
# Fstore_DET = SharedArray{Float64}(C.NPart,S.FullSteps,nsamples);
#
# Pstore_STOC = SharedArray{Float64}(C.NPart,S.FullSteps,nsamples);
# Qstore_STOC = SharedArray{Float64}(C.NPart,S.FullSteps,nsamples);
# Fstore_STOC = SharedArray{Float64}(C.NPart,S.FullSteps,nsamples);
#
# for i in 1:nsamples
#     (Pstore_DET[:,:,i],Qstore_DET[:,:,i],Fstore_DET[:,:,i],Kmean,Hmean) = LangevinFullTraj(P0,Q0,C,S);
#     (Pstore_STOC[:,:,i],Qstore_STOC[:,:,i],Fstore_STOC[:,:,i]) = LangConstrainedTraj(Pstore_DET[:,end,i],Qstore_DET[:,end,i],C,S);
# end
#
#
# # Write data to file
# open("1_Qstore_DET.csv","w+") do x
#     writedlm(x,Qstore_DET[:,:,1]);
# end
# open("1_Pstore_DET.csv","w+") do x
#     writedlm(x,Pstore_DET[:,:,1]);
# end
# open("1_Fstore_DET.csv","w+") do x
#     writedlm(x,Fstore_DET[:,:,1]);
# end
#
# open("1_Qstore_STOC.csv","w+") do x
#     writedlm(x,Qstore_STOC[:,:,1]);
# end
# open("1_Pstore_STOC.csv","w+") do x
#     writedlm(x,Pstore_STOC[:,:,1]);
# end
# open("1_Fstore_STOC.csv","w+") do x
#     writedlm(x,Fstore_STOC[:,:,1]);
# end
