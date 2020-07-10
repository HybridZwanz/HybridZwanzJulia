# NB: to run this script from the terminal is not necessary to give the following command "julia -p n scritpt.jl"
#     because I explicitly added in the script the command-line option "-p n".

using Distributed

@time begin

# Set up parallel computing
# addprocs(Sys.CPU_THREADS);
addprocs(4);

@everywhere include("LJChain.jl")
@everywhere include("CGStats.jl")
using Main.LJChain,Main.CGStats,SharedArrays,Parameters,Statistics,LinearAlgebra,StatsBase,DelimitedFiles

####################################################################################################
# General system

# C = Chain();     # system C ≣ [Patt=[0,1,0],NBeads=10,ConstrainedBeads=2,SEps=1,LEps=100,SM=10,LM=1,Dist=1,Beta=1,Gamma=1]
# S = Simulator(); # [dt=1e-4,Stochdt=1e-4,FullSteps=Int(1e5),OrthoSteps=1,ApproxSteps=Int(1e5),OrthoSamples=Int(1e5),MacroSamples=50,LagStep=400,NumLags=30]

# (RR,PP,Veff,Fmean,Fstd) = VEffConstrainedLangDyn_NVT(C,S);


#####################################################################################################

S = Simulator(dt=1e-3, FullSteps=Int(1e6), OrthoSamples=Int(1e6), check=Int(5e3), MacroSamples=100);

######################################################################################################
# system A
# ( [0,1,0], SM (mass of particles 0), LM (mass of particles 1), SEps=1.0, LEps=100.0, γ=1, β=1 )

# C = Chain(SM=1, LM=10, NBeads = 100);
# (RR,PP,Veff,Fmean,Fstd) = VEffConstrainedLangDyn_NVT(C,S);
#
# # Write data to file
# open("./ABCD_Systems/Outcomes/RR_sA.csv","w+") do x
#    writedlm(x,RR);
# end
# open("./ABCD_Systems/Outcomes/PP_sA.csv","w+") do x
#    writedlm(x,PP);
# end
# open("./ABCD_Systems/Outcomes/Veff_sA.csv","w+") do x
#    writedlm(x,Veff);
# end
# open("./ABCD_Systems/Outcomes/Fmean_sA.csv","w+") do x
#    writedlm(x,Fmean);
# end
# open("./ABCD_Systems/Outcomes/Fstd_sA.csv","w+") do x
#    writedlm(x,Fstd);
# end

#################################################################################
# system B
# ( [0,1,0], SM (mass of particles 0), LM (mass of particles 1), γ=1, β=1 )

# C = Chain(SM=1, LM=1, SEps=10, LEps=.1, NBeads = 10);
# (RR, PP, Veff, Fmean, Fstd, FdataMF, FdataFF) = VEffConstrainedLangDyn_NVT(C, S);
#
# # Write data to file
# open("./ABCD_Systems/Outcomes/RR_sB.csv","w+") do x
#    writedlm(x, RR);
# end
# open("./ABCD_Systems/Outcomes/PP_sB.csv","w+") do x
#    writedlm(x, PP);
# end
# open("./ABCD_Systems/Outcomes/Veff_sB.csv","w+") do x
#    writedlm(x, Veff);
# end
# open("./ABCD_Systems/Outcomes/Fmean_sB.csv","w+") do x
#    writedlm(x, Fmean);
# end
# open("./ABCD_Systems/Outcomes/Fstd_sB.csv","w+") do x
#    writedlm(x, Fstd);
# end
# open("./ABCD_Systems/Outcomes/Coeff_MF_SplineI_sB.csv","w+") do x
#    writedlm(x, FdataMF);
# end
# open("./ABCD_Systems/Outcomes/Coeff_FF_SplineI_sB.csv","w+") do x
#    writedlm(x, FdataFF);
# end


######################################################################################################
# system C
# ( [0,1,0], SM (mass of particles 0) = 10, LM (mass of particles 1) = 1, SEps=1.0, LEps=100.0, γ=1, β=1.096 )

C = Chain(SM=0.1, LM=1, Beta=1.096);
(RR,PP,Veff,Fmean,Fstd,FdataMF,FdataFF) = VEffConstrainedLangDyn_NVT(C,S);

# Write data to file
open("./ABCD_Systems/Outcomes/RR_sC.csv","w+") do x
   writedlm(x,RR);
end
open("./ABCD_Systems/Outcomes/PP_sC.csv","w+") do x
   writedlm(x,PP);
end
open("./ABCD_Systems/Outcomes/Veff_sC.csv","w+") do x
   writedlm(x,Veff);
end
open("./ABCD_Systems/Outcomes/Fmean_sC.csv","w+") do x
   writedlm(x,Fmean);
end
open("./ABCD_Systems/Outcomes/Fstd_sC.csv","w+") do x
   writedlm(x,Fstd);
end
open("./ABCD_Systems/Outcomes/Coeff_MF_SplineI_sC_M0.1m1.csv","w+") do x
   writedlm(x, FdataMF);
end
open("./ABCD_Systems/Outcomes/Coeff_FF_SplineI_sC_M0.1m1.csv","w+") do x
   writedlm(x, FdataFF);
end



#################################################################################
# system D
# ( [0,1,0], SM (mass of particles 0) = 10, LM (mass of particles 1) = 1, γ=1, β=1 )

# C = Chain(SM=20, LM=1, SEps=10, LEps=.1, NBeads = 10);
# (RR,PP,Veff,Fmean,Fstd) = VEffConstrainedLangDyn_NVT(C,S);
#
# # Write data to file
# open("./ABCD_Systems/Outcomes/RR_sD.csv","w+") do x
#    writedlm(x,RR);
# end
# open("./ABCD_Systems/Outcomes/PP_sD.csv","w+") do x
#    writedlm(x,PP);
# end
# open("./ABCD_Systems/Outcomes/Veff_sD.csv","w+") do x
#    writedlm(x,Veff);
# end
# open("./ABCD_Systems/Outcomes/Fmean_sD.csv","w+") do x
#    writedlm(x,Fmean);
# end
# open("./ABCD_Systems/Outcomes/Fstd_sD.csv","w+") do x
#    writedlm(x,Fstd);
# end


#################################################################################
end # @time
