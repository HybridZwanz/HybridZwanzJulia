include("../LJChain.include")
jl("../CGStats.jl")

using Main.LJChain,Main.CGStats,Plots,Statistics,DelimitedFiles

###################################################################################################
# Set up simulation parameters

# to change the number of samples to compute, i.e. OrthoSamples Value
# Oin = Int(1e3);
# Ofin = Int(1e6);
OrthoSamplesTot = [exp10(x) for x in 2:7];
nsamples = length(OrthoSamplesTot);

# to change the time step (dt=1e-4)
dtt = 1e-4;

# to change the gap between samples, i.e. OrthoSteps value (OrthoSteps=1)
Os = Int(1e1);

# to change the number of FullSteps to compute full dynamics (FullSteps=Int(1e5))
Fs = Int(1e5);


C = Chain(Beta=1,Gamma=1);
S = Simulator(dt=dtt, FullSteps=Fs);    # NB: if you change the parameters, add new information about the simulated system to the titles of the plots

###################################################################################################
###################################################################################################
# to compute the orthogonal dynamic with the constrained dynamic by varying
# the number of samples

# Initial condition
(P0,Q0,F0,H0,PCG0,QCG0) = InitConds(C);

(Pstore,Qstore,Fstore) = FullTraj(P0,Q0,C,S);
Portho = Pstore[:,end];
Qortho = Qstore[:,end];

# to compute the distance between centres of mass of beads
RR = CGDists(Qortho,C)[1]; # NB: the function CgDists(Q,C) returns an array

println("Computing orthogonal dynamics through constrained dynamics ...");

# OrthoSamplesTot = Vector(Oin:Oin:Ofin);
diff = circshift(circshift(OrthoSamplesTot,-1)-OrthoSamplesTot,1);
diff[1] = OrthoSamplesTot[1];
StM     = zeros(nsamples,2);
VarStM  = zeros(nsamples);

for i = 1:nsamples
	S = Simulator(OrthoSamples=diff[i],OrthoSteps=Os,dt=dtt);
	# S = Simulator(OrthoSamples=OrthoSamplesTot[i],OrthoSteps=Os,dt=dtt);

	println("  Start sample ",i,": Orthosample=",OrthoSamplesTot[i],".");
	(StressMoments,Portho,Qortho) = ConstrainedSample(Portho,Qortho,C,S);

	if i == 1
		StM[1,:] = StressMoments;
	else
		StM[i,1] = (OrthoSamplesTot[i-1]*StM[i-1,1]+diff[i]*StressMoments[1,1])/OrthoSamplesTot[i];
		StM[i,2] = (OrthoSamplesTot[i-1]*StM[i-1,2]+diff[i]*StressMoments[1,2])/OrthoSamplesTot[i];
	end
	VarStM[i]  = StM[i,2]-(StM[i,1])^2;

	global Portho = Portho;
	global Qortho = Qortho;
end

println(" ... done!");
println("Distance between centre of mass of beads, RR = ", RR);

# Write data to file
dataDet = zeros(nsamples,5);
dataDet[1,1] = RR;
dataDet[:,2] = OrthoSamplesTot;
dataDet[:,3] = StM[:,1];
dataDet[:,4] = StM[:,2];
dataDet[:,5] = VarStM;
open("./OrthogonalDynamics_DetAndStoc/data_sC_DET.csv","w+") do x
   writedlm(x,dataDet);
end

# s1 = scatter(OrthoSamplesTot,VarStM,
#          legend = false,
#          markershape = :auto,
#          markercolor = :red,
#          markersize = 2.5,
# 		   tickfontsize = 5,
# 		   guidefontsize = 8,
#          xscale = :log10,
#          xlabel = "OrthoSamples",
#          yscale = :log10,
#          ylabel = "VarStress",
#          title  = "VarStress vs OrthoSamples
# 		   (constraint dynamics)",
# 		   titlefontsize = 7);
# savefig("../Outcomes/OrthogonalDynamics_DetAndStoc/plotVarStressVSOrthoSamples_1.pdf");

#########################################################################################################################
#########################################################################################################################
# to compute the orthogonal dynamics with the Langevin constrained dynamics by varying
# the number of samples

# Initial condition
P0 = zeros(C.NPart);
Q0 = Vector(1.0:1.0:30.0);

(Pstore,Qstore,Fstore,Kmean,Hmean) = LangevinFullTraj(P0,Q0,C,S);
Portho = Pstore[:,end];
Qortho = Qstore[:,end];

# to compute the distance between centres of mass of beads
RR = CGDists(Qortho,C)[1]; # NB: the function CgDists(Q,C) returns an array

println("Computing orthogonal dynamics through Langevin constraint dynamics ...");

# OrthoSamplesTot = Vector(Oin:Oin:Ofin);
diff = circshift(circshift(OrthoSamplesTot,-1)-OrthoSamplesTot,1);
diff[1] = OrthoSamplesTot[1];
StM     = zeros(nsamples,2);
VarStM  = zeros(nsamples);

for i = 1:nsamples

    S = Simulator(OrthoSamples=diff[i],OrthoSteps=Os,dt=dtt);
	# S = Simulator(OrthoSamples=OrthoSamplesTot[i],OrthoSteps=Os,dt=dtt);

	println("  Start sample ",i,": Orthosample=",OrthoSamplesTot[i],".");
	(StressMoments,Portho,Qortho) = ConstrSampleLangDyn(Portho,Qortho,C,S);

	if i == 1
		StM[1,:] = StressMoments;
	else
		StM[i,1] = (OrthoSamplesTot[i-1]*StM[i-1,1]+diff[i]*StressMoments[1,1])/OrthoSamplesTot[i];
		StM[i,2] = (OrthoSamplesTot[i-1]*StM[i-1,2]+diff[i]*StressMoments[1,2])/OrthoSamplesTot[i];
	end
	VarStM[i]  = StM[i,2]-(StM[i,1])^2;

	global Portho = Portho;
	global Qortho = Qortho;
end

println(" ... done!");
println("Distance between centre of mass of beads, RR = ", RR);

# Write data to file
dataStoc = zeros(nsamples,5);
dataStoc[1,1] = RR;
dataStoc[:,2] = OrthoSamplesTot;
dataStoc[:,3] = StM[:,1];
dataStoc[:,4] = StM[:,2];
dataStoc[:,5] = VarStM;
open("./OrthogonalDynamics_DetAndStoc/data_sC_STOC.csv","w+") do x
   writedlm(x,dataStoc);
end

# s1 = scatter(OrthoSamplesTot,StM[:,1],
#          legend = false,
#          markershape = :auto,
#          markercolor = :red,
#          markersize = 2.5,
# 		   tickfontsize = 5,
# 		   guidefontsize = 8,
#          xscale = :log10,
#          xlabel = "OrthoSamples",
#          yscale = :log10,
#          ylabel = "First Stress Moment",
#          title  = "First stress moment vs OrthoSamples
# 		   (constraint Langevin dynamics)",
# 		   titlefontsize = 7)
#
# s2 = scatter(OrthoSamplesTot,VarStM,
#          legend = false,
#          markershape = :auto,
#          markercolor = :red,
#          markersize = 2.5,
# 		   tickfontsize = 5,
# 		   guidefontsize = 8,
#          xscale = :log10,
#          xlabel = "OrthoSamples",
#          yscale = :log10,
#          ylabel = "VarStress",
#          title  = "VarStress vs OrthoSamples
# 		   (constraint Langevin dynamics)",
# 		   titlefontsize = 7);

#################################################################################################
# # to plot the two graphs together
#
# l = @layout [a b];
# plot(s1, s2, layout = l);
# savefig("../Outcomes/OrthogonalDynamics_DetAndStoc/plotFirstStressMoment_VS_OrthoSamples_2.pdf");

################################################################################################
# Autocorrelation analysis to compute the autocorrelation function and the
# corrisponding integrated autocorrelation time τ

# Kmax = 600;                             # index autocorrelation function
# n    = length(OrthoSamplesTot);
# M    = zeros(2,1);
# P    = zeros(n-1,1);
# C    = zeros(Kmax,2);
# tau  = zeros(Kmax,2);
#
# for i in 1:2
# 	M[i,1] = (1/length(OrthoSamplesTot))*sum(StressMoments[1:end,i]);
#     for k in 1:Kmax
#         for r in 1:n-k
#             P[r,1] = (StressMoments[r+k,i]-M[i,1])*(StressMoments[r,i]-M[i,1]);
#         end
#         C[k,i] = sum(P[:,1])/(n-k);
#         tau[k,i] = (1/2)+(sum(C[1:k,i])/var(StressMoments[:,i]));
#     end
# end
#
# plot(1:1:Kmax,C[:,1]);   # plot the autocorrelation of the StressMoments[:,i] where i=1,2
# plot(1:1:Kmax,C[:,2]);
# plot(1:1:Kmax,tau[:,1]); # plot the autocorrelation time of the StressMoments[:,i] where i=1,2
# plot(1:1:Kmax,tau[:,2]);
