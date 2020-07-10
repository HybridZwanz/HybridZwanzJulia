################################################################################
# Module CGStats
################################################################################
module CGStats

# Include relevant modules
using Main.LJChain, Distributed, SharedArrays, Parameters, Statistics, Dates, SparseArrays, LinearAlgebra, Printf, DelimitedFiles

# Parameter object
export Simulator
# Functions to output raw trajectory data
export FullTraj_NVE, ConstrainedTraj_NVE, DetCGTraj_NVE,
       LangFullTraj_NVT, LangConstrainedTraj_NVT, DetCGTraj_NVT,
	   ApproximateTraj_MMZD
# High-level function to extract data
export VEffConstrained_NVE, VEffConstrainedLangDyn_NVT,
       FullAutocorr, ConstrainedAutocorr, FullAutocorrLangDyn, ConstrainedAutocorrLangDyn, ApproximateAutocorr,
	   FullAutocorrPar, ConstrainedAutocorrPar,  FullAutocorrLangDynPar, ConstrainedAutocorrLangDynPar, ApproximateAutocorrPar
# Orthogonal dynamics sampling
export ConstrainedSample, ConstrSampleLangDyn
# Functions to extract approximate dynamics
export InterpolationAlgorithm_Kriging, ApproxMeanForce, ApproxMeanStresses, ApproxFlucStep, ComputeVEff

#######################################
# Simulator type, encoding parameters #
#######################################
@with_kw struct Simulator
  dt::Float64		  = 1e-4
  Stochdt::Float64	  = 1e-4
  FullSteps::Int	  = Int(1e5)
  OrthoSteps::Int 	  = 1
  ApproxSteps::Int	  = Int(1e5)
  FullTime::Float64	  = FullSteps *dt
  ApproxTime::Float64 = ApproxSteps *Stochdt
  LagStep			  = 400       # spacing of lag times for autocorrelation
  NumLags			  = 30        # number of lag times to compute
  MacroSamples::Int	  = 50
  OrthoSamples::Int	  = Int(1e5)
  Moments::Int		  = 2
  check::Int          = Int(5e3)  # to print the number of iteration in a loop

  #NB: the below coeffcients must be checked and entered every time a new look-up table of Mean and Fluctuating Force is generated.
  #NB: the following values were generated for the look-up tables in the folder Julia_Code/Store_Trajectory for the system C (M-m-M) with masses M=10 and M=1
  #    consisting of 30 particles.

  # Coeff. of line (y = m*x + q) in the ApproxMeanStresses function
  m_MF::Float64       = 537.59;   # gradient of the line
  q_MF::Float64       = -1567.1;  #  y-intercept of the line

  # Coeff. of line (y = m*x + q) in the ApproxFlucStep function
  m_FF::Float64       = -135.21;  # gradient of the line
  q_FF::Float64       = 403.79;   #  y-intercept of the line

  # Coeff. of exponential function (y = A*exp(-B*x)) in the ApproxFlucStep function
  coeffEXP_A::Float64 = 3.3199e6;
  coeffEXP_B::Float64 = 4.8887;
end

#############################################
# Functions to generate raw trajectory data #
#############################################
"""
    (Pstore,Qstore,Fstore) = FullTraj_NVE(P0,Q0,C::Chain,S::Simulator)

    Outputs full trajectory data in matrices (Pstore,Qstore,Fstore).
"""
function FullTraj_NVE(P0,Q0,C::Chain,S::Simulator)
  # Report to user.
  println(Dates.now());
  println("Running full simulation and storing trajectory.");
  println("  NPart, number of particles simulated = ",C.NPart);
  println("  dt, time step size = ",S.dt," seconds.");
  println("  Total simulation time = ",S.FullTime," seconds.");

  # Set up.
  P = P0; Q = Q0; F = Forces(Q,C);
  Pstore = zeros(Float64,C.NPart,S.FullSteps);
  Qstore = zeros(Float64,C.NPart,S.FullSteps);
  Fstore = zeros(Float64,C.NPart,S.FullSteps);
  Pstore[:,1] = P0;
  Qstore[:,1] = Q0;
  Fstore[:,1] = F;

  # Do simulation.
  for i = 2:S.FullSteps
    (P,Q,F) = LeapfrogStep(P,Q,F,C,S.dt);
    Pstore[:,i] = P;
    Qstore[:,i] = Q;
    Fstore[:,i] = F;
  end

  return Pstore,Qstore,Fstore
end

"""
    (Pstore,Qstore,Fstore,Kmean,Hmean) = LangFullTraj_NVT(P0,Q0,C::Chain,S::Simulator)

    Outputs full trajectory data in matrices (Pstore,Qstore,Fstore) and the mean values of
	the kinetic energy and the total energy.
"""
function LangFullTraj_NVT(P0,Q0,C::Chain,S::Simulator)
  # Report to user.
  println(Dates.now());
  println("Running full simulation of Langevin dynamics, storing trajectory and computing");
  println("the mean values per particle of the kinetic energy and the total energy.");
  println("  NPart, number of particles simulated = ",C.NPart);
  println("  Beta, temperature parameter = ",C.Beta);
  println("  Gamma, friction coefficient = ",C.Gamma);
  println("  dt, time step size = ",S.dt);
  println("  Total simulation time = ",S.FullTime);

  # Set up.
  P = P0; Q = Q0; F = Forces(Q,C);
  K = KineticEnergy(P,C); H = Potential(Q,C) + K;

  r = Int(S.FullSteps/100) + 1;
  Pstore = zeros(Float64,C.NPart,r);
  Qstore = zeros(Float64,C.NPart,r);
  Fstore = zeros(Float64,C.NPart,r);
  Pstore[:,1] = P0;
  Qstore[:,1] = Q0;
  Fstore[:,1] = F;

  # Do simulation.
  j = 2;
  for i = 2:S.FullSteps
	  if i % S.check == 0
		  @printf "    ITERATION %0.3e OF %0.3e\n" i S.FullSteps
	  end
	  (P,Q,F) = LangDynIntegrators(P,Q,F,C,S.dt);
	  if i % 100 == 0
		  Pstore[:,j] = P;
	      Qstore[:,j] = Q;
	      Fstore[:,j] = F;
		  j += 1;
	  end
	  K += KineticEnergy(P,C);
	  H += Potential(Q,C) + K;
  end

  K /= S.FullSteps;
  H /= S.FullSteps;
  Kmean = sum(K)/C.NPart;
  Hmean = sum(H)/C.NPart;

  return Pstore,Qstore,Fstore,Kmean,Hmean
end

"""
    (Pstore,Qstore,Fstore) = ConstrainedTraj_NVE(P0,Q0,C::Chain,S::Simulator)

    Outputs full trajectory data in matrices (Pstore,Qstore,Fstore).
"""
function ConstrainedTraj_NVE(P0,Q0,C::Chain,S::Simulator)
  # Report to user.
  println(now());
  println("Running full simulation and storing trajectory.");
  println("  NPart, number of particles simulated = ",C.NPart);
  println("  dt, time step size = ",S.dt," seconds.");
  println("  Total simulation time = ",S.FullTime," seconds.");

  # Set up.
  P = P0; Q = Q0; F = Forces(Q,C);
  Pstore = zeros(Float64,C.NPart,S.FullSteps);
  Qstore = zeros(Float64,C.NPart,S.FullSteps);
  Fstore = zeros(Float64,C.NPart,S.FullSteps);
  Pstore[:,1] = P0;
  Qstore[:,1] = Q0;
  Fstore[:,1] = F;

  # Do simulation.
  for i = 2:S.FullSteps
    (P,Q,F) = ConstrainedStep(P,Q,F,C,S.dt);
    Pstore[:,i] = P;
    Qstore[:,i] = Q;
    Fstore[:,i] = F;
  end

  return Pstore,Qstore,Fstore
end

"""
    (Pstore,Qstore,Fstore) = LangConstrainedTraj_NVT(P0,Q0,C::Chain,S::Simulator)

    Outputs full trajectory data in matrices (Pstore,Qstore,Fstore).
"""
function LangConstrainedTraj_NVT(P0,Q0,C::Chain,S::Simulator)
  # Report to user.
  println(Dates.now());
  println("Running full simulation of constrained Langevin dynamics and storing trajectory");
  println("  NPart, number of particles simulated = ",C.NPart);
  println("  Beta, temperature parameter = ",C.Beta);
  println("  Gamma, friction coefficient = ",C.Gamma);
  println("  dt, time step size = ",S.dt," seconds.");
  println("  Total simulation time = ",S.FullTime," seconds.");

  # Set up.
  P = C.MomProj*P0; Q = Q0; F = Forces(Q,C);

  r = Int(S.FullSteps/100) + 1;
  Pstore = zeros(Float64,C.NPart,r);
  Qstore = zeros(Float64,C.NPart,r);
  Fstore = zeros(Float64,C.NPart,r);
  Pstore[:,1] = P;
  Qstore[:,1] = Q;
  Fstore[:,1] = F;

  # Do simulation.
  j = 2;
  for i = 2:S.FullSteps
	  if i % S.check == 0
		  @printf "    ITERATION %0.3e OF %0.3e\n" i S.FullSteps
	  end
	  (P,Q,F) = LangDynConstrainedStep(P,Q,F,C,S.dt);
	  if i % 100 == 0
		  Pstore[:,j] = P;
		  Qstore[:,j] = Q;
		  Fstore[:,j] = F;
		  j += 1;
	  end
  end

  return Pstore,Qstore,Fstore
end

"""
    (Pstore,Qstore,Fstore) = ApproximateTraj_MMZD(P0,Q0,FmeanData,FstdData,C::Chain,S::Simulator)

	Outputs full trajectory data for approximate Langevin dynamics, using approximate potentials generated by the data encoded
	in FmeanData and FstdData. Where this data are the coefficients for the cubic spline interpolation. These coefficients are
	generated from kriging interpolation data, generated separately to eliminate noise from the original Mean Force and Fluctuating
	force data.
"""

function ApproximateTraj_MMZD(P0,Q0,FmeanData::Array{Float64,2},FstdData::Array{Float64,2},C::Chain,S::Simulator)
 # Report to user.
 println(Dates.now());
 println("Running full simulation of MMZD and storing trajectory.");
 println("  NPart, number of particles simulated = ",C.NPart);
 println("  Beta, temperature parameter = ",C.Beta);
 println("  Gamma, friction coefficient = ",C.Gamma);
 println("  dt, time step size = ",S.Stochdt);
 println("  Total simulation time = ",S.ApproxTime);

  # Set up initial conditions
  P = P0; Q = Q0;

  # Generate mean force function
  MeanForce = X -> ApproxMeanForce(X,FmeanData,C,S)[1];
  F = MeanForce(Q0);

  # Set up approximate fluctuation generator
  FlucStep = (X,Y) -> ApproxFlucStep(X,Y,FstdData,C,S)[1];

  # Set up storage
  r = Int(S.ApproxSteps/100) + 1;
  Pstore = zeros(r,length(P0));
  Qstore = zeros(r,length(P0));
  Fstore = zeros(r,length(P0));
  Pstore[1,:] = P;
  Qstore[1,:] = Q;
  Fstore[1,:] = F;

  # Do simulation.
  j =2;
  for i = 2:S.ApproxSteps
	if i % S.check == 0
		@printf "    ITERATION %0.3e OF %0.3e\n" i S.ApproxSteps
	end
    (P,Q,F) = ApproxDynStep(P,Q,F,MeanForce,FlucStep,C,S.Stochdt);
	if i % 100 == 0
		Pstore[j,:] = P;
	    Qstore[j,:] = Q;
	    Fstore[j,:] = F;
		# if j % 10 == 0
		# 	open("Store_Trajectory/Outcomes/Q_MMZD_1.csv", "w+") do x
		# 		# @printf(x, "%s\n", map(y -> @sprintf(%f, y), Q) )
		# 		writedlm(x, Qstore);
		# 	end
		# 	open("Store_Trajectory/Outcomes/P_MMZD_1.csv", "w+") do x
		# 		writedlm(x, Pstore);
		# 	end
		# end
		j += 1;
	end
  end

  return Pstore, Qstore, Fstore
end

"""
    (Pstore,Qstore,Fstore) = DetCGTraj_NVE(P0,Q0,FmeanData,C::Chain,S::Simulator)

    Outputs full trajectory data for Deterministi CG dynamics in the Microcanonical ensemble, using approximate potential generated
    by the data encoded in FmeanData. Where Fmeandata are the coefficients for the cubic spline interpolation. These coefficients
	are generated from kriging interpolation data, generated separately to eliminate noise from the original Mean Force data.
	We do all this to speed up the DCGD, because the cubic spline interpolation is much faster than the Kriging interpolation.
"""
function DetCGTraj_NVE(P0,Q0,FmeanData::Array{Float64,2},C::Chain,S::Simulator)
 # Report to user.
 println(Dates.now());
 println("Running full simulation of DCGD in NVE ensemble and storing trajectory.");
 println("  NPart, number of particles simulated = ",C.NPart);
 println("  dt, time step size = ",S.dt);
 println("  Total simulation time = ",S.FullTime);

  # Set up initial conditions
  P = P0; Q = Q0;

  # Generate mean force function
  MeanForce = X -> ApproxMeanForce(X,FmeanData,C,S)[1];
  F = MeanForce(Q0);

  # Set up storage
  r = Int(S.FullSteps/100) + 1;
  Pstore = zeros(r,length(P0));
  Qstore = zeros(r,length(P0));
  Fstore = zeros(r,length(P0));
  Pstore[1,:] = P;
  Qstore[1,:] = Q;
  Fstore[1,:] = F;

  # Do simulation.
  j = 2;
  for i = 2:S.FullSteps
	  if i % S.check == 0
		  @printf "    ITERATION %0.2e OF %0.2e\n" i S.FullSteps
	  end
	  (P,Q,F) = DetCGDynNVEStep(P,Q,F,MeanForce,C,S.dt);
	  if i % 100 == 0
		  Pstore[j,:] = P;
		  Qstore[j,:] = Q;
		  Fstore[j,:] = F;
		  j += 1;
	  end
  end

  return Pstore,Qstore,Fstore
end

"""
    (Pstore,Qstore,Fstore) = DetCGTraj_NVT(P0,Q0,FmeanData,C::Chain,S::Simulator)

    Outputs full trajectory data for Deterministi CG dynamics in the canonical ensemble, using approximate potential generated
    by the data encoded in FmeanData. Where Fmeandata are the coefficients for the cubic spline interpolation. These coefficients
	are generated from kriging interpolation data, generated separately to eliminate noise from the original Mean Force data.
	We do all this to speed up the DCGD, because the cubic spline interpolation is much faster than the Kriging interpolation.
"""
function DetCGTraj_NVT(P0,Q0,FmeanData::Array{Float64,2},C::Chain,S::Simulator)
 # Report to user.
 println(Dates.now());
 println("Running full simulation of DCGD in NVT ensemble and storing trajectory.");
 println("  NPart, number of particles simulated = ",C.NPart);
 println("  Beta, temperature parameter = ",C.Beta);
 println("  Gamma, friction coefficient = ",C.Gamma);
 println("  dt, time step size = ",S.dt);
 println("  Total simulation time = ",S.FullTime);

  # Set up initial conditions
  P = P0; Q = Q0;

  # Generate mean force function
  MeanForce = X -> ApproxMeanForce(X,FmeanData,C,S)[1];
  F = MeanForce(Q0);

  # Set up storage
  r = Int(S.FullSteps/100) + 1;
  Pstore = zeros(r,length(P0));
  Qstore = zeros(r,length(P0));
  Fstore = zeros(r,length(P0));
  Pstore[1,:] = P;
  Qstore[1,:] = Q;
  Fstore[1,:] = F;

  # Do simulation.
  j = 2;
  for i = 2:S.FullSteps
	  if i % S.check == 0
		  @printf "    ITERATION %0.2e OF %0.2e\n" i S.FullSteps
	  end
	  (P,Q,F) = DetCGDynNVTStep(P,Q,F,MeanForce,C,S.dt);
	  if i % 100 == 0
		  Pstore[j,:] = P;
		  Qstore[j,:] = Q;
		  Fstore[j,:] = F;
		  j += 1;
	  end
  end

  return Pstore,Qstore,Fstore
end

################################
# Orthogonal dynamics sampling #
################################
function ConstrainedSample(P::Array{Float64,1},Q::Array{Float64,1},C::Chain,S::Simulator,Obs=(P,Q,C) -> 0)

  # Mean Forces
  StressMoments = zeros(C.ConstrainedBeads-1,S.Moments);

  #PCG = C.PhiConstr*P;
  #QCG = C.PsiConstr*Q;
  #E0  = Potential(Q,C);
  F   = Forces(Q,C);
  #H0  = E0 + KineticEnergy(P,C);

  s = 0;
  ObsMean = 0;
  while s < S.OrthoSamples
    # Run simulation to generate point with same QCG. (Constrained dynamics)
    for i = 1:S.OrthoSteps
      (P,Q,F) = ConstrainedStep(P,Q,F,C,S.dt);
    end

    # Record successful sample.
    s+=1;
    SigmaCG = C.PhiFConstr*Stresses(Q,C);
    for j = 1:S.Moments
      StressMoments[:,j] += SigmaCG.^j;
    end
    # ObsMean += Obs(P,Q,C); ## TODO FIX THIS CALL
  end

  # ObsMean /= S.OrthoSamples;
  # Efinal = Potential(Q,C);
  # Hfinal = Efinal + KineticEnergy(P,C);

  #println("        Error in momentum constraint = ",maximum(abs.(PCG-C.PhiConstr*P)));
  #println("        Error in position constraint = ",maximum(abs.(QCG-C.PsiConstr*Q)));
  #println("        Relative error in energy     = ",abs(Hfinal/H0-1));
  #if Obs!=0
  #  println("        Value of observable          = ",ObsMean);
  #end

  StressMoments ./= S.OrthoSamples
  return StressMoments,P,Q
end

function ConstrSampleLangDyn(P::Array{Float64,1},Q::Array{Float64,1},C::Chain,S::Simulator,Obs=(P,Q,C) -> 0)

  # Mean Forces
  StressMoments = zeros(C.ConstrainedBeads-1,S.Moments);

  #PCG = C.PhiConstr*P;
  #QCG = C.PsiConstr*Q;
  #E0  = Potential(Q,C);
  F   = Forces(Q,C);
  #H0  = E0 + KineticEnergy(P,C);

  s = 0;
  ObsMean = 0;
  while s < S.OrthoSamples
    # Run simulation to generate point with same QCG. (Constrained dynamics)
    for i = 1:S.OrthoSteps
      (P,Q,F) = LangDynConstrainedStep(P,Q,F,C,S.dt);
    end

    # Record successful sample.
    s+=1;
    SigmaCG = C.PhiFConstr*Stresses(Q,C);
    for j = 1:S.Moments
      StressMoments[:,j] += SigmaCG.^j;
    end
    # ObsMean += Obs(P,Q,C); ## TODO FIX THIS CALL
  end

  # ObsMean /= S.OrthoSamples;
  # Efinal = Potential(Q,C);
  # Hfinal = Efinal + KineticEnergy(P,C);

  # println("        Error in momentum constraint = ",maximum(abs.(PCG-C.PhiConstr*P)));
  # println("        Error in position constraint = ",maximum(abs.(QCG-C.PsiConstr*Q)));
  # println("        Relative error in energy     = ",abs(Hfinal/H0-1));
  # if Obs!=0
   # println("        Value of observable          = ",ObsMean);
  # end

  StressMoments ./= S.OrthoSamples
  return StressMoments,P,Q
end

function ConstrainedAutocorr(P::Array{Float64,1},Q::Array{Float64,1},C::Chain,S::Simulator)

  # Mean Forces
  StressStats 	= zeros(S.NumLags,C.ConstrainedBeads-1);
  StressStorage = zeros(S.NumLags,C.ConstrainedBeads-1);
  StressMean   	= zeros(C.ConstrainedBeads-1);
  StressSqMean 	= zeros(C.ConstrainedBeads-1);
  MomStats    	= zeros(S.NumLags,C.ConstrainedBeads);
  MomStorage 	= zeros(S.NumLags,C.ConstrainedBeads);
  MomMean	    = zeros(C.ConstrainedBeads);
  MomSqMean 	= zeros(C.ConstrainedBeads);

  E0  = Potential(Q,C);
  F   = Forces(Q,C);
  H0  = E0 + KineticEnergy(P,C);

  # first run the dynamics for one window
    m = 0;
    while m < S.NumLags
        for i = 1:S.LagStep
          (P,Q,F) = ConstrainedStep(P,Q,F,C,S.dt);
        end

        m+=1;

        SigmaCG = C.PhiFConstr*Stresses(Q,C);
        StressStorage = circshift(StressStorage,(1,0));
        StressStorage[1,:] = SigmaCG;

        PCG = C.PhiConstr*P;
        MomStorage = circshift(MomStorage,(1,0));
        MomStorage[1,:] = PCG;
    end

  # Now start sampling
  s=0;

  while s < S.OrthoSamples
    # Advance to next sample time
    for i = 1:S.LagStep
      (P,Q,F) = ConstrainedStep(P,Q,F,C,S.dt);
    end

    # Record successful sample.
    s+=1;

    # Update stress storage
    SigmaCG = C.PhiFConstr*Stresses(Q,C);
    StressStorage = circshift(StressStorage,(1,0));
    StressStorage[1,:] = SigmaCG;

    # Update moment estimators
    StressMean += SigmaCG;
    StressSqMean += SigmaCG.^2;

    # Update momentum storage
    PCG = C.PhiConstr*P;
    MomStorage = circshift(MomStorage,(1,0));
    MomStorage[1,:] = PCG;

    # Update moment estimators
    MomMean += PCG;
    MomSqMean += PCG.^2;

    # Add to estimator of expectation of product
    StressStats += StressStorage.*kron(SigmaCG',ones(S.NumLags));
    MomStats += MomStorage.*kron(PCG',ones(S.NumLags));
  end

  StressMean /= S.OrthoSamples;
  StressSqMean /= S.OrthoSamples;
  MomMean /= S.OrthoSamples;
  MomSqMean /= S.OrthoSamples;

  StressStats /= S.OrthoSamples;
  MomStats /= S.OrthoSamples;

  # Efinal = Potential(Q,C);
  # Hfinal = Efinal + KineticEnergy(P,C);

  # println("        Error in momentum constraint = ",maximum(abs(PCG-C.PhiConstr*P)));
  # println("        Error in position constraint = ",maximum(abs(QCG-C.PsiConstr*Q)));
  # println("        Relative error in energy     = ",abs(Hfinal/H0-1));
  # if Obs!=0
   # println("        Value of observable          = ",ObsMean);
  # end

  return StressStats,StressMean,StressSqMean,MomStats,MomMean,MomSqMean
end

function ConstrainedAutocorrLangDyn(P::Array{Float64,1},Q::Array{Float64,1},C::Chain,S::Simulator)

  # Mean Forces
  StressStats 	= zeros(S.NumLags,C.ConstrainedBeads-1);
  StressStorage = zeros(S.NumLags,C.ConstrainedBeads-1);
  StressMean   	= zeros(C.ConstrainedBeads-1);
  StressSqMean 	= zeros(C.ConstrainedBeads-1);
  MomStats    	= zeros(S.NumLags,C.ConstrainedBeads);
  MomStorage 	= zeros(S.NumLags,C.ConstrainedBeads);
  MomMean	    = zeros(C.ConstrainedBeads);
  MomSqMean 	= zeros(C.ConstrainedBeads);

  # E0  = Potential(Q,C);
  F   = Forces(Q,C);
  # H0  = E0 + KineticEnergy(P,C);

  # first run the dynamics for one window
    m = 0;
    while m < S.NumLags
        for i = 1:S.LagStep
          (P,Q,F) = LangDynConstrainedStep(P,Q,F,C,S.dt);
        end

        m+=1;

        SigmaCG = C.PhiFConstr*Stresses(Q,C);
        StressStorage = circshift(StressStorage,(1,0));
        StressStorage[1,:] = SigmaCG;

        PCG = C.PhiConstr*P;
        MomStorage = circshift(MomStorage,(1,0));
        MomStorage[1,:] = PCG;
    end

  # Now start sampling
  s=0;

  while s < S.OrthoSamples
    # Advance to next sample time
    for i = 1:S.LagStep
      (P,Q,F) = LangDynConstrainedStep(P,Q,F,C,S.dt);
    end

    # Record successful sample.
    s+=1;

    # Update stress storage
    SigmaCG = C.PhiFConstr*Stresses(Q,C);
    StressStorage = circshift(StressStorage,(1,0));
    StressStorage[1,:] = SigmaCG;

    # Update moment estimators
    StressMean += SigmaCG;
    StressSqMean += SigmaCG.^2;

    # Update momentum storage
    PCG = C.PhiConstr*P;
    MomStorage = circshift(MomStorage,(1,0));
    MomStorage[1,:] = PCG;

    # Update moment estimators
    MomMean += PCG;
    MomSqMean += PCG.^2;

    # Add to estimator of expectation of product
    StressStats += StressStorage.*kron(SigmaCG',ones(S.NumLags));
    MomStats += MomStorage.*kron(PCG',ones(S.NumLags));
  end

  StressMean /= S.OrthoSamples;
  StressSqMean /= S.OrthoSamples;
  MomMean /= S.OrthoSamples;
  MomSqMean /= S.OrthoSamples;

  StressStats /= S.OrthoSamples;
  MomStats /= S.OrthoSamples;

  # Efinal = Potential(Q,C);
  # Hfinal = Efinal + KineticEnergy(P,C);

  # println("        Error in momentum constraint = ",maximum(abs(PCG-C.PhiConstr*P)));
  # println("        Error in position constraint = ",maximum(abs(QCG-C.PsiConstr*Q)));
  # println("        Relative error in energy     = ",abs(Hfinal/H0-1));
  # if Obs!=0
   # println("        Value of observable          = ",ObsMean);
  # end

  return StressStats,StressMean,StressSqMean,MomStats,MomMean,MomSqMean
end

###################################
# Full Trajectory Autocorrelation #
###################################
function FullAutocorr(P::Array{Float64,1},Q::Array{Float64,1},C::Chain,S::Simulator,MeanForce)

  # Mean Forces
  StressStats 	= zeros(S.NumLags,C.NBeads);
  StressStorage = zeros(S.NumLags,C.NBeads);
  StressMean   	= zeros(C.NBeads);
  StressSqMean 	= zeros(C.NBeads);
  MomStats    	= zeros(S.NumLags,C.NBeads);
  MomStorage 	= zeros(S.NumLags,C.NBeads);
  MomMean	= zeros(C.NBeads);
  MomSqMean = zeros(C.NBeads);

  E0  = Potential(Q,C);
  F   = Forces(Q,C);
  H0  = E0 + KineticEnergy(P,C);

  # first run the dynamics for one window
    m = 0;
    while m < S.NumLags
        for i = 1:S.LagStep
          (P,Q,F) = LeapfrogStep(P,Q,F,C,S.dt);
        end

        m+=1;

        FCG = MeanForce(C.Psi*Q);
	    SigmaCG = C.PhiF*Stresses(Q,C)-FCG;
	    StressStorage = circshift(StressStorage,(1,0));
        StressStorage[1,:] = SigmaCG;

        PCG = C.Phi*P;
        MomStorage = circshift(MomStorage,(1,0));
        MomStorage[1,:] = PCG;
    end


    # Now start sampling
    s=0;

    while s < S.OrthoSamples
      # Advance to next sample time
      for i = 1:S.LagStep
        (P,Q,F) = LeapfrogStep(P,Q,F,C,S.dt);
      end

      # Record successful sample.
      s+=1;

      # Update stress storage
      FCG = MeanForce(C.Psi*Q);
      SigmaCG = C.PhiF*Stresses(Q,C)-FCG;
      StressStorage = circshift(StressStorage,(1,0));
      StressStorage[1,:] = SigmaCG;

      # Update stress estimators
      StressMean += SigmaCG;
      StressSqMean += SigmaCG.^2

      # Update momentum storage
      PCG = C.Phi*P;
      MomStorage = circshift(MomStorage,(1,0));
      MomStorage[1,:] = PCG;

      # Update moment estimators
      MomMean += PCG;
      MomSqMean += PCG.^2

      # Add to estimator of expectation of product
      StressStats += StressStorage.*kron(SigmaCG',ones(S.NumLags));
      MomStats += MomStorage.*kron(PCG',ones(S.NumLags));
    end

    StressMean /= S.OrthoSamples;
    StressSqMean /= S.OrthoSamples;
    MomMean /= S.OrthoSamples;
    MomSqMean /= S.OrthoSamples;

    StressStats /= S.OrthoSamples;
    MomStats /= S.OrthoSamples;

    # Efinal = Potential(Q,C);
    # Hfinal = Efinal + KineticEnergy(P,C);

    # println("        Error in momentum constraint = ",maximum(abs(PCG-C.PhiConstr*P)));
    # println("        Error in position constraint = ",maximum(abs(QCG-C.PsiConstr*Q)));
    # println("        Relative error in energy     = ",abs(Hfinal/H0-1));
    # if Obs!=0
     # println("        Value of observable          = ",ObsMean);
    # end

    return StressStats,StressMean,StressSqMean,MomStats,MomMean,MomSqMean
end

function FullAutocorrLangDyn(P::Array{Float64,1},Q::Array{Float64,1},C::Chain,S::Simulator,MeanForce)

  # Mean Forces
  StressStats 	= zeros(S.NumLags,C.NBeads);
  StressStorage = zeros(S.NumLags,C.NBeads);
  StressMean   	= zeros(C.NBeads);
  StressSqMean 	= zeros(C.NBeads);
  MomStats    	= zeros(S.NumLags,C.NBeads);
  MomStorage 	= zeros(S.NumLags,C.NBeads);
  MomMean	= zeros(C.NBeads);
  MomSqMean = zeros(C.NBeads);

  # E0  = Potential(Q,C);
  F   = Forces(Q,C);
  # H0  = E0 + KineticEnergy(P,C);

  # first run the dynamics for one window
    m = 0;
    while m < S.NumLags
        for i = 1:S.LagStep
          (P,Q,F) = LangDynIntegrators(P,Q,F,C,S.dt);
        end

        m+=1;

        FCG = MeanForce(C.Psi*Q);
        SigmaCG = C.PhiF*Stresses(Q,C)-FCG;
		# SigmaCG = C.PhiF*Stresses(Q,C);
        StressStorage = circshift(StressStorage,(1,0));
        StressStorage[1,:] = SigmaCG;

        PCG = C.Phi*P;
        MomStorage = circshift(MomStorage,(1,0));
        MomStorage[1,:] = PCG;
    end


    ### Now start sampling
    s=0;

    while s < S.OrthoSamples
      # Advance to next sample time
      for i = 1:S.LagStep
        (P,Q,F) = LangDynIntegrators(P,Q,F,C,S.dt);
      end

      # Record successful sample.
      s+=1;

      # Update STRESS storage
      FCG = MeanForce(C.Psi*Q);
      SigmaCG = C.PhiF*Stresses(Q,C)-FCG;
	  # SigmaCG = C.PhiF*Stresses(Q,C);
      StressStorage = circshift(StressStorage,(1,0));
      StressStorage[1,:] = SigmaCG;

      # Update STRESS estimators
      StressMean += SigmaCG;
      StressSqMean += SigmaCG.^2;

      # Update MOMENTUM storage
      PCG = C.Phi*P;
      MomStorage = circshift(MomStorage,(1,0));
      MomStorage[1,:] = PCG;

      # Update MOMENTUM estimators
      MomMean += PCG;
      MomSqMean += PCG.^2;

      # Add to estimator of expectation of product
      StressStats += StressStorage.*kron(SigmaCG',ones(S.NumLags));
      MomStats += MomStorage.*kron(PCG',ones(S.NumLags));
    end

    StressMean /= S.OrthoSamples;
    StressSqMean /= S.OrthoSamples;
    MomMean /= S.OrthoSamples;
    MomSqMean /= S.OrthoSamples;

    StressStats /= S.OrthoSamples;
    MomStats /= S.OrthoSamples;

    # Efinal = Potential(Q,C);
    # Hfinal = Efinal + KineticEnergy(P,C);

    # println("        Error in momentum constraint = ",maximum(abs(PCG-C.PhiConstr*P)));
    # println("        Error in position constraint = ",maximum(abs(QCG-C.PsiConstr*Q)));
    # println("        Relative error in energy     = ",abs(Hfinal/H0-1));
    # if Obs!=0
     # println("        Value of observable          = ",ObsMean);
    # end

    return StressStats,StressMean,StressSqMean,MomStats,MomMean,MomSqMean
end

function ApproximateAutocorr(P::Array{Float64,1},Q::Array{Float64,1},
	                           C::Chain,S::Simulator,FmeanData::Array{Float64,2},FstdData::Array{Float64,2})
  # Generate mean force function
  MeanForce = X -> ApproxMeanForce(X,FmeanData,C,S)[1];
  F = MeanForce(Q);

  # Set up approximate fluctuation generator
  FlucStep = (X,Y) -> ApproxFlucStep(X,Y,FstdData,C,S)[1];

  # Mean Forces
  StressStats 	= zeros(S.NumLags,C.NBeads);
  StressStorage = zeros(S.NumLags,C.NBeads);
  StressMean   	= zeros(C.NBeads);
  StressSqMean 	= zeros(C.NBeads);
  MomStats    	= zeros(S.NumLags,C.NBeads);
  MomStorage 	= zeros(S.NumLags,C.NBeads);
  MomMean       = zeros(C.NBeads);
  MomSqMean 	= zeros(C.NBeads);



  # first run the dynamics for one window
    m = 0;
    while m < S.NumLags
        for i = 1:S.LagStep
          (P,Q,F) = ApproxDynStep(P,Q,F,MeanForce,FlucStep,C,S.Stochdt);
        end

        m+=1;

        SigmaCG = circshift(F,-1)-F;
        StressStorage = circshift(StressStorage,(1,0));
        StressStorage[1,:] = SigmaCG;

        PCG = P;
        MomStorage = circshift(MomStorage,(1,0));
        MomStorage[1,:] = PCG;
    end

  # Now start sampling
  s=0;

  while s < S.OrthoSamples
    # Advance to next sample time
    for i = 1:S.LagStep
      (P,Q,F) = ApproxDynStep(P,Q,F,MeanForce,FlucStep,C,S.Stochdt);
    end

    # Record successful sample.
    s+=1;

    # Update stress storage
    SigmaCG = circshift(F,-1)-F;
    StressStorage = circshift(StressStorage,(1,0));
    StressStorage[1,:] = SigmaCG;

    # Update moment estimators
    StressMean += SigmaCG;
    StressSqMean += SigmaCG.^2;

    # Update momentum storage
    PCG = P;
    MomStorage = circshift(MomStorage,(1,0));
    MomStorage[1,:] = PCG;

    # Update moment estimators
    MomMean += PCG;
    MomSqMean += PCG.^2

    # Add to estimator of expectation of product
    StressStats += StressStorage.*kron(SigmaCG',ones(S.NumLags));
    MomStats += MomStorage.*kron(PCG',ones(S.NumLags));
  end

  StressMean /= S.OrthoSamples;
  StressSqMean /= S.OrthoSamples;
  MomMean /= S.OrthoSamples;
  MomSqMean /= S.OrthoSamples;

  StressStats /= S.OrthoSamples;
  MomStats /= S.OrthoSamples;

  # Efinal = Potential(Q,C);
  # Hfinal = Efinal + KineticEnergy(P,C);

  # println("        Error in momentum constraint = ",maximum(abs(PCG-C.PhiConstr*P)));
  # println("        Error in position constraint = ",maximum(abs(QCG-C.PsiConstr*Q)));
  # println("        Relative error in energy     = ",abs(Hfinal/H0-1));
  # if Obs!=0
	  # println("        Value of observable          = ",ObsMean);
  # end

  return StressStats,StressMean,StressSqMean,MomStats,MomMean,MomSqMean
end

####################################
# Functions to compute observables #
####################################
"""
    Observables(P,Q,C::Chain)

"""
function Observables(P,Q,C::Chain)

  Obs = (circshift(Q,-1) - Q) - C.Dist;
  Obs[end] += C.Period;

  return [mean(Obs),dot(Obs,Obs)/C.NPart]
end

"""
    CGObservables(P,Q,C::Chain)

"""
function CGObservables(P,Q,C::Chain)

  Obs = (circshift(C.Psi*Q,-1) - C.Psi*Q) - C.BeadsDist;
  Obs[end] += C.Period;

  return [mean(Obs),dot(Obs,Obs)/C.NBeads]
end

###################################
# High-level simulation functions #
###################################
function VEffConstrained_NVE(C::Chain,S::Simulator)
  println("Computing coarse graining through constrained dynamics ...");
  if S.dt > 0.1*sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt));
    warn("Timestep = ",S.dt,
      " not ≪  shortest timescale = ",sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt)));
  end

  # Generate initial condition
  (P,Q,F,H,PCG,QCG) = InitConds(C);

  # Set up storage
  RSamples 	= SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads-1);
  PSamples 	= SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads-1);
  MomentSamples = SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads-1,S.Moments);

  # Perform sampling
  np = nprocs();
  ss = 0;

  @sync begin
    while ss < S.MacroSamples
      for p=1:np
        if p != myid() || np == 1
          if ss == S.MacroSamples
            break;
          end
          println("    Starting macrosample ", ss+1);
          for i = 1:S.FullSteps
            (P,Q,F) = LeapfrogStep(P,Q,F,C,S.dt);
          end

          ss+=1;

          @async begin
            ii = ss;
            PP = P;
            QQ = Q;
            # Get distances between CG particles
            println("      Sample ",ii,"  over fine-grained variables...");
            RSamples[ii,:] = CGDists(QQ,C);
            PSamples[ii,:] = C.PsiFConstr*PP;
            # Print max: if it gets too large, the integrator becomes unstable.
            println("        Maximum distance = ",maximum(RSamples[ii,:]));
            MomentSamples[ii,:,:] = remotecall_fetch(ConstrainedSample,p,PP,QQ,C,S,CGStats.Observables)[1];
			# NB:the function "ConstrainedSample" returns a Tuple and "remotecall_fetch()[1]" returns the first value of this Tuple.
            println("      Sample ",ii," complete!");
          end
        end
      end
    end
  end
  println(" ...ComputeVEff done!");

  return ComputeVEff(RSamples,PSamples,MomentSamples,C)
end

function VEffConstrainedLangDyn_NVT(C::Chain,S::Simulator)
  println("Computing coarse graining through constrained Langevin dynamics ...")
  if S.dt > 0.1*sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt));
    warn("Timestep = ",S.dt,
      " not ≪  shortest timescale = ",sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt)));
  end

  # Generate initial condition
  (P,Q,F,H,PCG,QCG) = InitConds(C);

  # Set up storage
  RSamples 	= SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads-1);
  PSamples 	= SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads-1);
  MomentSamples = SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads-1,S.Moments);

  # Perform sampling
  np = nprocs();
  ss = 0;

  @sync begin
    while ss < S.MacroSamples
      for p=1:np
        if p != myid() || np == 1
          if ss == S.MacroSamples
            break
          end
          println("    Starting macrosample ", ss+1);
          for i = 1:S.FullSteps
			  if i % S.check == 0
				  @printf "      ITERATION %0.3e OF %0.3e\n" i S.FullSteps
			  end
			  (P,Q,F) = LangDynIntegrators(P,Q,F,C,S.dt);
          end

          ss+=1;

          @async begin
            ii = ss;
            PP = P;
            QQ = Q;
            # Get distances between CG particles
            println("      Sample ",ii,"  over fine-grained variables...");
            RSamples[ii,:] = CGDists(QQ,C);
            PSamples[ii,:] = C.PsiFConstr*PP;
            # Print max: if it gets too large, the integrator becomes unstable.
            println("        Maximum distance = ",maximum(RSamples[ii,:]));
            MomentSamples[ii,:,:] = remotecall_fetch(ConstrSampleLangDyn,p,PP,QQ,C,S,CGStats.Observables)[1];
			# NB:the function "ConstrSampleLangDyn" returns a Tuple and "remotecall_fetch()[1]" returns the first value of this Tuple.
            println("      Sample ",ii," complete!");
          end
        end
      end
    end
  end
  println(" ...ComputeVEff done!")

  return ComputeVEff(RSamples,PSamples,MomentSamples,C)
end

function ConstrainedAutocorrPar(C::Chain,S::Simulator)
  if S.dt > 0.1*sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt));
    warn("Timestep = ",S.dt,
      " not ≪  shortest timescale = ",sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt)));
  end

  # Generate initial condition
  (P,Q,F,H,PCG,QCG) = InitConds(C);

  println("Equilibration ");
  # Run in to equilibrate
  for i = 1:S.FullSteps
    P,Q,F = LeapfrogStep(P,Q,F,C,S.dt)
  end
  println("Equilibration done");


  # Set up storage
  StressStats 	= SharedArray{Float64}(S.MacroSamples,S.NumLags,C.ConstrainedBeads-1);
  StressMean   	= SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads-1);
  StressSqMean 	= SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads-1);
  MomStats  = SharedArray{Float64}(S.MacroSamples,S.NumLags,C.ConstrainedBeads);
  MomMean	= SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads);
  MomSqMean = SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads);


  # Perform sampling
  np = nprocs()
  ss = 0

  println("Macrosampling ");
  @sync begin
    while ss < S.MacroSamples
      for p=1:np
        if p != myid() || np == 1
          if ss == S.MacroSamples
            break
          end
          println("    Starting macrosample ", ss+1);
          for i = 1:S.FullSteps
            (P,Q,F) = LeapfrogStep(P,Q,F,C,S.dt)
          end

          ss+=1;

          @async begin
            ii = ss;
            PP = P;
            QQ = Q;
            # Get distances between CG particles
            println("      Sample ",ii,"  over fine-grained variables...");
            # RSamples[ii,:] = CGDists(QQ,C);
            # PSamples[ii,:] = C.PsiFConstr*PP;
            # Print max: if it gets too large, the integrator becomes unstable.
            # println("        Maximum distance = ",maximum(RSamples[ii,:]));
            (StressStats[ii,:,:],StressMean[ii,:],StressSqMean[ii,:],MomStats[ii,:,:],MomMean[ii,:,:],MomSqMean[ii,:]) = remotecall_fetch(ConstrainedAutocorr,p,PP,
            					    QQ,C,S);
            println("      Sample ",ii," complete!");
          end
        end
      end
    end
  end

  # Now compute autocorrelations
  StressAutocorr = zeros(S.MacroSamples,S.NumLags,C.ConstrainedBeads-1)
  for i=1:S.NumLags
      StressAutocorr[:,i,:] = (StressStats[:,i,:].-StressMean[:,:].^2)./(StressSqMean[:,:].-StressMean[:,:].^2);
  end
  StressAutocorr = mean(mean(StressAutocorr,dims=1),dims=3);
  StressAutocorr = reshape(StressAutocorr,length(StressAutocorr));

  MomAutocorr = zeros(S.MacroSamples,S.NumLags,C.ConstrainedBeads)
  for i=1:S.NumLags
      MomAutocorr[:,i,:] = (MomStats[:,i,:].-MomMean[:,:].^2)./(MomSqMean[:,:].-MomMean[:,:].^2);
  end
  MomAutocorr = mean(mean(MomAutocorr,dims=1),dims=3);
  MomAutocorr = reshape(MomAutocorr,length(MomAutocorr));

  println(" ConstrainedAutocorrPar Done!")

  return StressAutocorr, MomAutocorr
end

function ConstrainedAutocorrLangDynPar(C::Chain,S::Simulator)
  if S.dt > 0.1*sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt));
    warn("Timestep = ",S.dt,
      " not ≪ shortest timescale = ",sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt)));
  end

  # Generate initial condition
  (P,Q,F,H,PCG,QCG) = InitConds(C);

  println("Equilibration ");
  # Run in to equilibrate
  for i = 1:S.FullSteps
    (P,Q,F) = LangDynIntegrators(P,Q,F,C,S.dt)
  end
  println("Equilibration done");


  # Set up storage
  StressStats 	= SharedArray{Float64}(S.MacroSamples,S.NumLags,C.ConstrainedBeads-1);
  StressMean   	= SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads-1);
  StressSqMean 	= SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads-1);
  MomStats  = SharedArray{Float64}(S.MacroSamples,S.NumLags,C.ConstrainedBeads);
  MomMean	= SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads);
  MomSqMean = SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads);


  # Perform sampling
  np = nprocs()
  ss = 0

  println("Macrosampling ");
  @sync begin
    while ss < S.MacroSamples
      for p=1:np
        if p != myid() || np == 1
          if ss == S.MacroSamples
            break
          end
          println("    Starting macrosample ", ss+1);
          for i = 1:S.FullSteps
            (P,Q,F) = LangDynIntegrators(P,Q,F,C,S.dt)
          end

          ss+=1;

          @async begin
            ii = ss;
            PP = P;
            QQ = Q;
            # Get distances between CG particles
            println("      Sample ",ii,"  over fine-grained variables...");
            # RSamples[ii,:] = CGDists(QQ,C);
            # PSamples[ii,:] = C.PsiFConstr*PP;
            # Print max: if it gets too large, the integrator becomes unstable.
            # println("        Maximum distance = ",maximum(RSamples[ii,:]));
            (StressStats[ii,:,:],StressMean[ii,:],StressSqMean[ii,:],MomStats[ii,:,:],MomMean[ii,:,:],MomSqMean[ii,:]) =
			              remotecall_fetch(ConstrainedAutocorrLangDyn,p,PP,QQ,C,S);
            println("      Sample ",ii," complete!");
          end
        end
      end
    end
  end

  # Now compute autocorrelations
  StressAutocorr = zeros(S.MacroSamples,S.NumLags,C.ConstrainedBeads-1)
  for i=1:S.NumLags
      StressAutocorr[:,i,:] = (StressStats[:,i,:].-StressMean[:,:].^2)./(StressSqMean[:,:].-StressMean[:,:].^2);
  end
  StressAutocorr = mean(mean(StressAutocorr,dims=1),dims=3);
  StressAutocorr = reshape(StressAutocorr,length(StressAutocorr));

  MomAutocorr = zeros(S.MacroSamples,S.NumLags,C.ConstrainedBeads)
  for i=1:S.NumLags
      MomAutocorr[:,i,:] = (MomStats[:,i,:].-MomMean[:,:].^2)./(MomSqMean[:,:].-MomMean[:,:].^2);
  end
  MomAutocorr = mean(mean(MomAutocorr,dims=1),dims=3);
  MomAutocorr = reshape(MomAutocorr,length(MomAutocorr));

  println(" ConstrainedAutocorrParLangDyn Done!")

  return StressAutocorr, MomAutocorr
end

function FullAutocorrPar(C::Chain,S::Simulator,dQData::Array{Float64,1},FmeanData::Array{Float64,1})
  if S.dt > 0.1*sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt));
    warn("Timestep = ",S.dt,
      " not ≪  shortest timescale = ",sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt)));
  end

  # Generate initial condition
  (P,Q,F,H,PCG,QCG) = InitConds(C);

  # Generate mean force function
  MeanForce = X -> ApproxMeanStresses(X,FmeanData,C,S)[1];

  println("Equilibration ");
  # Run in to equilibrate
  for i = 1:S.FullSteps
    P,Q,F = LeapfrogStep(P,Q,F,C,S.dt)
  end
  println("Equilibration done");


  # Set up storage
  StressStats 	= SharedArray{Float64}(S.MacroSamples,S.NumLags,C.NBeads);
  StressMean   	= SharedArray{Float64}(S.MacroSamples,C.NBeads);
  StressSqMean  = SharedArray{Float64}(S.MacroSamples,C.NBeads);
  MomStats    	= SharedArray{Float64}(S.MacroSamples,S.NumLags,C.NBeads);
  MomMean	= SharedArray{Float64}(S.MacroSamples,C.NBeads);
  MomSqMean	= SharedArray{Float64}(S.MacroSamples,C.NBeads);


  # Perform sampling
  np = nprocs()
  ss = 0

  println("Macrosampling ");
  @sync begin
    while ss < S.MacroSamples
      for p=1:np
        if p != myid() || np == 1
          if ss == S.MacroSamples
            break
          end
          println("    Starting macrosample ", ss+1);
          for i = 1:S.FullSteps
            (P,Q,F) = LeapfrogStep(P,Q,F,C,S.dt)
          end

          ss+=1;

          @async begin
            ii = ss;
            PP = P;
            QQ = Q;
            # Get distances between CG particles
            println("      Sample ",ii,"  over fine-grained variables...");
            # RSamples[ii,:] = CGDists(QQ,C);
            # PSamples[ii,:] = C.PsiFConstr*PP;
            # Print max: if it gets too large, the integrator becomes unstable.
            # println("        Maximum distance = ",maximum(RSamples[ii,:]));
            (StressStats[ii,:,:],StressMean[ii,:],StressSqMean[ii,:],MomStats[ii,:,:],MomMean[ii,:],MomSqMean[ii,:]) = remotecall_fetch(FullAutocorr,p,PP,QQ,C,S,MeanForce);
            println("      Sample ",ii," complete!");
          end
        end
      end
    end
  end


  # Compute autocorrelations
  StressAutocorr = zeros(S.MacroSamples,S.NumLags,C.NBeads)
  for i=1:S.NumLags
      StressAutocorr[:,i,:] = (StressStats[:,i,:].-StressMean[:,:].^2)./(StressSqMean[:,:].-StressMean[:,:].^2);
  end
  StressAutocorr = mean(mean(StressAutocorr,dims=1),dims=3);
  StressAutocorr = reshape(StressAutocorr,length(StressAutocorr));

  MomAutocorr = zeros(S.MacroSamples,S.NumLags,C.NBeads)
  for i=1:S.NumLags
      MomAutocorr[:,i,:] = (MomStats[:,i,:].-MomMean[:,:].^2)./(MomSqMean[:,:].-MomMean[:,:].^2);
  end
  MomAutocorr = mean(mean(MomAutocorr,dims=1),dims=3);
  MomAutocorr = reshape(MomAutocorr,length(MomAutocorr));

  println(" FullAutocorrPar Done!")

  return StressAutocorr, MomAutocorr
end

function FullAutocorrLangDynPar(C::Chain,S::Simulator,FmeanData::Array{Float64,2})
  if S.dt > 0.1*sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt));
    warn("Timestep = ",S.dt,
      " not ≪ shortest timescale = ",sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt)));
  end

  # Generate initial condition
  (P,Q,F,H,PCG,QCG) = InitConds(C);

  # Generate mean force function
  MeanForce = X -> ApproxMeanStresses(X,FmeanData,C,S)[1];

  println("Equilibration ");
  # Run in to equilibrate
  for i = 1:S.FullSteps
    (P,Q,F) = LangDynIntegrators(P,Q,F,C,S.dt)
  end
  println("Equilibration done");


  # Set up storage
  StressStats 	= SharedArray{Float64}(S.MacroSamples,S.NumLags,C.NBeads);
  StressMean   	= SharedArray{Float64}(S.MacroSamples,C.NBeads);
  StressSqMean  = SharedArray{Float64}(S.MacroSamples,C.NBeads);
  MomStats    	= SharedArray{Float64}(S.MacroSamples,S.NumLags,C.NBeads);
  MomMean	= SharedArray{Float64}(S.MacroSamples,C.NBeads);
  MomSqMean	= SharedArray{Float64}(S.MacroSamples,C.NBeads);


  # Perform sampling
  np = nprocs()
  ss = 0

  println("Macrosampling ");
  @sync begin
    while ss < S.MacroSamples
      for p=1:np
        if p != myid() || np == 1
          if ss == S.MacroSamples
            break
          end
          println("    Starting macrosample ", ss+1);
          for i = 1:S.FullSteps
            (P,Q,F) = LangDynIntegrators(P,Q,F,C,S.dt)
          end

          ss+=1;

          @async begin
            ii = ss;
            PP = P;
            QQ = Q;
            # Get distances between CG particles
            println("      Sample ",ii,"  over fine-grained variables...");
            # RSamples[ii,:] = CGDists(QQ,C);
            # PSamples[ii,:] = C.PsiFConstr*PP;
            # Print max: if it gets too large, the integrator becomes unstable.
            # println("        Maximum distance = ",maximum(RSamples[ii,:]));
            (StressStats[ii,:,:],StressMean[ii,:],StressSqMean[ii,:],MomStats[ii,:,:],MomMean[ii,:],MomSqMean[ii,:]) =
			                                                                         remotecall_fetch(FullAutocorrLangDyn,p,PP,QQ,C,S,MeanForce);
            println("      Sample ",ii," complete!");
          end
        end
      end
    end
  end


  # Compute autocorrelations
  StressAutocorr = zeros(S.MacroSamples,S.NumLags,C.NBeads)
  for i=1:S.NumLags
      StressAutocorr[:,i,:] = (StressStats[:,i,:].-StressMean[:,:].^2)./(StressSqMean[:,:].-StressMean[:,:].^2);
  end
  StressAutocorr = mean(mean(StressAutocorr,dims=1),dims=3);
  StressAutocorr = reshape(StressAutocorr,length(StressAutocorr));

  MomAutocorr = zeros(S.MacroSamples,S.NumLags,C.NBeads)
  for i=1:S.NumLags
      MomAutocorr[:,i,:] = (MomStats[:,i,:].-MomMean[:,:].^2)./(MomSqMean[:,:].-MomMean[:,:].^2);
  end
  MomAutocorr = mean(mean(MomAutocorr,dims=1),dims=3);
  MomAutocorr = reshape(MomAutocorr,length(MomAutocorr));

  println(" FullAutocorrParLangDyn Done!")



  return StressAutocorr, MomAutocorr
end

function ApproximateAutocorrPar(C::Chain,S::Simulator,FmeanData::Array{Float64,2},FstdData::Array{Float64,2})
  if S.dt > 0.1*sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt));
    warn("Timestep = ",S.dt,
      " not ≪  shortest timescale = ",sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt)));
  end

  # Set up storage
  StressStats 	= SharedArray{Float64}(S.MacroSamples,S.NumLags,C.NBeads);
  StressMean   	= SharedArray{Float64}(S.MacroSamples,C.NBeads);
  StressSqMean  = SharedArray{Float64}(S.MacroSamples,C.NBeads);
  MomStats    	= SharedArray{Float64}(S.MacroSamples,S.NumLags,C.NBeads);
  MomMean	= SharedArray{Float64}(S.MacroSamples,C.NBeads);
  MomSqMean	= SharedArray{Float64}(S.MacroSamples,C.NBeads);


  # Perform sampling
  np = nprocs()
  ss = 0

  println("Macrosampling ");
  @sync begin
    while ss < S.MacroSamples
      for p=1:np
        if p != myid() || np == 1
          if ss == S.MacroSamples
            break
          end
          println("    Starting macrosample ", ss+1);
          # Generate initial condition
          P = randn(C.NBeads);
          Q = [0.0:C.BeadsDist:C.BeadsDist*(C.NBeads-1);];
          ss+=1;

          @async begin
            ii = ss;
            PP = P;
            QQ = Q;
            # Get distances between CG particles
            println("      Sample ",ii,"  over fine-grained variables...");
			(StressStats[ii,:,:],StressMean[ii,:],StressSqMean[ii,:],MomStats[ii,:,:],MomMean[ii,:],MomSqMean[ii,:]) = remotecall_fetch(ApproximateAutocorr,p,PP,QQ,C,S,FmeanData,FstdData);
            println("      Sample ",ii," complete!");
          end
        end
      end
    end
  end


  # Compute autocorrelations
  StressAutocorr = zeros(S.MacroSamples,S.NumLags,C.NBeads)
  for i=1:S.NumLags
      StressAutocorr[:,i,:] = (StressStats[:,i,:].-StressMean[:,:].^2)./(StressSqMean[:,:].-StressMean[:,:].^2);
  end
  StressAutocorr = mean(mean(StressAutocorr,dims=1),dims=3);
  StressAutocorr = reshape(StressAutocorr,length(StressAutocorr));

  MomAutocorr = zeros(S.MacroSamples,S.NumLags,C.NBeads)
  for i=1:S.NumLags
      MomAutocorr[:,i,:] = (MomStats[:,i,:].-MomMean[:,:].^2)./(MomSqMean[:,:].-MomMean[:,:].^2);
  end
  MomAutocorr = mean(mean(MomAutocorr,dims=1),dims=3);
  MomAutocorr = reshape(MomAutocorr,length(MomAutocorr));

  println(" ApproximateAutocorrPar Done!")

  return StressAutocorr, MomAutocorr
end

###################
# Data processing #
###################
function ComputeVEff(R,P,M,C)
  l = length(R)

  RR = reshape(R[:,:],l)
  i = sortperm(RR)
  RR = RR[i];
  dR = circshift(RR,-1)-RR;

  PP = reshape(P[:,:],l)
  PP = PP[i];

  Fm = reshape(M[:,:,1],l);
  Fm = Fm[i];

  V = cumsum(0.5*(Fm+circshift(Fm,-1)).*dR); #trapezoidal rule
  V = circshift(V,1);
  V[1] = 0.0;

  Fstd = reshape(sqrt.(M[:,:,2]-M[:,:,1].^2),l)
  Fstd = Fstd[i];

  ## Kriging-type interpolation

  # N = 500;
  # Rmin = 1.5;
  # Rmax = 4.5;
  # range = Rmax - Rmin;
  # delta = range / N;
  Rmin = 2.802; #NB: with these values the Kriging-type interpolation algorithm does not give problems
  Rmax = 3.504;
  delta = 0.006;

  dQ = [Rmin:delta:Rmax;]

  FmINT = InterpolationAlgorithm_Kriging(dQ,RR,Fm,C)[1];
  FmINT_D = (circshift(FmINT,-1) - circshift(FmINT,1)) ./(2*delta);

  FstdINT = InterpolationAlgorithm_Kriging(dQ,RR,Fstd,C)[1];
  FstdINT_D = (circshift(FstdINT,-1) - circshift(FstdINT,1)) ./(2*delta);


  ## to generate the COEFFICIENTS for the cubic spline interpolation from the kriging-type interpolation of the MEAN FORCE DATA
  ## NB: the grid points dQ *must* be equally spaced

  Fmean = FmINT[2:end-1];
  Fmean_d = FmINT_D[2:end-1];

  Fmean_coeff1 = (2. * (Fmean - circshift(Fmean, -1)) + (Fmean_d + circshift(Fmean_d, -1)) * delta) / delta;
  Fmean_coeff2 = circshift(Fmean, -1) - Fmean + (-Fmean_d - Fmean_coeff1) * delta;

  # Generate look-up table for Mean Force
  FdataMF = zeros(length(Fmean) - 1, 5);
  FdataMF[:, 1] = dQ[2:end - 2];
  FdataMF[:, 2] = Fmean[1:end - 1];
  FdataMF[:, 3] = Fmean_d[1:end - 1];
  FdataMF[:, 4] = Fmean_coeff1[1:end - 1];
  FdataMF[:, 5] = Fmean_coeff2[1:end - 1];

  ## to generate the COEFFICIENTS for the cubic spline interpolation from the kriging-type interpolation of the FLUCTUATING FORCE DATA
  ## NB: the grid points dQ *must* be equally spaced

  FF = FstdINT[2:end-1];
  FF_d = FstdINT_D[2:end-1];

  FF_coeff1 = (2. * (FF - circshift(FF, -1)) + (FF_d + circshift(FF_d, -1)) * delta) / delta;
  FF_coeff2 = circshift(FF, -1) - FF + (-FF_d - FF_coeff1) * delta;

  # Generate look-up table for Fluctuating Force
  FdataFF = zeros(length(FF) - 1, 5);
  FdataFF[:, 1] = dQ[2:end - 2];
  FdataFF[:, 2] = FF[1:end - 1];
  FdataFF[:, 3] = FF_d[1:end - 1];
  FdataFF[:, 4] = FF_coeff1[1:end - 1];
  FdataFF[:, 5] = FF_coeff2[1:end - 1];


  return RR, PP, V, Fm, Fstd, FdataMF, FdataFF
end

##################################
# Potential generation functions #
##################################
"""
    Sigma, dQ = InterpolationAlgorithm_Kriging(Q::Array{Float64,1},dQvals::Array{Float64,1},Fvals::Array{Float64,1},C::Chain)

	Generates approximate mean and fluctuating stresses based on precomputed data dQvals, Fvals and Chain parameters via a Kriging-type interpolation algorithm.
"""
function InterpolationAlgorithm_Kriging(dQ::Array{Float64,1},dQvals::Array{Float64,1},Fvals::Array{Float64,1},C::Chain,σ::Float64=0.1,N::Int=2)

  Sigma = zeros(size(dQ));

  for i=1:length(dQ)

	  W = spdiagm(0 => exp.(-(dQvals .- dQ[i]) .^2 ./ σ^2));

	  # Form Jacobian matrix
      J = ones(N+1, length(dQvals));
      for j = 1:N
		  J[j+1, :] = (dQ[i] .- dQvals).^j;
	  end
	  K = J * W;
	  c = (K * J') \ (K * Fvals);

	  Sigma[i] = c[1];
  end

  return Sigma, dQ
end

"""
    Sigma, dQ = ApproxMeanStresses(Q::Array{Float64,1},Fdata::Array{Float64,2},C::Chain,S::Simulator)

	Generates approximate mean stresses based on precomputed data dQvals, Fvals and Chain parameters via a Cubic Spline interpolation algorithm.
"""
function ApproxMeanStresses(Q::Array{Float64,1},Fdata::Array{Float64,2},C::Chain,S::Simulator,σ::Float64=0.1,N::Int=2)
	if sum(Fdata) == 0
		dQ = circshift(Q, -1) - Q;
		dQ[C.NBeads] += C.Period;

		return zeros(size(Q)), dQ
	else
		Sigma = zeros(size(Q));

  		### First compute inter-bead distances.
  		dQ = circshift(Q, -1) - Q;
  		dQ[C.NBeads] += C.Period;

  		Rmin = Fdata[1, 1];
  		Rmax = Fdata[end, 1];
  		dx = Fdata[2, 1] - Fdata[1, 1];
  		inv_sqr_delta = 1. / (dx * dx);

  		### Next, for each inter-bead distance, compute the stress
  		for i=1:length(dQ)
    		if dQ[i] < Rmin
        		Sigma[i] = S.m_MF * dQ[i] + S.q_MF;
    		elseif dQ[i] > Rmax
        		Sigma[i] = 0.;
    		else
        		index = trunc(Int, (dQ[i] - Rmin) / dx) + 1;
        		distance = dQ[i] - (Rmin + (index - 1) * dx);
        		Sigma[i] = Fdata[index, 2] + distance * (Fdata[index, 3] + distance * (Fdata[index, 5] + distance * Fdata[index, 4]) * inv_sqr_delta);
    		end
  		end

  		return Sigma, dQ
	end
end

"""
    F, Sigma, dQ = ApproxMeanForce(Q::Array{Float64,1},Fvals::Array{Float64,2},C::Chain,S::Simulator)

	Generates approximate mean forces F based on precomputed data dQvals, Fvals and Chain parameters via a Cubic Spline interpolation algorithm.
"""
function ApproxMeanForce(Q::Array{Float64,1},Fvals::Array{Float64,2},C::Chain,S::Simulator,σ::Float64=0.1,N::Int=2)

  F = zeros(size(Q));
  (Sigma, dQ) = ApproxMeanStresses(Q, Fvals, C, S, σ, N);
  F = Sigma - circshift(Sigma, 1);

  return F, Sigma, dQ
end

"""
    Thermal_bath_step, Gvec, dQ = ApproxFlucStep(P::Array{Float64,1},Q::Array{Float64,1},StressSqvals::Array{Float64,2},C::Chain,S::Simulator)

	Generates thermal bath intergration step and approximate fluctuating forces Gvec based on precomputed data StressSqvals and Chain parameters
	via a Kriging-type interpolation algorithm.
"""
function ApproxFlucStep(P::Array{Float64,1},Q::Array{Float64,1},StressSqvals::Array{Float64,2},C::Chain,S::Simulator,σ::Float64=0.1,N::Int=2)

	P = P .- mean(P);

    Gvec = zeros(size(Q));

    ### First compute inter-bead distances.
    dQ = circshift(Q, -1) - Q;
    dQ[C.NBeads] += C.Period;

    Rmin = StressSqvals[1, 1];
    Rmax = StressSqvals[end, 1];
    dx = StressSqvals[2, 1] - StressSqvals[1, 1];
    inv_sqr_delta = 1. / (dx * dx);

    ### Next, for each inter-bead distance, compute the stress
    for i=1:length(dQ)
        if dQ[i] < Rmin
            Gvec[i] = S.m_FF * dQ[i] + S.q_FF;
        elseif dQ[i] > Rmax
            Gvec[i] = S.coeffEXP_A * exp(-S.coeffEXP_B * dQ[i]);
        else
            index = trunc(Int, (dQ[i] - Rmin) / dx) + 1;
            distance = dQ[i] - (Rmin + (index - 1) * dx);
            Gvec[i] = StressSqvals[index, 2] + distance * (StressSqvals[index, 3] +
			          distance * (StressSqvals[index, 5] + distance * StressSqvals[index, 4]) * inv_sqr_delta);
        end
    end

    #println("dQ", dQ);

    ### Set diffusivity to zero if negative.
    Gvec = max.(Gvec, 0);

    ### Form memory matrix
    diagonal = Gvec + circshift(Gvec, 1);
    off_diagonal = -Gvec;
    Gamma = diagm(diagonal) + circshift(diagm(off_diagonal), -1) + circshift(diagm(off_diagonal), 1);

    ### Perform eigenvalue decomposition
	# Computes the eigenvalue decomposition, returning an Eigen factorization object
	# which contains the eigenvalues in F.values and the eigenvectors in the columns of the matrix F.vectors.
	F = eigen(Gamma);

    values = max.(real.(F.values), 0.);
    vectors = real.(F.vectors);

    ### Calculate exact Ornstein-Uhlenbeck step using decomposition
    exp_part = exp.(-values * C.Beta * S.Stochdt / (2.0 * C.BeadMass));
    c1 = diagm(0 => exp_part);

    # we want to avoid numerical issues and thus enforce the positiveness of the sqrt argument
    #sqrt_part = max.(1.0 .- exp_part.^2, 0);
    sqrt_part = max.(1.0 .- exp_part.^2, 0.) .* (C.BeadMass / C.Beta);
    c3 = diagm(0 => sqrt.(sqrt_part));

    return vectors * c1 * inv(vectors) * P + vectors * c3 * inv(vectors) * randn(C.NBeads), Gvec, dQ
end

end # module CGStats
