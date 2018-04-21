__precompile__()


################################################################################
# Module LJChain
################################################################################

"""
    LJChain is a module providing physics functions and a parameter object for a
    1D chain of 2 species of particle.
"""
module LJChain

using Parameters,StatsBase

# Physics functions
export Chain, LennardJones, LennardJonesPrime, Potential, Forces, Stresses,
  KineticEnergy
export InitConds, CGDists
# Time stepping functions
export LeapfrogStep, ConstrainedStep, ApproxDynStep

"""
    Chain(Patt,NBeads,SM,LM,SEps,LEps,Dist,Beta,Gamma)

Constructs an abstract type containing parameters necessary for simulation
of a periodic chain of particles of 2 types, with NBeads repeating units with
pattern Patt.
"""
@with_kw immutable Chain
  # Main parameters
  Patt::Array{Int,1} 		= [0,1,0] # Pattern of particles
  NBeads::Int			= 10	  # Number of beads
  ConstrainedBeads::Int 	= 2	  # Number of beads to constrain
  SM::Float64 			= 10.0	  # Mass of 0 particles
  LM::Float64 			= 1.0	  # Mass of 1 particles
  SEps::Float64			= 1.0     # Scaling of LJ potential for 0 particles
  LEps::Float64 		= 100.0   # Scaling of LJ potential for 1 particles
  Dist::Float64 		= 1.0     # Equilibrium distance for all LJ potentials
  Beta::Float64 		= 1.0	  # Inverse temperature
  Gamma::Float64		= 1.0	  # Friction coefficient for Langevin dynamics

  # The following are inherited from the previous parameters: changing them may break the code.
  NPart::Int 			= length(Patt)*NBeads
  Period::Float64 		= NPart*Dist
  BeadsDist 			= length(Patt)*Dist
  LTotPatt::Array{Int,1} 	= repmat(Patt,NBeads)
  STotPatt::Array{Int,1}	= 1.-LTotPatt
  # Mass parameters
  MassPatt::Array{Float64,1} 	= LM*Patt+SM*(1.-Patt)
  MassTotPatt::Array{Float64,1} = repmat(MassPatt,NBeads)
  MInv::Array{Float64,2}	= diagm(1./MassTotPatt)
  MHalf::Array{Float64,2}	= diagm(sqrt.(MassTotPatt))
  BeadMass::Float64		= LM*sum(Patt.==1)+SM*sum(Patt.==0)
  EpsPatt::Array{Float64,1} 	= (sqrt.(LEps)*LTotPatt+sqrt.(SEps)*STotPatt).*circshift(sqrt.(LEps)*LTotPatt+sqrt.(SEps)*STotPatt,-1)
  # Full CG Projections
  Phi::Array{Float64,2}		= kron(eye(NBeads),ones(Float64,length(Patt)))'
  Psi 				= kron(eye(NBeads),MassPatt./BeadMass)'
  Xi 				= MInv\(eye(NPart)-Phi'Psi)*MInv;
  Zeta 				= eye(NPart)-Psi'Phi;

  # Constrained versions of CG Projections
  PhiConstr::Array{Float64,2}	= kron(eye(NBeads,ConstrainedBeads),ones(Float64,length(Patt)))'
  PsiConstr 			= kron(eye(NBeads,ConstrainedBeads),MassPatt./BeadMass)'
  XiConstr 			= MInv\(eye(NPart)-PhiConstr'PsiConstr)*MInv;
  ZetaConstr 			= eye(NPart)-PsiConstr'PhiConstr;

  # Relative CG projections between identical beads (full and constrained)
  PhiF::Array{Float64,2} 	= kron(eye(NBeads),[zeros(Float64,length(Patt)-1);1.0])'  # Forces between beads
  PhiFConstr::Array{Float64,2} 	= kron(eye(NBeads,ConstrainedBeads-1),[zeros(Float64,length(Patt)-1);1.0])'  # Forces between constrained beads
  PsiF			        = (circshift(Phi,(0,length(Patt)))-Phi)/2.;  # Relative momenta
  PsiFConstr			= (circshift(PhiConstr[1:ConstrainedBeads-1,:],(0,length(Patt)))-PhiConstr[1:ConstrainedBeads-1,:])/2.;  # Relative momenta
end


"""
    U = LennardJones(Q,Dist,EpsPatt)

Lennard-Jones 6-12 potential with minimum at Dist and well depth EpsPatt.
"""
function LennardJones(Q::Array{Float64,1},Dist::Float64,
						EpsPatt::Array{Float64,1})
  U = 4.*EpsPatt.*((Dist./Q).^12-2*(Dist./Q).^6);
  return U
end

"""
    G = LennardJonesPrime(Q,Dist,EpsPatt)

Derivative of Lennard-Jones 6-12 potential with minimum at Dist and well depth
EpsPatt.
"""
function LennardJonesPrime(Q::Array{Float64,1},Dist::Float64,EpsPatt::Array{Float64,1})
  G = -48.*EpsPatt.*(Dist^12./(Q.^13)-Dist^6./(Q.^7));
  return G
end

"""
    E = Potential(Q,C::Chain)

Computes internal energy of chain.
"""
function Potential(Q,C::Chain)
  # Right-hand finite difference vector
  dQ 	 = circshift(Q,-1)-Q;
  dQ[C.NPart] += C.Period;

  # Get bond energies
  BondEnergies = LennardJones(dQ,C.Dist,C.EpsPatt)

  # Compute total potential energyf
  E = sum(BondEnergies);

  return E
end

"""
    F = Forces(Q,C::Chain)

Computes forces for a periodic chain described by the parameters encoded
in C.
"""
function Forces(Q,C::Chain)

  Sigma = Stresses(Q,C);

  F = Sigma-circshift(Sigma,1);

  return F
end

"""
    Sigma = Stresses(Q,C::Chain)

Computes stresses for a periodic chain described by the parameters encoded
in C.
"""
function Stresses(Q,C::Chain)
  # Get number of particles
  N = C.NPart;

  # Right-hand finite difference vector
  dQ = circshift(Q,-1)-Q;
  dQ[N] = dQ[N]+C.Period;

  Sigma = LennardJonesPrime(dQ,C.Dist,C.EpsPatt)

  return Sigma
end

"""
    KE = KineticEnergy(P,C::Chain)

Computes KE of chain with momenta P.
"""
function KineticEnergy(P::Array{Float64,1},C::Chain)
  return 0.5*dot(P,C.MInv*P)
end


# Time-stepping functions
"""
    (PNext,QNext,FNext) = LeapfrogStep(P,Q,F,C::Chain,dt)

Computes forward step of Leapfrog scheme for one-dimensional arrays of momentum P,
positions Q, forces F and chain parameters C, with timestep dt.
"""
function LeapfrogStep(P::Array{Float64,1},Q::Array{Float64,1},F::Array{Float64,1},C::Chain,dt::Float64)
  PNext = P + 0.5.*F.*dt 	   # B
  QNext = Q + dt.*C.MInv*PNext;    # A
  FNext = Forces(QNext,C);         # Force evaluation
  PNext = PNext + 0.5.*FNext.*dt;  # B
  return PNext,QNext,FNext
end

"""
    (PNext,QNext,FNext) = ConstrainedStep(P,Q,F,C::Chain,dt)

Computes leapfrog step of dynamics with timestep dt, constrained to maintain CG
positions and momenta.
"""
function ConstrainedStep(P,Q,F,C,dt)
  PNext = P + 0.5.*C.ZetaConstr*F.*dt	 	# B
  QNext = Q + C.MInv*C.ZetaConstr*PNext.*dt;	# A
  FNext = Forces(QNext,C);			# Force evaluation
  PNext = PNext + 0.5*C.ZetaConstr*FNext.*dt;	# B
  return PNext,QNext,FNext
end

"""
    (PNext,QNext,FNext) = ApproxDynStep(P,Q,FGamma,C::Chain,dt)

Computes forward step of Langevin dynamics with timestep dt using Leimkuhler-Matthews
"BAOAB" integrator, adapted to non-constant diffusion; diffusion coefficient
is evaluated after performing `half' of the deterministic part of the timestep.
"""
function ApproxDynStep(P::Array{Float64,1},Q::Array{Float64,1},F::Array{Float64,1},MeanForce,FlucStep,C::Chain,dt)
  PNext = P+0.5*dt*F;			# B
  QNext = Q+0.5*dt/C.BeadMass*PNext;  	# A
  PNext = FlucStep(PNext,QNext);
  QNext = QNext+0.5.*dt/C.BeadMass*PNext;         	 # A
  FNext = MeanForce(QNext);	  			 # Force evaluation
  PNext = PNext+0.5*dt*FNext;                		 # B
  return PNext,QNext,FNext
end

# Initial condition generation

"""
    (P0,Q0,F0,H0,PCG0,QCG0) = InitConds(C::Chain)

Generates initial momenta P0, positions Q0, forces F0, Hamiltonian H0, CG momenta PCG0 and
CG positions QCG0. P0 is random based on Boltzmann distribution, Q0 lie at equilibrium
positions.
"""
function InitConds(C::Chain)
  # Fix the initial condition for CG variables.
  P0 = sqrt.(C.MassTotPatt./C.Beta).*randn(C.NPart);
  # Set total momentum to 0 to avoid drift.
  P0 -= mean(P0);
  P0 *= sqrt.((C.NPart/C.Beta) / KineticEnergy(P0,C));
  Q0 = convert(Array{Float64,1},linspace(0.0,(C.NPart-1)*C.Dist,C.NPart));
  Q0 = Q0+0.001*C.Beta*randn(C.NPart)

  # Fix reference energy
  E0 = Potential(Q0,C);
  H0 = KineticEnergy(P0,C) + E0;

  # Compute initial forces
  F0 = Forces(Q0,C);

  # Extract IC for CG variables
  PCG0 = C.Phi*P0;
  QCG0 = C.Psi*Q0;

  return P0,Q0,F0,H0,PCG0,QCG0
end

"""
    R = CGDists(Q::Array{Float64,1},C::Chain)

Returns distances between centres of mass of beads in a vector R.
"""
function CGDists(Q::Array{Float64,1},C::Chain)
  R = circshift(C.PsiConstr*Q,-1)-C.PsiConstr*Q;
  if C.ConstrainedBeads == C.NBeads
    R[C.NBeads] += C.Period;
  end
  return R[1:end-1]
end

end # module LJChain

################################################################################
# Module CGStats
################################################################################
module CGStats

# Include relevant modules
using LJChain, Parameters, StatsBase

# Parameter object
export Simulator
# Functions to output raw trajectory data
export FullTraj, ConstrainedTraj, ApproximateTraj
# High-level function to extract data
export VEffConstrained, ConstrainedAutocorr, FullAutocorr, ApproximateAutocorr, ConstrainedAutocorrPar, FullAutocorrPar, ApproximateAutocorrPar
# Orthogonal dynamics sampling
export ConstrainedSample
# Functions to extract approximate dynamics
export ApproxMeanForce, ApproxMeanStresses, ApproxFlucStep, ComputeVEff

#######################################
# Simulator type, encoding parameters #
#######################################
@with_kw immutable Simulator
  dt::Float64		= 1e-3
  Stochdt::Float64	= 1e-3
  FullSteps::Int	= Int(1e5)
  OrthoSteps::Int 	= 1
  ApproxSteps::Int	= Int(1e5)
  FullTime::Float64	= FullSteps*dt
  ApproxTime::Float64 	= ApproxSteps*dt
  LagStep		= 400
  NumLags		= 30
  MacroSamples::Int	= 50
  OrthoSamples::Int	= Int(1e5)
  Moments::Int		= 2
end

#############################################
# Functions to generate raw trajectory data #
#############################################
"""
    (Pstore,Qstore,Fstore) = FullTraj(P0,Q0,C::Chain,S::Simulator)

    Outputs full trajectory data in matrices (Pstore,Qstore,Fstore).
"""
function FullTraj(P0,Q0,C::Chain,S::Simulator)
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
    (P,Q,F) = LeapfrogStep(P,Q,F,C,S.dt);
    Pstore[:,i] = P;
    Qstore[:,i] = Q;
    Fstore[:,i] = F;
  end

  return Pstore,Qstore,Fstore
end

"""
    (Pstore,Qstore,Fstore) = ConstrainedTraj(P0,Q0,C::Chain,S::Simulator)

    Outputs full trajectory data in matrices (Pstore,Qstore,Fstore).
"""
function ConstrainedTraj(P0,Q0,C::Chain,S::Simulator)
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
    (Pstore,Qstore,Fstore) = ApproximateTraj(P0,Q0,dQData,FmeanData,FstdData,
    							S::Simulator,C::Chain)

    Outputs full trajectory data for approximate Langevin dynamics, using
    approximate potentials generated from data encoded in 3 vectors
    (dQData,FmeanData,FstdData).
"""
function ApproximateTraj(P0,Q0,dQData::Array{Float64,1},FmeanData::Array{Float64,1},FstdData::Array{Float64,1},C::Chain,S::Simulator)
  # Set up initial conditions
  P = P0; Q = Q0;

  # Generate mean force function
  MeanForce = X -> ApproxMeanForce(X,dQData,FmeanData,C);
  F = MeanForce(Q0);

  # Set up approximate fluctuation generator
  FlucStep = (X,Y) -> ApproxFlucStep(X,Y,dQData,FstdData,C,S.Stochdt)

  # Set up storage
  Pstore = zeros(S.ApproxSteps,length(P0));
  Qstore = zeros(S.ApproxSteps,length(P0));
  Fstore = zeros(S.ApproxSteps,length(P0));
  Pstore[1,:] = P;
  Qstore[1,:] = Q;
  Fstore[1,:] = F;

  # Do simulation.
  for i = 2:S.ApproxSteps
    (P,Q,F) = ApproxDynStep(P,Q,F,MeanForce,FlucStep,C,S.Stochdt);
    Pstore[i,:] = P;
    Qstore[i,:] = Q;
    Fstore[i,:] = F;
  end

  return Pstore,Qstore,Fstore
end

################################
# Orthogonal dynamics sampling #
################################
function ConstrainedSample(P::Array{Float64,1},Q::Array{Float64,1},C::Chain,
						S::Simulator,Obs=(P,Q,C) -> 0)
  # Mean Forces
  StressMoments = zeros(C.ConstrainedBeads-1,S.Moments);

  PCG = C.PhiConstr*P;
  QCG = C.PsiConstr*Q;
  E0  = Potential(Q,C);
  F   = Forces(Q,C);
  H0  = E0 + KineticEnergy(P,C);

  s = 0;
  ObsMean = [0,0];
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
    ObsMean += Obs(P,Q,C);
  end

  ObsMean /= S.OrthoSamples;
  Efinal = Potential(Q,C);
  Hfinal = Efinal + KineticEnergy(P,C);

  #println("        Error in momentum constraint = ",maximum(abs(PCG-C.PhiConstr*P)));
  #println("        Error in position constraint = ",maximum(abs(QCG-C.PsiConstr*Q)));
  #println("        Relative error in energy     = ",abs(Hfinal/H0-1));
  #if Obs!=0
  #  println("        Value of observable          = ",ObsMean);
  #end

  StressMoments ./= S.OrthoSamples
  return StressMoments
end

function ConstrainedAutocorr(P::Array{Float64,1},Q::Array{Float64,1},
							C::Chain,S::Simulator)
  # Mean Forces
  StressStats 	= zeros(S.NumLags,C.ConstrainedBeads-1);
  StressStorage = zeros(S.NumLags,C.ConstrainedBeads-1);
  StressMean   	= zeros(C.ConstrainedBeads-1);
  StressSqMean 	= zeros(C.ConstrainedBeads-1);
  MomStats    	= zeros(S.NumLags,C.ConstrainedBeads);
  MomStorage 	= zeros(S.NumLags,C.ConstrainedBeads);
  MomMean	= zeros(C.ConstrainedBeads);
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
    StressSqMean += SigmaCG.^2

    # Update momentum storage
    PCG = C.PhiConstr*P;
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

#   Efinal = Potential(Q,C);
#   Hfinal = Efinal + KineticEnergy(P,C);

  #println("        Error in momentum constraint = ",maximum(abs(PCG-C.PhiConstr*P)));
  #println("        Error in position constraint = ",maximum(abs(QCG-C.PsiConstr*Q)));
  #println("        Relative error in energy     = ",abs(Hfinal/H0-1));
  #if Obs!=0
  #  println("        Value of observable          = ",ObsMean);
  #end

  return StressStats,StressMean,StressSqMean,MomStats,MomMean,MomSqMean
end

###################################
# Full Trajectory Autocorrelation #
###################################

function FullAutocorr(P::Array{Float64,1},Q::Array{Float64,1},
							C::Chain,S::Simulator,MeanForce)
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

    #   Efinal = Potential(Q,C);
    #   Hfinal = Efinal + KineticEnergy(P,C);

    #println("        Error in momentum constraint = ",maximum(abs(PCG-C.PhiConstr*P)));
    #println("        Error in position constraint = ",maximum(abs(QCG-C.PsiConstr*Q)));
    #println("        Relative error in energy     = ",abs(Hfinal/H0-1));
    #if Obs!=0
    #  println("        Value of observable          = ",ObsMean);
    #end

    return StressStats,StressMean,StressSqMean,MomStats,MomMean,MomSqMean
end

function ApproximateAutocorr(P::Array{Float64,1},Q::Array{Float64,1},
							C::Chain,S::Simulator,dQData::Array{Float64,1},FmeanData::Array{Float64,1},FstdData::Array{Float64,1})
  # Generate mean force function
  MeanForce = X -> ApproxMeanForce(X,dQData,FmeanData,C);
  F = MeanForce(Q);

  # Set up approximate fluctuation generator
  FlucStep = (X,Y) -> ApproxFlucStep(X,Y,dQData,FstdData,C,S.Stochdt)

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
    StressSqMean += SigmaCG.^2

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

#   Efinal = Potential(Q,C);
#   Hfinal = Efinal + KineticEnergy(P,C);

  #println("        Error in momentum constraint = ",maximum(abs(PCG-C.PhiConstr*P)));
  #println("        Error in position constraint = ",maximum(abs(QCG-C.PsiConstr*Q)));
  #println("        Relative error in energy     = ",abs(Hfinal/H0-1));
  #if Obs!=0
  #  println("        Value of observable          = ",ObsMean);
  #end

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
function VEffConstrained(C::Chain,S::Simulator)
  println("Computing coarse graining through constrained dynamics ...")
  if S.dt > 0.1*sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt));
    warn("Timestep = ",S.dt,
      " not ≪  shortest timescale = ",sqrt.(minimum(C.MassPatt)/maximum(C.EpsPatt)));
  end

  # Generate initial condition
  (P,Q,F,H,PCG,QCG) = InitConds(C);

  # Set up storage
  RSamples 	= SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads-1);
  PSamples 	= SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads-1);
  MomentSamples = SharedArray{Float64}(S.MacroSamples,C.ConstrainedBeads-1,
  								S.Moments);

  # Perform sampling
  np = nprocs()
  ss = 0

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
            RSamples[ii,:] = CGDists(QQ,C);
            PSamples[ii,:] = C.PsiFConstr*PP;
            # Print max: if it gets too large, the integrator becomes unstable.
            println("        Maximum distance = ",maximum(RSamples[ii,:]));
            MomentSamples[ii,:,:] = remotecall_fetch(ConstrainedSample,p,PP,
            					    QQ,C,S,CGStats.Observables);
            println("      Sample ",ii," complete!");
          end
        end
      end
    end
  end
  println(" ...ComputeVEff done!")

  return ComputeVEff(RSamples,PSamples,MomentSamples)
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
#            RSamples[ii,:] = CGDists(QQ,C);
#            PSamples[ii,:] = C.PsiFConstr*PP;
            # Print max: if it gets too large, the integrator becomes unstable.
#            println("        Maximum distance = ",maximum(RSamples[ii,:]));
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
  StressAutocorr = mean(mean(StressAutocorr,1),3);
  StressAutocorr = reshape(StressAutocorr,length(StressAutocorr));

  MomAutocorr = zeros(S.MacroSamples,S.NumLags,C.ConstrainedBeads)
  for i=1:S.NumLags
      MomAutocorr[:,i,:] = (MomStats[:,i,:].-MomMean[:,:].^2)./(MomSqMean[:,:].-MomMean[:,:].^2);
  end
  MomAutocorr = mean(mean(MomAutocorr,1),3);
  MomAutocorr = reshape(MomAutocorr,length(MomAutocorr));

  println(" ConstrainedAutocorrPar Done!")

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
  MeanForce = X -> ApproxMeanStresses(X,dQData,FmeanData,C);

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
            #RSamples[ii,:] = CGDists(QQ,C);
            #PSamples[ii,:] = C.PsiFConstr*PP;
            # Print max: if it gets too large, the integrator becomes unstable.
            #println("        Maximum distance = ",maximum(RSamples[ii,:]));
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
  StressAutocorr = mean(mean(StressAutocorr,1),3);
  StressAutocorr = reshape(StressAutocorr,length(StressAutocorr));

  MomAutocorr = zeros(S.MacroSamples,S.NumLags,C.NBeads)
  for i=1:S.NumLags
      MomAutocorr[:,i,:] = (MomStats[:,i,:].-MomMean[:,:].^2)./(MomSqMean[:,:].-MomMean[:,:].^2);
  end
  MomAutocorr = mean(mean(MomAutocorr,1),3);
  MomAutocorr = reshape(MomAutocorr,length(MomAutocorr));

  println(" FullAutocorrPar Done!")



  return StressAutocorr, MomAutocorr
end

function ApproximateAutocorrPar(C::Chain,S::Simulator,dQData::Array{Float64,1},FmeanData::Array{Float64,1},FstdData::Array{Float64,1})
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
            (StressStats[ii,:,:],StressMean[ii,:],StressSqMean[ii,:],MomStats[ii,:,:],MomMean[ii,:],MomSqMean[ii,:]) = remotecall_fetch(ApproximateAutocorr,p,PP,QQ,C,S,dQData,FmeanData,FstdData);
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
  StressAutocorr = mean(mean(StressAutocorr,1),3); #TODO: Check averaging
  StressAutocorr = reshape(StressAutocorr,length(StressAutocorr));

  MomAutocorr = zeros(S.MacroSamples,S.NumLags,C.NBeads)
  for i=1:S.NumLags
      MomAutocorr[:,i,:] = (MomStats[:,i,:].-MomMean[:,:].^2)./(MomSqMean[:,:].-MomMean[:,:].^2);
  end
  MomAutocorr = mean(mean(MomAutocorr,1),3);
  MomAutocorr = reshape(MomAutocorr,length(MomAutocorr));

  println(" ApproximateAutocorrPar Done!")

  return StressAutocorr, MomAutocorr
end




###################
# Data processing #
###################
function ComputeVEff(R,P,M)
  l = length(R)
  RR = reshape(R[:,:],l)
  i = sortperm(RR)
  RR = RR[i];
  dR = circshift(RR,-1)-RR;

  PP = reshape(P[:,:],l)
  PP = PP[i];

  Fm = reshape(M[:,:,1],l);
  Fm = Fm[i];

  V = cumsum(0.5*(Fm+circshift(Fm,-1)).*dR);
  V = circshift(V,1);
  V[1] = 0.0;

  Fstd = reshape(sqrt.(M[:,:,2]-M[:,:,1].^2),l)
  Fstd = Fstd[i];
  return RR,PP,V,Fm,Fstd
end

##################################
# Potential generation functions #
##################################
"""
    F = ApproxMeanStresses(Q::Array{Float64,1},dQvals::Array{Float64,1},Fvals::Array{Float64,1},C::Chain)

Generates approximate mean stresses based on precomputed data dQvals, Fvals and Chain parameters.
"""
function ApproxMeanStresses(Q::Array{Float64,1},dQvals::Array{Float64,1},Fvals::Array{Float64,1},
				C::Chain,σ::Float64=0.1,N::Int=2)
  Sigma = zeros(size(Q));

  # First compute inter-bead distances.
  dQ = circshift(Q,-1)-Q;
  dQ[C.NBeads] += C.Period;

  # Next, for each inter-bead distance, compute the stress
  for i=1:length(dQ)
    # Define (diagonal) weight matrix
    W = spdiagm(exp.(-(dQvals-dQ[i]).^2/σ^2));

    # Form Jacobian matrix
    J = ones(N+1,length(dQvals));
    for j=1:N
      J[j+1,:] = (dQ[i]-dQvals).^j;
    end

    K = J*W
    c = (K*J')\(K*Fvals);

    Sigma[i] = c[1];
  end

  return Sigma
end

"""
    F = ApproxMeanForce(Q::Array{Float64,1},dQvals::Array{Float64,1},Fvals::Array{Float64,1},C::Chain)

Generates approximate mean forces based on precomputed data dQvals, Fvals and Chain parameters.
"""
function ApproxMeanForce(Q::Array{Float64,1},dQvals::Array{Float64,1},Fvals::Array{Float64,1},
				C::Chain,σ::Float64=0.1,N::Int=2)

  F = zeros(size(Q));
  Sigma = ApproxMeanStresses(Q,dQvals,Fvals,C,σ,N)
  F = Sigma-circshift(Sigma,1);

  return F
end


function ApproxFlucStep(P::Array{Float64,1},Q::Array{Float64,1},dQvals::Array{Float64,1},
	StressSqvals::Array{Float64,1},C::Chain,dt::Float64,σ::Float64=0.1,N::Int=2)

  # Set up storage for diffusivity
  Gvec = zeros(size(Q));

  # Compute inter-bead distances
  dQ = circshift(Q,-1)-Q;
  dQ[C.NBeads] += C.Period;

  # For each inter-bead distance, compute diffusivity
  for i=1:length(dQ)
    # Define (diagonal) weight matrix
    W = spdiagm(exp.(-(dQvals-dQ[i]).^2/σ^2));

    # Form Jacobian matrix
    J = ones(N+1,length(dQvals));
    for j=1:N
      J[j+1,:] = (dQ[i]-dQvals).^j;
    end

    # Least squares approximation of interpolant coefficients
    K = J*W
    coeffs = (K*J')\(K*StressSqvals);

    # Extract value of interpolant
    Gvec[i] = coeffs[1];
  end

  # Set diffusivity to zero if negative.
  Gvec = max.(Gvec,0);

  # Form memory matrix
  Gamma = diagm(Gvec)-circshift(diagm(Gvec),-1);
  Gamma = (Gamma-circshift(Gamma,1));

  # Perform eigenvalue decomposition
  D,V = eig(Gamma);
  D[1] = 0;

  # Calculate exact Ornstein-Uhlenbeck step using decomposition
  c1 = exp.(-D/(2*C.BeadMass)*dt);
  c3 = sqrt(C.BeadMass)*sqrt.(1-c1.^2).*randn(C.NBeads);

  # Generate new momenta and return
  return V*(c1.*(V'*P)+c3);
end

end # module CGStats
