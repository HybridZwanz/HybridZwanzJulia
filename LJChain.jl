__precompile__()

################################################################################
# Module LJChain
################################################################################

"""
    LJChain is a module providing physics functions and a parameter object for a
    1D chain of 2 species of particle.
"""
module LJChain

using Parameters,Statistics,LinearAlgebra

# Physics functions
export Chain, LennardJones, LennardJonesPrime, Potential, Forces, Stresses, KineticEnergy
export InitConds, CGDists
# Time stepping functions
export LeapfrogStep, ConstrainedStep, DetCGDynNVEStep,
       LangDynIntegrators, LangDynConstrainedStep, OrnUhlenbeckPart, DetCGDynNVTStep, OrnUhlenbeckPartDCG,
	   ApproxDynStep

"""
    Chain(Patt,NBeads,SM,LM,SEps,LEps,Dist,Beta,Gamma)

Constructs an abstract type containing parameters necessary for simulation
of a periodic chain of particles of 2 types, with NBeads repeating units with
pattern Patt.
"""
@with_kw struct Chain
  # Main parameters
  Patt::Array{Int,1} 	= [0,1,0] # Pattern of particles
  NBeads::Int			= 10	  # Number of beads
  ConstrainedBeads::Int = 2	  	  # Number of beads to constrain
  SM::Float64 			= 10.0	  # Mass of 0 particles
  LM::Float64 			= 1.0	  # Mass of 1 particles
  SEps::Float64			= 1.0     # Scaling of LJ potential for 0 particles
  LEps::Float64 		= 100.0   # Scaling of LJ potential for 1 particles
  Dist::Float64 		= 1.0     # Equilibrium distance for all LJ potentials
  Beta::Float64 		= 1.0	  # Inverse temperature
  Gamma::Float64        = 1.0     # Friction parameter

  # The following are inherited from the previous parameters: changing them may break the code.
  NPart::Int 				= length(Patt) *NBeads
  Period::Float64 			= NPart *Dist
  BeadsDist 				= length(Patt) *Dist
  LTotPatt::Array{Int,1} 	= repeat(Patt, NBeads)
  STotPatt::Array{Int,1}	= 1 .-LTotPatt

  # Mass parameters
  MassPatt::Array{Float64,1} 	= LM*Patt + SM*(1 .-Patt)
  MassTotPatt::Array{Float64,1} = repeat(MassPatt, NBeads)
  MInv::Array{Float64,2}		= Diagonal(1 ./MassTotPatt)
  MHalf::Array{Float64,2}		= Diagonal(sqrt.(MassTotPatt))
  BeadMass::Float64				= LM *sum(Patt.==1) + SM *sum(Patt.==0)
  EpsPatt::Array{Float64,1} 	= (sqrt.(LEps) *LTotPatt + sqrt.(SEps) *STotPatt) .*circshift( sqrt.(LEps) *LTotPatt + sqrt.(SEps) *STotPatt, -1)

  # Friction coefficient for Langevin dynamics
  GammaTilde::Array{Float64,2}  = Diagonal(ones(NPart)) *MInv *Gamma
  GammaTildeCG::Array{Float64,2}  = Diagonal(ones(NBeads)) /BeadMass *Gamma

  # Full CG Projections
  Phi::Array{Float64,2}	= kron(Diagonal(ones(NBeads)), ones(Float64, length(Patt))')
  Psi 				    = kron(Diagonal(ones(NBeads)), MassPatt' ./BeadMass)

  # Constrained versions of CG Projections
  PhiConstr::Array{Float64,2}	= kron(Diagonal(ones(NBeads))[1:ConstrainedBeads, :], ones(Float64,length(Patt))')
  PsiConstr 					= kron(Diagonal(ones(NBeads))[1:ConstrainedBeads, :], MassPatt' ./BeadMass)
  MomProj 						= I - PsiConstr'*PhiConstr # this is correct up to a constant pinv(Matrix(PsiConstr))*PhiConstr
  NetMomProj					= I - ones(NPart) *pinv(ones(NPart))
  NetMomProjCG					= I - ones(NBeads) *pinv(ones(NBeads))

  # Relative CG projections between identical beads (full and constrained)
  PhiF::Array{Float64,2} 		= kron(Diagonal(ones(NBeads)), [zeros(Float64, length(Patt)-1) ;1.0])'  # Forces between beads
  PhiFConstr::Array{Float64,2} 	= kron(Diagonal(ones(NBeads))[1:ConstrainedBeads-1, :], [zeros(Float64, length(Patt)-1); 1.0]')  # Forces between constrained beads
  PsiF			        		= (circshift(Phi, (0,length(Patt))) -Phi) /2.;  # Relative momenta
  PsiFConstr					= (circshift(PhiConstr[1:ConstrainedBeads-1, :], (0,length(Patt))) -PhiConstr[1:ConstrainedBeads-1, :]) /2.;  # Relative momenta

  # to compute the CdM of the chain
  vectorCdM = BeadMass * ones(NBeads) / (NBeads * BeadMass) # NB: unless at the moment
end


"""
    U = LennardJones(Q,Dist,EpsPatt)

Lennard-Jones 6-12 potential with minimum at Dist and well depth EpsPatt.
"""
function LennardJones(Q::Array{Float64,1},Dist::Float64,
						EpsPatt::Array{Float64,1})
  U = EpsPatt .*((Dist./Q).^12-2*(Dist./Q).^6);
  return U
end

"""
    G = LennardJonesPrime(Q,Dist,EpsPatt)

Derivative of Lennard-Jones 6-12 potential with minimum at Dist and well depth
EpsPatt.
"""
function LennardJonesPrime(Q::Array{Float64,1},Dist::Float64,EpsPatt::Array{Float64,1})
  G = -12 .*EpsPatt.*(Dist^12 ./(Q.^13)-Dist^6 ./(Q.^7));
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
  PNext = P + 0.5.*F.*dt 	       # B
  QNext = Q + dt.*C.MInv*PNext;    # A
  FNext = Forces(QNext,C);         # Force evaluation
  PNext = PNext + 0.5.*FNext.*dt;  # B
  return PNext,QNext,FNext
end

"""
    (PNext,QNext,FNext) = DetCGDynNVEStep(P,Q,F,MeanForce,C::Chain,dt)

Computes forward step of Deterministic CG Dynamics in the Microcanonical ensemble with timestep dt.
"""
function DetCGDynNVEStep(P::Array{Float64,1},Q::Array{Float64,1},F::Array{Float64,1},MeanForce,C::Chain,dt::Float64)
  PNext = P + 0.5.*F.*dt				# B
  QNext = Q + dt./C.BeadMass.*PNext;    # A
  FNext = MeanForce(QNext);             # Force evaluation
  PNext = PNext + 0.5.*FNext.*dt;       # B
  return PNext,QNext,FNext
end

"""
    (PNext,QNext,FNext) = ConstrainedStep(P,Q,F,C::Chain,dt)

Computes leapfrog step of dynamics with timestep dt, constrained to maintain CG
positions and momenta.
"""
function ConstrainedStep(P,Q,F,C,dt)
  PNext = P + 0.5.*C.MomProj*F.*dt;	 	      # B
  QNext = Q + C.MInv*C.MomProj*PNext.*dt;     # A
  FNext = Forces(QNext,C);			          # Force evaluation
  PNext = PNext + 0.5*C.MomProj*FNext.*dt;	  # B
  return PNext,QNext,FNext
end

"""
    (PNext,QNext,FNext) = LangDynIntegrators(P,Q,F,C::Chain,dt)

Computes forward step of Langevin dynamics with timestep dt using Leimkuhler-Matthews
"BAOAB" integrator; diffusion coefficient is evaluated after performing `half' of the deterministic part of the timestep.
"""
function LangDynIntegrators(P::Array{Float64,1},Q::Array{Float64,1},F::Array{Float64,1},C::Chain,dt::Float64)
  PNext = P + 0.5.*F.*dt 	                      		    # B
  QNext = Q + 0.5.*dt.*C.MInv*PNext;                		# A
  # PNext = C.NetMomProj*OrnUhlenbeckPart(PNext,C,dt)[1];     # O
  PNext = OrnUhlenbeckPart(PNext,C,dt)[1];                  # O
  QNext = QNext+0.5.*dt*C.MInv*PNext;            		    # A
  FNext = Forces(QNext,C);                       		    # Force evaluation
  PNext = PNext + 0.5.*FNext.*dt;                		    # B
  return PNext,QNext,FNext
end

"""
    (PNext,QNext,FNext) = LangDynConstrainedStep(P,Q,F,C::Chain,dt)

Computes forward step of constrained Langevin dynamics with timestep dt using Leimkuhler-Matthews "g-BAOAB" integrator;
diffusion coefficient is evaluated after performing `half' of the deterministic part of the timestep.
"""
function LangDynConstrainedStep(P,Q,F,C,dt)
	PNext = P + 0.5.*C.MomProj*F.*dt;                         # B
	QNext = Q + 0.5.*C.MInv*C.MomProj*PNext.*dt;              # A
	PNext = C.NetMomProj*OrnUhlenbeckPart(PNext,C,dt)[2];     # O
	QNext = QNext + 0.5*C.MInv*C.MomProj*PNext.*dt;           # A
	FNext = Forces(QNext,C);                                  # Force evaluation
	PNext = PNext + 0.5.*C.MomProj*FNext.*dt;                 # B
	return PNext,QNext,FNext
end

function OrnUhlenbeckPart(P::Array{Float64,1},C::Chain,dt::Float64)
 	# Calculate exact Ornstein-Uhlenbeck step using decomposition
	c1 = exp(-C.GammaTilde*dt);
	c3 = sqrt(1/C.Beta)*sqrt(Diagonal(ones(C.NPart))-c1^2);

	# Generate new momenta and return
	return c1*P+c3*C.MHalf*randn(C.NPart), c1*P+c3*C.MomProj*C.MHalf*randn(C.NPart);
end

"""
    (PNext,QNext,FNext) = ApproxDynStep(P,Q,FGamma,C::Chain,dt)

Computes forward step of Langevin dynamics with timestep dt using Leimkuhler-Matthews
"BAOAB" integrator, adapted to non-constant diffusion; diffusion coefficient
is evaluated after performing `half' of the deterministic part of the timestep.
"""
function ApproxDynStep(P::Array{Float64,1},Q::Array{Float64,1},F::Array{Float64,1},MeanForce,FlucStep,C::Chain,dt)
    # println("  P (input) = ",P);
    # println("  Q (input) = ",Q);
    # println("  F (input) = ",F);
    PNext = P + 0.5.*F.*dt;			              # B
    QNext = Q + 0.5.*dt/C.BeadMass.*PNext;        # A
    # println("     PNExt = ",PNext);
    # println("     QNext = ",QNext);
    PNext = FlucStep(PNext,QNext);
    QNext = QNext + 0.5.*dt/C.BeadMass.*PNext;    # A
    # println("   Q = ",QNext);
    FNext = MeanForce(QNext);	  			      # Force evaluation
    # println("  F (output) = ",FNext);
    # println(" ...end");
    PNext = PNext + 0.5.*FNext.*dt;               # B
    return PNext,QNext,FNext
end

"""
    (PNext,QNext,FNext) = DetCGDynNVTStep(P,Q,F,MeanForce,C::Chain,dt)

Computes forward step of Deterministic CG Dynamics in the Canonical ensemble with timestep dt using Leimkuhler-Matthews
"BAOAB" integrator; diffusion coefficient is evaluated after performing `half' of the deterministic part of the timestep.
"""
function DetCGDynNVTStep(P::Array{Float64,1},Q::Array{Float64,1},F::Array{Float64,1},MeanForce,C::Chain,dt::Float64)
  PNext = P + 0.5.*F.*dt 	                      		    # B
  QNext = Q + 0.5.*dt./C.BeadMass.*PNext;                   # A
  # PNext = C.NetMomProjCG*OrnUhlenbeckPartDCG(PNext,C,dt);   # O
  PNext = OrnUhlenbeckPartDCG(PNext,C,dt);                  # O
  QNext = QNext+0.5.*dt/C.BeadMass.*PNext;            		# A
  FNext = MeanForce(QNext);                       		    # Force evaluation
  PNext = PNext + 0.5.*FNext.*dt;                		    # B
  return PNext,QNext,FNext
end

function OrnUhlenbeckPartDCG(P::Array{Float64,1},C::Chain,dt::Float64)
 	# Calculate exact Ornstein-Uhlenbeck step using decomposition
	c1 = exp(-C.GammaTildeCG*dt);
	c3 = sqrt(1/C.Beta)*sqrt(Diagonal(ones(C.NBeads))-c1^2);

	# Generate new momenta and return
	return c1*P+c3*sqrt(C.BeadMass)*randn(C.NBeads)
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
  P0 = sqrt.(C.MassTotPatt ./ C.Beta).*randn(C.NPart);
  # Set total momentum to 0 to avoid drift.
  P0 .-= mean(P0);
  P0 *= sqrt.((C.NPart / C.Beta) / KineticEnergy(P0, C));
  Q0 = convert(Array{Float64,1}, range(0.0, stop = (C.NPart-1) * C.Dist, length = C.NPart));
  Q0 = Q0 + 0.001 * C.Beta * randn(C.NPart)

  # Fix reference energy
  E0 = Potential(Q0, C);
  H0 = KineticEnergy(P0, C) + E0;

  # Compute initial forces
  F0 = Forces(Q0, C);

  # Extract IC for CG variables
  PCG0 = C.Phi * P0;
  QCG0 = C.Psi * Q0;

  return P0, Q0, F0, H0, PCG0, QCG0
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
