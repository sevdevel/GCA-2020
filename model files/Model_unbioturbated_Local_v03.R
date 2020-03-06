###############################################################################
# Model of early diagenesis
# Type: Carbon,sulfur and iron cycling  
# Author: Filip Meysman (filip.meysman@nioz.nl)
#         Sebastiaan van de Velde (sevdevel@vub.ac.be) 
# Affiliation: NIOZ, Korringaweg 7, 4401 NT, Yerseke 
#              VUB, Pleinlaan 2, 1050 BE, Brussel
###############################################################################

# public packages
require(AquaEnv)
require(rootSolve)
require(deSolve)
require(shape)
require(ReacTran)
require(gsw)
require(oce)
require(seacarb)
require(marelac)

# non-public packages
require(diagenesis)

#=============================================================================
# Definition of specific functions
#=============================================================================

#-----------------------------------------------------------------------------
# Saturation function
#-----------------------------------------------------------------------------

FSAT <- function(C,K,n) (C/K)^n/((C/K)^n + 1)

#-----------------------------------------------------------------------------
# Function: tortuosity 
#-----------------------------------------------------------------------------

tortuosity <- function(por) 1-2*log(por)

#-----------------------------------------------------------------------------
# Function: init.diffusion.coefficient 
#-----------------------------------------------------------------------------

init.diffusion.coefficient <- function (species,S,TC,P,conv.fac,grid,tort.grid)
{
  Dmol <- conv.fac*diffcoeff(S=S,t=TC,P=P,species=species)[[species]]
  D.grid <- setup.prop.1D(value=Dmol,grid=grid)
  D.grid$mid <- D.grid$mid/tort.grid$mid
  D.grid$int <- D.grid$int/tort.grid$int
  return(D.grid)
}

#=============================================================================
# Definition of objects: species, reactions, elements, simulations
#=============================================================================

# SL = list of species incorporated in the model 
SL <- list()

# chemical species incorporated in the model 

SL <- initialise.list.species(SL,"1D",name="OC.f",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="OC.s",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="OC.r",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeOOH.f",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeOOH.a",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeS",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeS2",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="X.Fe",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="X.NH4",phase="solid",type="chemical")

SL <- initialise.list.species(SL,"1D",name="HCO3",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="NH4",phase="solute",type="chemical")

SL <- initialise.list.species(SL,"1D",name="O2",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="SO4",phase="solute",type="chemical")

SL <- initialise.list.species(SL,"1D",name="Fe",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="HS",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="CH4",phase="solute",type="chemical")

SL <- initialise.list.species(SL,"1D",name="H",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="H2O",phase="solute",type="chemical")

# Composite species incorporated in the model 

SL <- initialise.list.species(SL,"1D",name="OC",phase="solid",type="composite")
SL <- initialise.list.species(SL,"1D",name="FeOOH",phase="solid",type="composite")

# Artificial species incorporated in the model 

SL <- initialise.list.species(SL,"1D",name="Omega.FeS",phase="solute",type="artificial")

# Initialise species matrix 

spec.mat <- initialise.species.matrix(SL)
spec.mat["OC",c("OC.f","OC.s","OC.r")] <- c(1,1,1)
spec.mat["FeOOH",c("FeOOH.f","FeOOH.a")] <- c(1,1)

spec.table <- as.data.frame(spec.mat)
remove(spec.mat)


# List of reactions included 

RL <- list()

#---------------------
# mineralization terms
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="Cmin.f",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="Cmin.s",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="Cmin.r",type="kinetic")

RL <- initialise.list.reaction(RL,"1D",name="Cmin.ftot",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="Cmin.stot",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="Cmin.rtot",type="kinetic")

#---------------------
# individual mineralization terms
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="AR.f",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="FR.f",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SR.f",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="MG.f",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="AR.s",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="FR.s",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SR.s",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="MG.s",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="AR.r",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="FR.r",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SR.r",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="MG.r",type="kinetic")

#---------------------
# aerobic oxidation terms
# excluding reactions with sulfide minerals
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="AMO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="FIO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SIO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="CSO",type="kinetic")

#---------------------
# anaerobic oxdiation terms
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="CSFO.f",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="CSFO.a",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SMO",type="kinetic")

#---------------------
# sulfide mineral reactions
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="ISP",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="ISD",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="ISO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="PyP",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="PyO",type="kinetic")

#---------------------
# iron oxide ageing
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="IOA",type="kinetic")

# List of elements included 

EL <- list() # create empty element object list

EL <- initialise.list.element(EL,"1D",name="C")
EL <- initialise.list.element(EL,"1D",name="O")
EL <- initialise.list.element(EL,"1D",name="H")
EL <- initialise.list.element(EL,"1D",name="S")
EL <- initialise.list.element(EL,"1D",name="Fe")
EL <- initialise.list.element(EL,"1D",name="N")
EL <- initialise.list.element(EL,"1D",name="e")

# Initialise element matrix 

elt.mat <- initialise.element.matrix(SL, EL)

elt.mat["OC.f",c("C","H","N")] <- c(1,25/11,1/11)
elt.mat["OC.s",c("C","H","N")] <- c(1,25/11,1/11)
elt.mat["OC.r",c("C","H","N")] <- c(1,25/11,1/11)

elt.mat["FeOOH.f",c("Fe","O","H")] <- c(1,2,1)
elt.mat["FeOOH.a",c("Fe","O","H")] <- c(1,2,1)
elt.mat["FeS",c("S","Fe")] <- c(1,1)
elt.mat["FeS2",c("S","Fe")] <- c(2,1)

elt.mat["HCO3",c("H","C","O","e")] <- c(1,1,3,1)
elt.mat["NH4",c("N","H","e")] <- c(1,4,-1)
elt.mat["X.NH4",c("N","H","e")] <- c(1,4,-1)

elt.mat["O2",c("O")] <- c(2)
elt.mat["SO4",c("S","O","e")] <- c(1,4,2)

elt.mat["Fe",c("Fe","e")] <- c(1,-2)
elt.mat["X.Fe",c("Fe","e")] <- c(1,-2)
elt.mat["HS",c("H","S","e")] <- c(1,1,1)
elt.mat["CH4",c("C","H")] <- c(1,4)

elt.mat["H",c("H","e")] <- c(1,-1)
elt.mat["H2O",c("H","O")] <- c(2,1)

elt.table <- as.data.frame(elt.mat)
remove(elt.mat)

################################################################################
# initialise.parameters function
# This function initialises a default set of parameters (baseline simulation),
# but this function is also used when updating the parameter set.
################################################################################

initialise.parameters <- function(PL=NULL)
  
  # Model parameters are intialised on the first call (PL = NULL), and updated 
  # when a parmeter list is supplied (PL not NULL)
  
  # begin initialise.parameters
{
  
  #=============================================================================
  # Physical units 
  #=============================================================================
  
  if (is.null(PL$units$M)) PL$units$M <- "umol" # Mass unit
  if (is.null(PL$units$T)) PL$units$T <- "yr"   # Time unit
  if (is.null(PL$units$L)) PL$units$L <- "cm"   # Length unit
  
  # Limiting concentration: numerical limit on consumption in kinetic reactions
  PL$C.lim <- 1E-3
  
  #=============================================================================
  # Model domain and grid definition
  #=============================================================================
  
  if (is.null (PL$L)) PL$L <- 30   # depth of sediment domain [cm]
  if (is.null (PL$N)) PL$N <- 400  # number of grid layers
  
  # Call "setup.grid.1D" from ReacTran
  if (is.null (PL$grid.type)) PL$grid.type <- "uniform"  # number of grid layers
  if (is.null (PL$dx.1)) PL$dx.1 <- PL$L/(10*PL$N)
  #if (PL$grid.type == "uniform")
  #{
    # even grid - all depth layers have same thickness 
    #PL$grid <- setup.grid.1D(x.up=0,L=PL$L,N=PL$N)
  #} else {
    # uneven grid - higher resolution near the sediment water interface
    PL$grid <- setup.grid.1D(x.up=0,L=PL$L,N=PL$N,dx.1=PL$dx.1)
  #}
  
  #=============================================================================
  # Sediment parameters: 
  #=============================================================================
  
  # Environmental parameters
  
  if (is.null(PL$S))  PL$S  <- 32    # salinity
  if (is.null(PL$TC)) PL$TC <- 10    # temperature [deg C]
  if (is.null(PL$P))  PL$P  <- 1.013 # pressure [bar]
  
  # Porosity profile 
  
  if (is.null(PL$por.0)) PL$por.0 <- 0.96     # porosity at the sediment-water interface
  if (is.null(PL$por.inf)) PL$por.inf  <- 0.78   # asymptotic porosity at depth
  if (is.null(PL$por.x.att)) PL$por.x.att <- 9    # attenuation depth [cm]
    
  PL$por.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$por.0,y.inf=PL$por.inf,x.L=0,x.att=PL$por.x.att)
  PL$svf.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=(1-PL$por.0),y.inf=(1-PL$por.inf),x.L=0,x.att=PL$por.x.att)
  
  
  # Initialisation tortuosity 

  PL$tort.grid <- setup.prop.1D(value=0,grid=PL$grid)
  PL$tort.grid$mid <- tortuosity(PL$por.grid$mid)
  PL$tort.grid$int <- tortuosity(PL$por.grid$int)

  # Fixed profiles: pH
  
  if (is.null(PL$pH)) PL$pH      <- 6.75
  PL$pH.grid <-  setup.prop.1D(value=PL$pH,grid=PL$grid)
  
  #=============================================================================
  # Diffusion coefficients 
  # Uses routine 'diffcoeff' from 'marelac' package
  #=============================================================================
  
  c.fac <- 10000*(3600*24*365.25) # conversion from [m2 s-1] to [cm2 yr-1]
  
  # Diffusion coefficients from marelac package
  PL$D.HCO3.grid <- init.diffusion.coefficient(species="HCO3",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.NH4.grid <- init.diffusion.coefficient(species="NH4",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  
  PL$D.O2.grid <- init.diffusion.coefficient(species="O2",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.SO4.grid <- init.diffusion.coefficient(species="SO4",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  
  PL$D.Fe.grid <- init.diffusion.coefficient(species="Fe",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.HS.grid <- init.diffusion.coefficient(species="HS",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.CH4.grid <- init.diffusion.coefficient(species="CH4",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  
  PL$D.H2O.grid <- init.diffusion.coefficient(species="H2O",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.H.grid <- init.diffusion.coefficient(species="H",S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  
  #=============================================================================
  # Transport parameters
  #=============================================================================
  
  # Diffusive boundary layer
  
  if (is.null (PL$x.DBL)) PL$x.DBL <- 0 # thickness of DBL [cm]
  
  # Advective velocities 
  
  if (is.null(PL$rho.sed)) PL$rho.sed <- 2.0  # density solid sediment [g cm-3]
  if (is.null(PL$u.inf)) PL$u.inf <- 0.1/PL$rho.sed/(1-PL$por.inf)  # sedimentation velocity pore water [cm yr-1]
  if (is.null(PL$v.inf)) PL$v.inf <- 0.1/PL$rho.sed/(1-PL$por.inf)  # sedimentation velocity solids [cm yr-1]
  
  PL$SedFlux <- 0.1  # sedimentation flux [g cm-2 yr-1]
  
  PL$v.grid <- setup.compaction.1D(v.inf=PL$v.inf,por.0=PL$por.0,por.inf=PL$por.inf,por.grid=PL$por.grid)$v # advection velocity solid
  PL$u.grid <- setup.compaction.1D(v.inf=PL$v.inf,por.0=PL$por.0,por.inf=PL$por.inf,por.grid=PL$por.grid)$u # advection velocity pore water
  
  # Bioturbation profile
  
  if (is.null(PL$x.L))  PL$x.L    <- 0
  if (is.null(PL$x.att))PL$x.att  <- 0
  if (is.null(PL$y.inf))PL$y.inf  <- 0
  if (is.null(PL$Db.0)) PL$Db.0   <- 0     # biodiffusion coefficient Db at SWI [cm2 yr-1]
  PL$Db.grid <- setup.prop.1D(func=p.sig,y.0=PL$Db.0,y.inf=PL$y.inf,x.att=PL$x.att,x.L=PL$x.L, grid = PL$grid)
  
  # Irrigation profile
  
  if (is.null(PL$irr.0)) PL$irr.0 <- 0     # irrigation rate at SWI [yr-1]
  if (is.null(PL$x.irr)) PL$x.irr <- 3     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.at)) PL$irr.att <- 2
  if (is.null(PL$irr.inf)) PL$irr.inf <- 0
  
  if (is.null(PL$irr.Fe))   PL$irr.Fe   <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.HS))   PL$irr.HS   <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.NH4))  PL$irr.NH4  <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.N2))   PL$irr.N2   <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.HCO3)) PL$irr.HCO3 <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.SO4))  PL$irr.SO4  <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.O2))   PL$irr.O2   <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.H2))   PL$irr.H2   <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.NO3))  PL$irr.NO3  <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.CH4))  PL$irr.CH4  <- 1     # irrigation rate at SWI [yr-1]
  
  #PL$irr.grid <- setup.prop.1D(func=p.exp,y.0=PL$irr.0,y.inf=0,x.L=0,x.att=PL$x.irr,grid=PL$grid)
  PL$irr.grid <- setup.prop.1D(func=p.sig,y.0=PL$irr.0,y.inf=PL$irr.inf,x.att=PL$irr.att,x.L=PL$x.irr, grid = PL$grid)
  
  #=============================================================================
  # Reaction parameters
  #=============================================================================
  
  # Decay constants organic matter
  
  if (is.null(PL$k.f))         PL$k.f   <- 10     # decay constant of organic matter [yr-1]
  if (is.null(PL$k.s))         PL$k.s   <- 0.04      # decay constant of organic matter [yr-1]
  if (is.null(PL$priming.factor)) PL$priming.factor <- 0     # stimulation factor of the presence of reactive OM
  if (is.null(PL$priming.ref))    PL$priming.ref    <- 1000  # saturation factor of the presence of reactive OM
  if (is.null(PL$aerobic.factor)) PL$aerobic.factor <- 1     # stimulation factor aerobic respiration
  
  if (is.null(PL$meth.factor)) PL$meth.factor <- 1  # difference factor between aerobic mineralization and methanogenesis
  if (is.null(PL$k.r))         PL$k.r     <- 0      # decay constant of organic matter [yr-1] (Meysman et al., 2003)
    
  if (is.null(PL$K_O2))   PL$K_O2     <- 0.001         # Monod constant O2 consumption [umol cm-3 or mM] (Meysman et al., 2003)
  if (is.null(PL$K_NO3))  PL$K_NO3    <- 0.001         # Monod constant NO3 consumption [umol cm-3 or mM] (Meysman et al., 2003)
  if (is.null(PL$K_FeOOH))PL$K_FeOOH <- 0.4*PL$rho.sed  # Monod constant FeIII reduction [umol cm-3 or mM] (this paper -> 2 fractions, avoid too much overlap)
  if (is.null(PL$K_SO4))  PL$K_SO4   <- 0.9           # Monod constant SO4 reduction [umol cm-3 or mM] (Meysman et al., 2003)
  
  if (is.null(PL$R_CN))   PL$R_CN    <- 35       # from data
  
  # aerobic oxidation reactions
  # excluding sulfide mineral reactions
  
  if (is.null(PL$k.FIO)) PL$k.FIO <- 1E+07  # kinetic constant Ferrous iron oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.SIO)) PL$k.SIO <- 1E+07  # kinetic constant Ferrous iron oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.CSO)) PL$k.CSO  <- 1E+07    # kinetic constant sulfide oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.AMO)) PL$k.AMO <- 1E+04  # kinetic constant aerobic methane oxidation [µmol-1 cm3 yr-1] (Meysman et al., 2003)

  # anaerobic oxidation reactions
  
  if (is.null(PL$k.CSFOf))PL$k.CSFOf <- 494     # kinetic constant sulfide oxidation by fresh FeIII [umol-1 cm3 yr-1] (Poulton et al., 2004)
  if (is.null(PL$k.CSFOa))PL$k.CSFOa <- 3.6     # kinetic constant sulfide oxidation by aged FeIII [umol-1 cm3 yr-1] (Poulton et al., 2004)
  if (is.null(PL$k.SMO))  PL$k.SMO   <- 1E+04     # kinetic constant Sulfate methane oxidation [µmol-1 cm3 yr-1] (-)
  
  # reactions involving sulfide minerals
  
  if (is.null(PL$k.ISO))  PL$k.ISO  <- 1E+07    # kinetic constant iron sulfide oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996)
  if (is.null(PL$k.PyO))  PL$k.PyO  <- 9.47    # kinetic constant Pyrite oxidation [umol-1 cm3 yr-1] (Berg et al., 2003)
  if (is.null(PL$k.PyP))  PL$k.PyP  <- 0.05    # kinetic constant Pyrite precipitation [umol-1 cm3 yr-1] (this paper)
    
  if (is.null(PL$k.ISP)) PL$k.ISP <- 1E+04   # kinetic constant FeS precipitation [umol-1 cm3 yr-1] (Meysman et al., 2015)
  if (is.null(PL$k.ISD)) PL$k.ISD <- 3       # kinetic constant FeS dissolution [yr-1] (Meysman et al., 2015)
  if (is.null(PL$n.ISP)) PL$n.ISP <- 1       # kinetic exponent FeS precipitation [] (Meysman et al., 2003)
  if (is.null(PL$n.ISD)) PL$n.ISD <- 1       # kinetic exponent FeS dissolution [] (Meysman et al., 2003)
  
  if (is.null(PL$K_IS))  PL$K_IS  <- (10^-PL$pH)*1000*3.16    # Saturation constant [umol-1 cm3] (Rickard, 2006)
  
  # iron oxide ageing and adsorption reactions 
  
  if (is.null(PL$k.IOA)) PL$k.IOA <- 0.57  # kinetic constant iron oxide ageing [yr-1] (Berg et al., 2003)
  
  #if (is.null(PL$k.FIS)) PL$k.FIS <- 1    # kinetic constant Ferrous iron sorption  [yr-1] (Berg et al., 2003)
  if (is.null(PL$K.FIS)) PL$K.FIS <- PL$rho.sed*268   # Equilibrium constant Ferrous iron sorption [] (Berg et al., 2003)
  if (is.null(PL$K.AMS)) PL$K.AMS <- PL$rho.sed*1.6   # Equilibrium constant ammonium sorption [] (Berg et al., 2003; Mackin and Aller, 1984)
  
  #=============================================================================
  # Boundary conditions
  #=============================================================================
  
  # Upper boundary fluxes

  F.OC.f <- 15 # mmol m-2 d-1
  if (is.null(PL$F.OC.f)) PL$F.OC.f <- F.OC.f*365.25/10 # conversion to umol cm-2 yr-1 
  
  F.OC.s <- 22 # mmol m-2 d-1
  if (is.null(PL$F.OC.s)) PL$F.OC.s <- F.OC.s*365.25/10 # conversion to umol cm-2 yr-1 
  
  F.OC.r <- 7 # mmol m-2 d-1
  if (is.null(PL$F.OC.r)) PL$F.OC.r <- F.OC.r*365.25/10 # conversion to umol cm-2 yr-1 
  
  F.FeS <- 0 # mmol m-2 d-1
  if (is.null(PL$F.FeS)) PL$F.FeS <- F.FeS*365.25/10 # conversion to umol cm-2 yr-1 
  
  F.FeS2 <- 0 # mmol m-2 d-1
  if (is.null(PL$F.FeS2)) PL$F.FeS2 <- F.FeS2*365.25/10 # conversion to umol cm-2 yr-1 
  
  # 1% Fe content of mineral deposition 
  F.FeOOH.f <- 16/36.525#1E06*0.01*PL$SedFlux/(55.847) # mmol m-2 d-1
  if (is.null(PL$F.FeOOH.f)) PL$F.FeOOH.f <- F.FeOOH.f*365.25/10 # conversion to umol cm-2 yr-1 
  F.FeOOH.a <- 16/36.525 # mmol m-2 d-1
  if (is.null(PL$F.FeOOH.a)) PL$F.FeOOH.a <- F.FeOOH.a*365.25/10 # conversion to umol cm-2 yr-1 
  
  F.X.Fe <- 0 # mmol m-2 d-1
  if (is.null(PL$F.X.Fe)) PL$F.X.Fe <- F.X.Fe*365.25/10 # conversion to umol cm-2 yr-1 
  F.X.NH4 <- 0 # mmol m-2 d-1
  if (is.null(PL$F.X.NH4)) PL$F.X.NH4 <- F.X.NH4*365.25/10 # conversion to umol cm-2 yr-1 
  
  # Upper boundary concentrations

  if (is.null(PL$HCO3.ow)) PL$HCO3.ow <- 3.9 # SumCO2 concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$NH4.ow)) PL$NH4.ow   <- 0.01
  
  if (is.null(PL$O2.ow)) PL$O2.ow   <- 0.28 # O2 concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$SO4.ow)) PL$SO4.ow <- 20 #140.0*PL$S/((32.065+4*15.999)*1.80655)
  
  if (is.null(PL$Fe.ow)) PL$Fe.ow   <- 0
  if (is.null(PL$HS.ow)) PL$HS.ow   <- 0   # SumH2S concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$CH4.ow)) PL$CH4.ow <- 0
  
  if (is.null(PL$H2O.ow)) PL$H2O.ow <- 0
  if (is.null(PL$H.ow)) PL$H.ow     <- 0 # Proton concentration bottom water, without equilibrium with H2O
    
  # Boundary condition lower boundary 
  # All species: no gradient
  
  if (is.null(PL$HCO3.ds)) PL$HCO3.ds <- NA
  if (is.null(PL$NH4.ds)) PL$NH4.ds <- NA
  
  if (is.null(PL$O2.ds)) PL$O2.ds   <- NA
  if (is.null(PL$SO4.ds)) PL$SO4.ds <- NA
    
  if (is.null(PL$Fe.ds)) PL$Fe.ds   <- NA
  if (is.null(PL$HS.ds)) PL$HS.ds   <- NA
  if (is.null(PL$CH4.ds)) PL$CH4.ds <- NA
  
  if (is.null(PL$H2O.ds)) PL$H2O.ds <- NA
  if (is.null(PL$H.ds)) PL$H.ds     <- NA
  
  #=============================================================================
  # Flags 
  #=============================================================================
  
  if (is.null(PL$simulation.type)) PL$simulation.type <- "direct.steady.state"
  #simulation.type <- "dynamic.steady.state"
  #simulation.type <- "time.dependent"
  
  ##############################################################################
  return(PL)
  # end initialise.parameters
}
##############################################################################

#=============================================================================
# Model formulation
#=============================================================================

CSFe.model <- function (t,state,parameters,summary.call=FALSE) 
{
  
  with(as.list(c(parameters)),{
    
    #---------------------------------------------------------------------------
    # Initialisation of 
    # SV : matrix of state variables
    # OL : Object List
    #---------------------------------------------------------------------------
    
    SV <- matrix(nrow=N,ncol=N.var,data=state)
    for (i in 1:N.var) OL[[var.names[i]]]$C <- SV[,i]
    
    #---------------------------------------------------------------------------
    OL <- within(OL,{
    #---------------------------------------------------------------------------
      
    #=========================================================================
    # Initialisation of concentration depth profiles
    #=========================================================================
    #-------------------------------------------------------------------------
    # Initialisation of species with fixed concentration depth profiles 
    #-------------------------------------------------------------------------
    
    if (!("OC.f" %in% var.names)) OC.f$C <- OC.f.fix
    if (!("OC.s" %in% var.names)) OC.s$C <- OC.s.fix
    if (!("OC.r" %in% var.names)) OC.r$C <- OC.r.fix
    
    X.Fe$C <- K.FIS*Fe$C
    X.NH4$C <- K.AMS*NH4$C
    
    #-------------------------------------------------------------------------
    # Upper and lower boundary conditions
    #-------------------------------------------------------------------------

    HCO3$C.up <- HCO3.ow
    NH4$C.up  <- NH4.ow
    
    O2$C.up   <- O2.ow
    SO4$C.up  <- SO4.ow
    
    Fe$C.up   <- Fe.ow
    HS$C.up   <- HS.ow
    CH4$C.up  <- CH4.ow
    
    H2O$C.up  <- H2O.ow
    H$C.up    <- H.ow
    
    if (is.na(HCO3.ds)) HCO3$C.down <- HCO3$C[N] else HCO3$C.down <- HCO3.ds 
    if (is.na(NH4.ds))  NH4$C.down <- NH4$C[N] else NH4$C.down <- NH4.ds 
    
    if (is.na(O2.ds))   O2$C.down  <- O2$C[N] else O2$C.down <- O2.ds 
    if (is.na(SO4.ds))  SO4$C.down <- SO4$C[N] else SO4$C.down <- SO4.ds 
    
    if (is.na(Fe.ds))   Fe$C.down  <- Fe$C[N] else Fe$C.down <- Fe.ds 
    if (is.na(HS.ds))   HS$C.down  <- HS$C[N] else HS$C.down <- HS.ds 
    if (is.na(CH4.ds))  CH4$C.down <- CH4$C[N] else CH4$C.down <- CH4.ds 
    
    if (is.na(H2O.ds))  H2O$C.down <- H2O$C[N] else H2O$C.down <- H2O.ds 
    if (is.na(H.ds))    H$C.down   <- H$C[N] else H$C.down <- H.ds 
    
    
    #=========================================================================
    # TRANSPORT terms using tran.1D from ReacTran
    #=========================================================================
    
    if (summary.call) {ind <- c("tran","C.up","C.down","dif.flux","adv.flux","flux","flux.up","flux.down")} else {ind <- c("tran","flux.up","flux.down")}
    
    # Transport terms: Solids [umol cm-3 solid yr-1]
    
    OC.f[ind]  <- tran.1D(C=OC.f$C,flux.up=F.OC.f,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    OC.s[ind]  <- tran.1D(C=OC.s$C,flux.up=F.OC.s,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    OC.r[ind]  <- tran.1D(C=OC.r$C,flux.up=F.OC.r,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    FeOOH.f[ind]<- tran.1D(C=FeOOH.f$C,flux.up=F.FeOOH.f,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    FeOOH.a[ind]<- tran.1D(C=FeOOH.a$C,flux.up=F.FeOOH.a,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    FeS[ind]   <- tran.1D(C=FeS$C,flux.up=F.FeS,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    FeS2[ind]  <- tran.1D(C=FeS2$C,flux.up=F.FeS2,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    X.Fe[ind]  <- tran.1D(C=X.Fe$C,flux.up=F.X.Fe,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    X.NH4[ind] <- tran.1D(C=X.NH4$C,flux.up=F.X.NH4,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    # Transport terms: Solutes [umol cm-3 pore water yr-1]

    HCO3[ind] <- tran.1D(C=HCO3$C,C.up=HCO3$C.up,C.down=HCO3$C.down,v=u.grid,D=D.HCO3.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    NH4[ind]  <- tran.1D(C=NH4$C,C.up=NH4$C.up,C.down=NH4$C.down,v=u.grid,D=D.NH4.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    O2[ind]  <- tran.1D(C=O2$C,C.up=O2$C.up,C.down=O2$C.down,v=u.grid,D=D.O2.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    #N2[ind]  <- tran.1D(C=N2$C,C.up=N2$C.up,C.down=N2$C.down,v=u.grid,D=D.N2.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    #NO3[ind] <- tran.1D(C=NO3$C,C.up=NO3$C.up,C.down=NO3$C.down,v=u.grid,D=D.NO3.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    SO4[ind] <- tran.1D(C=SO4$C,C.up=SO4$C.up,C.down=SO4$C.down,v=u.grid,D=D.SO4.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    Fe[ind]  <- tran.1D(C=Fe$C,C.up=Fe$C.up,C.down=Fe$C.down,v=u.grid,D=D.Fe.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    HS[ind]  <- tran.1D(C=HS$C,C.up=HS$C.up,C.down=HS$C.down,v=u.grid,D=D.HS.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    CH4[ind] <- tran.1D(C=CH4$C,C.up=CH4$C.up,C.down=CH4$C.down,v=u.grid,D=D.CH4.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    H2O[ind] <- tran.1D(C=H2O$C,C.up=H2O$C.up,C.down=H2O$C.down,v=u.grid,D=D.H2O.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    H[ind]   <- tran.1D(C=H$C,C.up=H$C.up,C.down=H$C.down,v=u.grid,D=D.H.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
        
    #-----------------------------------------------------------------------------
    # Saturation state calculation
    #-----------------------------------------------------------------------------
    
    # Saturation state FeS 
    
    Omega.FeS$C      <- Fe$C*HS$C/(K_IS)#*FSAT(Fe$C,C.lim/100,5)*FSAT(HS$C,C.lim/100,5)
    Omega.FeS$C.up   <- Fe.ow*HS$C.up/K_IS
    Omega.FeS$C.down <- Fe$C.down*HS$C.down/K_IS
    
    #=========================================================================
    # Rates of Kinetic reactions 
    # Expressed per volume of bulk sediment [umol cm-3 yr-1]
    #=========================================================================
    
    # Mineralisation reactions
    
#    Cmin.f$R <- svf.grid$mid*k.f*OC.f$C*FSAT(OC.f$C,C.lim,5)
#    Cmin.s$R <- svf.grid$mid*k.s*OC.s$C*FSAT(OC.s$C,C.lim,5)

    Cmin.f$R <- k.f*svf.grid$mid*OC.f$C*FSAT(OC.f$C,C.lim,5)
    Cmin.s$R <- k.s*(1+priming.factor*(1-exp(-OC.f$C*(OC.f$C>0)/priming.ref)))*svf.grid$mid*OC.s$C*FSAT(OC.s$C,C.lim,5) 
    Cmin.r$R <- k.r*(1+priming.factor*(1-exp(-OC.f$C*(OC.f$C>0)/priming.ref)))*svf.grid$mid*OC.r$C*FSAT(OC.r$C,C.lim,5)

    O2.lim <- O2$C/(O2$C+K_O2)*FSAT(O2$C,C.lim,5) #*(O2$C>sqrt(.Machine$double.eps))
    O2.inh <- K_O2/(O2$C+K_O2)
    FeOOH.lim <- FeOOH.f$C/(FeOOH.f$C+K_FeOOH)*FSAT(FeOOH.f$C,C.lim,5)
    FeOOH.inh <- K_FeOOH/(FeOOH.f$C+K_FeOOH)
    SO4.lim <- SO4$C/(SO4$C+K_SO4)*FSAT(SO4$C,C.lim,5)
    SO4.inh <- K_SO4/(SO4$C+K_SO4)
    
    AR.f$R <- O2.lim*Cmin.f$R  # Aerobic respiration
    FR.f$R <- FeOOH.lim*O2.inh*Cmin.f$R # FeIII reduction
    SR.f$R <- SO4.lim*O2.inh*FeOOH.inh*Cmin.f$R # Sulfate reduction
    MG.f$R <- SO4.inh*O2.inh*FeOOH.inh*Cmin.f$R # methanogenesis

    AR.s$R <- aerobic.factor*O2.lim*Cmin.s$R  # Aerobic respiration
    FR.s$R <- FeOOH.lim*O2.inh*Cmin.s$R # FeIII reduction
    SR.s$R <- SO4.lim*O2.inh*FeOOH.inh*Cmin.s$R # Sulfate reduction
    MG.s$R <- SO4.inh*O2.inh*FeOOH.inh*Cmin.s$R # methanogesis

    AR.r$R <- aerobic.factor*O2.lim*Cmin.r$R  # Aerobic respiration
    FR.r$R <- FeOOH.lim*O2.inh*Cmin.r$R # FeIII reduction
    SR.r$R <- SO4.lim*O2.inh*FeOOH.inh*Cmin.r$R # Sulfate reduction
    MG.r$R <- SO4.inh*O2.inh*FeOOH.inh*Cmin.r$R # methanogesis

    #a.f <- aerobic.factor*O2.lim/(aerobic.factor*O2.lim + FeOOH.lim*O2.inh + SO4.lim*FeOOH.inh*O2.inh + SO4.inh*FeOOH.inh*O2.inh)
    #f.f <- (FeOOH.lim*O2.inh)/(aerobic.factor*O2.lim + FeOOH.lim*O2.inh + SO4.lim*FeOOH.inh*O2.inh + SO4.inh*FeOOH.inh*O2.inh)
    #s.f <- (SO4.lim*FeOOH.inh*O2.inh)/(aerobic.factor*O2.lim + FeOOH.lim*O2.inh + SO4.lim*FeOOH.inh*O2.inh + SO4.inh*FeOOH.inh*O2.inh)
    #m.f <- (SO4.inh*FeOOH.inh*O2.inh)/(aerobic.factor*O2.lim + FeOOH.lim*O2.inh + SO4.lim*FeOOH.inh*O2.inh + SO4.inh*FeOOH.inh*O2.inh)
     
    #AR.f$R <- a.f*Cmin.f$R # Aerobic respiration
    #FR.f$R <- f.f*Cmin.f$R # FeIII reduction
    #SR.f$R <- s.f*Cmin.f$R # Sulfate reduction
    #MG.f$R <- m.f*Cmin.f$R # methanogenesis
    
    #AR.s$R <- a.f*Cmin.s$R # Aerobic respiration
    #FR.s$R <- f.f*Cmin.s$R # FeIII reduction
    #SR.s$R <- s.f*Cmin.s$R # Sulfate reduction
    #MG.s$R <- m.f*Cmin.s$R # methanogenesis
    
    #AR.r$R <- a.f*Cmin.r$R # Aerobic respiration
    #FR.r$R <- f.f*Cmin.r$R # FeIII reduction
    #SR.r$R <- s.f*Cmin.r$R # Sulfate reduction
    #MG.r$R <- m.f*Cmin.r$R # methanogenesis

    Cmin.ftot$R <- AR.f$R + FR.f$R + SR.f$R + MG.f$R
    Cmin.stot$R <- AR.s$R + FR.s$R + SR.s$R + MG.s$R
    Cmin.rtot$R <- AR.r$R + FR.r$R + SR.r$R + MG.r$R

    # Aerobic oxidation reactions
    
    #AAO$R  <- por.grid$mid*k.AAO*O2$C*FSAT(O2$C,C.lim/100,5)*NH4$C*FSAT(NH4$C,C.lim/100,5) # ammonium oxidation
    #SAO$R  <- svf.grid$mid*k.SAO*O2$C*FSAT(O2$C,C.lim/100,5)*X.NH4$C*FSAT(X.NH4$C,C.lim/100,5) # Sorbed ammonium oxidation
    
    FIO$R  <- por.grid$mid*k.FIO*O2$C*FSAT(O2$C,C.lim,5)*Fe$C*FSAT(Fe$C,C.lim/100,5) # iron re-oxidation
    SIO$R  <- svf.grid$mid*k.SIO*O2$C*FSAT(O2$C,C.lim,5)*X.Fe$C*FSAT(X.Fe$C,C.lim,5) # sorbed iron re-oxidation
    CSO$R  <- por.grid$mid*k.CSO*O2$C*FSAT(O2$C,C.lim/100,5)*HS$C*FSAT(HS$C,C.lim/100,5) # Sulfide re-oxidation
    AMO$R  <- por.grid$mid*k.AMO*O2$C*FSAT(O2$C,C.lim/100,5)*CH4$C*FSAT(CH4$C,C.lim/100,5) # aerobic methane oxidation
 
    # Anaerobic oxidation reactions
    
    CSFO.f$R <- svf.grid$mid*k.CSFOf*FeOOH.f$C*FSAT(FeOOH.f$C,C.lim,5)*HS$C*FSAT(HS$C,C.lim/100,5)# Sulfide re-oxidation by FeIII
    CSFO.a$R <- svf.grid$mid*k.CSFOa*FeOOH.a$C*FSAT(FeOOH.a$C,C.lim,5)*HS$C*FSAT(HS$C,C.lim/100,5)# Sulfide re-oxidation by FeIII
    SMO$R  <- por.grid$mid*k.SMO*SO4$C*FSAT(SO4$C,C.lim/100,5)*CH4$C*FSAT(CH4$C,C.lim/100,5) # Sulfate methane oxidation
    
    # Reactions involving sulfide minerals
   
    ISO$R  <- svf.grid$mid*k.ISO*O2$C*FSAT(O2$C,C.lim,5)*FeS$C*FSAT(FeS$C,C.lim,5)*(FeS$C>0)*(O2$C>0) # iron sulfide re-oxidation
    PyO$R  <- svf.grid$mid*k.PyO*O2$C*FSAT(O2$C,C.lim,5)*FeS2$C*FSAT(FeS2$C,C.lim,5)*(FeS2$C>0)*(O2$C>0) # Pyrite re-oxidation
    PyP$R  <- svf.grid$mid*k.PyP*FeS$C*FSAT(FeS$C,C.lim,5)*HS$C*FSAT(HS$C,C.lim,5)*FSAT(SO4$C,C.lim,5)*(FeS$C>0)*(HS$C>0)*(SO4$C>0) # sorbed iron re-oxidation
    ISP$R  <- svf.grid$mid*k.ISP*((Omega.FeS$C-1)^n.ISP)*(Omega.FeS$C>1)*(Fe$C>0)*(HS$C>0)              #iron sulfide precipitation
    ISD$R  <- svf.grid$mid*k.ISD*FeS$C*FSAT(FeS$C,C.lim,5)*(1-Omega.FeS$C)^n.ISD*(Omega.FeS$C<1)*(FeS$C>0) #iron sulfide dissolution

    # Iron oxide ageing

    IOA$R <- svf.grid$mid*k.IOA*FeOOH.f$C*FSAT(FeOOH.f$C,C.lim,5) # Iron oxide ageing
    #FIS$R  <- svf.grid$mid*k.FIS*(K.FIS*Fe$C - X.Fe$C)*(Fe$C>0)*(X.Fe$C>0) # ferric iron sorption

    #=========================================================================
    # Consumption terms: kinetic reactions
    # Solutes [umol cm-3 pore water yr-1] Solids [umol cm-3 solid yr-1] 
    #=========================================================================
    
    OC.f$kin.reac  <- - (AR.f$R + FR.f$R + SR.f$R + MG.f$R)/svf.grid$mid
    OC.s$kin.reac  <- - (AR.s$R + FR.s$R + SR.s$R + MG.s$R)/svf.grid$mid
    OC.r$kin.reac  <- - (AR.r$R + FR.r$R + SR.r$R + MG.r$R)/svf.grid$mid
    FeOOH.f$kin.reac <- (- 6*FR.f$R - 6*FR.s$R - 6*FR.r$R + 4*FIO$R + 4*SIO$R + ISO$R - 8*CSFO.f$R + PyO$R - IOA$R)/svf.grid$mid  
    FeOOH.a$kin.reac <- (- 8*CSFO.a$R + IOA$R)/svf.grid$mid    
    FeS$kin.reac     <- (- ISO$R - PyP$R + 8*CSFO.f$R + 8*CSFO.a$R)/svf.grid$mid
    FeS2$kin.reac    <- (PyP$R - PyO$R)/svf.grid$mid
    X.Fe$kin.reac    <- (- 4*SIO$R)/svf.grid$mid
    
    HCO3$kin.reac<- (AR.f$R + FR.f$R + SR.f$R + AR.s$R + FR.s$R + SR.s$R + AR.r$R + FR.r$R + SR.r$R + 1/4*MG.f$R + 1/4*MG.s$R + 1/4*MG.r$R + AMO$R + SMO$R)/por.grid$mid
    NH4$kin.reac <-  (+ 1/R_CN*AR.f$R + 1/R_CN*FR.f$R + 1/R_CN*SR.f$R + 1/R_CN*MG.f$R + 1/R_CN*AR.s$R + 1/R_CN*FR.s$R + 1/R_CN*SR.s$R + 1/R_CN*MG.s$R + 1/R_CN*AR.r$R + 1/R_CN*FR.r$R + 1/R_CN*SR.r$R + 1/R_CN*MG.r$R)/por.grid$mid
    
    O2$kin.reac  <- (- 1.5*AR.f$R - 1.5*AR.s$R - 1.5*AR.r$R - FIO$R - SIO$R - 2*CSO$R - 9/4*ISO$R - 15/4*PyO$R - 2*AMO$R)/por.grid$mid
    SO4$kin.reac <- (- 0.75*SR.f$R - 0.75*SR.s$R - 0.75*SR.r$R  + CSFO.f$R  + CSFO.a$R + CSO$R + ISO$R - 0.25*PyP$R + 2*PyO$R - SMO$R)/por.grid$mid
    
    Fe$kin.reac  <- (+ 6*FR.f$R + 6*FR.s$R + 6*FR.r$R - 4*FIO$R)/por.grid$mid
    HS$kin.reac  <- (+ 0.75*SR.f$R + 0.75*SR.s$R + 0.75*SR.r$R - 9*CSFO.f$R - 9*CSFO.a$R - CSO$R - 0.75*PyP$R + SMO$R)/por.grid$mid #- PyP$R 
    CH4$kin.reac <- (+ 3/4*MG.f$R + 3/4*MG.s$R + 3/4*MG.r$R - AMO$R - SMO$R)/por.grid$mid
    
    H2O$kin.reac <- (9*FR.f$R + 9*FR.s$R + 9*FR.r$R - 1/4*MG.f$R - 1/4*MG.s$R - 1/4*MG.r$R - 6*FIO$R - 6*SIO$R + 12*CSFO.f$R + 12*CSFO.a$R - 3/2*ISO$R -5/4*PyP$R - 5/2*PyO$R + 2*AMO$R + 2*SMO$R)/por.grid$mid
    H$kin.reac   <- (+ (R_CN - 1)/R_CN*AR.f$R + (R_CN - 1)/R_CN*AR.s$R + (R_CN - 1)/R_CN*AR.r$R - (11*R_CN + 1)/R_CN*FR.f$R - (11*R_CN + 1)/R_CN*FR.s$R - (11*R_CN + 1)/R_CN*FR.r$R + (2*R_CN-4)/(4*R_CN)*SR.f$R + (2*R_CN-4)/(4*R_CN)*SR.s$R + (2*R_CN-4)/(4*R_CN)*SR.r$R - (2*R_CN + 4)/(4*R_CN)*MG.f$R - (2*R_CN + 4)/(4*R_CN)*MG.s$R - (2*R_CN + 4)/(4*R_CN)*MG.r$R + 8*FIO$R + 8*SIO$R + CSO$R - 7*CSFO.f$R - 7*CSFO.a$R + 2*ISO$R + PyP$R + 4*PyO$R - SMO$R)/por.grid$mid
    
    
    #=========================================================================
    # Equilibrium reactions
    # Solutes [umol cm-3 pore water yr-1] Solids [umol cm-3 solid yr-1] 
    #=========================================================================

    FeS$eq.reac  <- (ISP$R - ISD$R)/svf.grid$mid
    Fe$eq.reac   <- (- ISP$R + ISD$R)/por.grid$mid
    HS$eq.reac   <- (- ISP$R + ISD$R)/por.grid$mid
    H$eq.reac    <- (ISP$R - ISD$R)/por.grid$mid

    #=========================================================================
    # Irrigation 
    # Solutes [umol cm-3 pore water yr-1] Solids [umol cm-3 solid yr-1] 
    #=========================================================================

    HCO3$irr <- irr.grid$mid*irr.HCO3*(HCO3.ow - HCO3$C)
    NH4$irr  <- irr.grid$mid*irr.NH4*(NH4.ow - NH4$C)    

    O2$irr   <- irr.grid$mid*irr.O2*(O2.ow - O2$C)
    SO4$irr  <- irr.grid$mid*irr.SO4*(SO4.ow - SO4$C)

    Fe$irr   <- irr.grid$mid*irr.Fe*(Fe.ow - Fe$C)
    HS$irr   <- irr.grid$mid*irr.HS*(HS.ow - HS$C)
    CH4$irr  <- irr.grid$mid*irr.CH4*(CH4.ow - CH4$C)
    
    H2O$irr  <- irr.grid$mid*(H2O.ow - H2O$C)
    
    #=========================================================================
    # Rate of changes composite species 
    #=========================================================================
    
    OC$tran     <- OC.f$tran + OC.s$tran + OC.r$tran
    OC$kin.reac <- OC.f$kin.reac + OC.s$kin.reac + OC.r$kin.reac
    OC$irr      <- OC.f$irr + OC.s$irr + OC.r$irr
    OC$eq.reac  <- OC.f$eq.reac + OC.s$eq.reac + OC.r$eq.reac

    FeOOH$tran     <- FeOOH.f$tran + FeOOH.a$tran
    FeOOH$kin.reac <- FeOOH.f$kin.reac + FeOOH.a$kin.reac
    FeOOH$irr      <- FeOOH.f$irr + FeOOH.a$irr
    FeOOH$eq.reac  <- FeOOH.f$eq.reac + FeOOH.a$eq.reac

    #=========================================================================
    # Rate of changes of all species
    #=========================================================================
    
  # Fe$irr <- alfa*(PL$Fe.ow - Fe$C)

    OC.f$ddt  <- OC.f$tran  + OC.f$kin.reac  + OC.f$eq.reac  + OC.f$irr
    OC.s$ddt  <- OC.s$tran  + OC.s$kin.reac  + OC.s$eq.reac  + OC.s$irr
    OC.r$ddt  <- OC.r$tran  + OC.r$kin.reac  + OC.r$eq.reac  + OC.r$irr
    FeOOH.f$ddt<- FeOOH.f$tran + FeOOH.f$kin.reac + FeOOH.f$eq.reac + FeOOH.f$irr
    FeOOH.a$ddt<- FeOOH.a$tran + FeOOH.a$kin.reac + FeOOH.a$eq.reac + FeOOH.a$irr
    FeS$ddt   <- FeS$tran   + FeS$kin.reac   + FeS$eq.reac   + FeS$irr
    FeS2$ddt  <- FeS2$tran  + FeS2$kin.reac  + FeS2$eq.reac  + FeS2$irr
    
    HCO3$ddt  <- HCO3$tran  + HCO3$kin.reac  + HCO3$eq.reac  + HCO3$irr
    NH4$ddt   <- 1/(por.grid$mid+svf.grid$mid*K.AMS)*(por.grid$mid*(NH4$tran   + NH4$kin.reac   + NH4$eq.reac   + NH4$irr) + svf.grid$mid*(X.NH4$tran   + X.NH4$kin.reac   + X.NH4$eq.reac   + X.NH4$irr))
 
    O2$ddt    <- O2$tran    + O2$kin.reac    + O2$eq.reac    + O2$irr
    SO4$ddt   <- SO4$tran   + SO4$kin.reac   + SO4$eq.reac   + SO4$irr

    Fe$ddt    <- 1/(por.grid$mid+svf.grid$mid*K.FIS)*(por.grid$mid*(Fe$tran + Fe$kin.reac + Fe$eq.reac + Fe$irr) + svf.grid$mid*(X.Fe$tran + X.Fe$kin.reac + X.Fe$eq.reac + X.Fe$irr))
    HS$ddt    <- HS$tran    + HS$kin.reac    + HS$eq.reac    + HS$irr
    CH4$ddt   <- CH4$tran    + CH4$kin.reac    + CH4$eq.reac    + CH4$irr
    
    H2O$ddt   <- H2O$tran + H2O$kin.reac + H2O$eq.reac + H2O$irr
    H$ddt     <- H$tran + H$kin.reac + H$eq.reac + H$irr
    
    OC$ddt <- OC$tran + OC$kin.reac + OC$eq.reac + OC$irr
    FeOOH$ddt <- FeOOH$tran + FeOOH$kin.reac + FeOOH$eq.reac + FeOOH$irr

    #=========================================================================
    }) # end within call
  #=========================================================================
  
    #-----------------------------------------------------------------------------
    # Assemble matrix with rate of changes (RC) 
    #-----------------------------------------------------------------------------
    
    RC <- 0*SV
    for (i in 1:N.var) RC[,i] <- OL[[var.names[i]]]$ddt
    
    #-----------------------------------------------------------------------------
    # Monitoring of simulation progession during dynamic simulation 
    #-----------------------------------------------------------------------------
    
    if ((simulation.type %in% c("dynamic.steady.state","time.dependent")) && (t > sim.info$tc)){
      
      # sim.info$fc = keeps track of number of function calls
      # sim.info$tc = keeps track of incremental simulation time
      
      sim.info$fc <<- sim.info$fc + 1  
      sim.info$tc <<- t*1
      
      # Calculate for each state variable 
      #  ACC = absolute concentration change 
      #  RCC = relative concentration change 
      
      SV.ref <- rep(1E-03,N.var)
      ACC.t <- ACC(RC)
      RCC.t <- RCC(RC,SV,SV.ref)
      
      # Output of simulation progress to the screen 
      screenprint.simulation.progress(sim.info$fc,t,SV,ACC.t,RCC.t,var.names)
      
      # Store simulation progress results 
      sim.info$progress <<- rbind(sim.info$progress,c(t,ACC.t,RCC.t))
      points(log10(t),log10(sqrt(sum(RCC.t^2)/length(RCC.t))),pch=16,col="red")
      
    }
  
  #===========================================================================
  # return statement 
  #===========================================================================
  
  if (!summary.call) {return(list(RC))}
  else {
    
    # Initialisation of lists
    
    SL <- OL[names(SL)]
    RL <- OL[names(RL)]
    EL <- EL
    
    # Updating species information
    
    SL <- updating.chemical.species(SL,PL) 
    SL <- updating.composite.species(spec.table,SL) 
    
    # Updating reaction information
    
    for (i in 1:length(RL)) RL[[i]]$R.int <- sum(RL[[i]]$R*grid$dx) 
    
    # Updating element information
    
    EL <- updating.elements(elt.table,EL,SL) 
    
    return(list(SL=SL,RL=RL,EL=EL,PL=parameters))
    
   }
   
   #-----------------------------------------------------------------------------
 }) # end with call
}
#-----------------------------------------------------------------------------

################################################################################
# Initialisation of state variable matrix
################################################################################

initialise.state.variables <- function(PL)
{
  
  PL$var.names <- c("OC.f","OC.s","OC.r","FeOOH.f","FeOOH.a","FeS","FeS2","HCO3","NH4","O2","SO4","Fe","HS","CH4","H2O","H")
  PL$non.var.names <- NULL
  
  PL$N <- PL$grid$N
  
  PL$N.var <- length(PL$var.names)
  PL$N.non.var <- length(PL$non.var.names)
  
  return(PL)
}

initialise.species.arrays <- function(PL,SL)
{
  
  for (i in 1:length(SL))
  {
    SL[[i]]$C <- rep(0,PL$N)
    SL[[i]]$tran <- rep(0,PL$N)
    SL[[i]]$kin.reac <- rep(0,PL$N)
    SL[[i]]$eq.reac <- rep(0,PL$N)
    SL[[i]]$irr <- rep(0,PL$N)
  }
  return(SL)
}

# Function that creates a rough guess for the initial state

crude.initial.state <- function(PL)
{
  with(PL,{
    
    # Initialisation State Variable matrix (SV) 
    
    SV.guess <- matrix(nrow=N,ncol=N.var,dimnames=list(1:N,var.names))
    
    if ("OC.f" %in% PL$var.names) SV.guess[,"OC.f"] <- rep(0,length.out=N)
    if ("OC.s" %in% PL$var.names) SV.guess[,"OC.s"] <- rep(0,length.out=N)
    if ("OC.r" %in% PL$var.names) SV.guess[,"OC.r"] <- rep(0,length.out=N)
    
    if ("FeOOH.f" %in% PL$var.names)  SV.guess[,"FeOOH.f"] <- rep(0,length.out=N)
    if ("FeOOH.a" %in% PL$var.names)  SV.guess[,"FeOOH.a"] <- rep(0,length.out=N)
    if ("FeS" %in% PL$var.names)  SV.guess[,"FeS"] <- rep(0,length.out=N)
    if ("FeS2" %in% PL$var.names)  SV.guess[,"FeS2"] <- rep(0,length.out=N)
    
    if ("HCO3" %in% PL$var.names) SV.guess[,"HCO3"] <- rep(HCO3.ow,length.out=N)
    if ("NH4" %in% PL$var.names) SV.guess[,"NH4"] <- rep(NH4.ow,length.out=N)
    
    if ("O2" %in% PL$var.names) SV.guess[,"O2"] <- rep(O2.ow,length.out=N)
    if ("SO4" %in% PL$var.names) SV.guess[,"SO4"] <- rep(SO4.ow,length.out=N)
    
    if ("Fe" %in% PL$var.names) SV.guess[,"Fe"] <- rep(Fe.ow,length.out=N)
    if ("HS" %in% PL$var.names) SV.guess[,"HS"] <- rep(HS.ow,length.out=N)
    if ("CH4" %in% PL$var.names) SV.guess[,"CH4"] <- rep(CH4.ow,length.out=N)
    
    if ("H2O" %in% PL$var.names) SV.guess[,"H2O"] <- rep(H2O.ow,length.out=N)
    if ("H" %in% PL$var.names) SV.guess[,"H"] <- rep(H.ow,length.out=N)
    
    return(SV.guess)
    
  })
}

################################################################################
# Test trial model runs
################################################################################

model <- CSFe.model
sim.info <- list(tc= 0,fc=0)

# Initialisation of object list (OL)  

PL <- initialise.parameters(PL=list())
PL$simulation.type <- "direct.steady.state"
SL <- initialise.species.arrays(PL,SL) # default action
OL <- c(SL,RL) # default action

# Test trial 

PL <- initialise.state.variables(PL)
SV <- crude.initial.state(PL)

system.time({
  result <- model(t=0,state=SV,parameters=PL,summary.call=FALSE)
  result <- model(t=0,state=SV,parameters=PL,summary.call=TRUE)
})

remove(result)
remove(SV)

