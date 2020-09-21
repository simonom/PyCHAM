'''implement the fully-moving size structure'''
# following Section 13.5.2 of Jacobson (2005) (Fundamentals of Atmospheric Modeling),
# the full-moving size structure revalues the size of particles per size bin
# following gas-particle partitioning

import scipy.constants as si
import numpy as np

def fullmov(num_sb, n0, num_comp, Cp, MV, Vol0): # define module

	# inputs: -----------------------------------------------------
	# num_sb - number of size bins excluding wall (if present)
	# n0 - particle number concentration per size bin
	# num_comp - number of components
	# Cp - particle-phase concentration (molecules/cc (air))
	# MV - molar volume per component (um3/mol)
	# Vol0 - initial volume per size bin at bin centre (um3)
	# -------------------------------------------------------------

	NA = si.Avogadro # Avogadro's number (molecules/mol)
	ish = n0[:, 0]>1.0e-10 # size bins containing particles


	nmolC = np.zeros((num_comp, ish.sum())) # empty array for molar concentration
	# number of moles of each component in a single particle (mol/cc (air))
	nmolC[:,:] = ((Cp[:, ish]/(NA*n0[ish, 0])))

	Vnew = np.zeros((num_sb)) # empty array for new volumes
	# new volume of single particle per size bin (um3), including volume of water
	Vnew[ish] = np.sum(nmolC*MV*1.0e12, axis=0)
	# if no particles in a size bin, assign starting volume (um3)
	Vnew[n0[:, 0]<=1.0e-10] = Vol0[0::][n0[:, 0]<=1.0e-10]
	# new radii per size bin (um)
	x = ((3.0*Vnew)/(4.0*np.pi))**(1.0/3.0)

	return(Vnew, x)
