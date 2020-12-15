'''implement the fully-moving size structure'''
# following Section 13.5.2 of Jacobson (2005) (Fundamentals of Atmospheric Modeling),
# the full-moving size structure revalues the size of particles per size bin
# following gas-particle partitioning

import scipy.constants as si
import numpy as np

def fullmov(num_sb, n0, num_comp, Cp, MV, Vol0, Vbou, rbou): # define module

	# inputs: -----------------------------------------------------
	# num_sb - number of size bins excluding wall (if present)
	# n0 - particle number concentration per size bin
	# num_comp - number of components
	# Cp - particle-phase component concentrations (molecules/cc (air))
	# MV - molar volume per component (um3/mol)
	# Vol0 - initial volume per size bin at bin centre (um3)
	# Vbou - volume bounds between size bins (um3)
	# rbou - radius bounds (um)
	# -------------------------------------------------------------

	# particle-phase concentrations	(molecules/cc (air))
	Cp = np.transpose(Cp.reshape(num_sb, num_comp))
	NA = si.Avogadro # Avogadro's number (molecules/mol)
	ish = n0[:, 0]>0. # size bins containing particles


	nmolC = np.zeros((num_comp, ish.sum())) # empty array for molar concentration
	# number of moles of each component in a single particle (mol/cc (air))
	nmolC[:, :] = ((Cp[:, ish]/(NA*n0[ish, 0])))

	Vnew = np.zeros((num_sb)) # empty array for new volumes
	MVrep = np.repeat(MV, num_sb, axis=1)
	# new volume of single particle per size bin (um3), including volume of water
	Vnew[ish] = np.sum(nmolC*MVrep[:, ish], axis=0)
	
	# if no particles in a size bin, assign starting volume (um3)
	Vnew[n0[:, 0]<=1.0e-10] = Vol0[0::][n0[:, 0]<=1.0e-10]

	# flatten particle-phase concentrations (molecules/cc (air))
	Cp = np.ravel(np.transpose(Cp))

	Vnew_ord = sorted(Vnew) # arrange size bins in ascending size order
	Cpn = np.zeros((len(Cp)))
	nn = np.zeros((n0.shape[0], n0.shape[1]))
	# arrange particle-phase concentrations and particle number concentrations similarly
	for i in range(len(Vnew)):
		ind = Vnew_ord.index(Vnew[i]) # new index
		Cpn[(ind)*num_comp:(ind+1)*num_comp] = Cp[(i)*num_comp:(i+1)*num_comp]
		nn[ind] = n0[i]
	
	Vnew_ord = np.array((Vnew_ord)) # transform list to numpy array

	# new volume bounds (um3) between size bin	
	Vbou[1:-1] = Vnew_ord[0:-1]+(Vnew_ord[1::]-Vnew_ord[0:-1])/2.

	# new radii per size bin (um)
	x = ((3.*Vnew_ord)/(4.*np.pi))**(1./3.)
		
	# new radius bounds (um)
	rbou = ((3.*Vbou)/(4.*np.pi))**(1./3.)


	return(Vnew_ord, x, Cpn, nn, Vbou, rbou)
