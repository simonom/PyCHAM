'''module to allow complete evaporation of the smallest size bin when conditions set in Vchange_check.py met'''
# rather than reduce the time step to unpractical length, allow the smallest size bin
# to completely evaporate using this module

def compl_evap(ytest, n0, Vnew, Vol0, nc, sbi):

	# inputs: ----------------------------------------
	# ytest - concentration (molecules/cc (air)) of components
	# n0 - number concentration of particles per size bin (# particles/cc (air))
	# Vnew - new estimates of single particle volume (um3)
	# Vol0 - default volume per size bin (um3), volume at centre of size bin bounds
	# nc - number of components
	# sbi - index of size bins in question
	# ------------------------------------------------

	# allow complete evaporation of the smallest size bin
	# components move to gas-phase (molecules/cc (air))
	ytest[0:nc] += ytest[(sbi+1)*nc:(sbi+2)*nc]
	# effectively zero concentration in particle-phase (molecules/cc (air))
	ytest[(sbi+1)*nc:(sbi+2)*nc] = 1.0e-40
	# remove particle concentration from size bin (# particles/cc (air))
	n0[sbi] = 1.0e-40
	# with no particles in smalles size bin now, default to volume at mid-point between 
	# size bin bounds (um3)
	Vnew[0] = Vol0[0]
	
	return(ytest, n0, Vnew)
