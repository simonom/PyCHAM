'''module to prepare PyCHAM for partitioning coefficient calculation (particle and wall)'''

import numpy as np
import scipy.constants as si

def kimt_prep(y_mw, TEMP, num_speci, testf, Cw, kgwt):
	
	# ------------------------------------------------------------------
	# inputs:
	# y_mw - molecular weight of components (g/mol) (num_speci,1)
	# TEMP - temperature (K)
	# num_speci - number of components
	# testf - flag for whether in normal mode (0) or testing mode (1)
	# Cw - effective absorbing mass of wall (g/m3 (air))
	# kgwt - mass transfer coefficient for vapor-wall partitioning (m3/g.s)
	# -----------------------------------------------------------------
	
	if testf == 1: # if in testing mode (for test_front.py)
		return(0,0,0,0,0,0,0) # return dummies
	
	surfT = 72.0 # assume surface tension of water (g/s2==mN/m==dyn/cm) for all particles
	# molecular diffusion coeffficient of each species in air (m2/s) 
	# Equation from JL Schnoor, Environmental Modelling: fate and transport of pollutants 
	# in water, air and soil, 1996, ISBN : 0471124362, page 331.
	# Note that no reference given in this text book for this equation, but results for
	# Xe and Kr are within 10 % of those calculated independently in Table 16.1 of 
	# Jacobson (2005)
	# Scale by 1e-4 to convert from cm2/s to m2/s
	DStar_org = (1.9E0*(y_mw**(-2.0E0/3.0E0)))*1e-4
	
	# mean thermal speed of each molecule (m/s) (11.151 Jacobson 2005)
	# note that we need the weight of one molecule, which is why y_mw is divided by
	# Avogadro's constant, and we need it in kg, which is why we multiply by 1e-3
	therm_sp = (np.power((8.0E0*si.k*TEMP)/(np.pi*(y_mw/si.N_A)*1.0E-3), 0.5E0))
	
	# mean free path (m) for each species (16.23 of Jacobson 2005)
	# molecular weight of air (28.966 g/mol taken from table 16.1 Jacobson 2005)
	mfp = (((64.0*DStar_org)/(5*np.pi*therm_sp))*(28.966/(28.966+y_mw))).reshape(-1, 1)

	# accommodation coefficient of particles
	accom_coeff = np.ones((num_speci, 1))*1.0e0

	# convert Cw (effective absorbing mass of wall) from g/m3 (air) to 
	# molecules/cc (air), assuming a molecular weight of 200g/mol
	Cw = ((Cw*1.0e-6)/200.0)*si.N_A

	# convert kgwt from m3/g.s to cm3/molecules.s
	kgwt = (kgwt*1.0e6*200.0)/si.N_A
	return DStar_org, mfp, accom_coeff, therm_sp, surfT, Cw, kgwt