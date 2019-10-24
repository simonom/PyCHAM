'''module to calculate density of particles (kg/m3)'''

import numpy as np

def part_prop(y, num_speci, num_sb, NA, y_mw, y_dens, n0):

	# inputs:
	
	# y - concentration of species in particle size bins (molecules/cc (air)), excluding
	# gas and wall
	# ---------------------------------------------------------------------------

	y_asmat = (y.reshape(num_sb, num_speci))
	y_asmat = y_asmat.transpose() # species in rows and size bins in columns
	
	# convert molecules/cc (air) to moles/cc (air) and then multiply by molecular weight  
	# to get g/cc (air) of all species in size bins
	y_mass_array = (y_asmat/NA)*y_mw # mass conc. ind. species
	
	mass_fracs = np.zeros((num_speci, num_sb)) # empty matrix for mass fractions
	ish = n0>0.0 # size bins with particles
	# mass fractions of each species including core
	mass_fracs[:, ish] = y_mass_array[: , ish]/((np.sum(y_mass_array, 0))[ish])
	# total density of particles per size bin (g/cm3)
	tot_rho = np.zeros((num_sb))
	tot_rho[ish] = 1.0E-3/((np.sum(mass_fracs/y_dens, 0))[ish])
	# average molecular weight of particles (g/mol)
	avMW = (mass_fracs*y_mw).sum(0)
	
	return tot_rho, ish, avMW