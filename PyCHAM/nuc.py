'''module to implement nucleation in PyCHAM'''

import numpy as np
import scipy.constants as si

def nuc(sumt, new_part_sum1, n0, y, MW, rho, num_speci, V1stbin, x1stbin, new_partr, MV, 
		nucv1, nucv2, nucv3, nuc_comp):


	# ---------------------------------------------------------------
	# input:
	
	# sumt - time through simulation (s)
	# new_part_sum - number concentration of newly nucleated particles already 
	# added (# particles/cc (air))
	# n0 - original number concentration of particles (# particles/cc (air))
	# y - molecular concentration of components in the gas and particle phase
	# (molecules/cc (air))
	# MW - molecular weight (g/mol)
	# rho - density (g/cc)
	# num_speci - number of species
	# V1stbin - particles volume in smallest size bin(s) (um3)	
	# x1stbin - particles radius in smallest size bin(s) (um)
	# new_partr - radius of two ELVOC molecules together in a newly nucleating 
	# particle (cm)
	# MV - molar volume (cc/mol)
	# nucv1/v2/v3 - parameter values for nucleation equation
	# nuc_comp - index of the nucleating component
	# ---------------------------------------------------------------
	
	# total number of particles that should have been produced through nucleation 
	# by now given by integration:

	# Gompertz function (for asymmetric sigmoid curve): https://en.wikipedia.org/wiki/Gompertz_function
	nsum1 = nucv1*np.exp(nucv2*(np.exp(-sumt/nucv3)))

	# new number of particles this time step
	new_part1 = (nsum1)-new_part_sum1
	if (new_part1 <1. ): # reality
		new_part1 = 0.

	n0[0] += new_part1 # add to smallest size bin
	new_part_sum1 += new_part1
	
	# volume concentration of new particles (cc/cc (air))
	new_vol1 = new_part1*((4.0/3.0)*np.pi*(new_partr)**3.0)
	
	# molecular volume of nucleating component
	Vpermolec = (MV[nuc_comp])/si.N_A # molecular volume (cc/molecule)
	# concentration of nucleating component this represents (molecules/cc (air))
	nuc_conc1 = new_vol1/Vpermolec
	
	# remove from gas-phase (molecules/cc (air)), note commented out as violation of
	# mass conservation assumed negligible
# 	y[nuc_comp] -= nuc_conc1
		
	# addition to particle-phase (molecules/cc (air))
	y[num_speci*1+nuc_comp] += nuc_conc1
	
	# average volume of single particles in first bin now (scale by 1.0e12 to convert 
	# from cm3 to um3)
	if (n0[0]>0.):
		y1stbin = y[num_speci:num_speci*2]
		V1stbin = np.array((((y1stbin/(si.N_A*n0[0]))*MV[:]*1.0e12).sum()))
		# average radius of first bin particles (um)
		x1stbin = ((3.0*V1stbin)/(4.0*np.pi))**(1.0/3.0)
	
	
	return(n0, y, x1stbin, V1stbin, new_part_sum1)
