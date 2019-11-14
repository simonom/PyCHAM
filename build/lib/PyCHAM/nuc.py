'''module to implement nucleation in PyCHAM'''

import numpy as np
import scipy.constants as si
import ipdb

def nuc(sumt, new_part_sum1, n0, y, MW, rho, num_speci, x, new_partr, t, MV, 
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
	# x - particles radius in smallest size bin(s) (um)
	# new_partr - radius of two ELVOC molecules together in a newly nucleating 
	# particle (cm)
	# t - time step (s)
	# MV - molar volume (cc/mol)
	# nucv1/v2/v3 - parameter values for nucleation equation
	# nuc_comp - index of the nucleating component
	# ---------------------------------------------------------------
	# output:
	
	# n0 - end of time step particle number concentration per size bin
	# (# particle/cc (air))
	# new_part_sum - number concentration of new particles already added 
	# (# particles/cc (air))
	# ---------------------------------------------------------------
	
	# total number of particles that should have been produced through nucleation 
	# by now given by integration:

	# Gompertz function (for asymmetric sigmoid curve): https://en.wikipedia.org/wiki/Gompertz_function
	nsum1 = nucv1*np.exp(nucv2*(np.exp(-sumt/nucv3)))

	# new number of particles this time step
	new_part1 = (nsum1)-new_part_sum1
	if new_part1<1.0:
		new_part1 = 0.0
	n0[0] += new_part1
	new_part_sum1 += new_part1
	
	# volume concentration of new particles (cc/cc (air))
	new_vol1 = new_part1*((4.0/3.0)*np.pi*(new_partr)**3.0)
	
	# molecular volume of ELVOC
	Vpermolec = (MV[nuc_comp, 0])/si.N_A # molecular volume (cc/molecule)
	# concentration of ELVOC this represents (molecules/cc (air))
	ELVOC_conc1 = new_vol1/Vpermolec
	 		
	# if insufficient ELVOC present
	if y[nuc_comp]<(ELVOC_conc1):
	
		# reverse above changes
		n0[0] -= new_part1
		
		new_part_sum1 -= new_part1
		
		ELVOC_conc1 = 0.0
		
		print('insufficient ELVOC for nucleation')
	
	
	# remove from gas-phase (molecules/cc (air))
	y[nuc_comp] -= ELVOC_conc1
		
	# add to particle-phase (molecules/cc (air))
	y[num_speci*1+nuc_comp] += ELVOC_conc1
	
	# average volume of single particles in first bin now (scale by 1.0e12 to convert 
	# from cm3 to um3)
	if n0[0]>0.0:
		# only use non-zero elements (zero ones given by 1.0e-40)
		y1stbin = (y[num_speci:num_speci*2][y[num_speci:num_speci*2]!=1.0e-40])
		mv = MV[y[num_speci:num_speci*2]!=1.0e-40]
		V1stbin = ((y1stbin/(si.N_A*n0[0]))*mv[:,0]*1.0e12).sum()
		# average radius of first bin particles (um)
		x1stbin = ((3.0*V1stbin)/(4.0*np.pi))**(1.0/3.0)
	else:
		x1stbin = x
		V1stbin = (4.0/3.0)*np.pi*(x**3.0)
	
	# radii and volumes	
	xret = np.array((x1stbin))
	Vret = np.array((V1stbin))
	
	return(n0, y, xret, Vret, new_part_sum1)