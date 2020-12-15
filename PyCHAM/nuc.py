'''module to implement nucleation in PyCHAM'''

import numpy as np
import scipy.constants as si

def nuc(sumt, new_part_sum1, n0, y, MW, rho, num_comp, Varr, x, new_partr, MV, 
		nucv1, nucv2, nucv3, nuc_comp, siz_str, rbou, Vbou, num_sb):


	# ---------------------------------------------------------------
	# input:
	
	# sumt - time through simulation (s)
	# new_part_sum1 - number concentration of newly nucleated particles already 
	# added (# particles/cc (air))
	# n0 - original number concentration of particles (# particles/cc (air))
	# y - molecular concentration of components in the gas and particle phase
	# (molecules/cc (air))
	# MW - molecular weight (g/mol)
	# rho - density (g/cc)
	# num_comp - number of components
	# Varr - particles volumes per size bin(s) (um3)	
	# x - particles radius per size bin(s) (um)
	# new_partr - radius of newly nucleated particle (cm)
	# MV - molar volume (cc/mol)
	# nucv1/v2/v3 - parameter values for nucleation equation
	# nuc_comp - index of the nucleating component
	# siz_str - the size structure
	# rbou - radius bounds (um)
	# Vbou - volume bounds (um3)
	# num_sb - number of size bins, excluding wall
	# ---------------------------------------------------------------
	
	# total number of particles that should have been produced through nucleation 
	# by now given by integration:

	# Gompertz function (for asymmetric sigmoid curve): https://en.wikipedia.org/wiki/Gompertz_function
	nsum1 = nucv1*np.exp(nucv2*(np.exp(-sumt/nucv3)))

	# new number of particles this time step
	new_part1 = (nsum1)-new_part_sum1
	
	# note, setting condition to less than zero may cause very small concentrations of particles 
	# to appear very far from the peak nucleation time
	if (new_part1 < 1.):
		new_part1 = 0.
		
	else: # if new particles are being added
		if (siz_str == 1): # if full-moving being used
			# if there is geometric space for newly formed particles to form a new size bin
			if (rbou[1]/2. > new_partr*1.e4): # scale new_partr to convert from cm to um
				# form new size bin for newly formed particles
				n0[-2] += n0[-1]
				n0 = np.concatenate((np.zeros((1, 1)), n0[0:-1]), axis = 0)
				y[-num_comp*3:-num_comp*2] += y[-num_comp*2:-num_comp*1]
				# remove final size bin
				y = np.concatenate((y[0:-num_comp*2], y[-num_comp::]))
				rbou = np.concatenate((rbou[0:-2], rbou[-1].reshape(1)))
				# add new size bin
				y = np.concatenate((y[0:num_comp], np.zeros((num_comp)), y[num_comp:]))
				rbou = np.concatenate((rbou[0].reshape(1), (rbou[1]/2.).reshape(1), rbou[1::]))
				Vbou = (4./3.)*np.pi*rbou**3.
				
	# identify size bin newly formed particles should enter
	sbi = sum(new_partr*1.e4 > rbou[1::]) # scale new_partr to convert from cm to um
	
	n0[sbi] += new_part1 # add to this size bin
	new_part_sum1 += new_part1
	
	# volume concentration of new particles (cc/cc (air))
	new_vol1 = new_part1*((4./3.)*np.pi*(new_partr)**3.)
	
	# molecular volume of nucleating component
	Vpermolec = (MV[nuc_comp])/si.N_A # molecular volume (cc/molecule)
	# concentration of nucleating component this represents (molecules/cc (air))
	nuc_conc1 = new_vol1/Vpermolec
	
	# remove from gas-phase (molecules/cc (air)), note commented out as violation of
	# mass conservation assumed negligible
# 	y[nuc_comp] -= nuc_conc1
		
	# addition to particle-phase (molecules/cc (air))
	y[num_comp*(1+sbi)+nuc_comp] += nuc_conc1
	
	# average volume of single particles now (scale MV by 1.e12 to convert 
	# from cm3 to um3), note, consider all size bins since if full-moving used,
	# the smallest and largest size bin will be affected
	ypsb = y[num_comp:num_comp*(num_sb+1)].reshape(num_sb, num_comp, order = 'C')
	
	ish = (n0[:, 0] > 0) # index of size bins containing particles
	if (sum(ish) > 0): # if any bins contain particles
		Varr[ish] = ((ypsb[ish, :]/(si.N_A*n0[ish, :]))*(MV[:].reshape(1, -1)*1.e12)).sum(axis = 1)
		# average radius of particles (um)
		x = ((3.*Varr)/(4.*np.pi))**(1./3.)
		x = x.reshape(-1) # ensure x remains as 1D array
	
	return(n0, y, x, Varr, new_part_sum1, rbou, Vbou)
