'''module to track particle number size distribution using moving centre size structure (p. 416 of Jacobson 2000)'''

import numpy as np
import scipy.constants as si
import matplotlib.pyplot as plt
from v_check import Vchange_check as Vchange_check

def mov_cen_main(n0, s0, sbn, nc, MW, x, Vol0, t, tmax, C0, MV, 
				Psat, ic_red, res, solv_time, wall_on):


	# input:---------------------------------------------------------
	
	# n0 - particle number concentration per size bin before integration
	# (# particle/cc (air))
	# s0 - volume bounds per size bin (um3) (1st dim.)
	# before a time step over which gas phase reaction and 
	# gas-particle partitioning have occurred (molecules/cc (air))
	# sbn - number of size bins
	# nc - number of components
	# MW - molar weight of components (g/mol)
	# x - original particle size bin radii (um)
	# Vol0 - original volume of size bins (um3) (excluding wall, i.e. volume at mid-point
	# 		between size bin bounds)
	# t - total integration time (s)
	# tmax - maximum integration time (s)
	# tinc_count - count on number of time steps since time interval last required 
	# decreasing
	# C0 - original concentrations (molecules/cc (air))
	# MV - molar volume (cc/mol)
	# Psat - saturation vapour pressures (molecules/cm3 (air))
	# ic_red - flag for time step reduction due to changing initial conditions
	# res - results from every adaptive time step of ode solver, with time steps in rows
	# 		and estimated concentrations of components in all phases in columns 
	#		(molecules/cc (air)) 
	# solv_time - times at which integration solved (s)
	# wall_on - flag for whether wall turned on
	# ---------------------------------------------------------------
	
	NA = si.Avogadro # Avogadro's number (molecules/mol)

	sbn -= wall_on # exclude wall from size bin number
	
	# get new volumes of single particles per size bin and
	# check whether volume change is acceptable
	(redt, t, ic_red, Vnew, tsi) = Vchange_check(res, MV, s0, sbn, NA, 
					n0, nc, solv_time, t, ic_red, Vol0, Psat)
	
	if redt == 1: # repeat integration with new smaller time step
		return(n0, Vol0, C0, x, redt, t, ic_red)
	
	y = np.zeros((nc*(sbn+1+wall_on))) # empty array for holding new concentrations
	# gas and wall concentrations (molecules/cc (air))
	y[0:nc] = res[0:nc]
	if (wall_on == 1):
		y[-nc::] = res[-nc::]
		
	# empty array for holding new particle number concentrations
	N_perbin = np.zeros((sbn, 1))
	
	# if volume condition met, then redistribute particles and components based on
	# the new volume
	for sbi in range(sbn): # size bin loop
		
		sbi_new = sum(Vnew[sbi]>s0[1::]) # index of size bin these particles fit now
		
		# add number concentration (# particles/cc (air))
		N_perbin[sbi_new, 0] += n0[sbi]
		# add components (molecules/cc (air))
		y[((sbi_new+1)*nc):((sbi_new+2)*nc)] += res[((sbi+1)*nc):((sbi+2)*nc)]
	
	# reshape particle-phase concentrations into components in rows and size bins in
	# columns
	num_molec_new = (y[nc:nc*(sbn+1)].reshape(sbn, nc))
	# need to find new volumes of single particles (um3)
	# total volume of components 
	# ((um3 (all particles)/cc (air))/(particle number/cc (air))) 
	# calculation is:
	# divide number of molecules/cc (air) by Na to get moles/cc(air), then 
	# multiply by um3/mol (MV*1.e12) to get ug3 (of each component)/cc (air),
	# then sum volume of components per size bin to get ug3 (all particles)/cc (air)
	ish = N_perbin[:, 0]> 0.0
	Vsing = np.zeros((sbn))
	
	for sbi in range(sbn): # size bin loop
		if (N_perbin[sbi, 0] > 0.): # single particle volume (um3) if particles present
			Vsing[sbi] = np.sum(((y[nc*(sbi+1):nc*(sbi+2)]/(NA*N_perbin[sbi, :]))*(MV[:, 0]*1.e12)), 0)
	
	Vsing[N_perbin[:, 0]<1.e-20] = Vol0[N_perbin[:, 0]<1.e-20] # assume volume at size bin centre
	
	rad = ((3.*Vsing)/(4.*np.pi))**(1./3.) # new radius per size bin (um)
		   
	return(N_perbin, Vsing, y, rad, redt, t, ic_red)
