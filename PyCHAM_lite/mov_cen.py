##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2022 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
#                                                                                       											 #
#    All Rights Reserved.                                                                									 #
#    This file is part of PyCHAM                                                         									 #
#                                                                                        											 #
#    PyCHAM is free software: you can redistribute it and/or modify it under              						 #
#    the terms of the GNU General Public License as published by the Free Software       					 #
#    Foundation, either version 3 of the License, or (at your option) any later          						 #
#    version.                                                                            										 #
#                                                                                        											 #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT                						 #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       			 #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              				 #
#    details.                                                                            										 #
#                                                                                        											 #
#    You should have received a copy of the GNU General Public License along with        					 #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                 							 #
#                                                                                        											 #
##########################################################################################
'''module to track particle number size distribution using moving centre size structure (p. 416 of Jacobson 2000)'''

import numpy as np
import scipy.constants as si
import matplotlib.pyplot as plt
from v_check import Vchange_check as Vchange_check

def mov_cen_main(n0, s0, sbn, nc, MW, x, Vol0, t, C0, MV, 
			ic_red, res, solv_time, self):


	# input:---------------------------------------------------------
	
	# n0 - particle number concentration per size bin before integration
	# (# particle/cm3 (air))
	# s0 - volume bounds per size bin (um3) (1st dim.)
	# before a time step over which gas phase reaction and 
	# gas-particle partitioning have occurred (# molecules/cm3 (air))
	# sbn - number of size bins (including any walls)
	# nc - number of components
	# MW - molar weight of components (g/mol)
	# x - original particle size bin radii (um)
	# Vol0 - original volume of size bins (um3) (excluding wall, i.e. volume at mid-point
	# 		between size bin bounds)
	# t - total integration time (s)
	# tinc_count - count on number of time steps since time interval last required 
	# decreasing
	# C0 - original concentrations (# molecules/cm3 (air))
	# MV - molar volume of all components (cm3/mol)
	# self.Psat - saturation vapour pressures (# molecules/cm3 (air))
	# ic_red - flag for time step reduction due to changing initial conditions
	# res - resulting concentrations from ode solver (# molecules/cm3 (air)) 
	# solv_time - times at which integration solved (s)
	# self.wall_on - flag for whether wall turned on
	# self - reference to PyCHAM
	# ---------------------------------------------------------------
	
	NA = si.Avogadro # Avogadro's number (# molecules/mol)

	sbn -= self.wall_on # exclude wall from particle size bin number
	
	# get new volumes of single particles per size bin and
	# check whether volume change is acceptable
	(redt, t, ic_red, Vnew, tsi) = Vchange_check(res, s0, sbn, NA, 
					n0, nc, solv_time, t, ic_red, Vol0, MV, self)
	
	if (redt == 1): # repeat integration with new smaller time step
		return(n0, Vol0, C0, x, redt, t, ic_red)
	
	y = np.zeros((nc*(sbn+1+self.wall_on))) # empty array for holding new concentrations
	# gas concentrations (# molecules/cm3 (air))
	y[0:nc] = res[0:nc]
	if (self.wall_on > 0): # wall concentrations (# molecules/cm3 (air))
		y[-nc*self.wall_on::] = res[-nc*self.wall_on::]
		
	# empty array for holding new particle number concentrations
	N_perbin = np.zeros((sbn, 1))
	
	# if volume condition met, then redistribute particles and components based on
	# the new volume
	for sbi in range(sbn): # particle size bin loop
		
		sbi_new = sum(Vnew[sbi]>s0[1::]) # index of size bin these particles fit now
		
		# add number concentration (# particles/cm3 (air))
		N_perbin[sbi_new, 0] += n0[sbi]
		# add components (# molecules/cm3 (air))
		y[((sbi_new+1)*nc):((sbi_new+2)*nc)] += res[((sbi+1)*nc):((sbi+2)*nc)]
	
	# reshape particle-phase concentrations into components in rows and size bins in
	# columns
	num_molec_new = (y[nc:nc*(sbn+1)].reshape(sbn, nc))
	
	# need to find new volumes of single particles (um3)
	# total volume of components 
	# ((um3 (all particles)/cm3 (air))/(# particle number/cm3 (air))) 
	# calculation is:
	# divide number of # molecules/cm3 (air) by Na to get moles/cm3 (air), then 
	# multiply by um3/mol (MV*1.e12) to get ug3 (of each component)/cm3 (air),
	# then sum volume of components per size bin to get ug3 (all particles)/cm3 (air)
	ish = N_perbin[:, 0]> 0.
	Vsing = np.zeros((sbn))
	
	for sbi in range(sbn): # size bin loop
		if (N_perbin[sbi, 0] > 0.): # single particle volume (um3) if particles present
			Vsing[sbi] = np.sum(((y[nc*(sbi+1):nc*(sbi+2)]/(NA*N_perbin[sbi, :]))*(MV[:, 0]*1.e12)), 0)
	
	Vsing[N_perbin[:, 0]<1.e-20] = Vol0[N_perbin[:, 0]<1.e-20] # assume volume at size bin centre
	
	rad = ((3.*Vsing)/(4.*np.pi))**(1./3.) # new radius per size bin (um)
		   
	return(N_perbin, Vsing, y, rad, redt, t, ic_red)
