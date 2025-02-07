########################################################################
#								       #
# Copyright (C) 2018-2025					       #
# Simon O'Meara : simon.omeara@manchester.ac.uk			       #
#								       #
# All Rights Reserved.                                                 #
# This file is part of PyCHAM                                          #
#                                                                      #
# PyCHAM is free software: you can redistribute it and/or modify it    #
# under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or    #
# (at  your option) any later version.                                 #
#                                                                      #
# PyCHAM is distributed in the hope that it will be useful, but        #
# WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# General Public License for more details.                             #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with PyCHAM.  If not, see <http://www.gnu.org/licenses/>.      #
#                                                                      #
########################################################################
'''module to set up particle phase'''
# using the user-defined or default values, the initial number size 
# distribution is determined here

import numpy as np
import scipy.constants as si
import part_nsd # calculating number size distributions
from pp_water_equil import pp_water_equil

def pp_intro(y, num_comp, TEMP, H2Oi,
		mfp, accom_coeff, y_mw, surfT, 
		siz_str, num_asb, lowersize, uppersize, 
		testf, std, therm_sp, act_coeff, Press,
		seed_mw, R_gas, self):
	
	# inputs -----------------------------------
	# TEMP - temperature (K) in chamber at start of experiment
	# y_mw - molecular weight (g/mol) of components (num_comp, 1)
	# surfT - surface tension of particles (g/s2 == mN/m == dyn/cm)
	# siz_str - the size structure
	# num_asb - number of size bins (excluding wall)
	# lowersize - lowest size bin radius bound (um)
	# uppersize - largest size bin radius bound (um)
	# self.pmode - whether particle number concentrations given as 
	#	modes or explicitly
	# self.pconc - starting particle concentration (# particle/cm3 (air))
	#	 - if scalar then gets split between size bins in 
	#	Size_distributions call, or if an array, elements 
	# 	are allocated to corresponding size bins
	# self.pconct - time(s) through experiment that particle number 
	#	concentrations correspond to (s)
	# self.nuc_comp - name of the nucleating component
	# testf - test flag to say whether in normal mode (0) or test 
	#	mode for front.py (1) or test mode for pp_intro.py
	# std - geometric standard deviation of the particle number 
	#	concentration (dimensionless)
	# self.mean_rad - either the mean radius (um) of particles in 
	#	lognormal number-size distribution, 
	#	or mean radius of particles where 
	#	just one size bin present (in which case pconc is also
	#	scalar)
	# self.y_dens - liquid density of components (kg/m3) 
	#	(num_comp, 1)
	# self.Psat - saturation vapour pressure of components 
	# 	(# molecules/cm3 (air))
	# self.seedx - mole ratio of non-water components comprising seed 
	# particles
	# self.comp_namelist - names of components noted in chemical 
	# scheme file
	# act_coeff - activity coefficient of components
	# self.wall_on - whether or not to consider wall
	# Press - pressure inside chamber
	# self.seedi - index of seed components
	# self.pcont - whether particle injections are instantaneous 
	#	or continuous (flag)
	# seed_mw - molecular weight of seed components (g/mol)
	# R_gas - the universal gas constant (cm3.Pa/K.mol == 
	#	kg.m2.s-2.K-1.mol-1)
	# self.Vwat_inc - flag for whether (1) or not (0) the number size 
	# 	distribution of seed particles includes the volume of water
	# seed_eq_wat - flag for whether (1) or not (0) to allow water
	# equilibrium with seed particles prior to experiment start
	# self - reference to program
	# ------------------------------------------
	
	if (testf == 1): # in test mode
		return(0,0,0,0,0,0,0,0,0,0,0,0) # return dummies
	
	# index of initial information
	i = (np.where(self.pconct[0, :] == 0))[0]
	
	# if seed particle present at start
	if (sum(self.pconct[0, :] == 0) > 0):
	
		# get the starting seed mole fraction
		seedx_now = np.squeeze(self.seedx[:, :, i]).reshape(
			self.seedx.shape[0], self.seedx.shape[1])

	else: # if no seed particle present at start
		seedx_now = np.zeros((self.seedx.shape[0], self.seedx.shape[1]))

	if (i.size == 0): # if no initial information provide fillers
		# note that this line changed from pconcn = np.zeros((1)) 
		# on 27/04/2023 so that 
		# correct number of size bins given even if particle 
		# number is zero at start of simulation
		pconcn = np.zeros((num_asb))
		stdn = [1.e20]
		mean_radn = [-1.e6]
		
	else: # obtain initial information
		
		pconcn = self.pconc[:, i]

		try: # if std already an array
			if (std[:, i[0]].ndim > 1):
				stdn = std[:, i[0]]
			if (std[:, i[0]].ndim == 1):
				stdn = std[:, i[0]]
			if (std[:, i[0]].ndim == 0):
				stdn = [std[:, i[0]]]
		except:
			stdn = [std[:, i[0]]]
		
		try: # if mean_rad already an array
			if (self.mean_rad[:, i[0]].ndim > 1):
				mean_radn = self.mean_rad[:, i[0]]
			if (self.mean_rad[:, i[0]].ndim == 1):
				mean_radn = self.mean_rad[:, i[0]]
			if (self.mean_rad[:, i[0]].ndim == 0):
				mean_radn = [self.mean_rad[:, i[0]]]
		except:
			mean_radn = [self.mean_rad[:, i[0]]]
	
	# index of nucleating component
	if (len(self.nuc_comp) > 0):
		nuc_compi = self.comp_namelist.index(self.nuc_comp[0])
		self.nuc_comp = np.empty(1, dtype=int)
		self.nuc_comp[0] = nuc_compi

	NA = si.Avogadro # Avogadro's number (# molecules/mol)
	
	# if in test mode, let user know what's happening
	if (testf == 2):
		print('calling size_distr.lognormal via part_nsd')

	if (num_asb == 0): # create dummy variables if no size bins
	
		N_perbin = np.zeros((1, 1))
		x = np.zeros((1, 1))
		Varr = np.zeros((1, 1))
		Vbou = np.zeros((1, 1))
		rad0 = np.zeros((1, 1))
		Vol0 = np.zeros((1, 1))
		rbou = np.zeros((1, 1))
		upper_bin_rad_amp = 1.e6
		# empty array for concentration of components on wall due to 
		# particle deposition to wall (# molecules/cm3)
		self.C_p2w  = 0.

	# if size bins present, get the particle number size distribution from inputs
	else:
		[N_perbin, x, rbou, Vbou, Varr, 
		upper_bin_rad_amp] = part_nsd.part_nsd(lowersize, 
		num_asb, uppersize, mean_radn, stdn, pconcn, testf, self)
		
		# if injection of particles at start of experiment is 
		# continuous, then even for the 
		# start of the experiment, this will be dealt with in cham_up
		if (i.size > 0 and num_asb > 0 and self.pcont[0, 0] == 1):
			N_perbin[:, :] = 0

	# if the concentration of seed particles has a dimension
	# for time, but the mean radius has only one element
	# along the time dimension, then tile the mean radius
	# over all the times represented by number concentration
	if (self.pconc.shape[1]>1 and self.mean_rad.shape[1] == 1):
		self.mean_rad = np.tile(x.reshape(-1, 1), (1, self.pconc.shape[1]))

	# set first volume and radius bound to zero, thereby allowing shrinkage to zero in the 
	# smallest bin
	# remember initial first radius bound for saving
	rbou00 = rbou[0]
	Vbou[0] = 0.
	rbou[0] = 0. # this reversed in saving.py back to rbou00
	
	if (num_asb > 0):
		# remember the radii (um) and volumes (um3) at size bin 
		# centre now
		rad0 = np.zeros((len(x)))
		rad0[:] = x[:]
		Vol0 = np.zeros((len(Varr)))
		Vol0[:] = Varr[:]
		# empty array for concentration of components on wall due
		# to particle deposition to wall (# molecules/cm3)
		self.C_p2w = np.zeros((num_asb*num_comp))
	
	if (self.wall_on > 0):
		num_sb = num_asb+self.wall_on  # account for wall
	else:
		num_sb = num_asb

	# remember wall concentrations (set in init_conc)
	y_w = y[num_comp::]

	# include particle-phase concentrations of components to y 
	# (# molecules/cm3 (air))
	y = np.append(y[0:num_comp], np.zeros(((num_sb-
		self.wall_on)*num_comp)))

	# include wall concentrations (# molecules/cm3)
	y = np.append(y, y_w)	

	# molar volume of components (multiply self.y_dens by 1e-3 to 
	# convert from kg/m3 to g/cm3 and give
	# MV in units cm3/mol)
	MV = (y_mw/(self.y_dens*1.e-3)).reshape(num_comp, 1)

	# number of size bins, excluding wall
	num_aasb = num_sb-self.wall_on

	# account for particle-phase concentration of components 
	# contained in seed particles
	if (sum(pconcn) > 0.):
		
		# volume concentration of new seed particles (um3/cm3 (air)),
		# note that for N_perbin, size bins are in rows (as set in size_distr.py)
		Vperbin = Varr*N_perbin[:, 0]
		
		# concentrations of components in new seed particles 
		# (# molecules/cm3 (air))
		yn = np.zeros((num_comp*(num_aasb)))
		
		yn = pp_water_equil(y[self.H2Oi], yn, seedx_now, num_aasb, y_mw, 
			R_gas, TEMP, surfT, act_coeff, Vperbin, x, num_comp, self)
		
		# include particle-phase concentations in y 
		# (molecules/cm3)
		y[num_comp:(num_comp*(num_aasb+1))] = yn[:]
		
	# if wall present and water to be equilibrated between wall and
	# gas
	if (self.wall_on > 0 and self.Vwat_inc == 2):
		# loop through walls
		for walli in range(self.wall_on):
			y[num_comp*(num_asb+1+walli)+H2Oi] = ((y[H2Oi]*self.Cw[walli, H2Oi])/
			(self.Psat[num_asb+walli, H2Oi]*act_coeff[num_asb+walli, H2Oi])) 
			
	# mass concentration of particles (scale self.y_dens by 1e-3 
	# to convert from kg/m3 to g/cm3)
	if (num_aasb > 0): # with particles
		mass_conc = 0. # start cumulation
		# as size bin now account for wall too
		for i in range(num_aasb): 
			mass_conc += (sum((self.y_dens[:, 0]*1.0e-3)*
			((y[num_comp*(i+1):num_comp*(i+2)]/si.N_A)*
			MV[:, 0])))
			mass_conc -= ((self.y_dens[int(H2Oi), 0]*1.e-3)
			*((y[num_comp*(i+1)+int(H2Oi)]/si.N_A)*
			MV[int(H2Oi), 0]))
		mass_conc = mass_conc*1.e12 # convert from g/cm3 (air) to ug/m3 (air)

	# start counter on number concentration of newly 
	# nucleated particles (# particles/cm3(air))
	np_sum = 0.

	try: # in case called from autorun, in which a finisher simulation is setup
		if (self.param_const['sim_type'] == 'finisher'):
			# particle and wall concentrations (# molecules/cm3)
			y[num_comp::] = self.param_const['ynow'][num_comp::]
			N_perbin[:, 0] = self.param_const['Nnow'][:]
			
	except: # not called from finisher simulation
		N_perbin[:] = N_perbin[:]

	# if starting concentrations of non-gas components supplied
	# from end of previous simulation 
	if hasattr(self, 'y0_other_phase'): 
		y[num_comp::] = self.y0_other_phase
		N_perbin = self.N_perbin0_prev_sim
		x = self.x0_prev_sim
		Varr = self.Varr0_prev_sim

	# prepare matrix to hold dissociation constants of components (columns)
	# with respect to water over all size bins (rows)
	self.diss_wrtw = np.zeros((num_asb, num_comp))
	self.diss_wrtw[:, :] = self.noncore_diss_wrtw
	self.diss_wrtw[:, self.seedi] = self.core_diss_wrtw
	self.diss_wrtw[:, H2Oi] = 1. # water with respect to itself
	
	return(y, N_perbin, x, Varr, Vbou, rad0, Vol0, rbou, MV, num_sb,
		 rbou00, upper_bin_rad_amp, np_sum)
