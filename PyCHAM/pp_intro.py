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
'''module to set up particle phase'''
# using the user-defined or default values, the initial number size distribution is determined here

import numpy as np
import scipy.constants as si
import part_nsd # calculating number size distributions

def pp_intro(y, num_comp, Pybel_objects, TEMP, H2Oi,
		mfp, accom_coeff, y_mw, surfT, 
		siz_str, num_asb, lowersize, uppersize, pmode, pconc, 
		pconct, nuc_comp, testf, std, mean_rad, therm_sp,
		y_dens, core_diss, space_mode, seedx, 
		act_coeff, partit_cutoff, Press,
		pcont, seed_mw, R_gas, Vwat_inc, seed_eq_wat, self):
	
	# inputs -----------------------------------
	# TEMP - temperature (K) in chamber at start of experiment
	# y_mw - molecular weight (g/mol) of components (num_comp, 1)
	# surfT - surface tension of particles (g/s2 == mN/m == dyn/cm)
	# siz_str - the size structure
	# num_asb - number of size bins (excluding wall)
	# lowersize - lowest size bin radius bound (um)
	# uppersize - largest size bin radius bound (um)
	# pmode - whether particle number concentrations given as modes or explicitly
	# pconc - starting particle concentration (# particle/cc (air)) - if scalar then
	# gets split between size bins in Size_distributions call, or if an array, elements 
	# are allocated to corresponding size bins
	# pconct - time(s) through experiment that particle number concentrations 
	#	correspond to (s)
	# nuc_comp - name of the nucleating component
	# testf - test flag to say whether in normal mode (0) or test mode for front.py (1)
	#       or test mode for pp_intro.py
	# std - geometric standard deviation of the particle number concentration 
	# 		(dimensionless)
	# mean_rad - either the mean radius (um) of particles in lognormal number-size 
	#			distribution (in which case pconc should be scalar), or mean radius of
	#			particles where just one size bin present (in which case pconc is also
	#			scalar)
	# y_dens - liquid density of components (kg/m3) (num_comp, 1)
	# self.Psat - saturation vapour pressure of components (# molecules/cm3 (air))
	# core_diss - core dissociation constant
	# space_mode - string specifying whether to space size bins logarithmically or 
	# linearly
	# seedx - mole ratio of non-water components comprising seed particles
	# self.comp_namelist - names of components noted in chemical scheme file
	# act_coeff - activity coefficient of components
	# self.wall_on - whether or not to consider wall
	# partit_cutoff - product of vapour pressure and activity coefficient
	#		at which gas-particle partitioning assumed zero (Pa)
	# Press - pressure inside chamber
	# self.seedi - index of seed components
	# pcont - whether particle injections are instantaneous 
	#	or continuous (flag)
	# seed_mw - molecular weight of seed components (g/mol)
	# R_gas - the universal gas constant (cc.Pa/K.mol == kg.m2.s-2.K-1.mol-1)
	# Vwat_inc - flag for whether (1) or not (0) the number size 
	# distribution of seed particles includes the volume of water
	# seed_eq_wat - flag for whether (1) or not (0) to allow water
	# equilibrium with seed particles prior to experiment start
	# self - reference to program
	# ------------------------------------------
	
	if (testf == 1): # in test mode
		return(0,0,0,0,0,0,0,0,0,0,0,0) # return dummies
	
	# isolate the starting number size distribution information
	i = (np.where(pconct[0, :] == 0))[0] # index of initial information
	
	if (i.size == 0): # if no initial information provide fillers
		pconcn = np.zeros((1))
		stdn = [1.e20]
		mean_radn = [-1.e6]
		
	else: # obtain initial information
		
		pconcn = pconc[:, i]

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
			if (mean_rad[:, i[0]].ndim > 1):
				mean_radn = mean_rad[:, i[0]]
			if (mean_rad[:, i[0]].ndim == 1):
				mean_radn = mean_rad[:, i[0]]
			if (mean_rad[:, i[0]].ndim == 0):
				mean_radn = [mean_rad[:, i[0]]]
		except:
			mean_radn = [mean_rad[:, i[0]]]
	
	# index of nucleating component
	if len(nuc_comp)>0:
		nuc_compi = self.comp_namelist.index(nuc_comp[0])
		nuc_comp = np.empty(1, dtype=int)
		nuc_comp[0] = nuc_compi

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

	else: # get the particle number size distribution from inputs
	
		[N_perbin, x, rbou, Vbou, Varr, 
		upper_bin_rad_amp] = part_nsd.part_nsd(lowersize, 
		num_asb, uppersize, mean_radn, stdn, pmode, pconcn, space_mode, testf, self)
		
		# if injection of particles at start of experiment is continuous, then even for the 
		# start of the experiment, this will be dealt with in cham_up
		if (i.size > 0 and num_asb > 0 and pcont[0, 0] == 1):
			N_perbin[:, :] = 0
		
	# set first volume and radius bound to zero, thereby allowing shrinkage to zero in the 
	# smallest bin
	# remember initial first radius bound for saving
	rbou00 = rbou[0]
	Vbou[0] = 0.
	rbou[0] = 0. # this reversed in saving.py back to rbou00
	
	if (num_asb > 0):
		# remember the radii (um) and volumes (um3) at size bin centre now
		rad0 = np.zeros((len(x)))
		rad0[:] = x[:]
		Vol0 = np.zeros((len(Varr)))
		Vol0[:] = Varr[:]
		# empty array for concentration of components on wall due to 
		# particle deposition to wall (# molecules/cm3)
		self.C_p2w = np.zeros((num_asb*num_comp))
	
	if (self.wall_on > 0):
		num_sb = num_asb+self.wall_on  # account for wall
	else:
		num_sb = num_asb

	# remember wall concentrations (set in init_conc)
	y_w = y[num_comp::]

	# include particle-phase concentrations of components to y (# molecules/cm3 (air))
	y = np.append(y[0:num_comp], np.zeros(((num_sb-self.wall_on)*num_comp)))

	# include wall concentrations (# molecules/cm3)
	y = np.append(y, y_w)	

	# molar volume of components (multiply y_dens by 1e-3 to convert from kg/m3 to g/cm3 and give
	# MV in units cm3/mol)
	MV = (y_mw/(y_dens*1.e-3)).reshape(num_comp, 1)

	# number of size bins, excluding wall
	num_aasb = num_sb-self.wall_on
	# account for particle-phase concentration of components contained in seed particles
	if (sum(pconcn) > 0.):
		
		# check whether water to be equilibrated with seed particles prior to experiment start
		if (seed_eq_wat == 1 or Vwat_inc == 1): # if yes, water is to be equilibrated
		
			# check whether the stated initial number size distribution included the 
			# volume of water
			if (Vwat_inc == 1): # if number size distribution does include volume of water
				
				avMW0 = np.ones((num_aasb)) # first guess of average molecular weight

				# average molecular weight of seed in each size bin (g/mol)
				avMW1 = np.ones((num_aasb))*(np.sum(y_mw[self.seedi])/len(self.seedi))
				
				avMW = avMW1 # for calculating Kelvin factor
				lcnt = 1 # loop count
				
				while (np.max((avMW1-avMW0)/avMW1) > 0.1):

					# calculate Kelvin effect factor for the provided number size distribution
					# kelvin factor for each size bin (excluding wall), eq. 16.33 Jacobson et al. (2005)
					# note that seed_mw has units g/mol, surfT (g/s2==mN/m==dyn/cm), R_gas is multiplied by 
					# 1e7 for units g cm2/s2.mol.K, 
					# TEMP is K, x (radius) is multiplied by 1e-4 to give cm from um and 
					# y_dens multiplied by  by 1e-3 to convert from kg/m3 to g/cc
					kelv = np.exp((2.e0*avMW*surfT)/(R_gas*1.e7*TEMP*x*1.e-4*(sum(y_dens[self.seedi[:], 0])/len(self.seedi)*1.e-3)))
				
					# equilibrium mole fraction of water per size bin
					# from the ode solver equation for vapour-particle partitioning of water
					xwat = y[H2Oi]/(self.Psat[0:num_aasb, H2Oi]*kelv*act_coeff[0:num_aasb, H2Oi])
					
					# allow for mole fraction of water in mole fraction of non-water seed components
					# for all size bins
					seedx = seedx*(1./sum(seedx)) # ensure the non-water mole fractions sum to one
					seedxn = seedx*(1.-xwat)
					
					# average molar volume of seed components (cm3/mol) for all size bins
					avMV = (sum(seedxn*MV[self.seedi[:], 0])+xwat*MV[H2Oi, 0])
					
					# total molecular concentration of seed components including water (# molecules/cm3) per size bin,
					# note that volume multiplied by 1e-12 to convert from um3 to cm3
					tmc = (((Varr*1.e-12)*N_perbin[:, 0])/avMV)*NA
					
					# concentration of particle-phase seed components in all size bin
					for ci in range(len(self.seedi)): # loop through indices of seed components
						
						# non-water (# molecules/cm3)
						y[num_comp+self.seedi[ci]:(num_comp*(num_aasb)+self.seedi[ci])+1:num_comp] = tmc*(seedxn[ci, :])
						
						if (ci == len(self.seedi)-1): # reached water component
							# water (# molecules/cm3)
							y[num_comp+H2Oi:(num_comp*(num_aasb)+H2Oi)+1:num_comp] = tmc*xwat

		
				
					# for average molecular weight, first get fraction of each component in each size bin
					avMW = y[num_comp:num_comp*(num_aasb+1)].reshape(num_comp, num_aasb, order='F')
					avMW = avMW/np.sum(avMW, axis=0)
					# average molecular weight of seed in each size bin (g/mol)
					if (lcnt % 2 != 0): # if on even count
						avMW0 = np.sum(avMW*y_mw.reshape(-1, 1), axis = 0)
						avMW = avMW0 # for calculating Kelvin factor
						
					else: # if on odd count
						avMW1 = np.sum(avMW*y_mw.reshape(-1, 1), axis = 0)
						avMW = avMW1 # for calculating Kelvin factor
					
					lcnt += 1 # loop count
				
			if (Vwat_inc == 0): # if number size distribution does not include volume of water
			
				# calculate Kelvin effect factor for the provided number size distribution
				# kelvin factor for each size bin (excluding wall), eq. 16.33 Jacobson et al. (2005)
				# note that seed_mw has units g/mol, surfT (g/s2==mN/m==dyn/cm), R_gas is multiplied by 
				# 1e7 for units g cm2/s2.mol.K, 
				# TEMP is K, x (radius) is multiplied by 1e-4 to give cm from um and 
				# y_dens multiplied by  by 1e-3 to convert from kg/m3 to g/cc
				kelv = np.exp((2.e0*sum(seed_mw)/len(seed_mw)*surfT)/(R_gas*1.e7*TEMP*x*1.e-4*(sum(y_dens[self.seedi[:], 0])/len(self.seedi)*1.e-3)))
				
				# equilibrium mole fraction of water from the ode solver 
				# equation for vapour-particle partitioning
				xwat = y[H2Oi]/(self.Psat[:, H2Oi]*kelv*act_coeff[:, H2Oi])
				
				seedx = seedx*(1./sum(seedx)) # ensure the non-water mole fractions sum to one
				
				# average molar volume of dry (no water) seed components (cc/mol) for all size bins
				avMV = (sum(seedx*MV[self.seedi[:], 0]))
				
				# total molecular concentration of dry (no water) seed components 
				# including water (molecules/cm3) per size bin, 
				# note that volume multiplied by 1e-12 to convert from um3 to cm3
				tmc = (((Varr*1.e-12)*N_perbin[:, 0])/avMV)*NA
				
				# concentration of particle-phase seed components in all size bin
				for ci in range(len(self.seedi)): # loop through indices of seed components
					# non-water (# molecules/cm3)
					y[num_comp+self.seedi[ci]:(num_comp*(num_aasb)+self.seedi[ci])+1:num_comp] = tmc*(seedx[ci, :])
						
					if (ci == len(self.seedi)-1): # reached final non-water seed component
						# water (# molecules/cm3)
						y[num_comp+H2Oi:(num_comp*(num_aasb)+H2Oi)+1:num_comp] = (-xwat*tmc)/(xwat-1.)
					
		# if water not to be equilibrated with seed particles prior to experiment start
		if (seed_eq_wat == 0 and Vwat_inc == 0):			
	
			# initial seed particles will not contain water and any water vapour in chamber
			# will try to equilibrate with particles through the ODE solver		
			
			# mole-fraction weighted molar volume (cm3/mole) (average molar volume of particles)
			mfwMV = sum(MV[self.seedi[:], 0]*(seedx/sum(seedx)))
			# convert to molecular volume (cm3/molecule)
			mfwMV = mfwMV/NA
			# total molecular concentration of seed components per size bin (molecules/cm3)
			# note that volume multiplied by 1e-12 to convert from um3 to cm3
			ytot = ((Varr*1.e-12)*N_perbin[:, 0])/mfwMV
		
			for ci in range(len(self.seedi)): # loop through indices of seed components
		
				# concentration of this component in all size bins (molecules/cc (air)):
				y[num_comp+self.seedi[ci]:(num_comp*(num_aasb)+self.seedi[ci])+1:num_comp] = ytot*(seedx[ci]/sum(seedx))
				
				
	# mass concentration of particles (scale y_dens by 1e-3 to convert from kg/m3
	# to g/cm3)
	if (num_aasb > 0): # with particles
		mass_conc = 0. # start cumulation
		for i in range(num_aasb): # as size bin now account for wall too
			mass_conc += sum((y_dens[:, 0]*1.0e-3)*((y[num_comp*(i+1):num_comp*(i+2)]/si.N_A)*MV[:, 0]))
			mass_conc -= (y_dens[int(H2Oi), 0]*1.e-3)*((y[num_comp*(i+1)+int(H2Oi)]/si.N_A)*MV[int(H2Oi), 0])
		mass_conc = mass_conc*1.e12 # convert from g/cc (air) to ug/m3 (air)

	# start counter on number concentration of newly nucleated particles (# particles/cm3(air))
	np_sum = 0.
	
	return(y, N_perbin, x, Varr, Vbou, rad0, Vol0, rbou, MV, num_sb, nuc_comp, rbou00, 
			upper_bin_rad_amp, np_sum)
