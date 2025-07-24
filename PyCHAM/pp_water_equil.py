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
'''get the concentrations of particle-phase components'''
# provides the molecules concentration of particle-phase components
# (molecules/cm3)

import numpy as np
import scipy.constants as si

def pp_water_equil(H2Ogc, yn, seedx_now, num_asb, y_mm, R_gas, 
		TEMP, surfT, act_coeff, Vperbin, x, num_comp, self):
	
	# inputs: ---------------------------------------------------------
	# H2Ogc - concentration of water in gas-phase (molecules/cm3)
	# yp - container for concentration of particle-phase component
	# concentrations (molecules/cm3)
	# seedx_now - mole fraction of particle-phase components (may or 
	#	may not include water), components in rows, size bins in columns
	# num_asb - number of particle size bins
	# y_mm - molar mass of components (g/mol)
	# R_gas - the universal gas constant (cm3.Pa/K.mol == 
	#	kg.m2.s-2.K-1.mol-1)
	# TEMP - temperature (K) in experiment
	# surfT - surface tension of particles (g/s2 == mN/m == dyn/cm)
	# act_coeff - activity coefficient of components
	# Vperbin - volume concentration of new seed particles 
	#	(um^3/cm^3 (air))
	# x - particle radii (um)
	# num_comp - number of components
	# self - parent PyCHAM object
	# -----------------------------------------------------------------
	
	NA = si.Avogadro # Avogadro's number (# molecules/mol)

	# check whether water to be equilibrated with seed particles 
	# prior to experiment start
	if (self.seed_eq_wat == 1 or self.Vwat_inc > 0): 
		
		# check whether the stated initial number size 
		# distribution included the 
		# volume of water
		# if number size distribution DOES include volume of 
		# water
		if (self.Vwat_inc > 0):
			
			# hold seed indices in new array
			seedi_here = np.zeros((len(self.seedi))).astype('int')
			seedi_here[:] = self.seedi[:].astype('int')


			# ensure water index is included in indices of
			# seed material for new seed particles
			if self.H2Oi not in seedi_here:

				seedi_here = np.concatenate((seedi_here, 
				np.array((self.H2Oi)).reshape(1)), axis=0)
				
				seedx_now = np.concatenate((seedx_now, 
				np.zeros((num_asb)).reshape(1, -1)), axis=0)

				# indices of non-water seed components
				seedi_nw = (seedi_here != self.H2Oi)
				# indices of water seed component
				seedi_w = (seedi_here == self.H2Oi)
				# ensure mole fractions of non-water seed
				# components in 
				# each size bin sum to 1
				seedx_now = seedx_now*(1./np.sum(seedx_now, axis=0))

			# if water index is already included in indices
			# of seed material for new seed particles, 
			# then zero its mole fraction to begin the 
			# equilibration process
			else:
		
				seedx_now[seedi_here==self.H2Oi, :] = 0.
				# indices of non-water seed components
				seedi_nw = (seedi_here != self.H2Oi)
				# indices of water seed component
				seedi_w = (seedi_here == self.H2Oi)
				# raise sum of mole fractions of non-water
				# components to one
				seedx_now[seedi_nw, :] = 1./np.sum(
				seedx_now[seedi_nw, :], axis=0) 
				
			# first guess of average molar mass of new seed particles 
			# per size bin (g/mol)	
			avMM0 = np.ones((num_asb))

			lcnt = 1 # loop count

			# prepare density of seed components, including water,
			# * by 1.e-3 to convert from kg/m3
			# to g/cm3
			dens_seed = self.y_dens[seedi_here, 0].reshape(-1, 1)*1.e-3

			# average molar mass of seed in each size
			# bin (g/mol)
			avMM1 = np.sum(y_mm[seedi_here]*seedx_now)

			avMM = avMM1 # for calculating Kelvin factor
	
			# maximum difference when comparing across
			# iterations
			while (np.max((avMM1-avMM0)/avMM1) > 0.1):
				
				# average density of seed components 
				# and water (g/cm3)
				av_dens = np.sum(seedx_now*dens_seed, axis=0)

				# calculate Kelvin effect factor for the 
				# provided number size distribution
				# kelvin factor for each size bin 
				# (excluding wall), eq. 16.33 Jacobson 
				# et al. (2005)
				# note that seed_mw has units g/mol, 
				# surfT (g/s2==mN/m==dyn/cm), R_gas is 
				# multiplied by 
				# 1e7 for units g cm2/s2.mol.K, 
				# TEMP is K, x (radius) is multiplied by 
				# 1e-4 to give cm from um and
				kelv = (np.exp((2.e0*avMM*surfT)/
					(R_gas*1.e7*TEMP*(x*1.e-4)*av_dens)))
				
				# equilibrium mole fraction of water per 
				# size bin
				# from the ode solver equation for 
				# vapour-particle partitioning of 
				# water
				xwat = H2Ogc/(self.Psat[0:(num_asb), 
					self.H2Oi]*kelv*act_coeff[0:(num_asb),
				 	self.H2Oi])

				# allow for mole fraction of water in 
				# mole fraction of non-water 
				# seed components
				# for all size bins
				seedxn = np.zeros((seedx_now.shape[0], seedx_now.shape[1]))
				seedxn[seedi_w, :] = xwat
				seedxn[seedi_nw, :] = seedx_now[seedi_nw, :]*(1.-xwat)
				
				# average molar mass of seed particles
				# including water (g/mol) per size bin
				av_MM = np.sum(seedxn*y_mm[seedi_here], axis=0)

				# average liquid-phase density of seed particles
				# including water (g/cm3) per size bin
				av_dens = np.sum(seedxn*dens_seed, axis=0)

				# molar volume averaged over seed 
				# components (including water) (cm3/mol) 
				# per size bin
				avMV = av_MM/av_dens

				# total molecular concentration of seed 
				# components including water 
				# (# molecules/cm3) per size bin, note 
				# that volume multiplied by 
				# 1e-12 to convert from um3 to cm3
				tmc = ((Vperbin*1.e-12)/avMV)*NA

				seed_cnt = 0 # count on seed components

				# concentration of particle-phase seed 
				# components in all size bin
				# loop through indices of seed 
				# components
				for ci in seedi_here: 
					# concentration per size bins (# molecules/cm3)
					yn[ci:
					(num_comp*(num_asb-1)+ci)+1:
					num_comp] = tmc*(seedxn[seed_cnt, :])

					seed_cnt += 1 # count on seed components

				# prepare for checking convergence
				
				# average molar mass of seed in 
				# each size bin (g/mol)
				if (lcnt % 2 != 0): # if on even count
					avMM0 = av_MM
						
				else: # if on odd count
					avMM1 = av_MM
				
				lcnt += 1 # loop count
		
		# if number size distribution does 
		# not include volume of water, i.e. the volume
		# represented by the distribution represents
		# non-water components, and we want to equilibrate
		# water with these particles
		if (self.Vwat_inc == 0):
			
			# hold seed indices in new array
			seedi_here = np.zeros((len(self.seedi))).astype('int')
			seedi_here[:] = self.seedi[:].astype('int')

			# ensure water index is included in indices of
			# seed material for new seed particles
			if self.H2Oi not in seedi_here:
				seedi_here = np.concatenate((seedi_here, 
				np.array((self.H2Oi)).reshape(1)), axis=0)
				seedx_now = np.concatenate((seedx_now, 
				np.zeros((num_asb)).reshape(1, -1)), axis=0)

			# if water index is already included in indices
			# of seed material for new seed particles, 
			# then zero its mole fraction to begin the 
			# equilibration process
			else:
		
				seedx_now[seedi_here==self.H2Oi, :] = 0.

			# indices of non-water seed components
			seedi_nw = (seedi_here != self.H2Oi)
			# indices of water seed component
			seedi_w = (seedi_here == self.H2Oi)
			# ensure mole fractions of non-water seed
			# components in 
			# each size bin sum to 1
			seedx_now = seedx_now*(1./np.sum(seedx_now, axis=0))

			# molar mass of seed 
			# components (g/mol) spread across 
			# size bins
			seed_mm = y_mm[seedi_here] 
			
			# density of seed components, including water,
			# * by 1.e-3 to convert from kg/m3
			# to g/cm3
			dens_seed = self.y_dens[seedi_here, 0].reshape(-1, 1)*1.e-3
	
			# average molar mass of 
			# dry (no water) seed components 
			# (cm3/mol) per size bins
			av_MM = np.sum(seedx_now*seed_mm, axis=0)

			# average liquid-phase density of seed particles
			# including water (g/cm3) per size bin
			av_dens = np.sum(seedx_now*dens_seed, axis=0)

			# molar volume averaged over seed 
			# components (including water) (cm3/mol) 
			# per size bin
			avMV = np.sum(av_MM/av_dens, axis=0)				

			# calculate Kelvin effect factor 
			# for the provided number size 
			# distribution
			# kelvin factor for each size bin 
			#(excluding wall), eq. 16.33 
			# Jacobson et al. (2005)
			# note that avMW has units 
			# g/mol, surfT (g/s2==mN/m==dyn/cm), 
			# R_gas is multiplied by 
			# 1e7 for units g cm2/s2.mol.K, 
			# TEMP is K, x (radius) is multiplied 
			# by 1e-4 to give cm from um and 
			# y_dens multiplied by  by 1e-3 to 
			# convert from kg/m3 to g/cm3
			kelv = np.exp((2.e0*av_MM*surfT)/(R_gas*1.e7*TEMP*(x)*1.e-4*av_dens))
								
			# equilibrium mole fraction of water from the ode solver 
			# equation for vapour-particle partitioning per size bin, note that
			# self.Psat has particle size bins and walls in rows
			xwat = H2Ogc/(
			self.Psat[0:num_asb, self.H2Oi]*kelv*act_coeff[0:num_asb, self.H2Oi])
				
			# total molecular concentration of dry (no water) seed 
			# components 
			# (molecules/cm3) per size bin, 
			# note that volume multiplied by 1e-12 to convert 
			# from um3 to cm3 (this is no water seed components
			# as the condition above says that volume size 
			# distribution exclused volume of water)
			tmc = ((Vperbin*1.e-12)/avMV)*NA

			# the molecular concentration of water per size bin
			# (molecules/cm3), note that (xwat/(1.-wat)) gives 
			# the fraction of tmc that gives the molecular 
			# concentration of water
			wmc = tmc*(xwat/(1.-xwat))
			seed_cnt = 0 # count on seed components

			# concentration of particle-phase seed 
			# components in all size bin
			# loop through indices of seed 
			# components
			for ci in seedi_here: 
				
				if (ci != self.H2Oi): # if non-water
				
					# concentration per size bins (# molecules/cm3)
					yn[ci:
					(num_comp*(num_asb-1)+ci)+1:
					num_comp] = tmc*(seedx_now[seed_cnt, :])

				else: # if water
					
					# concentration per size bins (# molecules/cm3)
					yn[ci:
					(num_comp*(num_asb-1)+ci)+1:
					num_comp] = wmc

				seed_cnt += 1 # count on seed components

	# if water not to be equilibrated with seed 
	# particles prior to experiment start (possibly because
	# particle-phase mole fraction of water already supplied by
	# user)
	if (self.seed_eq_wat == 0):			
	
		# initial seed particles will not be equilibrated with 
		# water here, therefore water
		# will try to equilibrate between the gas, 
		# particle and any other surface 
		# through the ODE solver		
		
		# hold seed indices in new array
		seedi_here = np.zeros((len(self.seedi))).astype('int')
		seedi_here[:] = self.seedi[:].astype('int')

		# ensure where no mole fractions are present (e.g. because a size bin
		# has no new particles, seedx_now is 0 rather than nan)
		seedx_zeros_indx = np.sum(seedx_now, axis=0) == 0.
	
		# index of size bins where mole fractions are not zero
		seedx_nzeros_indx = np.sum(seedx_now, axis=0) != 0.

		# ensure seedx sums to 1 per size bin (columns) across
		# components (rows)
		seedx_now[:, seedx_nzeros_indx] = (seedx_now[:, seedx_nzeros_indx]/
			np.sum(seedx_now[:, seedx_nzeros_indx], axis=0))
		
		# zero any nan values
		seedx_now[:, seedx_zeros_indx] = 0.

		# molar mass of seed 
		# components (g/mol) (rows)
		seed_mm = y_mm[seedi_here] 
			
		# density of seed components (rows), including water,
		# * by 1.e-3 to convert from kg/m3
		# to g/cm3
		dens_seed = self.y_dens[seedi_here, 0].reshape(-1, 1)*1.e-3
	
		# average molar mass of 
		# dry (no water) seed components (rows)
		# (cm^3/mol) per size bin (columns)
		av_MM = (np.sum(seedx_now*seed_mm, axis=0)).reshape(num_asb)

		# average liquid-phase density of seed particles
		# including water (g/cm^3) per size bin
		av_dens = (np.sum(seedx_now*dens_seed, axis=0)).reshape(num_asb)

		# molar volume averaged over seed 
		# components (including water) (cm^3/mol) 
		# per size bin (columns)
		avMV = np.zeros((num_asb))
		avMV[seedx_nzeros_indx] = np.sum(av_MM[seedx_nzeros_indx]/
			av_dens[seedx_nzeros_indx], axis=0)

		# total molecular concentration of seed 
		# components 
		# (molecules/cm^3) per size bin, 
		# note that volume multiplied by 1e-12 to convert 
		# from um^3 to cm^3, multiply by NA to convert
		# moles/cm^3 to molecules/cm^3
		tmc = ((Vperbin*1.e-12)/avMV)*NA

		# convert any nans to 0
		tmc[seedx_zeros_indx] = 0.
		
		seed_cnt = 0 # count on seed components

		# concentration of particle-phase seed 
		# components in all size bin
		# loop through indices of seed 
		# components
		for ci in seedi_here: 
					
			# concentration per size bins (# molecules/cm3)
			yn[ci:(num_comp*(num_asb-1)+ci)+1:
			num_comp] = tmc*(seedx_now[seed_cnt, :])

			seed_cnt += 1 # count on seed components
	
	return(yn)