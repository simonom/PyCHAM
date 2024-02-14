########################################################################
#								       #
# Copyright (C) 2018-2024					       #
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
'''module to set up particle phase part of box model, calls on 
init_water_partit to initiate water partitioning with seed particles 
and wall'''
# gets called when there is an injection of seed particle during an 
# experiment

import numpy as np
from init_water_partit import init_water_partit
import scipy.constants as si
from scipy import stats # import the scipy.stats module

def pp_dursim(y, N_perbin0, mean_rad, pmode, pconc, seedx, lowersize, uppersize, num_comp, 
		num_sb, MV, rad0, radn, std, H2Oi, rbou, y_mw, surfT, TEMP, 
		act_coeff, pcontf, H2Ogc, self):
	
			
	# inputs -----------------------------------
	# y - concentrations of components in particle phase (# molecules/cm3 (air))
	# N_perbin0 - starting number concentration of particles (# particles/cm3 (air))
	# mean_rad - mean radius of seed particles at this time (um)
	# pmode - whether particle number size distribution stated by mode or explicitly
	# pconc - number concentration of seed particles (# particles/cm3 (air))
	# self.seedi - index of seed material
	# seedx - mole ratio of non-water components comprising seed particles
	# lowersize - smallest radius bound (um)
	# uppersize - greatest radius bound (um)
	# num_comp - number of components
	# num_sb - number of size bins (excluding wall)
	# MV - molar volume of all components (cm3/mol) (shape: number of components, 1)
	# rad0 - original radius at size bin centres (um)
	# radn - current radius at size bin centres (um)
	# std - standard deviation for lognormal size distribution calculation (dimensionless)
	# self.y_dens  - density of components (kg/m3)
	# H2Oi - index of water
	# rbou - radius bounds per size bin (um)
	# y_mw - molecular weight of components (g/mol)
	# surfT - surface tension of particles (g/s2 == mN/m == dyn/cm)
	# TEMP - chamber temperature (K)
	# self.Psat - saturation vapour pressure of components (# molecules/cm3 (air))	
	# act_coeff - activity coefficient of components
	# self.seed_eq_wat - whether seed particles to be equilibrated with water prior to ODE solver
	# self.Vwat_inc - whether suppled seed particle volume contains equilibrated water
	# pcontf - flag for whether injection of particles is continuous or instantaneous
	# H2Ogc - gas-phase concentration of water (# molecules/cm3)
	# self - reference to program
	# ------------------------------------------
	
	# if mean radius not stated explicitly calculate from size ranges (um)
	if (any(mean_rad == -1.e6) and num_sb > 0):
		if (lowersize > 0.):
			mean_rad = 10**((np.log10(lowersize)+np.log10(uppersize))/2.)
		else:
			mean_rad = 10**(np.log10(uppersize)-1)
	
	R_gas = si.R # ideal gas constant (kg.m2.s-2.K-1.mol-1)
	NA = si.Avogadro # Avogadro's number (molecules/mol)
	
	# prepare for new concentration of particles (# particles/cm3)
	N_perbin = np.zeros((N_perbin0.shape[0], N_perbin0.shape[1])) 

	if (num_sb == 1):
		
		N_perbin[:, :] = N_perbin0[:, :] + np.array((pconc)) # (# particles/cm3 (air))
		pconc_new = pconc

	# number concentration stated per size bin in multi size bin simulation
	if (pmode == 1 and num_sb > 1):
		N_perbin[:, :] = N_perbin0[:, :] + np.array((pconc)).reshape(-1, 1) # (# particles/cm3 (air))
		pconc_new = pconc

	# total number concentration per mode stated in multi size bin simulation
	if (pmode == 0 and num_sb > 1):
		
		for i in range(len(pconc)): # loop through modes
			# set scale and standard deviation input for lognormal probability distribution 
			# function, following guidance here: 
			# http://all-geo.org/volcan01010/2013/09/how-to-use-lognormal-distributions-in-python/
			if np.isscalar(mean_rad):
				scale = np.exp(np.log(mean_rad))
			else:
				scale = np.exp(np.log(mean_rad[i]))
			if np.isscalar(std):
				std_now = np.log(std)
			else:
				std_now = np.log(std[i])
			loc = 0. # no shift
		
			# number fraction-size distribution - enforce high resolution to ensure size
			# distribution of seed particles fully captured
			if (lowersize > 0.): # if lowermost size bin bound useful
				hires = 10**(np.linspace(np.log10(lowersize), np.log10(uppersize), int(num_sb*1.e2)))
			else: # enforce lowermost size bin radius bound of 1 nm (1e-3 um)
				hires = 10**(np.linspace(np.log10(1.e-3), np.log10(uppersize), int(num_sb*1.e2)))
			pdf_output = stats.lognorm.pdf(hires, std_now, loc, scale)
			pdf_out = np.interp(rad0, hires, pdf_output)	
			# number concentration of seed in all size bins (# particle/cm3 (air))
			pconc_new = (pdf_out/sum(pdf_out))*pconc[i]
				
			# number concentration realism
			pconc_new[pconc_new < 1.e-2] = 0.
			
			N_perbin[:, 0] = N_perbin0[:, 0] +  pconc_new # (# particles/cm3 (air))
					
	# volume concentration of new seed particles (um3/cm3 (air))
	Vperbin = ((pconc_new*(4./3.)*np.pi*(rad0)**3.))

	# concentrations of components in new seed particles (# molecules/cm3 (air))
	yn = np.zeros((num_comp*(num_sb)))
	
	# account for particle-phase concentration of components contained in seed particles --------------------
	
	# check whether water to be equilibrated with seed particles prior to experiment start
	if (self.seed_eq_wat == 1 or self.Vwat_inc == 1): # if yes, water is to be equilibrated
		
		# check whether the stated initial number size distribution included the 
		# volume of water
		if (self.Vwat_inc == 1): # if number size distribution does include volume of water

			avMW0 = np.ones((num_sb)) # first guess of average molecular weight	

			# average molecular weight of seed in each size bin (g/mol)
			avMW1 = np.ones((num_sb))*(np.sum(y_mw[self.seedi])/len(self.seedi))

			avMW = avMW1 # for calculating Kelvin factor
			lcnt = 1 # loop count
	
			while (np.max((avMW1-avMW0)/avMW1) > 0.1):

				# calculate Kelvin effect factor for the provided number size distribution
				# kelvin factor for each size bin (excluding wall), eq. 16.33 Jacobson et al. (2005)
				# note that seed_mw has units g/mol, surfT (g/s2==mN/m==dyn/cm), R_gas is multiplied by 
				# 1e7 for units g cm2/s2.mol.K, 
				# TEMP is K, x (radius) is multiplied by 1e-4 to give cm from um and 
				# y_dens multiplied by  by 1e-3 to convert from kg/m3 to g/cm3
				kelv = np.exp((2.e0*avMW*surfT)/(R_gas*1.e7*TEMP*(radn)*1.e-4*(sum(self.y_dens[self.seedi[:], 0]*1.e-3)/len(self.seedi)*1.e-3)))
				
				# equilibrium mole fraction of water per size bin
				# from the ode solver equation for vapour-particle partitioning of 
				# water
				xwat = H2Ogc/(self.Psat[0:(num_sb), H2Oi]*kelv*act_coeff[0:(num_sb),
				 H2Oi])
				
				# allow for mole fraction of water in mole fraction of non-water 
				# seed components
				# for all size bins
				# ensure the non-water mole fractions sum to one
				seedx = seedx*(1./sum(seedx))
				seedxn = seedx*(1.-xwat)
				
				# average molar volume of seed components (cm3/mol) 
				# for all size bins
				avMV = (sum(seedxn*MV[self.seedi[:], 0].reshape(-1, 1))
						+xwat*MV[H2Oi, 0])
				# total molecular concentration of seed components including water 
				# (# molecules/cm3) per size bin, note that volume multiplied by 
				# 1e-12 to convert from um3 to cm3
				tmc = ((Vperbin*1.e-12)/avMV)*NA

				# concentration of particle-phase seed components in all size bin
				# loop through indices of seed components
				for ci in range(len(self.seedi)): 
					
					# non-water (# molecules/cm3)
					yn[self.seedi[ci]:(num_comp*(num_sb-1)+self.seedi[ci])+1:num_comp] = tmc*(seedxn[ci, :])
					
					if (ci == len(self.seedi)-1): # reached water component
						# water (# molecules/cm3)
						yn[H2Oi:(num_comp*(num_sb-1)+H2Oi)+1:num_comp] = tmc*xwat

				# for average molecular weight in new particles, first get fraction of 
				# each component in each size bin
				avMW = yn.reshape(num_comp, num_sb, order='F')
				
				# size bins with no particles
				ish = np.where(np.sum(avMW, axis=0) == 0.)[0][:]
				avMW = avMW/np.sum(avMW, axis=0)
				# fill the empty size bin values
				avMW[:, ish] = 1./num_comp
				
				# average molecular weight of seed in each size bin (g/mol)
				if (lcnt % 2 != 0): # if on even count
					avMW0 = np.sum(avMW*y_mw.reshape(-1, 1), axis = 0)
					avMW = avMW0 # for calculating Kelvin factor
						
				else: # if on odd count
					avMW1 = np.sum(avMW*y_mw.reshape(-1, 1), axis = 0)
					avMW = avMW1 # for calculating Kelvin factor
				
				lcnt += 1 # loop count
		

		if (self.Vwat_inc == 0): # if number size distribution does not include volume of water
			
				seed_mw = y_mw[self.seedi] # molecular weight of seed components (g/mol)

				seedx = seedx*(1./sum(seedx)) # ensure the non-water mole fractions sum to one

				# average molecular weight of dry (no water) seed components (cm3/mol) for all size bins
				avMW = (sum(seedx*seed_mw))				

				# average molar volume of dry (no water) seed components (cm3/mol) for all size bins
				avMV = (sum(seedx*MV[self.seedi[:], 0]))

				# calculate Kelvin effect factor for the provided number size distribution
				# kelvin factor for each size bin (excluding wall), eq. 16.33 Jacobson et al. (2005)
				# note that seed_mw has units g/mol, surfT (g/s2==mN/m==dyn/cm), R_gas is multiplied by 
				# 1e7 for units g cm2/s2.mol.K, 
				# TEMP is K, x (radius) is multiplied by 1e-4 to give cm from um and 
				# y_dens multiplied by  by 1e-3 to convert from kg/m3 to g/cm3
				kelv = np.exp((2.e0*avMW*surfT)/(R_gas*1.e7*TEMP*(radn)*1.e-4*(sum(self.y_dens[self.seedi[:], 0]*1.e-3)/len(self.seedi)*1.e-3)))
				
				# equilibrium mole fraction of water from the ode solver 
				# equation for vapour-particle partitioning
				xwat = H2Ogc/(self.Psat[:, H2Oi]*kelv*act_coeff[:, H2Oi])
				
				# total molecular concentration of dry (no water) seed components 
				# including water (molecules/cm3) per size bin, 
				# note that volume multiplied by 1e-12 to convert from um3 to cm3
				tmc = ((Vperbin*1.e-12)/avMV)*NA
				
				# concentration of particle-phase seed components in all size bin
				for ci in range(len(self.seedi)): # loop through indices of seed components

					# non-water (# molecules/cm3)
					yn[self.seedi[ci]:(num_comp*(num_sb-1)+self.seedi[ci])+1:num_comp] = tmc*(seedx[ci, :])
						
					if (ci == len(self.seedi)-1): # reached final non-water seed component
						# water (# molecules/cm3)
						yn[H2Oi:(num_comp*(num_sb-1)+H2Oi)+1:num_comp] = (xwat*tmc)/(1.-xwat)

	# if water not to be equilibrated with seed particles prior to experiment start
	if (self.seed_eq_wat == 0 and self.Vwat_inc == 0):			
	
		# initial seed particles will not contain water and any water vapour in chamber
		# will try to equilibrate with particles through the ODE solver		
		
		# ensure seedx sums to 1
		seedx = seedx/sum(seedx)

		# mole-fraction weighted molar volume (cm3/mole) (average molar volume of particles)
		mfwMV = sum(MV[self.seedi[:], 0]*seedx)
		# convert to molecular volume (cm3/molecule)
		mfwMV = mfwMV/NA
		# total molecular concentration of seed components per size bin (# molecules/cm3)
		# note that volume multiplied by 1e-12 to convert from um3 to cm3
		ytot = (Vperbin*1.e-12)/mfwMV
		
		for ci in range(len(self.seedi)): # loop through indices of seed components
		
			# concentration of this component in all size bins (# molecules/cm3 (air)):
			yn[self.seedi[ci]:(num_comp*(num_sb-1)+self.seedi[ci])+1:num_comp] = ytot*seedx[ci]

	# if instantaneous injection of particles
	# factor concentrations of components comprising 
	# new seed particles into existing concentration (# molecules/cm3)
	if (pcontf == 0):
		y += yn
		
	if (pcontf == 1): # if continuous injection of particles
		y += yn # adding to y array here rather than in ode_solv

	# loop through size bins to estimate new total volume concentrations (um3/cc (air))
	Vtot = np.zeros((num_sb))
	Varr = np.zeros((num_sb)) # volume concentration of single particles (um3/cc (air))
	mass_conc = 0. # mass concentration of particles (g/cm3)
	
	for i in range(num_sb): # size bin loop
		Vtot[i] = np.sum(y[num_comp*i:num_comp*(i+1)]/NA*(MV[:, 0]*1.e12))
		# new volume concentration of single particles (um3/cm3 (air))
		if (N_perbin[i, 0] > 0.):
			Varr[i] = Vtot[i]/N_perbin[i, 0]
		else:
			Varr[i] = (4./3.*np.pi)*rad0[i]**3.
		# multiply y_dens by 1e-9 to get ug/um3 (particle) from kg/m3, 
		# then multiplying ug/um3 (particle) by um3/cm3 (air) gives ug/cm3 (air)
		mass_conc += np.sum((self.y_dens[self.seedi, 0]*1.e-9)*Vtot[i])
	
	# new radius of single particles (um)
	radn = ((3./(4.*np.pi))*Varr)**(1./3.)
	

	return(y, N_perbin, radn, Varr)