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
'''module to set up particle phase part of box model'''
# gets called when there is an injection of seed particle during an 
# experiment

import numpy as np
from pp_water_equil import pp_water_equil
import scipy.constants as si
from scipy import stats # import the scipy.stats module

def pp_dursim(y, N_perbin0, mean_rad, pconc, lowersize, 
		uppersize, num_comp, 
		num_sb, MV, std, H2Oi, rbou, y_mm, surfT, 
		TEMP, act_coeff, pcontf, H2Ogc, 
		x, self):
	
			
	# inputs -----------------------------------
	# y - concentrations of components in particle phase 
	# 	(# molecules/cm3 (air))
	# N_perbin0 - starting number concentration of particles 
	# 	(# particles/cm3 (air))
	# mean_rad - mean radius of seed particles at this time (um)
	# self.pmode - whether particle number size distribution stated by 
	#	 mode or explicitly
	# pconc - number concentration of seed particles now
	# 	(# particles/cm3 (air))
	# self.seedi - index of seed material
	# self.seedx - mole ratio of non-water components comprising seed 
	# 	particles
	# lowersize - smallest radius bound (um)
	# uppersize - greatest radius bound (um)
	# num_comp - number of components
	# num_sb - number of size bins (excluding wall)
	# MV - molar volume of all components (cm3/mol) 
	#	(shape: number of components, 1)
	# std - standard deviation for lognormal size distribution 
	# 	calculation (dimensionless)
	# self.y_dens  - density of components (kg/m3)
	# H2Oi - index of water
	# rbou - radius bounds per size bin (um)
	# y_mm - molar mass of components (g/mol)
	# surfT - surface tension of particles (g/s2 == mN/m == dyn/cm)
	# TEMP - chamber temperature (K)
	# self.Psat - saturation vapour pressure of components 
	# 	(# molecules/cm3 (air))	
	# act_coeff - activity coefficient of components
	# self.seed_eq_wat - whether seed particles to be equilibrated 
	# 	with water prior to ODE solver
	# self.Vwat_inc - whether suppled seed particle volume 
	# contains equilibrated water
	# pcontf - flag for whether injection of particles is continuous 	
	#	or instantaneous
	# H2Ogc - gas-phase concentration of water (# molecules/cm3)
	# x - particle radius per size bin (um)
	# self - reference to PyCHAM object
	# ------------------------------------------
	
	# if mean radius not stated explicitly calculate from size 
	# ranges (um)
	if (sum(mean_rad[:] == -1.e6) > 0 and num_sb > 0):
		if (lowersize > 0.):
			radn = 10**((np.log10(lowersize)+
			np.log10(uppersize))/2.)
		else:
			radn = 10**(np.log10(uppersize)-1)

	# if mean radius is given (either per size bin)
	# or for number-size distribution modes
	else:
		radn = mean_rad[:]
	
	R_gas = si.R # ideal gas constant (kg.m2.s-2.K-1.mol-1)
	NA = si.Avogadro # Avogadro's number (molecules/mol)
	
	# prepare for new concentration of particles (# particles/cm3)
	N_perbin = np.zeros((N_perbin0.shape[0], N_perbin0.shape[1])) 

	if (num_sb == 1):
		
		# (# particles/cm3 (air))
		N_perbin[:, :] = N_perbin0[:, :] + np.array((pconc))
		pconc_new = pconc

	# number concentration stated per size bin in multi size bin 
	# simulation
	if (self.pmode == 1 and num_sb > 1):
		# if this is the new overall particle number 
		# concentration, e.g. from an observation 
		# file, as interpretted by obs_file_open
		if (self.pp_dil == 0):
			# (# particles/cm3 (air))
			N_perbin[:, :] = (np.array((
			pconc)).reshape(-1, 1))
		# if this new particle number represents
		# injection of new particle in addition to
		# the existing particles
		else:
			# (# particles/cm3 (air))
			N_perbin[:, :] = (N_perbin0[:, :] + 
			np.array((pconc)).reshape(-1, 1))

		pconc_new = pconc

	# total number concentration per mode stated in modal 
	# representation of multi size bin 
	# simulation
	if (self.pmode == 0 and num_sb > 1):
		
		# if this is the new overall particle number 
		# concentration, so inherently including any 
		# previously present particles, e.g. from an observation 
		# file, as interpreted by obs_file_open
		if (self.pp_dil == 0):
			N_perbin[:, :] = 0.
		# if this is adding to any previously present
		# particles
		N_perbin[:, 0] = N_perbin0[:, 0]
		
		for i in range(len(pconc)): # loop through modes
			# set scale and standard deviation input for 
			# lognormal probability distribution 
			# function, following guidance here: 
			# http://all-geo.org/volcan01010/2013/09/how-to-
			# use-lognormal-distributions-in-python/
			scale = np.exp(np.log(radn[i]))
			if np.isscalar(std):
				std_now = np.log(std)
			else:
				std_now = np.log(std[i])
			loc = 0. # no shift
		
			# number fraction-size distribution - enforce 
			# high resolution to ensure size
			# distribution of seed particles fully captured
			# if lowermost size bin bound useful
			if (lowersize > 0.): 
				hires = 10**(np.linspace(np.log10(
				lowersize), np.log10(uppersize), 
				int(num_sb*1.e2)))
			# enforce lowermost size bin radius bound of 
			# 1 nm (1e-3 um)
			else: 
				hires = 10**(np.linspace(
				np.log10(1.e-3), np.log10(uppersize), 
				int(num_sb*1.e2)))
			pdf_output = stats.lognorm.pdf(hires, std_now, 
			loc, scale)
			pdf_out = np.interp(x, hires, pdf_output)	
			# number concentration of seed in all size bins
			# (# particle/cm3 (air))
			pconc_new = (pdf_out/sum(pdf_out))*pconc[i]
				
			# number concentration realism
			pconc_new[pconc_new < 1.e-2] = 0.
			
			# include in number-size distribution 
			# array (# particles/cm3 (air))
			N_perbin[:, 0] += pconc_new
					
	# volume concentration of new seed particles (um3/cm3 (air))
	# per size bin summed across components
	Vperbin = pconc_new*((4./3.)*np.pi*(radn)**3.)

	# container for concentrations of components in new seed particles 
	# (# molecules/cm3 (air))
	yn = np.zeros((num_comp*(num_sb)))
	
	# account for particle-phase concentration of components 
	# contained in new particles --------
	# get supplied mole fractions of new seed material
	# components in rows, size bins in columns, note that
	# self.seedx_tcnt is set in cham_up
	seedx_now = np.squeeze(self.seedx[:, :, self.seedx_tcnt]).reshape(
		self.seedx.shape[0], self.seedx.shape[1])
	
	# consider water equilibrium between gas and particles
	yn = pp_water_equil(H2Ogc, yn, seedx_now, num_sb, y_mm, 
		R_gas, TEMP, surfT, act_coeff, Vperbin, radn, 
		num_comp, self)
	
	# if instantaneous injection of particles
	# factor concentrations of components comprising 
	# new seed particles into existing concentration (# molecules/cm3)
	if (pcontf == 0):
		# if this is the new overall particle component 
		# concentration, e.g. from an observation 
		# file, as interpretted by obs_file_open
		if (self.pp_dil == 0):
			y = yn
		# if this new particle represents
		# injection of new particle in addition to
		# the existing particles
		else:
			y += yn
		
	if (pcontf == 1): # if continuous injection of particles
		y += yn # adding to y array here rather than in ode_solv
	
	# loop through size bins to estimate new total 
	# volume concentrations (um3/cm3 (air))
	Vtot = np.zeros((num_sb))
	# volume concentration of single 
	# particles (um3/cm3 (air))
	Varr = np.zeros((num_sb))
	mass_conc = 0. # mass concentration of particles (g/cm3)
	
	# prepare to get new mean radius of particles per size bin (um)
	for i in range(num_sb): # size bin loop
		Vtot[i] = np.sum(y[num_comp*i:num_comp*(i+1)]/NA*(MV[:, 0]*1.e12))
		# new volume concentration of single particles (um3/cm3 (air))
		if (N_perbin[i, 0] > 0.):
			Varr[i] = Vtot[i]/N_perbin[i, 0]
		else:
			Varr[i] = (4./3.*np.pi)*x[i]**3.
		# multiply y_dens by 1e-9 to get ug/um3 (particle) from kg/m3, 
		# then multiplying ug/um3 (particle) by um3/cm3 (air) gives ug/cm3 (air)
		mass_conc += np.sum((self.y_dens[self.seedi, 0]*1.e-9)*Vtot[i])
	
	# new radius of single particles (um)
	radn = ((3./(4.*np.pi))*Varr)**(1./3.)

	return(y, N_perbin, radn, Varr)
