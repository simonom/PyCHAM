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
'''module to estimate the particle and wall partitioning coefficient'''
# the kimt_calc module is called at the start of the model loop time 
# interval to update the mass transfer coefficient of gases to particles 
# and gases to walls

import importlib
import numpy as np
from part_prop import part_prop
import scipy.constants as si

def kimt_calc(y, mfp, num_sb, num_comp, accom_coeff, y_mw, surfT, 
		R_gas, TEMP, NA, N_perbin, radius, therm_sp,
		H2Oi, act_coeff, caller, Press, 
		DStar_org, z_prt_coeff, chamSA, chamV, self):
	
	# inputs:-----------------------------------------------------
	
	# y - concentrations of components (# molecules/cm3 (air))
	# mfp - mean free path of gas molecules (m) (num_comp, 1)
	# accom_coeff - accommodation coefficients of components in each 
	# size bin including wall
	# y_mw - molecular weight of components (g/mol)
	# radius - particle radius (m)
	# surfT - surface tension (g/s2==mN/m==dyn/cm)
	# R_gas - the universal gas constant (cm3.Pa/K.mol)
	# TEMP - current temperature in chamber (K)
	# NA - Avogadro's constant (molecules/mol) 
	# self.y_dens - density of components (kg/m3)
	# N_perbin - number of particles per size bin (excluding wall)
	# self.Psat - liquid-phase saturation vapour pressures of 
	# components (# molecules/cm3 (air))	
	# therm_sp - thermal speed of components (m/s) (num_comp, 1)
	# H2Oi - water index (integer)
	# act_coeff - activity coefficient of components (dimensionless)
	# self.wall_on - marker for whether wall present
	# caller - marker for the calling function
	# self.partit_cutoff - the product of self.Psat and act_coeff above
	#	 which gas-particle partitioning assumed zero (Pa)
	# Press - air pressure (Pa)
	# DStar_org - gas-phase diffusion coefficients of components 
	#	(cm2/s)
	# z_prt_coeff - fraction of total gas-particle partitioning 
	#	coefficient 
	#	below which partitioning to a particle size bin is 
	#	treated as zero,
	#	e.g. because surface area of that size bin is tiny
	# chamSA - chamber surface area (m2)
	# chamV - chamber volume (m3)
	# self.kwf - gas-wall partitioning coefficient flag 
	#	(-1 means treat with Huang et al. 2018)
	# self - reference to program
	# ---------------------------------------------------
	
	if (num_sb == 0): # fillers
		kimt = np.zeros((num_sb, num_comp))
		kelv_fac = np.zeros((num_sb-self.wall_on, 1))

		return(kimt, kelv)

	if (num_sb > 0) and (self.wall_on > 0): # if wall present
		y_part = y[num_comp:-(num_comp*self.wall_on)]
	if (num_sb > 0) and (self.wall_on == 0): # if wall absent
		y_part = y[num_comp::]
	
	if (num_sb-self.wall_on > 0): # if particles present
		
		# density (g/cm3) and average molecular weight (g/mol) of 
		# particles (excluding wall)
		[tot_rho, ish, avMW] = part_prop(y_part, num_comp, 
			(num_sb-self.wall_on), NA, y_mw, N_perbin, self)
	
	
		# Knudsen number (dimensionless)
		Kn = np.repeat(mfp, 
		(num_sb-self.wall_on), 1)/np.repeat(radius, num_comp, 0)
	
		# update accommodation coefficients
		import accom_coeff_calc
		# import most recent version
		importlib.reload(accom_coeff_calc)
		accom_coeff_now = accom_coeff_calc.accom_coeff_func(
			accom_coeff, radius, self)
		
		# Non-continuum regime correction 
		# calculate a correction factor according to the
		# continuum versus non-continuum 
		# regimes
		# expression taken from Jacobson et al (2000), 
		# page 457, or Jacobson (2005), page 530 
		# (eq. 16.19).
		# They reference:
		# Fuchs and Sutugin 1971
		# Pruppacher and Klett 1997
		Inverse_Kn = Kn**-1.
		correct_1 = (1.33+0.71*Inverse_Kn)/(1.+Inverse_Kn)
		correct_2 = ((4.*(1.-accom_coeff_now[:, 
			0:num_sb-self.wall_on]))/
			(3.*accom_coeff_now[:, 0:num_sb-self.wall_on]))
		correct_3 = 1.e0+(correct_1+correct_2)*Kn
		correction = correct_3**-1.
		
		# kelvin factor for each size bin (excluding wall), 
		# eq. 16.33 Jacobson et al. (2005)
		# note that avMW has units g/mol, surfT 
		# (g/s2==mN/m==dyn/cm), R_gas is multiplied by 
		# 1e7 for units g cm2/s2.mol.K, 
		# TEMP is K, radius is multiplied by 1e2 to give cm 
		# and tot_rho is g/cm3
		kelv = np.zeros((num_sb-self.wall_on, 1))
		kelv[ish, 0] = np.exp((2.e0*avMW[ish]*surfT)/
			(R_gas*1.e7*TEMP*radius[0, ish]*
			1.e2*tot_rho[ish]))
	
		# ------------------------------------------------
		# gas-phase diffusion coefficient*Fuch-Sutugin 
		# correction (cm2/s)
		# eq. 5 Zaveri et al. (2008) 
		# (https://doi.org/10.1029/2007JD008782)
		kimt = (DStar_org)*correction
		
		# final partitioning coefficient (converting radius 
		# from m to cm)
		# eq. 16.2 of Jacobson (2005) and eq. 5 
		# Zaveri et al. (2008)
		# (https://doi.org/10.1029/2007JD008782)
		# components in rows and size bins in columns (/s)
		kimt = (4.*np.pi*(radius*1.e2)*N_perbin.reshape(1, 
			-1))*kimt
		
		# zero partitioning coefficient for particle size bins
		# with 
		# such little number concentration or radius that 
		# partitioning is relatively tiny
		if (z_prt_coeff > 0.):
			kimt[:, np.sum(kimt, axis = 0)/
				np.sum(np.sum(kimt)) < z_prt_coeff] = 0.
		
		# zero partitioning coefficient for any components 
		# with relatively tiny abundance 
		# in the gas and particle phase
		# when compared against non-water and non-seed 
		# components
		if (z_prt_coeff > 0. and num_sb-self.wall_on > 0):
			
			# get gas and particle phase concentrations
			if (self.wall_on > 0): # if wall present
				y_gp = np.zeros((len(y[0:
					-(num_comp*self.wall_on)])))
				y_gp[:] = y[0:-(num_comp*self.wall_on)]
			if (self.wall_on == 0): # if wall absent
				y_gp = np.zeros((len(y[:])))
				y_gp[:] = y[:]

			# rearrange so that components in rows and 
			# gas/size bins in columns
			y_gp = y_gp.reshape(num_comp, 
				num_sb-self.wall_on+1, order='F')
			# zero water
			y_gp[H2Oi, :] = 0.
			# zero seed
			y_gp[self.seedi, :] = 0.
			# sum over gas/size bins for each component
			y_gp = np.sum(y_gp, axis=1)
			y_gp_frac = y_gp/np.sum(y_gp)
			# if any components in tiny abundance but are present
			if (sum((y_gp_frac < z_prt_coeff)*(y_gp > 0)) >  0): 
				kimt[((y_gp_frac < z_prt_coeff)*(y_gp > 0)), :] = 0.
		
		# transpose kimt ready for multiplication inside ode solver, so size bins
		# in rows and components in columns
		kimt = np.transpose(kimt)
		
		# zero partitioning to particles for any components with low condensability
		# if a value provided (default is empty list)
		if (self.partit_cutoff):

			# convert partit_cutoff from Pa to molecules/cm3
			# (air), note README states
			# that just one value accepted for partit_cutoff 
			# input
			partit_cutoff_molcm = (self.partit_cutoff[0]*(
			NA/((R_gas*1.e6)*TEMP)))
			highVPi = ((self.Psat*act_coeff) > 
			partit_cutoff_molcm)
			# mask water to allow its partitioning
			highVPi[:, H2Oi] = 0
			kimt[highVPi[0:num_sb-self.wall_on, :]] = 0.
	
	else: # if no particles
		kimt = np.zeros((num_sb-self.wall_on, num_comp))
		kelv = np.zeros((num_sb-self.wall_on, 1))

	if (self.kwf == -1):
		# gas-wall partitioning coefficient (/s), from Huang et 
		# al. 2018 
		# (Eq. 2 and accompanying text), 
		# https://doi.org/10.1021/acs.est.7b05575
		# mass transport coefficient across gas-phase boundary 
		# layer, note the 
		# assumption of ke = 1./40. derived from fitting decay 
		# of ozone to observed
		# decay from MAC
		ve = ((2./np.pi)*((1./40.*DStar_org*1.e-4)**0.5))
		vc = 1.*therm_sp/4.

		# for Dekati Oxidation Flow Reactor
		# eddy diffusivity coefficient
		ke = 4.e-3+10**-2.25*chamV**0.74
		ve = ((2./np.pi)*((ke*DStar_org*1.e-4)**0.5))	
		vc = 0. # because glass not teflon

		if (vc != 0.):
			self.kw = (chamSA/chamV)*((1./ve + 1./vc)**-1)
		if (vc == 0.):
			self.kw = (chamSA/chamV)*((1./ve)**-1)

		# spread over wall bins
		self.kw = np.tile(self.kw.reshape(1, num_comp), 
			(self.wall_on, 1))
		
	
	# zero partitioning to walls for any components with low 
	# condensability
	# if a value provided (default is empty list)
	if (self.partit_cutoff):

		# convert partit_cutoff from Pa to molecules/cm3 (air), 
		# note README states
		# that just one value accepted for partit_cutoff input
		partit_cutoff_molcm = self.partit_cutoff[0]*(
		NA/((R_gas*1.e6)*TEMP))
		highVPi = (self.Psat*act_coeff) > partit_cutoff_molcm
		# mask water to allow its partitioning
		highVPi[:, H2Oi] = 0 
		self.kw[:, highVPi[-1, :]] = 0.


	# make any necessary adjustments for vapour pressure of 
	# components as a function of
	# wall, as defined by the user
	if 'self.P_wfunc' in locals():
		self.Psat[self.P_wfunc_wi, 
		self.P_wfunc_ci] = self.P_wfunc
	
	# concatenate kw onto kimt, ready for ode solver
	if (self.wall_on > 0):
		kimt = np.concatenate((kimt, self.kw), axis = 0)
	
	return(kimt, kelv)
