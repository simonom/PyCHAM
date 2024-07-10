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
'''update the chamber variables depending on time'''
# based on the user inputs and time through simulation,
# chamber variables are updated before calling the 
# ODE solver, so they are present at the start of the 
# integration interval

import numpy as np
import os
import pp_dursim
import volat_calc
import scipy.constants as si
from water_calc import water_calc

# define function
def cham_up(sumt, Pnow, 
	light_time_cnt, tnew, 
	new_part_sum1, update_count,
	injectt, gasinj_cnt, inj_indx, 
	Ct, seedt_cnt, num_comp, y0, y, 
	N_perbin0, corei, lowsize, 
	uppsize, num_sb, MV, radn, std, 
	H2Oi, rbou, infx_cnt, Cfactor, diff_vol, 
	DStar_org, tempt_cnt, RHt_cnt, nuci, y_mw, 
	temp_now, gpp_stab, t00, x, pcontf, Cinfl_now, surfT, 
	act_coeff, tot_in_res, Compti, self, vol_Comp, volP, ic_red):
	
	# inputs: ------------------------------------------------
	# sumt - cumulative time through simulation (s)
	# self.TEMP - temperature in chamber (K)
	# self.tempt - times that temperatures reached (s)
	# Pnow - pressure in chamber (Pa)
	# self.light_stat - status of lights
	# self.light_time - times that light attain status (s)
	# light_time_cnt - light status counter
	# self.light_ad - marker for whether to change time interval 
	#	in response to changing natural light intensity
	# tnew - time interval between chamber updates (s)
	# self.nuc_ad - flag for whether user wants time step adapted 
	# to nucleation
	# self.nucv1 - nucleation parameter one
	# self.nucv2 - nucleation parameter two
	# self.nucv3 - nucleation parameter three
	# new_part_sum1 - total number concentration of new 
	#	particles so far (#/cm3 (air))
	# self.update_stp - time interval between operator-split 
	#	updates (s)
	# update_count - count since operator-split last 
	#	updated (s)
	# self.lat - latitude (degrees)
	# self.lon - longitude (degrees)
	# self.dayOfYear - number of days since 31st December
	# self.photo_path - photochemistry parameter file
	# self.af_path - actinic flux file
	# injectt - time of instantaneous injections of 
	#	components (s)
	# gasinj_cnt - count on injection times of component(s)
	# inj_indx - index of components being instantaneously 
	# 	injected after experiment start
	# Ct - concentration(s) (ppb) of component(s) injected 
	#	instantaneously after experiment start
	# self.pmode - whether particle number size distributions stated 
	#	explicitly or by mode
	# self.pconc - concentration of injected particles 
	#	(# particles/cm3 (air))
	# self.pconct - times of particle injection (s)
	# seedt_cnt - count on injections of particles
	# num_comp - number of components
	# y0 - concentration of components prior to integration 
	# 	(# molecules/cm3 (air))
	# y - variable concentration of components prior to 
	#	integration (# molecules/cm3 (air))
	# N_perbin0 - concentration of particles 
	#	(# particles/cm3 (air)) at start of time interval
	# self.mean_rad - mean radius for particle number size 
	#	distribution (um)
	# corei - index of core component
	# self.seedx - mole ratio of non-water components 
	# comprising seed particles
	# self.seed_name - name(s) of component(s) comprising seed 
	#	particles
	# lowsize - lower size bin boundary (um)
	# uppsize - upper size bin boundary (um)
	# num_sb - number of size bins (including wall if turned on)
	# MV - molar volume of components (cm3/mol)
	# radn - current radius at size bin centres (um)
	# std - standard deviation for injected particle number size 
	#	distributions
	# self.y_dens - component densities (kg/m3)
	# H2Oi - index of water
	# rbou - size bin radius bounds (um)
	# self.con_infl_t - times for constant influxes (s)
	# infx_cnt - count on constant influx occurrences
	# self.Cinfl - influx rate for components with constant influx 
	#	(ppb/s)
	# self.wall_on - marker for whether wall is on
	# Cfactor - conversion factor from ppb to molecules/cm3 (air)
	# self.seedi - index of seed component(s)
	# diff_vol - diffusion volumes of components according to 
	#	Fuller et al. (1969)
	# DStar_org - gas-phase diffusion coefficients of components 
	#	(cm2/s)
	# self.RH - relative humidities (fraction 0-1)
	# self.RHt - times through experiment at which relative 
	# humidities reached (s)
	# tempt_cnt - count on temperatures
	# RHt_cnt - relative humidity counts
	# self.Pybel_objects - the pybel identifiers for components
	# nuci - index of nucleating component
	# y_mw - molar weight of components (g/mol)
	# temp_now - chamber temperature (K) prior to this update
	# self.Psat - saturation vapour pressures of components at the 
	#	current chamber temperature (# molecules/cm3)
	# gpp_stab - flag for whether to linearly interpolate any 
	#	change to chamber conditions (equals -1 if change 
	#	needed)
	# t00 - the initial integration step on the current 
	# 	integration step (s)
	# x - starting sizes of particles per size bin (um)
	# self.pcont - flags for whether particle injection instantaneous 
	# 	or continuous
	# pcontf - whether current state of particle injection is 
	#	continuous
	# Cinfl_now - influx rate of components with continuous influx 
	#	(ppb/s)
	# surfT - surface tension of particles (g/s2 == mN/m == dyn/cm)
	# act_coeff - activity coefficient of components
	# self.seed_eq_wat - whether seed particles to be equilibrated 
	#	with water prior to ODE solver
	# self.Vwat_inc - whether suppled seed particle volume 
	#	contains equilibrated water
	# tot_in_res - count on total injected concentration of 
	#	injected components (ug/m3)
	# Compti - index for total injection record for 
	#	instantaneously injected components
	# self.cont_inf_reci - index for total injection record for 
	#	continuously injected components
	# self.con_infl_indx - index for continuously injected 
	#	components from all components
	# ic_red - whether time step reduced
	# --------------------------------------------------------------

	# number of particle phase size bins
	num_asb = (num_sb-self.wall_on)

	# ensure N_perbin has a value
	N_perbin = np.zeros((N_perbin0.shape[0], N_perbin0.shape[1])) 
	N_perbin[:] = N_perbin0[:] # particle number concentration (# particles/cm3)

	# check on dilution factor setting --------------------------
	self.dil_fac_cnt =  sum(self.dil_fact <= (sumt))-1

	# update dilution factor
	self.dil_fac_now = self.dil_fac[self.dil_fac_cnt]
	# for water vapour
	# if water previously set to 0, e.g. below because of
	# fixing to observation, then keep at 0.
	if hasattr(self, 'dil_fac_H2O_now'):
		if (self.dil_fac_H2O_now == 0.):
			self.dil_fac_H2O_now = 0
		else:
			self.dil_fac_H2O_now = self.dil_fac[self.dil_fac_cnt]
	else:
		self.dil_fac_H2O_now = self.dil_fac[self.dil_fac_cnt]
		
	# check on change of light setting ------------------------------

	# begin by assuming no change to time interval required due to chamber 
	# condition/nucleation
	bc_red = 0
	
	if ((len(self.light_time)) > 0):
	
		# whether lights on (>1) or off (0) during this step
		self.light_stat_now = self.light_stat[int(sum(self.light_time<=sumt)-1)]
		
		# check whether changes occur at start of this time step
		if (sumt == self.light_time[light_time_cnt] and light_time_cnt>-1):
			
			if (light_time_cnt<(len(self.light_stat)-1)):
				light_time_cnt += 1 # keep count of light setting index
			else:
				light_time_cnt = -1 # reached end
			# reset flag for time step reduction due to chamber condition
			bc_red = 0
				
		
		# check whether light on/off changes during proposed integration time step
		if (sumt+tnew > self.light_time[light_time_cnt] and light_time_cnt != -1):
			# if yes, then reset integration time step so that next step coincides 
			# with change
			tnew = self.light_time[light_time_cnt]-sumt
			bc_red = 1 # flag for time step reduction due to boundary conditions
			
		# if reached final status of lights, then keep this status
		if (light_time_cnt == -1):
			self.light_stat_now = self.light_stat[light_time_cnt]
 
	# if lights are on during this step and lighting is natural, then check whether
	# proposed time step needs reducing to limit change to light intensity, if this
	# time interval adaption is requested
	# if using natural light
	cwd = os.getcwd() # address of current working directory

	if (self.light_stat_now >= 1 and self.photo_path == str(cwd + '/PyCHAM/photofiles/MCMv3.2') and self.af_path == 'no' and self.light_ad == 1):	
		# check time step required to limit change to rate of 
		# MCM photochemical equation number 6, 
		# which the unit test for
		# the zenith module shows to be most photosensitive (compared to
		# other photochemical equations)
		import zenith
		# photochemical rate now
		self.sumt = sumt
		(secxn, cosxn) = zenith.zenith(self)
		Jn =1.747e-1*cosxn**(0.155)*np.exp(-1.*0.125*secxn)
		
		# photochemical rate after proposed time step
		self.sumt += tnew # temporary changed (reversed below)
		(secxt, cosxt) = zenith.zenith(self)
		self.sumt -= tnew # reverse to temporary change above
		
		Jt =1.747e-1*cosxn**(0.155)*np.exp(-1.*0.125*secxn)
		
		# iteratively reduce proposed time interval until photochemical
		# rate changes by acceptable amount
		while (abs(Jt-Jn) > 5.e-3):
			tnew = tnew*0.9
			self.sumt += tnew # temporary changed (reversed below)
			# photochemical rate after proposed time step
			(secxt, cosxt) = zenith.zenith(self)
			self.sumt -= tnew # reverse to temporary change above
			Jt = 1.747e-1*cosxt**(0.155)*np.exp(-1.*0.125*secxt)
			bc_red = 1
			
	# check on updates to temperature (K) --------------------------------------	
	if (len(self.TEMP) > 1): # because a temperature must be given for experiment start
	
		# check whether changes occur at start of this time step
		if (sumt >= self.tempt[tempt_cnt] and tempt_cnt != -1):

			# new temperature (K)
			if (gpp_stab != -1): # if no linear interpolation required
			
				temp_nown = self.TEMP[tempt_cnt] # new temperature (K)
				if (tempt_cnt < (len(self.tempt)-1)):
					# keep count of temperature setting index
					tempt_cnt += 1
				else:
					tempt_cnt = -1 # reached end
				# reset flag for time step 
				# reduction due to boundary conditions
				bc_red = 0
			else:
				# new temperature (K)
				temp_nown = np.interp(tnew, [0, t00], 
				[temp_now, self.TEMP[tempt_cnt]])
				# reset flag for time step reduction 
				# due to boundary conditions
				bc_red = 1
			
			# update vapour pressure of water (log10(atm)),
			# but don't update gas-phase concentration of water, since
			# RH should be allowed to vary with temperature
			[_, Psat_water, _] = water_calc(temp_nown, self.RH[RHt_cnt], si.N_A)

			# update vapour pressures of all components (# molecules/cm3 and Pa), 
			# ignore density output
			[self, _] = volat_calc.volat_calc(0, temp_nown,
				 H2Oi, num_comp, Psat_water, vol_Comp, 
				volP, 0, corei, 
				0, 0.0, [], 1, nuci, self)
			
			# according to the ideal gas law, air pressure (Pa) inside chamber
			# is proportional to temperature, therefore pressure changes by 
			# the same factor 
			Pnow = Pnow*(temp_nown/temp_now)
			
			# update ppb to molecules/cm3 conversion factor concentrations
			# total number of molecules in 1 cc air using ideal gas law.  R has units cc.Pa/K.mol
			ntot = Pnow*(si.N_A/((si.R*1.e6)*temp_nown))
			# one billionth of number of molecules in chamber unit volume
			Cfactor = ntot*1.e-9 # ppb to molecules/cc conversion factor
			
			# dynamic viscosity of air (kg/m.s), eq. 4.54 of Jacobson 2005
			dyn_visc = 1.8325e-5*((416.16/(temp_nown+120.))*(temp_nown/296.16)**1.5)
	
			ma = 28.966e-3 # molecular weight of air (kg/mol) (Eq. 16.17 Jacobson 2005)
			
			# air density (kg/m3 (air)), ideal gas law
			rho_a =  (Pnow*ma)/((si.R)*temp_nown)
						
			# update mean free path and thermal speed
			# mean thermal speed of each molecule (m/s) (11.151 Jacobson 2005)
			# note that we need the weight of one molecule, which is why y_mw is divided by
			# Avogadro's constant, and we need it in kg, which is why we multiply by 1e-3
			therm_sp = ((8.*si.k*temp_nown)/(np.pi*(y_mw/si.N_A)*1.e-3))**0.5
			
			# mean free path (m) for each component (15.24 of Jacobson 2005)
			# molecular weight of air (28.966 g/mol taken from table 16.1 Jacobson 2005)
			mfp = (2.*dyn_visc/(rho_a*therm_sp)).reshape(-1, 1)
			
			# diffusion coefficient (m2/s) of components in gas phase (air), eq 4.1.4 of
			# the Taylor (1993) textbook 
			# Multicomponent Mass Transfer, ISBN: 0-471-57417-1, note diffusion 
			# volume for air (19.7) taken from Table 4.1 of Taylor (1993) and mw of 
			# air converted to g/mol from kg/mol.  This is a replication of the original method 			
			# from Fuller et al. (1969): doi.org/10.1021/j100845a020
			DStar_org = 1.013e-2*temp_nown**1.75*(((y_mw+ma*1.e3)/(y_mw*ma*1.e3))**0.5)/(Pnow*(diff_vol**(1./3.)+19.7**(1./3.))**2.)
			# convert to cm2/s
			DStar_org = DStar_org*1.e4
			
			temp_now = temp_nown # update current temperature (K)
		
		# check whether temperature changes during proposed integration time step
		if (sumt+tnew > self.tempt[tempt_cnt] and tempt_cnt != -1 and gpp_stab != -1):
			# if yes, then reset integration time step so that next step coincides 
			# with change
			tnew = self.tempt[tempt_cnt]-sumt
			bc_red = 1 # flag for time step reduction due to boundary conditions
		
	if (len(self.TEMP) == 1):
		temp_now = self.TEMP[0] # temperature constant if only one value given
	
	# check on instantaneous injection of components ---------------------------------------
	if (len(injectt) > 0 and gasinj_cnt > -1): # if any injections occur
	
		# check whether changes occur at start of this time step
		if (sumt >= injectt[gasinj_cnt] and gasinj_cnt != -1):
		
			if (gpp_stab != -1): # if no linear interpolation required
				Ct_gain = Ct[:, gasinj_cnt]
				
				if (gasinj_cnt < (Ct.shape[1]-1)):
					gasinj_cnt += 1 # update count on injections
				else:
					gasinj_cnt = -1 # reached end
				bc_red = 0 # reset flag for time step reduction due to boundary conditions
				
			else:
				# loop through components with instantaneous injection
				inj_cntn = 0 # keep count on components
				Ct_gain = np.zeros((len(inj_indx))) # empty results
				for inj_indxi in inj_indx:
					Ct_gain[inj_cntn] = np. interp(tnew, [0, t00], [y0[inj_indxi]/Cfactor, Ct[inj_cntn, gasinj_cnt]])
					inj_cntn += 1  # keep count on components
				bc_red = 1 # reset flag for time step reduction due to boundary conditions

			# record additional injection of components (ug/m3)
			tot_in_res[Compti] += (((Ct_gain*Cfactor-y[inj_indx])/si.N_A)*(y_mw[inj_indx].squeeze()))*1.e12

			# keep just the non-zero Ct gains, 
			# since injection to create a 0 abundance
			# is impossible
			nz_indx = (Ct_gain != 0.)

			# account for change in gas-phase concentration,
			# convert from ppb to molecules/cm3 (air)
			y[inj_indx[nz_indx]] = Ct_gain[nz_indx]*Cfactor
			
		# check whether changes occur during proposed integration time step
		# and that time step has not been forced to reduce due to unstable ode solver
		if (sumt+tnew > injectt[gasinj_cnt] and gasinj_cnt != -1 and gpp_stab != -1):
			# if yes, then reset integration time step so that next step coincides 
			# with change
			tnew = injectt[gasinj_cnt]-sumt
			bc_red = 1 # flag for time step reduction due to boundary conditions
	
	# check on instantaneous change in relative humidity -----------
	# if any injections occur
	if (len(self.RHt) > 0 and RHt_cnt > -1):
		
		# check whether changes occur at start of this time step
		if (sumt >= self.RHt[RHt_cnt] and RHt_cnt != -1):
			
			# if no linear interpolation required
			if (gpp_stab != -1):
				
				RHn = self.RH[RHt_cnt]
				
				if (RHt_cnt < (self.RHt.shape[0]-1)):
					
					if (self.RH[RHt_cnt] == self.RH[RHt_cnt+1]):
						
						# assume we don't want dilution of 
						# water vapour
						self.dil_fac_H2O_now = 0.
					RHt_cnt += 1 # update count on RH
					
				else:
					RHt_cnt = -1 # reached end
					# assume no dilution of water
					# vapour
					self.dil_fac_H2O_now = 0.
 
				# reset flag for time step reduction due 
				# to boundary conditions
				bc_red = 0				
			else:
				RHn = np. interp(tnew, [0, t00], [self.RH[RHt_cnt-1], 
					self.RH[RHt_cnt]])
				
				if (self.RH[RHt_cnt-1] == self.RH[RHt_cnt]):
					
					# assume we don't want dilution of 
					# water vapour
					self.dil_fac_H2O_now = 0.
				# reset flag for time step reduction due to boundary conditions
				bc_red = 1
		

			# update vapour pressure of water (log10(atm)), and change 
			# gas-phase concentration of water vapour since 
			# RH stays as stated in the RH and RHt model variables
			[y[H2Oi], _, _] = water_calc(temp_now, RHn, si.N_A)
				
		# check whether changes occur during next proposed integration time step
		# and that time step has not been forced to reduce due to unstable ode solvers
		if ((sumt+tnew > self.RHt[RHt_cnt]) and (RHt_cnt != -1) and gpp_stab != -1):
			# if yes, then reset integration time step so that next step coincides 
			# with change
			tnew = self.RHt[RHt_cnt]-sumt
			# flag for time step reduction 
			# due to boundary conditions
			bc_red = 1
	
	# get whether next/current injection of seed is instantaneous or continuous
	pcontf = self.pcont[0, seedt_cnt]
	
	# check on injection of particles --------------------------------------
	# filler for fraction of new seed particles injected so far
	pconcn_frac = 0.
	
	# if we are at start of simulation and there is 
	# continuous injection of particles from start, 
	# then set flag to solve continuous injection 
	# of particles. Note that pcontf set in 
	# ode_updater
	if (sumt == 0. and pcontf == 1):
		self.pcont_ongoing = 1
	
	# if influx occurs and we are not on the cham_up 
	# call from the rec module (which has tnew=0)
	# note, if this particle influx section called 
	# during the call to rec, then continuous influx is 
	# cancelled when called later
	if ((sum(self.pconct[0, :]) > 0) and (seedt_cnt > -1) and 
		(num_sb-self.wall_on > 0) and tnew > 0.):
		
		# in case particle influx repeated every 24 hours
		if (self.pconctf == 1 and sumt >= 24.*3.6e3):
			pinsumt = (sumt % (24.*3.6e3))
			if (sumt >= (24.*3.6e3) and 
			(sumt % (24.*3.6e3)) == 0):
				pinsumt = 0
			# if seedt_cnt already been reset to zero and time 
			# through simulation not yet caught up with reset 
			# time, then prevent influx at reset time 
			# occurring now
		else:
			pinsumt = sumt

		# if repeating every 24 hours
		if (seedt_cnt == 0 and self.pconctf == 1):
			# if on first 24 hours
			if (sumt > 0. and sumt < (24.*3.6e3)):
				pinsumt = -1
			# if already past 24 hours
			if (sumt % (24.*3.6e3) > self.pconct[0, seedt_cnt]):
				pinsumt = -1
		
		# check whether changes occur at start of this time step
		if (pinsumt >= self.pconct[0, seedt_cnt]):
			
			# get the current injection concentration
			# particles/cm3				
			pconcn = self.pconc[:, seedt_cnt]

			if (self.pmode == 0): # if in modal mode
				mean_radn = self.mean_rad[:, 
				seedt_cnt]
				stdn = std[:, seedt_cnt]
			# if number concentrations given per size bin
			else:
				mean_radn = self.mean_rad[:, seedt_cnt]
				stdn = std
			
			# if linear interpolation required and 
			# instantaneous injection of seed
			if (gpp_stab == -1 and self.pcont[0, seedt_cnt] == 0):
				# empty results array
				pconcn = np.zeros((self.pconc.shape[0])) 
				# loop through size bins for interpolation 
				# since interpolation 
				# is one dimensional
				for i in range(num_sb-self.wall_on):
					pconcn[i] = np.interp(tnew, [0, t00],
						[self.pconc[i, seedt_cnt-1], 
						self.pconc[i, seedt_cnt]])
				# remember the fraction of the number 
				# concentration added so far
				pconcn_frac = pconcn/self.pconc[:, seedt_cnt]
				# reset flag for time step reduction due to
				# boundary conditions
				bc_red = 1


			# account for instantaneous concentration in 
			# seed particles (continuous change dealt with below)
			if (pcontf == 0):
				
				# prepare to hold particle-phase concentrations (molecules/cm3)
				yp0 = np.zeros((num_comp*num_asb))
				yp0[:] = y0[num_comp:num_comp*(num_sb-
				self.wall_on+1)]

				# index for mole fraction of seed particles
				self.seedx_tcnt = seedt_cnt

				[y[num_comp:num_comp*(num_sb-self.wall_on+1)], 
				N_perbin, _,_] = pp_dursim.pp_dursim(
				yp0, N_perbin0, mean_radn, pconcn, lowsize, 
				uppsize, num_comp, (num_sb-self.wall_on), MV, 
				stdn, H2Oi, rbou, y_mw, surfT, 
				self.TEMP[tempt_cnt],
				act_coeff, pcontf, y[H2Oi], x, self)
					
				# turn off flag for ongoing injection of particles
				self.pcont_ongoing = 0
			
			# account for new continuous change in seed particles
			if (pcontf == 1):
					
				# turn on flag for ongoing injection of particles
				self.pcont_ongoing = 1
			
			# move count on particle injections up by one
			if (seedt_cnt < (self.pconct.shape[1]-1)):
				
				# if repeating every 24 hours
				if (self.pconctf == 1):
					if (self.pconct[0, seedt_cnt+1] >= 24.*3.6e3):
						seedt_cnt = 0 # reset if necessary
						self.seedx_tcnt = 0
					else:
						seedt_cnt += 1
						self.seedx_tcnt += 1

				if (self.pconctf == 0):
					seedt_cnt += 1
					self.seedx_tcnt += 1
					
			else:
				# if influx times treated explicitly
				if (self.pconctf == 0):
					seedt_cnt = -1 # reached end
					self.seedx_tcnt = -1
				# if repeating every 24 hours
				if (self.pconctf == 1):
					seedt_cnt = 0 # reset
					self.seedx_tcnt = 0
	
		# check whether changes occur during proposed 
		# integration time step
		# and that time step has not been forced to 
		# reduce due to unstable ode solvers
		if (pinsumt+tnew > self.pconct[0, seedt_cnt] 
		and seedt_cnt!=-1 and gpp_stab != -1): 

			# in case influxes repeated every 24 hours 
			# and seedt_cnt has been reset
			if (self.pconctf == 1):
				# in case seedt_cnt has been reset and 
				# new time interval is going into the 
				# next day	
				if (seedt_cnt == 0 and (pinsumt+tnew) >=
				 (24.*3.6e3)):
					if ((sumt+tnew) % (24.*3.6e3) > 
					self.pconct[0, seedt_cnt]):
						tnew -= ((sumt+tnew) % 
						(24.*3.6e3)-
						self.pconct[
						0, seedt_cnt])
						# time step change flag
						bc_red = 1
				# in case seedt_cnt hasn't been reset
				if (seedt_cnt != 0):
					if ((sumt+tnew) % (24.*3.6e3) >
					 self.pconct[0, seedt_cnt]):
						tnew = self.pconct[
						0, seedt_cnt]-(
						sumt % (24.*3.6e3))
			
			# if influxes treated explicitly
			else:
				# if yes, then reset integration 
				# time step 
				# so that next step coincides 
				# with change
				tnew = self.pconct[0, seedt_cnt]-sumt
				# flag for time step reduction 
				# due to boundary conditions
				bc_red = 1
	
	# account for ongoing continuous influx of particles
	# account for new continuous change in seed particles
	if (self.pcont_ongoing == 1):
		
		if (seedt_cnt != -1):
			
			# index for mole fraction of seed particles
			self.seedx_tcnt = seedt_cnt-1
			
			if (self.pmode == 0): # if in modal mode
				mean_radn = self.mean_rad[:, seedt_cnt-1]
				stdn = std[:, seedt_cnt-1]
			
			# if number concentration per size bin 
			# given explicitly
			if (self.pmode == 1): 				
				mean_radn = self.mean_rad[:, seedt_cnt-1]
				stdn = std		
		
			# injected seed particle number concentration 
			# integrated over proposed 
			# time step (# particles/cm3)
			pconcn = self.pconc[:, seedt_cnt-1]*tnew
		
		if (seedt_cnt == -1):
			pconcn = self.pconc[:, -1]
			mean_radn = self.mean_rad[:, -1]
			stdn = std[:, -1]		
		
			# injected seed particle number concentration 
			# integrated over proposed 
			# time step (# particles/cm3)
			pconcn = self.pconc[:, -1]*tnew
		
		# number of actual particle size bins
		num_asb = (num_sb-self.wall_on)

		# prepare to hold particle-phase concentrations (molecules/cm3)
		yp0 = np.zeros((num_comp*num_asb))
		yp0[:] = y0[num_comp:num_comp*(num_asb+1)]

		[y[num_comp:num_comp*(num_sb-self.wall_on+1)], N_perbin, _, 
		_] = pp_dursim.pp_dursim(yp0, 
		N_perbin0, mean_radn, pconcn, lowsize, 
		uppsize, num_comp, num_asb, MV, 
		stdn, H2Oi, rbou, y_mw, surfT, self.TEMP[tempt_cnt], act_coeff, 
		pcontf, y[H2Oi], x, self)
	
	# --------------------------------------------------------------------------
	
	# check on continuous influx of gas-phase components -------------
	if (len(self.con_infl_t) > 0): # if influx occurs

		# set relevant time if loop involved
		if (self.con_infl_tf == 0):
			ci_sumt = sumt
		if (self.con_infl_tf == 1):
			if (sumt >= 24.*3.6e3):
				ci_sumt = 0. + sumt % (24.*3.6e3)
			else:
				ci_sumt = sumt
		
		# in case influxes begin after simulation start 
		# create a zero array of correct shape
		# note that final condition (infx_cnt==0) means 
		# that this only activated if we're
		# really at the first supplied influx point 
		# (because cham_up in rec_prep could have)
		# moved infx_cnt up by 1
		if (sumt == 0. and self.con_infl_t[infx_cnt] != 0. and infx_cnt == 0):
			Cinfl_now = np.zeros((self.con_infl_C.shape[0], 1))
		
		# if the final input for influxes reached
		if (infx_cnt == -1):

			# influx of components now, convert from ppb/s 
			# to # molecules/cm3/s (air)
			if ('ppb' in self.abun_unit):
				Cinfl_now = np.float64(self.con_infl_C[:, infx_cnt]*self.Cfactor)
			if ('mol' in self.abun_unit):
				Cinfl_now = np.float64(self.con_infl_C[:, infx_cnt])

			# ensure correct shape
			Cinfl_now = np.array((Cinfl_now)).reshape(-1, 1)

			if (self.H2Oin == 1):
				if ('ppb' in self.abun_unit):
					# continuous influx rate of water now
					self.Cinfl_H2O_now = (self.con_infl_H2O[:, infx_cnt]*self.Cfactor)
				if ('mol' in self.abun_unit):
					self.Cinfl_H2O_now = self.con_infl_H2O[:, infx_cnt]

				# ensure correct shape
				self.Cinfl_H2O_now = np.float64(self.Cinfl_H2O_now.reshape(-1, 1))
			
			# record cumulative injection of components (ug/m3)
			tot_in_res[self.cont_inf_reci] += (((((Cinfl_now.squeeze())*(tnew))/si.N_A)*(y_mw[self.con_infl_indx].squeeze()))*1.e12).reshape(-1)
			
		# check whether changes occur at start of this time step
		if (ci_sumt == self.con_infl_t[infx_cnt] and (infx_cnt != -1)):
			
			# we need units of # molecules/cm3/s
			if ('ppb' in self.abun_unit):
				# influx of components now, convert from ppb/s 
				# to # molecules/cm3/s (air)
				Cinfl_now = np.float64(self.con_infl_C[:, infx_cnt]*Cfactor)

			if ('mol' in self.abun_unit):
				Cinfl_now = np.float64(self.con_infl_C[:, infx_cnt])
			
			# ensure correct shape
			Cinfl_now = np.array((Cinfl_now)).reshape(-1, 1)

			if (self.H2Oin == 1):
				if ('ppb' in self.abun_unit):
					# continuous influx rate of water now
					self.Cinfl_H2O_now = np.float64(
					self.con_infl_H2O[:, infx_cnt]*self.Cfactor)
				if ('mol' in self.abun_unit):
					self.Cinfl_H2O_now = np.float64(
					self.con_infl_H2O[:, infx_cnt])

				# ensure correct shape
				self.Cinfl_H2O_now = self.Cinfl_H2O_now.reshape(-1, 1)
			
			# record cumulative injection of components (ug/m3)
			tot_in_res[self.cont_inf_reci] += ((((((
			Cinfl_now.squeeze())*(tnew))/si.N_A)*(
			y_mw[self.con_infl_indx].squeeze()))*1.e12).reshape(-1))

			# update index counter for continuous influxes - 
			# used in integrator below
			if (infx_cnt < (self.con_infl_C.shape[1]-1)):
				infx_cnt += 1
			else:
				if (self.con_infl_tf == 1):
					infx_cnt = 0
				if (self.con_infl_tf == 0):
					infx_cnt = -1 # reached end

			# reset flag for time step reduction due to boundary conditions
			bc_red = 0
		
		# check whether changes occur during proposed 
		# integration time step
		if (ci_sumt+tnew > self.con_infl_t[infx_cnt] and (infx_cnt != -1)):
			if (self.con_infl_tf == 1 and infx_cnt == 0):
				if ((ci_sumt+tnew) % (24.*3.6e3) > self.con_infl_t[infx_cnt]):
					tnew = (24.*3.6e3+self.con_infl_t[infx_cnt])-ci_sumt	
				# flag for time step reduction to accommodate
				# changing boundary conditions
				bc_red = 1
			else:
				# if yes, then reset integration time step 
				# so that next step coincides 
				# with change
				tnew = self.con_infl_t[infx_cnt]-ci_sumt
			
				# flag for time step reduction due to 
				# boundary conditions
				bc_red = 1		
	else: # if no continuous influxes, provide filler
		Cinfl_now = np.zeros((1, 1))
	
	# check on nucleation ------------------------------------------
	# if automatic time step adaption to nucleation requested, 
	# check whether number of new particles
	# exceeds 10 % of total number formed during nucleation event.  	# Second part of condition is that
	# the specified nucleation event has not yet reached its defined 	# finishing particle number
	# concentration (# particles/cm3 (air))
	if ((self.nuc_ad == 1) and (new_part_sum1 < self.nucv1*0.9) 
		and ((num_sb-self.wall_on) > 0)):
	
		# the time step (s) needed to increase number 
		# concentration of nucleated particles by 10 %
		t_need = (0.1*self.nucv1+new_part_sum1)
		t_need = np.log(t_need/self.nucv1)
		t_need = np.log(t_need/self.nucv2)
		t_need = t_need*self.nucv3*-1.-sumt

		# if suggested time step exceeds this, 
		# then reduce to required time step	
		if (tnew > t_need):  
			tnew = t_need
			self.update_stp = t_need
			update_count = 0.
			bc_red = 1
			
	# nucleation check end -----------------------------------------
	
	# check on new vapour pressure of HOM-RO2+MCM-RO2 accretion 
	# products -------------
	
	# convert y into components in rows and phases in columns
	if ('RO2_POOL' in self.comp_namelist):
	
		y_mat = y.reshape(num_comp, num_sb+1, order='F')
		if self.wall_on > 0:
			y_mat = y_mat[:, 0:-self.wall_on]
	
		# get the average oxygen and carbon number of the gas- and particle-phase RO2
		Cnumav = sum(sum((self.Cnum[self.RO2_indices[:, 1], :])*(y_mat[self.RO2_indices[:, 1], :].sum(axis=1))/sum(sum(y_mat[self.RO2_indices[:, 1], :]))))

		Onumav = sum(sum((self.Onum[self.RO2_indices[:, 1], :])*(y_mat[self.RO2_indices[:, 1], :].sum(axis=1))/sum(sum(y_mat[self.RO2_indices[:, 1], :]))))
	
		# estimate vapour pressure (Pa) effect of the RO2 pool based on carbon and oxygen number
		RO2pool_effect_Pa = 10**(-0.12*Onumav + Cnumav*-0.22)*101325.
		
		# take effect on the HOM-RO2-MCM-RO2 accretion product, note that inside Psat_Pa_rec
		# is the estimated vapour pressure of the HOM-RO2 (Pa)
		self.Psat_Pa[:, self.RO2_POOL_APi] = self.Psat_Pa_rec[self.RO2_POOL_APi] + RO2pool_effect_Pa
		# convert to # molecules/cm3 (air) using ideal
		# gas law, R has units cm3.Pa/K.mol
		self.Psat[:, self.RO2_POOL_APi] = self.Psat_Pa[0, self.RO2_POOL_APi]*(si.N_A/((si.R*1.e6)*self.TEMP[tempt_cnt]))
		
		
	# end of check on new vapour pressure of HOM-RO2+MCM-RO2 accretion products ----------

	# in case time step reduction was actioned prior to call to 
	# cham_up, then keep the signal that this reduction occurred
	if (ic_red == 1 and bc_red == 0):
		bc_red = 1

	return(temp_now, Pnow, light_time_cnt, tnew, bc_red, update_count, 
		Cinfl_now, seedt_cnt, Cfactor, infx_cnt, gasinj_cnt, DStar_org, y, tempt_cnt, 
		RHt_cnt, N_perbin, x, pconcn_frac, pcontf, tot_in_res, self)
