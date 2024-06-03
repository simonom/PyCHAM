########################################################################							
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
'''module to prepare the recording matrices'''
# recording matrices are prepared and initial conditions stored

import numpy as np
import partit_var
import rrc_calc
import cham_up
import scipy.constants as si
import importlib
import act_coeff_update

# define function
def rec_prep(nrec_step, y, y0, 
	num_sb, num_comp, N_perbin, core_diss, mfp,
	accom_coeff, y_mw, surfT, R_gas, NA, 
	x, therm_sp, H2Oi, act_coeff, 
	sumt, Pnow, light_time_cnt, 
	Jlen, Cfactor, Vbou, tnew, 
	np_sum, update_count, injectt, gasinj_cnt, 
	inj_indx, Ct, seedt_cnt, corei, 
	lowsize, uppsize, radn, std, rbou, 
	infx_cnt, MV, diff_vol, 
	DStar_org, tempt_cnt, RHt_cnt, 
	nuci, t0, pcontf, NOi, HO2i, NO3i, z_prt_coeff,
	tot_in_res, Compti, 
	tot_in_res_indx, chamSA, chamV, wat_hist, self, vol_Comp, volP):
	
	# inputs: ------------------------------------------------------
	# nrec_step - number of steps to record on
	# y - initial concentrations (# molecules/cm3 (air))
	# y0 - component concentrations prior to integration 
	#	step (# molecules/cm3 (air))
	# self.rindx_g - indices of reactants
	# self.rstoi_g - stoichiometries of reactants
	# self.pindx_g - indices of products
	# self.pstoi_g - stoichiometries of products
	# self.nprod_g - number of products
	# self.nreac_g - number of reactants
	# num_sb - number of size bins
	# num_comp - number of components
	# N_perbin - particle concentration per size bin (#/cm3 (air))
	# core_diss - dissociation constant of seed
	# self.Psat - pure component saturation vapour pressure 
	#	(# molecules/cm3 (air))
	# mfp - mean free path (m)
	# accom_coeff - accommodation coefficient
	# y_mw - molecular weight (g/mol)
	# surfT - surface tension (g/s2)
	# R_gas - ideal gas constant (kg.m2.s-2.K-1.mol-1)
	# self.TEMP - temperature (K)
	# self.tempt - times temperatures achieved (s)
	# NA - Avogadro constants (molecules/mol)
	# self.y_dens - component densities (kg/m3)
	# x - particle radii (um)
	# therm_sp - thermal speed (m/s)
	# H2Oi - index for water
	# act_coeff - activity coefficient
	# self.RO2_indices - index of alkyl peroxy radicals in second column
	# sumt - cumulative time through simulation (s)
	# Pnow - chamber pressure (Pa)
	# self.light_stat - light status
	# self.light_time - times corresponding to light status
	# light_time_cnt - count on light status array
	# self.daytime - start time of simulation (s)
	# self.lat - latitude
	# self.lon - longitude
	# self.af_path - path to actinic flux file
	# self.dayOfYear - day of year
	# self.photo_path - path to photochemistry file
	# Jlen - length of photochemistry equations
	# self.kw - gas-wall partitioning coefficient (/s)
	# Cfactor - conversion factor (ppb/molecules/cc)
	# self.tf - transmission factor for natural sunlight
	# self.light_ad - marker for whether to adapt time interval to 
	#		changing natural light intensity
	# self.wall_on - marker for whether wall present
	# Vbou - volume boundary of particle size bins (um3)
	# tnew - the proposed integration time interval (s)
	# self.nuc_ad - marker for whether to adapt time step based on 
	#	nucleation
	# self.nucv1 - nucleation parameter one
	# self.nucv2 - nucleation parameter two
	# self.nucv3 - nucleation parameter three
	# np_sum - cumulative number concentration of new 
	#	particles so far (#/cm3 (air))
	# self.update_stp - time interval between operator-split processes
	# update_count - count on time since last operator-split (s)
	# injectt - time of instantaneous injection of components (s)
	# gasinj_cnt - count on injection times of components
	# inj_indx - index of component(s) being injected after 
	#	experiment start
	# Ct - concentration(s) (ppb) of component(s) injected 
	# instantaneously after experiment start
	# self.pmode - whether number size distributions expressed as modes or explicitly
	# self.pconc - concentration of injected particles (# particles/cm3 (air))
	# self.pconct - times of particle injection (s)
	# seedt_cnt - count on injection of seed particles
	# self.mean_rad - mean radius for particle number size 
	#	distribution (um)
	# corei - index of core component
	# self.seed_name - name(s) of component(s) comprising seed particles
	# self.seedx - mole ratio of components comprising seed particles
	# lowsize - lower size bin boundary (um)
	# uppsize - upper size bin boundary (um)
	# radn - current radius at size bin centres (um)	
	# std - standard deviation for injected particle number size 
	#	distributions
	# rbou - size bin radius bounds (um)
	# self.cont_infl_t - times for constant influxes (s)
	# infx_cnt - count on constant influx occurrences
	# self.con_infl_C - influx rates of components with constant influx 
	#	(ppb/s)
	# MV - molar volume (cc/mol)
	# self.partit_cutoff - the product of saturation vapour pressure and
	#		activity coefficient above which gas-particle
	#		partitioning assumed zero (Pa)
	# diff_vol - diffusion volume of components according to Fuller et al. (1969)
	# DStar_org - gas-phase diffusion coefficient of components (cm2/s)
	# self.seedi - index of seed component(s)
	# self.C_p2w - concentration of components on the wall due to 
	#	particle
	# deposition to wall (# molecules/cm3)
	# self.RH - relative humidities (fraction 0-1)
	# self.RHt - times through experiment at which relative 
	#	humidities reached (s)
	# tempt_cnt - count on chamber temperatures
	# RHt_cnt - chamber relative humidity counts
	# self.Pybel_objects - the pybel objects for components
	# nuci - index of nucleating component
	# self.nuc_comp - name of nucleating component
	# t0 - initial integration step (s)
	# self.pcont - flag for whether particle injection instantaneous or 
	#	continuous
	# pcontf - current status of particle injection (instantaneous 
	#	or continuous)
	# NOi - index of NO
	# HO2i - index of HO2
	# NO3i - index of NO3
	# z_prt_coeff - fraction of total gas-particle partitioning 
	#	coefficient below which partitioning to a particle 
	#	size bin is treated as zero,
	#	e.g. because surface area of that size bin is tiny
	# self.seed_eq_wat - whether seed particles to be equilibrated 
	#	with water prior to ODE solver
	# self.Vwat_inc - whether suppled seed particle volume contains equilibrated water
	# tot_in_res - record of total input of injected components (ug/m3)
	# Compti - index for total injection record for instantaneously injected components
	# self.tot_time - total experiment time (s)
	# self.cont_inf_reci - index of components with continuous influx in record
	# self.con_infl_indx - index of components with continuous influx in concentration array
	# tot_in_res_indx - index of components with recorded influx
	# chamSA - chamber surface area (m2)
	# chamV - chamber volume (m3)
	# tf_UVC - transmission factor for 254 nm wavelength light (0-1)
	# wat_hist - flag for water history with respect to particles
	# ----------------------------------------------------------------

	# update chamber variables
	
	# array to record cumulative influxes of influxing components
	self.tot_in_res_ft = np.zeros((nrec_step+1, len(tot_in_res)))
	# index of components being stored for influx
	self.tot_in_res_ft[0, :] = tot_in_res_indx
	self.tot_in_res_ft[1, :] = tot_in_res # influx at start
	
	# array to record carbon flux between reservoirs (influx to gas-phase 
	# due to emissions that exclude partitioning), net flux to/from 
	# particles, net flux to/from wall due to gas-wall partitioning,
	# loss to wall due to particle deposition on wall, net flux 
	# due to air exchange
	self.C_res_flux = np.zeros((nrec_step+1, 5))
	
	[temp_now, Pnow, light_time_cnt, tnew, ic_red, 
		update_count, Cinfl_now, seedt_cnt, Cfactor, infx_cnt, 
		gasinj_cnt, DStar_org, y, tempt_cnt, RHt_cnt, 
		N_perbin, x, pconcn_frac,  pcontf, tot_in_res, 
		self] = cham_up.cham_up(sumt, 
		Pnow, light_time_cnt, 0, 
		np_sum, update_count, 
		injectt, gasinj_cnt, inj_indx, Ct, 
		seedt_cnt, num_comp, y0, y, N_perbin, 
		corei, 
		lowsize, uppsize, num_sb, MV, radn, std, H2Oi, 
		rbou, infx_cnt, Cfactor, diff_vol, 
		DStar_org, tempt_cnt, RHt_cnt, nuci,
		y_mw, self.TEMP[0], 0, t0, x,  pcontf, 0., 
		surfT, act_coeff, tot_in_res, Compti, self, vol_Comp, 
		volP)
	
	# note that recording occurs after any instaneous changes-------
	# array to record time through simulation (s)
	trec = np.zeros((nrec_step))
	# array to record concentrations with time (# molecules/cm3)
	yrec = np.zeros((nrec_step, len(y)))
	# record initial concentrations (# molecules/cm3)
	yrec[0, :] = y[:] 
	# array to record conversion factor
	Cfactor_vst = np.zeros((nrec_step, 1))
	Cfactor_vst[0] = Cfactor

	# arrays to record particle radii and number size distributions	
	if ((num_sb-self.wall_on) > 0): # if particle size bins present
		x2 = np.zeros((nrec_step, (num_sb-self.wall_on)))
		Nres_dry = np.zeros((nrec_step, (num_sb-self.wall_on)))
		Nres_wet = np.zeros((nrec_step, (num_sb-self.wall_on)))
		rbou_rec = np.zeros((nrec_step, (num_sb-self.wall_on+1)))
		# concentration of components on the wall due to 
		# particle-wall loss, stacked by component first then by
		# size bin (# molecules/cM3)
		self.yrec_p2w = np.zeros((nrec_step, (num_sb-self.wall_on)*num_comp))
	else:
		x2 = 0.
		Nres_dry = 0.
		Nres_wet = 0.
		rbou_rec = 0.
		self.yrec_p2w = 0.
	
	if ((num_sb-self.wall_on) > 0 or self.wall_on == 1): # if particles or wall present
		
		# update partitioning variables
		[kimt, kelv_fac] = partit_var.kimt_calc(y, mfp, num_sb, 
		num_comp, accom_coeff, y_mw,   
		surfT, R_gas, temp_now, NA, N_perbin, 
		x.reshape(1, -1)*1.0e-6, therm_sp, H2Oi, act_coeff, 1, 
		Pnow, DStar_org, z_prt_coeff, chamSA, chamV, self)
		
	if (num_sb-self.wall_on) > 0: # if particles present
		# single particle radius (um) at size bin centre 
		# including contribution of water
		x2[0, :] = x

		# single particle size bin bounds including water (um)
		rbou_rec[0, :] = rbou
		
		# estimate particle number size distributions ----------------------------------
		Vnew = np.zeros((num_sb-self.wall_on))
		ish = N_perbin[:, 0] > 1.e-10
		
		if ish.sum()>0:
			# rearrange particle concentrations into size bins in rows, components in columns
			Cn = y[num_comp:num_comp*(num_sb-self.wall_on+1)].reshape(num_sb-self.wall_on, num_comp)
			# new volume of single particle per size bin (um3) excluding volume of water	
			Vnew[ish] = (np.sum((Cn[ish, :]/(si.N_A*N_perbin[ish]))*MV[:, 0]*1.e12, 1)-
					((Cn[ish, H2Oi]/(si.N_A*N_perbin[ish, 0]))*MV[H2Oi, 0]*1.e12))
			# loop through size bins to find number of particles in each 
			# (# particle/cm3 (air))
			for Ni in range(0, (num_sb-self.wall_on)):
				ish = (Vnew>=Vbou[Ni])*(Vnew<Vbou[Ni+1])
				Nres_dry[0, Ni] = N_perbin[ish, 0].sum()
			Nres_wet[0, :] = N_perbin[:, 0] # record with water
		
		# end of number size distribution part ----------------------------------------
		
		# concentration of components on the wall due to particle deposition to wall (# molecules/cm3)
		self.yrec_p2w[0, :] = self.C_p2w

	else: # fillers
		kimt = kelv_fac = 0.
	
	# update reaction rate coefficients
	[rrc, erf, err_mess] = rrc_calc.rrc_calc(y[H2Oi], temp_now, y, 
				Pnow, Jlen, y[NOi], y[HO2i], y[NO3i], 0., self)

	# chamber environmental conditions ----------------------------------
	# initiate the array for recording chamber temperature (K), pressure (Pa) 
	# and relative humidity (fraction (0-1))
	cham_env = np.zeros((nrec_step, 4))
	
	cham_env[0, 0] = temp_now # temperature (K)
	cham_env[0, 1] = Pnow # pressure (Pa)
	# relative humidity (fraction (0-1))
	cham_env[0, 2] = y[H2Oi]/self.Psat[0, H2Oi] 
	# transmission factor of light
	cham_env[0, 3] = self.tf[0]
	# --------------------------------------------------------------------------------

	# relative humidity now
	RH0 = cham_env[0, 2]

	# update particle-phase activity coefficients, note the output,
	# note that if ODE solver unstable, then y resets to y0 via
	# the cham_up module prior to this call
	[act_coeff, wat_hist, RHn, y, 
	dydt_erh_flag] = act_coeff_update.ac_up(y, H2Oi, RH0, temp_now, 
	wat_hist, act_coeff, num_comp, (num_sb-self.wall_on))

	# before solving ODEs for chemistry, gas-particle partitioning 
	# and gas-wall partitioning, 
	# estimate and record any change tendencies (# molecules/cm3/s)
	# resulting from these processes
	if (self.testf != 5 and len(self.dydt_vst) > 0):
		import dydt_rec
		importlib.reload(dydt_rec) # import most recent version
		dydt_cnt = 0 # index for row to record on
		self = dydt_rec.dydt_rec(y, rrc, dydt_cnt, num_sb, 
			num_comp, core_diss, kelv_fac, 
			kimt, act_coeff, dydt_erh_flag, H2Oi, wat_hist,
			self)
						
	return(trec, yrec, Cfactor_vst, Nres_dry, Nres_wet, x2, 
		seedt_cnt, rbou_rec, Cfactor, infx_cnt, temp_now, 
		cham_env, Pnow, cham_env[0, 2], Cinfl_now)
