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
'''updates integration constants and calls ODE solver'''
# constants for the integration of the relevant ODEs are updated
# at intervals and passed to the ODE solver to obtain solutions

import numpy as np
import math as math
import rrc_calc # function to update rate coefficients
import cham_up
import rec_prep
import partit_var
import rec
import jac_up
import mov_cen
import fullmov
import wallloss
import nuc
import coag
try:
	import ode_solv
except: # in case of a bad function
	import os
	if os.path.exists('ode_solv'): # remove any bad functions
		os.remove(ode_solv)
import ode_solv_wat
import dydt_rec
import importlib
import save
import time
import act_coeff_update
# providing error message if ODE solver produces 
# negative results below minimum integration time
import ode_brk_err_mess


def ode_updater(update_stp, 
	tot_time, save_stp, y, rindx, 
	pindx, rstoi, pstoi, nreac, nprod, jac_stoi, njac, 
	jac_den_indx, jac_indx, RO2_indx, RO_indx, H2Oi, temp, tempt, 
	Pnow, Jlen, con_infl_C, nrec_steps, 
	dydt_vst, siz_str, num_sb, num_comp, seed_name, seedx, 
	core_diss, Psat, mfp, therm_sp,
	accom_coeff, y_mw, surfT, R_gas, NA, y_dens, 
	x, Varr, act_coeff, Cw, kw, Cfactor, y_arr, 
	y_rind, uni_y_rind, y_pind, uni_y_pind, reac_col, prod_col, 
	rstoi_flat, pstoi_flat, rr_arr, rr_arr_p, rowvals, 
	colptrs, wall_on, jac_wall_indx, jac_part_indx, jac_extr_indx, Vbou,
	N_perbin, Vol0, rad0, np_sum, new_partr, nucv1, nucv2, 
	nucv3, nuci, nuc_comp, nuc_ad, RH, RHt, coag_on, inflectDp, pwl_xpre, 
	pwl_xpro, inflectk, chamR, McMurry_flag, p_char, e_field, 
	injectt, inj_indx, Ct, pmode, pconc, pconct, mean_rad, lowsize, 
	uppsize, std, rbou, const_infl_t, MV, rindx_aq, 
	pindx_aq, rstoi_aq, pstoi_aq, nreac_aq, nprod_aq, jac_stoi_aq, njac_aq, 
	jac_den_indx_aq, jac_indx_aq, y_arr_aq, y_rind_aq, 
	uni_y_rind_aq, y_pind_aq, uni_y_pind_aq, reac_col_aq, prod_col_aq, 
	rstoi_flat_aq, pstoi_flat_aq, rr_arr_aq, rr_arr_p_aq, eqn_num, 
	partit_cutoff, diff_vol, DStar_org, corei, ser_H2O, C_p2w, 
	sav_nam, comp_namelist, dydt_trak, space_mode, 
	rbou00, ub_rad_amp, indx_plot, comp0, rel_SMILES,
	Psat_Pa_rec, Psat_Pa, OC, wat_hist, Pybel_objects, pcont, dil_fac, NOi, 
	HO2i, NO3i, z_prt_coeff, seed_eq_wat, Vwat_inc, tot_in_res,
	Compti, cont_inf_reci, cont_inf_i, tot_in_res_indx, chamSA, 
	chamV, self):
	
	# inputs: ----------------------------------------------------
	# update_stp - interval at which to update integration 
	#		constants (s)
	# tot_time - total time to simulate (s)
	# save_stp - frequency to save results (s)
	# y - initial component concentrations (molecules/cc/s)
	# rindx - index of reactants per equation
	# pindx - index of products per equation
	# nreac - number of reactants per equation
	# nprod - number of products per equation
	# jac_stoi - stoichiometries relevant to Jacobian
	# njac - number of elements of Jacobian affected per equation
	# jac_den_indx - index of denominator components for Jacobian
	# jac_indx - index of Jacobian affected per equation
	# RO2_indx - index of components in peroxy radical list
	# RO_indx - index of components in alkoxy radical list
	# H2Oi - index of water
	# temp - temperature in chamber (K)
	# tempt - times that temperatures reached (s)
	# Pnow - pressure inside chamber (Pa)
	# self.light_stat - lights status
	# self.light_time - time light status attained (s)
	# self.daytime - time of day experiment starts (s)
	# self.lat - latitude of experiment (degrees)
	# self.lon - longitude of experiment (degrees)
	# self.af_path - path to actinic flux values
	# self.dayOfYear - number of days since 31st December experiment held on
	# self.photo_path - path to file with absorption cross-sections
	# 		and quantum yields
	# Jlen - number of photochemical reactions
	# con_infl_C - influx of components with constant concentration
	# (molecules/cc/s)
	# nrec_step - number of recording steps
	# dydt_vst - dictionary for holding change tendencies of specified
	#		components
	# siz_str - the size structure
	# num_sb - number of particle size bins
	# num_comp - number of components
	# self.seedi - index of components comprising seed material
	# seed_name - names of components comprising seed particles
	# seedx - mole ratio of components comprising seed material
	# core_diss - dissociation constant of seed
	# Psat - pure component saturation vapour pressure 
	# 	(# molecules/cm3 (air))
	# mfp - mean free path (m)
	# accom_coeff - accommodation coefficient
	# y_mw - molecular weight (g/mol)
	# surfT - surface tension (g/s2)
	# R_gas - ideal gas constant (kg.m2.s-2.K-1.mol-1)
	# NA - Avogadro's constant (molecules/mol)
	# y_dens - component densities (kg/m3)
	# x - particle radii (um)
	# Varr - particle volume (um3)
	# therm_sp - thermal speed (m/s)
	# act_coeff - activity coefficient
	# Cw - effective absorbing mass of wall (# molecules/cm3 (air))
	# kw - gas-wall mass transfer coefficient (/s)
	# Cfactor - conversion factor for concentrations (ppb/# molecules/cm3)
	# self.tf - transmission factor for natural sunlight
	# self.light_ad - marker for whether to adapt time interval for 
	#	changing natural light intensity
	# y_arr - index for arranging concentrations into matrix that 
	# 	allows reaction rate coefficient calculation
	# y_rind - index for the concentration array  that 
	# 	allows reaction rate coefficient calculation
	# uni_y_rind - unique index of reactants
	# y_pind - index for the concentration array  that 
	# 	allows product gains to be indexed
	# uni_y_pind - unique index of products
	# reac_col - column indices for sparse matrix of reaction 
	#		losses
	# prod_col - column indices for sparse matrix of production
	#		gains
	# rstoi_flat - 1D array of reactant stoichiometries per 
	#	equation
	# pstoi_flat - 1D array of product stoichiometries per
	#	equation
	# rr_arr - index for reaction rates to allow reactant loss
	# 	calculation
	# rr_arr_p - index for reaction rates to allow product gain
	#	calculation
	# rowvals - row indices of Jacobian elements
	# colptrs - indices of  rowvals corresponding to each column of
	# 	the Jacobian
	# wall_on - marker for whether wall partitioning turned on
	# jac_wall_indx - index of inputs to Jacobian from wall 
	# 	partitioning
	# jac_part_indx - index of inputs to Jacobian from particle
	#	partitioning
	# jac_extr_indx - index of inputs to Jacobian from extraction
	#	of chamber air
	# Vbou - volume boundary of size bins (um3)
	# N_perbin - number concentration of particles per size bin 
	#	(# particles/cm3 (air))
	# Vol0 - initial single particle volumes per size bin (um3)
	# rad0 - initial radius at particle centres (um)
	# np_sum - number concentration of newly nucleated particles 
	#		(#/cc (air))
	# new_partr - radius of newly nucleated particles (cm)
	# nucv1, v2, v3 - nucleation parameters
	# nuci - index of nucleating component
	# nuc_comp - the nucleating component
	# nuc_ad - marker for whether to reduce time step to allow 
	#	for accurate capture of nucleation
	# RH - relative humidities (fraction 0-1)
	# RHt - times through experiment at which relative humidities reached (s)
	# coag_on - whether coagulation to be modelled
	# inflectDp - particle diameter at which wall loss inflection occurs (m)
	# pwl_xpre - x value preceding inflection point
	# pwl_xpro - x value proceeding inflection point
	# inflectk - deposition rate at inflection (/s)
	# chamR - spherical-equivalent radius of chamber (m2)
	# McMurry_flag - marker for treament of particle deposition to walls
	# p_char - average number of charges per particle (/particle)
	# e_field - average electric field inside chamber (g.m/A.s3)
	# injectt - time of injection of components (s)
	# inj_indx - index of components being instantaneously injected after 
	#	experiment start
	# Ct - concentration(s) (ppb) of component(s) injected 
	#	instantaneously after experiment start
	# pmode - whether number size distributions expressed as modes or explicitly
	# pconc - concentration of injected particles (#/cc (air))
	# pconct - times of particle injection (s)
	# mean_rad - mean radius for particle number size 
	#	distribution (um)
	# lowsize - lower size bin boundary (um)
	# uppsize - upper size bin boundary (um)
	# std - standard deviation for injected particle number size 
	#	distributions
	# rbou - size bin radius bounds (um)
	# const_infl_t - times for constant influxes (s)
	# MV - molar volume of all components (cm3/mol)
	# rindx_aq - index of reactants for aqueous-phase 
	# pindx_aq - index of products for aqueous-phase
	# rstoi_aq - stoichiometry of reactants for aqueous-phase
	# pstoi_aq - stoichiometry of products for aqueous-phase
	# nreac_aq - number of reactants per aqueous-phase reaction
	# nprod_aq - number of products per aqueous-phase reaction
	# jac_stoi_aq - stoichiometry for Jacobian for aqueous-phase
	# njac_aq - number of Jacobian elements per aqueous-phase reaction 
	# jac_den_indx_aq - index of Jacobian denominators 
	# jac_indx_aq - index of  Jacobian for aqueous-phase
	# y_arr_aq - y indices for aqueous-phase
	# y_rind_aq - reactant indices for aqueous-phase 
	# uni_y_rind_aq - y indices for reactants for aqueous-phase
	# y_pind_aq - y indices for products for aqueous-phase
	# uni_y_pind_aq - y indices for products for aqueous-phase
	# reac_col_aq - columns of sparse matrix for aqueous-phase
	# prod_col_aq - columns of sparse matrix for aqueous-phase
	# rstoi_flat_aq - reactant stoichiometries for Jacobian for 
	# 	aqueous-phase
	# pstoi_flat_aq - product stoichiometries for Jacobian for
	#	aqueous-phase
	# rr_arr_aq - aqueous-phase reaction rate indices
	# rr_arr_p_aq - aqueous-phase reaction rate indices
	# eqn_num - number of reactions in gas- and aqueous-phase
	# partit_cutoff - the product of saturation vapour pressure
	#	and activity coefficient above which gas-particle
	#	partitioning assumed negligible
	# diff_vol - diffusion volumes of components according to 
	#		Fuller et al. (1969)
	# DStar_org - diffusion coefficient of components at initial temperature (cm2/s)
	# corei - index of core component
	# ser_H2O - whether to serialise the gas-particle partitioning of water
	# C_p2w - concentration of components on the wall due to particle
	# deposition to wall (# molecules/cm3)
	# the following inputs are used only for the saving module:
	# self.sch_name - path to chemical scheme
	# sav_nam - name of folder to save in
	# comp_namelist - chemical scheme name of components
	# dydt_trak - name of components to track change tendencies
	# space_mode - type of spacing used in particle size distribution
	# rbou00 - original particle size bin bounds
	# ub_rad_amp - amplificatin factor for upper bin size bound
	# indx_plot - indices of components to plot the gas-phase temporal profile of
	# comp0 - names of components to plot the gas-phase temporal profile of
	# self.inname - path to model variables file
	# rel_SMILES - SMILES strings of components in chemical scheme
	# Psat_Pa_rec - pure component saturation vapour pressures (Pa) at 298.15 K
	# Psat_Pa - pure component saturation vapour pressures (Pa) at starting temperature in chamber
	# OC - oxygen to carbon ratio of components
	# wat_hist - flag for history of particle-phase with respect to water partitioning,
	# 	where 0 is dry (therefore on the deliquescence curve) and 1 is wet 
	#	(therefore on the efflorescence curve)
	# Pybel_objects - the pybel objects for components
	# pcont - flag for whether particle injection continuous or instantaneous
	# dil_fac - chamber dilution factor (fraction of chamber/s)
	# NOi - NO index
	# HO2i - HO2 index
	# NO3i - NO3 index
	# z_prt_coeff - fraction of total gas-particle partitioning coefficient 
	#	below which partitioning to a particle size bin is treated as zero,
	#	e.g. because surface area of that size bin is tiny 
	# self.con_C_indx - index of components with constant 
	# 	gas-phase concentration
	# seed_eq_wat - whether seed particles to be equilibrated with water prior to ODE solver
	# Vwat_inc - whether suppled seed particle volume contains equilibrated water
	# tot_in_res - record of total input of injected components (ug/m3)
	# Compti - index for total injection record for instantaneously injected components
	# cont_inf_reci - index of components with continuous influx in record
	# cont_inf_i - index of components with continuous influx in concentration array
	# tot_in_res_indx - index of components with recorded influx
	# chamSA - chamber surface area (m2)
	# chamV - chamber volume (m3)
	# self.tf_UVC - transmission factor for 254 nm wavelength light
	# ------------------------------------------------------------

	# start timer
	st_time = time.time()

	step_no = 0 # track number of time steps
	sumt = 0. # track time through simulation (s)
	self.sumt = 0. # track time through simulation (s)
	# counters on updates
	light_time_cnt = 0 # light time status count
	gasinj_cnt = 0 # count on injection times of components
	if (pconct[0, 0] == 0. and len(pconct[0, :]) > 1 and pcont[0, 0] == 0): 
		seedt_cnt = 1 # count on injection times of particles
	else:
		seedt_cnt = 0
	
	# current status of whether injection of particles instantaneous or continuous,
	# if not stated assume instantaneous
	pcontf = 0
	if (pconct[0, 0] == 0 and pcont[0, 0] == 1):
		pcontf = 1
	infx_cnt = 0 # count on constant gas-phase influx occurrences
	infx_cnt0 = 0 # remember count at start of integration step
	tempt_cnt = 0 # count on chamber temperatures
	tempt_cnt0 = 0 # remember count at start of integration step
	RHt_cnt = 0 # count on chamber relative humidities
	RHt_cnt0 = 0 # remember count at start of integration step
	conPin_cnt = 0 # count on continuous influx of seed particles
	conPin_cnt0 = 0 # remember count at start of integration step
	# count on recording results, note starting on two because results at t=0 already stored
	save_cnt = 2
	
	# count on time since update to integration initial values/constants last called (s)
	update_count = 0.
	y0 = np.zeros((len(y))) # remember initial concentrations (molecules/cm3 (air))
	N_perbin0 = np.zeros((N_perbin.shape[0], N_perbin.shape[1])) # remember initial particle number concentrations (# particles/cm3)
	x0 = np.zeros((len(x)))# remember initial particle sizes (um)
	t0 = update_stp # remember initial integration step (s)
	# flag for changing integration time step due to changing initial values	
	ic_red = 0
	tnew = update_stp # the time to integrate over (s)
	# fraction of newly injected seed particles
	pconcn_frac = 0.
	comp_namelist_np = np.array(comp_namelist) # numpy array version of chemical scheme names
	save_cnt_chck = 1 # counting recording steps for time interval check
	
	# find out what to do with the gas-wall partitioning coefficient
	if (kw == -1):
		kwf = -1 # Huang et al. 2018 treatment
	else:
		# standard PyCHAM (GMD paper) treatment with same gas-wall 
		# partitioning coefficient for all components
		kwf = kw
		
	# prepare recording matrices, including recording of initial
	# conditions, note initial change tendencies not recorded 
	# in this call but are recorded below
	[trec, yrec, Cfactor_vst, Nres_dry, Nres_wet, x2, 
	seedt_cnt, rbou_rec, Cfactor, infx_cnt, 
	yrec_p2w, temp_now, cham_env, Pnow, Psat, 
	RHn, Cinfl_now, tot_in_res_ft] = rec_prep.rec_prep(nrec_steps, y, y0, rindx, 
	rstoi, pindx, pstoi, nprod, nreac, 
	num_sb, num_comp, N_perbin, core_diss, Psat, mfp,
	accom_coeff, y_mw, surfT, R_gas, temp, tempt, NA,
	y_dens*1.e-3, x, therm_sp, H2Oi, act_coeff,
	RO2_indx, sumt, Pnow, light_time_cnt, 
	Jlen, Cw, kw, Cfactor, 
	wall_on, Vbou, tnew, nuc_ad, nucv1, nucv2, nucv3, 
	np_sum, update_stp, update_count, injectt, gasinj_cnt, 
	inj_indx, Ct, pmode, pconc, pconct, seedt_cnt, mean_rad, corei, 
	seed_name, seedx, lowsize, uppsize, rad0, x, std, rbou, const_infl_t, 
	infx_cnt, con_infl_C, MV, partit_cutoff, diff_vol, DStar_org, 
	C_p2w, RH, RHt, tempt_cnt, RHt_cnt, Pybel_objects, nuci, 
	nuc_comp, t0, pcont, pcontf, NOi, HO2i, NO3i, z_prt_coeff,
	seed_eq_wat, Vwat_inc, tot_in_res, Compti, tot_time, cont_inf_reci, 
	cont_inf_i, tot_in_res_indx, chamSA, chamV, kwf, self)
	
	import ode_solv
	importlib.reload(ode_solv) # import most recent version
	importlib.reload(ode_solv_wat) # import most recent version
	importlib.reload(dydt_rec) # import most recent version

	while (tot_time-sumt) > (tot_time/1.e10):
		
		# remembering variables at the start of the integration step ------------------------------------------
		y0[:] = y[:] # remember initial concentrations (# molecules/cm3 (air))
		N_perbin0[:] = N_perbin[:] # remember initial particle number concentration (# particles/cm3)
		x0[:] = x[:] # remember initial particle sizes (um)
		temp_now0 = temp_now # remember temperature (K)
		wat_hist0 = wat_hist # remember water history flag at start of integration step
		RH0 = RHn # relative humidity at start of integration step
		Pnow0 = Pnow # pressure (Pa)
		
		# remember counts at start of integration step
		infx_cnt0 = infx_cnt
		tempt_cnt0 = tempt_cnt
		RHt_cnt0 = RHt_cnt
		seedt_cnt0 = seedt_cnt
		gasinj_cnt0 = gasinj_cnt
		light_time_cnt0 = light_time_cnt
		conPin_cnt0 = conPin_cnt
		
		# --------------------------------------------------------------------------------------------------------------------
		
		
		gpp_stab = 0 # flag for stability in gas-particle partitioning for solver while loop
		stab_red = 0 # flag for stability in gas-particle partitioning for time interval reset
		lin_int = 0 # flag to linearly interpolate changes to chamber
		t00 = tnew # remember the initial integration step for this integration step (s)
		save_cntf = 0 # flag for updating count on number of recordings
		
		while (gpp_stab != 1): # whilst ode solver flagged as unstable

			# if integration interval decreased, reset concentrations to those at start of interval
			if (gpp_stab == -1):
				y[:] = y0[:] # (# molecules/cm3)
			
			# update chamber variables
			[temp_now, Pnow, lightm, light_time_cnt, tnew, ic_red, update_stp, 
			update_count, Cinfl_now, seedt_cnt, Cfactor, infx_cnt, 
			gasinj_cnt, DStar_org, y, tempt_cnt, RHt_cnt, Psat, N_perbin, x,
			pconcn_frac,  pcontf, tot_in_res, Cinfl_nowp_indx, 
			Cinfl_nowp] = cham_up.cham_up(sumt, temp, tempt, 
			Pnow0, light_time_cnt0, 
			tnew, nuc_ad, nucv1, nucv2, nucv3, np_sum, 
			update_stp, update_count, 
			injectt, gasinj_cnt0, inj_indx, Ct, pmode, pconc, pconct, 
			seedt_cnt0, num_comp, y0, y, N_perbin0, mean_rad, corei, seedx, seed_name, 
			lowsize, uppsize, num_sb, MV, rad0, x0, std, y_dens, H2Oi, rbou, 
			const_infl_t, infx_cnt0, con_infl_C, wall_on, Cfactor, diff_vol, 
			DStar_org, RH, RHt, tempt_cnt0, RHt_cnt0, Pybel_objects, nuci, nuc_comp,
			y_mw, temp_now0, Psat, gpp_stab, t00, x0, pcont,  pcontf, Cinfl_now, surfT,
			act_coeff, seed_eq_wat, Vwat_inc, tot_in_res, Compti, tot_time, 
			cont_inf_reci, cont_inf_i, self)
			
			# aligning time interval with pre-requisites -------------------------
			# ensure end of time interval does not surpass recording time
			if ((sumt+tnew) > save_stp*save_cnt_chck):
				tnew = (save_stp*(save_cnt_chck))-sumt
				# temporarily set the update step for operator-split processes
				# to align with recording time step, this ensures that
				# recording and operator-split intervals don't fall out of sync
				update_stp = tnew
				update_count = 0.
				ic_red = 1
			
			# ensure update to operator-split processes interval not surpassed
			if (update_count+tnew > update_stp):
				tnew = (update_stp-update_count)
				ic_red = 1
			
			# ensure simulation end time not surpassed
			if (sumt+tnew > tot_time):
				tnew = (tot_time-sumt)
				ic_red = 1
			# ------------------------------------------------------------------		
			
			if ((num_sb-wall_on) > 0 or wall_on == 1): # if particles present
			
				# update partitioning variables
				[kimt, kelv_fac, kw] = partit_var.kimt_calc(y, mfp, num_sb, num_comp, accom_coeff, 
				y_mw, surfT, R_gas, temp_now, NA, y_dens, N_perbin, 
				x.reshape(1, -1)*1.e-6, Psat, therm_sp, H2Oi, act_coeff, wall_on, 1, partit_cutoff, 
				Pnow, DStar_org, z_prt_coeff, chamSA, chamV, kwf, self)
			
				# update particle-phase activity coefficients, note the output,
				# note that if ODE solver unstable, then y resets to y0 via
				# the cham_up module prior to this call
				[act_coeff, wat_hist, RHn, y, 
				dydt_erh_flag] = act_coeff_update.ac_up(y, H2Oi, RH0, temp_now, 
				wat_hist0, act_coeff, num_comp, (num_sb-wall_on))
				
			else: # fillers
				kimt = np.zeros((num_sb-wall_on, num_comp))
				kelv_fac = np.zeros((num_sb-wall_on, 1))
				dydt_erh_flag = 0
		
			# reaction rate coefficient
			[rrc, erf, err_mess] = rrc_calc.rrc_calc(RO2_indx, 
				y[H2Oi], temp_now, lightm, y, 
				Pnow, Jlen, y[NOi], y[HO2i], y[NO3i], 
				sumt, self)
			
			if (erf == 1): # if error message from reaction rate calculation
				yield(err_mess)
			
			# update Jacobian inputs based on particle-phase fractions of components
			[rowvalsn, colptrsn, jac_part_indxn, jac_mod_len, jac_part_hmf_indx, rw_indx, jac_wall_indxn, 
			jac_part_H2O_indx] = jac_up.jac_up(y[num_comp:num_comp*((num_sb-wall_on+1))], rowvals, 
			colptrs, (num_sb-wall_on), num_comp, jac_part_indx, H2Oi, y[H2Oi], jac_wall_indx, ser_H2O)
			
			
			# for change tendencies, do want to save from t=0
			# record any change tendencies of specified components
			if (len(dydt_vst) > 0 and save_cntf == 0):
				if (sumt == 0. or (sumt-(save_stp*(save_cnt-1)) > -1.e-10)):
					if (sumt == 0.):
						dydt_cnt = save_cnt-2
					if (sumt-(save_stp*(save_cnt-1)) > -1.e-10):
						dydt_cnt = save_cnt-1
					
					# before solving ODEs for chemistry, gas-particle partitioning and gas-wall partitioning, 
					# estimate and record any change tendencies (# molecules/cm3/s) resulting from these processes
					dydt_vst = dydt_rec.dydt_rec(y, rindx, rstoi, rrc, pindx, pstoi, nprod, dydt_cnt, 
						dydt_vst, nreac, num_sb, num_comp, pconc, core_diss, Psat, kelv_fac, 
						kimt, kw, Cw, act_coeff, dydt_erh_flag, H2Oi, wat_hist, self)
			
			# record output if on the first attempt at solving this time interval, 
			# note that recording here in this way means we include any
			# instantaneous changes at this time step without interpolation to
			# smaller instantaneous changes (when interpolation forced due to 
			# instability)
			if (save_cntf == 0 and (sumt-(save_stp*(save_cnt-1)) > -1.e-10)):
				
				[trec, yrec, Cfactor_vst, save_cntf, Nres_dry, Nres_wet, 
				x2, rbou_rec, yrec_p2w, cham_env, tot_in_res_ft] = rec.rec(save_cnt-1, trec, yrec, 
				Cfactor_vst, y, sumt, rindx, rstoi, rrc, pindx, pstoi, 
				nprod, nreac, num_sb, num_comp, N_perbin, core_diss, 
				Psat, kelv_fac, kimt, kw, Cw, act_coeff, Cfactor, Nres_dry, 
				Nres_wet, x2, x, MV, H2Oi, Vbou, rbou, wall_on, rbou_rec, 
				yrec_p2w, C_p2w, cham_env, temp_now, Pnow, tot_in_res, tot_in_res_ft, self)
				# prepare for recording next point
				save_cnt += 1

			if (ser_H2O == 1 and (num_sb-wall_on) > 0 and (sum(N_perbin) > 0)): # if water gas-particle partitioning serialised

				# if on the deliquescence curve rather than the 
				# efflorescence curve in terms of water gas-particle partitioning
				if (wat_hist == 1):
					
					# call on ode solver for water
					[y, res_t] = ode_solv_wat.ode_solv(y, tnew, rindx, pindx, rstoi, pstoi,
					nreac, nprod, rrc, jac_stoi, njac, jac_den_indx, jac_indx,
					Cinfl_now, y_arr, y_rind, uni_y_rind, y_pind, uni_y_pind, 
					reac_col, prod_col, rstoi_flat, 
					pstoi_flat, rr_arr, rr_arr_p, rowvalsn, colptrsn, num_comp, 
					num_sb, wall_on, Psat, Cw, act_coeff, kw, jac_wall_indxn,
					core_diss, kelv_fac, kimt, (num_sb-wall_on), 
					jac_part_indxn,
					rindx_aq, pindx_aq, rstoi_aq, pstoi_aq,
					nreac_aq, nprod_aq, jac_stoi_aq, njac_aq, jac_den_indx_aq, jac_indx_aq, 
					y_arr_aq, y_rind_aq, uni_y_rind_aq, y_pind_aq, uni_y_pind_aq, 
					reac_col_aq, prod_col_aq, rstoi_flat_aq, 
					pstoi_flat_aq, rr_arr_aq, rr_arr_p_aq, eqn_num, jac_mod_len, 
					jac_part_hmf_indx, rw_indx, N_perbin, jac_part_H2O_indx, H2Oi, self)
					
					if (any(y[H2Oi::num_comp] < 0.)): # check on stability of water partitioning
						
						# identify components with negative concentrations
						neg_comp_indx = y < 0.
						# transform into components in columns, locations in rows
						neg_comp_indx = neg_comp_indx.reshape(num_sb+1, num_comp)
						# get component indices with negative concentration
						neg_comp_indx = np.unique((np.where(neg_comp_indx == 1))[1])
						# get chemical scheme names of components with negative concentration
						neg_names = comp_namelist_np[neg_comp_indx]

						y_H2O = y[H2Oi::num_comp] # isolate just water concentrations (molecules/cm3)
						# sum the negative concentrations and convert to absolute value (molecules/cm3)
						neg_H2O = np.abs(sum(y_H2O[y_H2O<0.]))
						
						# allow a given fraction of water concentrations to be negative
						if (neg_H2O/sum(np.abs(y[H2Oi::num_comp])) > 0. ):
				
							gpp_stab = -1 # maintain unstable flag
							# tell user what's happening
							yield (str('Note: negative water concentration generated following call to ode_solv_wat module, the program assumes this is because of a change in relative humidity in chamber air, and will automatically half the integration time interval and linearly interpolate any change to chamber conditions supplied by the user.  To stop this the simulation must be cancelled using the Quit button in the PyCHAM graphical user interface.  Current update time interval is ' + str(tnew) + ' seconds'))
							
							if (tnew < 1.e-20): # if time step has decreased to unreasonably low and solver still unstable then break
								ode_brk_err_mess.ode_brk_err_mess(y0, neg_names, rindx, y_arr, 
									y_rind, rstoi, pstoi, rrc, nreac, nprod, wall_on, num_comp, 
									(num_sb-wall_on), Cw, Psat, act_coeff, kw, neg_comp_indx, 
									N_perbin, core_diss, kelv_fac, kimt, eqn_num, rindx_aq, y_rind_aq, y_arr_aq, 
									rstoi_aq, pstoi_aq, nreac_aq, nprod_aq, 0, H2Oi, y, self)

								yield (str('Error: negative concentrations generated following call to ode_solv_wat module, the program has assumed this is because of a change in chamber condition (e.g. injection of components), and has automatically halved the integration time interval and linearly interpolated any change to chamber conditions supplied by the user.  However, the integration time interval has now decreased to ' + str(tnew) + ' seconds, which is assumed too small to be useful, so the program has been stopped.  The components with negative concentrations are : ' + str(neg_names) + '.  The problem could be too stiff for the solver and the relevant fluxes (change tendencies) have been output to the file ODE_solver_break_relevant_fluxes.txt for your analysis of problem stiffness.  You could identify the maximum and minimum fluxes to gain indication of the components and/or processes making the problem stiff.  Therefafter you could modify the relevant model variables (supplied by the user) and the chemical scheme (supplied by the user).' ))
							# half the update and integration time step (s) if necessary
							tnew = tnew/2.
							stab_red = 1 # remember that time step temporarily reduced due to instability
							continue

						else: # if acceptable
							gpp_stab = 1 # change to stable flag
						
					else: # if solution stable, change stability flag to represent this
						gpp_stab = 1 # change to stable flag
					
				# zero partitioning of water to particles for integration without 
				# water gas-particle partitioning
				kimt[:, H2Oi] = 0.
			
			# model component concentration changes to get new concentrations
			# (# molecules/cm3 (air))
			[y, res_t] = ode_solv.ode_solv(y, tnew, rindx, pindx, rstoi, pstoi,
				nreac, nprod, rrc, jac_stoi, njac, jac_den_indx, jac_indx,
				Cinfl_now, y_arr, y_rind, uni_y_rind, y_pind, uni_y_pind, 
				reac_col, prod_col, rstoi_flat, 
				pstoi_flat, rr_arr, rr_arr_p, rowvalsn, colptrsn, num_comp, 
				num_sb, wall_on, Psat, Cw, act_coeff, kw, jac_wall_indxn,
				core_diss, kelv_fac, kimt, (num_sb-wall_on), 
				jac_part_indxn, jac_extr_indx,
				rindx_aq, pindx_aq, rstoi_aq, pstoi_aq,
				nreac_aq, nprod_aq, jac_stoi_aq, njac_aq, jac_den_indx_aq, jac_indx_aq, 
				y_arr_aq, y_rind_aq, uni_y_rind_aq, y_pind_aq, uni_y_pind_aq, 
				reac_col_aq, prod_col_aq, rstoi_flat_aq, 
				pstoi_flat_aq, rr_arr_aq, rr_arr_p_aq, eqn_num, jac_mod_len, 
				jac_part_hmf_indx, rw_indx, N_perbin, jac_part_H2O_indx, 
				H2Oi, dil_fac, RO2_indx[:, 1], comp_namelist, Psat_Pa, Cinfl_nowp_indx, 
				Cinfl_nowp, self)
			
			# if any components set to have constant gas-phase concentration
			if (any(self.con_C_indx)): # then keep constant
				y[self.con_C_indx] = y0[self.con_C_indx] # (# molecules/cm3)

			# if negative, suggests ODE solver instability, but could also be numerical limits, 
			# especially if concentrations are relatively close to zero, so allow some leeway
			if (any(y < -np.sum(np.abs(y))*1.e-50)):
			
				# identify components with negative concentrations
				neg_comp_indx = y < 0.
				# transform into components in columns, locations in rows
				neg_comp_indx = neg_comp_indx.reshape(num_sb+1, num_comp)
				# get component indices with negative concentration
				neg_comp_indx = np.unique((np.where(neg_comp_indx == 1))[1])
				# get chemical scheme names of components with negative concentration
				neg_names = comp_namelist_np[neg_comp_indx]				

				# loop through components with negative concentrations
				for ci in neg_comp_indx:
					all_c = y[ci::num_comp]
					# sum of absolute negative concentrations
					neg_c = np.abs(sum(all_c[all_c<0.]))
					# sum of absolute of all concentrations
					sum_c = sum(np.abs(all_c))
					
				gpp_stab = -1 # maintain unstable flag
				# tell user what's happening
				yield (str('Note: negative concentrations generated following call to ode_solv module, the program assumes this is because of a change in chamber condition (e.g. injection of components), and will automatically half the integration time interval and linearly interpolate any change to chamber conditions supplied by the user.  To stop this the simulation must be cancelled using the Quit button in the PyCHAM graphical user interface.  Current integration time interval is ' + str(tnew) + ' seconds'))
				#y_gp = np.zeros((len(y)))
				#y_gp[:] = y[:]
				#y_gp = y_gp.reshape(num_comp, num_sb-wall_on+1, order='F')
				#y_gpi = np.where(y_gp<0.)
				
				# in case more detailed look at negative issue wanted
				#for gpi in y_gpi[0]: # loop through negative components
				#	print('negative component name:', comp_namelist[gpi])
				#	print(y0[gpi::num_comp])
				#	print(y[gpi::num_comp])
				#	import ipdb; ipdb.set_trace()
				#print('negative component size bin:', np.unique(y_gpi[1][:]))
					
				if (tnew < 1.e-20): # if time step has decreased to unreasonably low and solver still unstable then break
					# estimate gas-phase reaction fluxes for all reactions and partitioning fluxes for troublesome components
					ode_brk_err_mess.ode_brk_err_mess(y0, neg_names, rindx, y_arr, 
						y_rind, rstoi, pstoi, rrc, nreac, nprod, wall_on, 
						num_comp, (num_sb-wall_on), Cw, Psat, act_coeff, kw, 
						neg_comp_indx, N_perbin, core_diss, kelv_fac, 
						kimt, eqn_num, rindx_aq, y_rind_aq, y_arr_aq, 
						rstoi_aq, pstoi_aq, nreac_aq, nprod_aq, 1, H2Oi, y, self)

					yield (str('Error: negative concentrations generated following call to ode_solv module, the program has assumed this is because of a change in chamber condition (e.g. injection of components), and has automatically halved the integration time interval and linearly interpolated any change to chamber conditions supplied by the user.  However, the integration time interval has now decreased to ' + str(tnew) + ' seconds, which is assumed too small to be useful, so the program has been stopped.  The components with negative concentrations are : ' + str(neg_names) + '.  The problem could be too stiff for the solver and the relevant fluxes (change tendencies) have been output to the file ODE_solver_break_relevant_fluxes.txt for your analysis of problem stiffness.  You could identify the maximum and minimum fluxes to gain indication of the components and/or processes making the problem stiff.  Therefafter you could modify the relevant model variables (supplied by the user) and the chemical scheme (supplied by the user).' ))
						
				# half the update and integration time step (s) if necessary	
				tnew = tnew/2.
				stab_red = 1 # remember that time step temporarily reduced due to instability
				
			else: # if solution stable, change stability flag to represent this
				
				# account for any partial addition of newly injected seed particles
				pconc[:, seedt_cnt] -= pconc[:, seedt_cnt]*pconcn_frac
				# reset fraction of newly injected seed particles
				pconcn_frac = 0.
				gpp_stab = 1 # change to stable flag
			
		# end of integration stability condition section ----------------------------
		step_no += 1 # track number of steps
		sumt += tnew # total time through simulation (s)
		self.sumt += tnew
		
		if (any(self.obs_comp_i)): # if any components constrained to observations
				# get observed concentrations now
				for ci in range(len(self.obs_comp_i)): # loop through components
					y[self.obs_comp_i[ci]] = np.interp(sumt, self.obs[:, 0], self.obs[:, ci+1])

		if (sumt%save_stp < 1.e-12 or (sumt%save_stp-save_stp) < 1.e-12): # get remainder
			save_cnt_chck += 1 # if need to move up count on recording
		
		# dilute chamber particle number following an integration time step, e.g. for flow-reactor -----
		# note that concentrations of components inside particles (and in the gas-phase) will have been
		# reduced due to dilution inside the ODE solver
		if (dil_fac > 0):
			N_perbin -= N_perbin*(dil_fac*tnew)
		
		if ((num_sb-wall_on) > 0): # if particle size bins present
			# update particle sizes
			if ((num_sb-wall_on) > 1) and (any(N_perbin > 1.e-10)): # if particles present
				
				if (siz_str == 0): # moving centre
					(N_perbin, Varr, y, x, redt, t, bc_red) = mov_cen.mov_cen_main(N_perbin, 
					Vbou, num_sb, num_comp, y_mw, x, Vol0, tnew, 
					update_stp, y0, MV, Psat[0, :], ic_red, y, res_t, wall_on)
				
				if (siz_str == 1): # full-moving
					(Varr, x, y[num_comp:(num_comp*(num_sb-wall_on+1))], 
					N_perbin, Vbou, rbou) = fullmov.fullmov((num_sb-wall_on), N_perbin,
 					num_comp, y[num_comp:(num_comp)*(num_sb-wall_on+1)], MV*1.e12, 
					Vol0, Vbou, rbou)
			
			update_count += tnew # time since operator-split processes last called (s)
			
			# if time met to implement operator-split processes
			if (update_count >= (update_stp*9.999999e-1)):
				if (any(N_perbin > 1.e-10)):
				
					# particle-phase concentration(s) (# molecules/cm3 (air))
					Cp = np.transpose(y[num_comp:(num_comp)*(num_sb-wall_on+1)].reshape(num_sb-wall_on, num_comp))
				
					# coagulation
					[N_perbin, y[num_comp:(num_comp)*(num_sb-wall_on+1)], x, Gi, eta_ai, 
						Varr, Vbou, rbou] = coag.coag(RH[RHt_cnt], temp_now, x*1.e-6, 
						(Varr*1.0e-18).reshape(1, -1), 
						y_mw.reshape(-1, 1), x*1.e-6, 
						Cp, (N_perbin).reshape(1, -1), update_count, 
						(Vbou*1.0e-18).reshape(1, -1), rbou,
						num_comp, 0, (np.squeeze(y_dens*1.e-3)), Vol0, rad0, Pnow, 0,
						Cp, (N_perbin).reshape(1, -1), (Varr*1.e-18).reshape(1, -1),
						coag_on, siz_str, wall_on)
					
					if ((McMurry_flag > -1) and (wall_on == 1)): #if particle loss to walls turned on
						# particle loss to walls
						[N_perbin, 
						y[num_comp:(num_comp)*(num_sb-wall_on+1)]] = wallloss.wallloss(
							N_perbin.reshape(-1, 1), 
							y[num_comp:(num_comp)*(num_sb-wall_on+1)], Gi, eta_ai,
 							x*2.0e-6, y_mw, 
							Varr*1.0e-18, (num_sb-wall_on), num_comp, temp_now, update_count, 
							inflectDp, pwl_xpre, pwl_xpro, inflectk, chamR, McMurry_flag, 
							0, p_char, e_field, (num_sb-wall_on), C_p2w)
			
				if (nucv1 > 0.): # nucleation
					
					[N_perbin, y, x, Varr, np_sum, rbou, Vbou] = nuc.nuc(sumt, np_sum, 
						N_perbin, y, y_mw.reshape(-1, 1), 
						np.squeeze(y_dens*1.0e-3),  
						num_comp, Varr, x, new_partr, MV, nucv1, nucv2, 
						nucv3, nuc_comp[0], siz_str, rbou, Vbou, (num_sb-wall_on))
				
				# reset count that tracks when next operator-split should be called (s)
				update_count = 0.
		
		# update the percentage time in the GUI progress bar
		yield (sumt/tot_time*100.)
		
		if (sumt >= (tot_time-tot_time/1.e10)): # record output at experiment end
			
			# update chamber variables, note this ensures that any changes made
			# coincidentally with the experiment end are captured
			[temp_now, Pnow, lightm, light_time_cnt, tnew, ic_red, update_stp, 
			update_count, Cinfl_now, seedt_cnt, Cfactor, infx_cnt, 
			gasinj_cnt, DStar_org, y, tempt_cnt, RHt_cnt, Psat, N_perbin, x,
			pconcn_frac,  pcontf, tot_in_res, Cinfl_nowp_indx, 
			Cinfl_nowp] = cham_up.cham_up(sumt, temp, tempt, 
			Pnow0, light_time_cnt0, 
			tnew, nuc_ad, nucv1, nucv2, nucv3, np_sum, 
			update_stp, update_count, 
			injectt, gasinj_cnt0, inj_indx, Ct, pmode, pconc, pconct, 
			seedt_cnt0, num_comp, y0, y, N_perbin0, mean_rad, corei, seedx, seed_name, 
			lowsize, uppsize, num_sb, MV, rad0, x0, std, y_dens, H2Oi, rbou, 
			const_infl_t, infx_cnt0, con_infl_C, wall_on, Cfactor, diff_vol, 
			DStar_org, RH, RHt, tempt_cnt0, RHt_cnt0, Pybel_objects, nuci, nuc_comp,
			y_mw, temp_now0, Psat, gpp_stab, t00, x0, pcont,  pcontf, Cinfl_now, surfT,
			act_coeff, seed_eq_wat, Vwat_inc, tot_in_res, Compti, tot_time, 
			cont_inf_reci, cont_inf_i, self)
			
			[trec, yrec, Cfactor_vst, save_cnt, Nres_dry, Nres_wet,
			x2, rbou_rec, yrec_p2w, cham_env, tot_in_res_ft] = rec.rec(save_cnt-1, 
			trec, yrec, Cfactor_vst, y, sumt, rindx, rstoi, rrc, pindx, pstoi, 
			nprod, nreac, num_sb, num_comp, N_perbin, core_diss, 
			Psat, kelv_fac, kimt, kw, Cw, act_coeff, Cfactor, Nres_dry, 
			Nres_wet, x2, x, MV, H2Oi, Vbou, rbou, wall_on, rbou_rec, 
			yrec_p2w, C_p2w, cham_env, temp_now, Pnow, tot_in_res, tot_in_res_ft, self)		
		
		# if time step was temporarily reduced, then reset
		if (ic_red == 1 or stab_red == 1):
			update_stp = t0
			tnew = update_stp
			ic_red = 0 # reset flag
			stab_red = 0 # reset flag
		
		# remember the gas-phase water concentration from previous 
		# integration step (# molecules/cm3)
		y_H2O0 = y[H2Oi]
	
	time_taken = time.time()-st_time

	# save results
	save.saving(yrec, Nres_dry, Nres_wet, trec, sav_nam, 
		dydt_vst, num_comp, Cfactor_vst, 0, 
		num_sb, comp_namelist, dydt_trak, y_mw, MV, time_taken, 
		seed_name, x2, rbou_rec, wall_on, space_mode, rbou00, ub_rad_amp, indx_plot, 
		comp0, yrec_p2w, rel_SMILES, Psat_Pa_rec, OC, H2Oi, 
		siz_str, cham_env, RO2_indx[:, 1], RO_indx, tot_in_res_ft, self)
	return()
