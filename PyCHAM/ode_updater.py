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
import importlib
import save
import time
import act_coeff_update
# providing error message if ODE solver produces 
# negative results below minimum integration time
import ode_brk_err_mess


def ode_updater(y, H2Oi, 
	Pnow, Jlen, nrec_steps, 
	siz_str, num_sb, num_comp, 
	core_diss, mfp, therm_sp,
	accom_coeff, y_mw, surfT, R_gas, NA, 
	x, Varr, act_coeff, Cfactor, rowvals, colptrs, Vbou,
	N_perbin, Vol0, rad0, np_sum, new_partr, 
	nuci, coag_on, inflectDp, pwl_xpre, 
	pwl_xpro, inflectk, chamR, McMurry_flag, p_char, e_field, 
	injectt, inj_indx, Ct, lowsize, 
	uppsize, std, rbou, MV, 
	diff_vol, DStar_org, corei, ser_H2O, 
	sav_nam, rbou00, ub_rad_amp, indx_plot,
	wat_hist, NOi, 
	HO2i, NO3i, z_prt_coeff, tot_in_res,
	Compti, tot_in_res_indx, chamSA, 
	chamV, tempt_cnt, self, vol_Comp, volP):
	
	# inputs: ----------------------------------------------------
	# self.update_stp - interval at which to update integration 
	#		constants (s)
	# self.tot_time - total time to simulate (s)
	# self.save_step - frequency to save results (s)
	# y - initial component concentrations (# molecules/cm3/s)
	# self.rindx_g - index of reactants per equation
	# self.pindx_g - index of products per equation
	# self.nreac_g - number of reactants per equation
	# self.nprod_g - number of products per equation
	# self.jac_stoi_g - stoichiometries relevant to Jacobian
	# self.njac_g - number of elements of Jacobian affected per 
	# equation
	# self.jac_den_indx_g - index of denominator components for 
	# Jacobian
	# self.jac_indx_g - index of Jacobian affected per equation
	# self.RO2_indx - index of components in alkyl peroxy radical 
	# list
	# self.RO_indx - index of components in alkoxy radical list
	# H2Oi - index of water
	# self.TEMP - temperature in chamber (K)
	# self.tempt - times that temperatures reached (s)
	# Pnow - pressure inside chamber (Pa)
	# self.light_stat - lights status
	# self.light_time - time light status attained (s)
	# self.daytime - time of day experiment starts (s)
	# self.lat - latitude of experiment (degrees)
	# self.lon - longitude of experiment (degrees)
	# self.af_path - path to actinic flux values
	# self.dayOfYear - number of days since 31st December 
	#	experiment held on
	# self.photo_path - path to file with absorption cross-sections
	# 		and quantum yields
	# Jlen - number of photochemical reactions
	# self.con_infl_C - influx of components with continuous 
	#	influx (# molecules/cm3/s)
	# nrec_step - number of recording steps
	# self.dydt_vst - dictionary for holding change tendencies of 
	#	specified components
	# siz_str - the size structure
	# num_sb - number of particle size bins
	# num_comp - number of components
	# self.seedi - index of components comprising seed material
	# self.seed_name - names of components comprising seed particles
	# self.seedx - mole ratio of components comprising seed material
	# core_diss - dissociation constant of seed
	# self.Psat - pure component saturation vapour pressure 
	# 	(# molecules/cm3 (air))
	# mfp - mean free path (m)
	# accom_coeff - accommodation coefficient
	# y_mw - molecular weight (g/mol)
	# surfT - surface tension (g/s2)
	# R_gas - ideal gas constant (kg.m2.s-2.K-1.mol-1)
	# NA - Avogadro's constant (molecules/mol)
	# self.y_dens - component densities (kg/m3)
	# x - particle radii (um)
	# Varr - particle volume (um3)
	# therm_sp - thermal speed (m/s)
	# act_coeff - activity coefficient
	# self.Cw - effective absorbing mass of wall (# molecules/cm3 
	#	(air))
	# self.kw - gas-wall mass transfer coefficient (/s)
	# Cfactor - conversion factor for concentrations (ppb/# 
	#	molecules/cm3)
	# self.tf - transmission factor for natural sunlight
	# self.light_ad - marker for whether to adapt time interval for 
	#	changing natural light intensity
	# self.y_arr_g - index for arranging concentrations into matrix 
	# that allows reaction rate coefficient calculation
	# self.y_rind_g - index for the concentration array  that 
	# 	allows reaction rate coefficient calculation
	# self.uni_y_rind_g - unique index of reactants
	# self.y_pind_g - index for the concentration array  that 
	# 	allows product gains to be indexed
	# self.uni_y_pind_g - unique index of products
	# self.reac_co_gl - column indices for sparse matrix of reaction 
	#		losses
	# self.prod_col_g - column indices for sparse matrix of 
	# production gains
	# self.rstoi_flat_g - 1D array of reactant stoichiometries per 
	#	equation
	# self.pstoi_flat_g - 1D array of product stoichiometries per
	#	equation
	# self.rr_arr_g - index for reaction rates to allow reactant loss
	# 	calculation
	# self.rr_arr_p_g - index for reaction rates to allow product gain
	#	calculation
	# rowvals - row indices of Jacobian elements
	# colptrs - indices of  rowvals corresponding to each column of
	# 	the Jacobian
	# self.wall_on - marker for whether wall partitioning turned on
	# self.jac_wall_indx - index of inputs to Jacobian from wall 
	# 	partitioning
	# self.jac_part_indx - index of inputs to Jacobian from particle
	#	partitioning
	# self.jac_extr_indx - index of inputs to Jacobian from 
	#	extraction of chamber air
	# Vbou - volume boundary of size bins (um3)
	# N_perbin - number concentration of particles per size bin 
	#	(# particles/cm3 (air))
	# Vol0 - initial single particle volumes per size bin (um3)
	# rad0 - initial radius at particle centres (um)
	# np_sum - number concentration of newly nucleated particles 
	#		(#/cc (air))
	# new_partr - radius of newly nucleated particles (cm)
	# self.nucv1, v2, v3 - nucleation parameters
	# nuci - index of nucleating component
	# self.nuc_comp - the nucleating component
	# self.nuc_ad - marker for whether to reduce time step to allow 
	#	for accurate capture of nucleation
	# self.RH - relative humidities (fraction 0-1)
	# self.RHt - times through experiment at which relative 
	#	humidities reached (s)
	# coag_on - whether coagulation to be modelled
	# inflectDp - particle diameter at which wall loss inflection 
	#	occurs (m)
	# pwl_xpre - x value preceding inflection point
	# pwl_xpro - x value proceeding inflection point
	# inflectk - deposition rate at inflection (/s)
	# chamR - spherical-equivalent radius of chamber (m2)
	# McMurry_flag - marker for treament of particle deposition to 
	#	walls
	# p_char - average number of charges per particle (/particle)
	# e_field - average electric field inside chamber (g.m/A.s3)
	# injectt - time of injection of components (s)
	# inj_indx - index of components being instantaneously injected 
	# after experiment start
	# Ct - concentration(s) (ppb) of component(s) injected 
	#	instantaneously after experiment start
	# self.pmode - whether number size distributions expressed as modes
	# 	or explicitly
	# self.pconc - concentration of injected particles (#/cm3 (air))
	# self.pconct - times of particle injection (s)
	# self.mean_rad - mean radius for particle number size 
	#	distribution (um)
	# lowsize - lower size bin boundary (um)
	# uppsize - upper size bin boundary (um)
	# std - standard deviation for injected particle number size 
	#	distributions
	# rbou - size bin radius bounds (um)
	# self.cont_infl_t - times for constant influxes (s)
	# MV - molar volume of all components (cm3/mol)
	# self.rindx_aq - index of reactants for aqueous-phase 
	# self.pindx_aq - index of products for aqueous-phase
	# self.rstoi_aq - stoichiometry of reactants for aqueous-phase
	# self.pstoi_aq - stoichiometry of products for aqueous-phase
	# self.nreac_aq - number of reactants per aqueous-phase reaction
	# self.nprod_aq - number of products per aqueous-phase reaction
	# self.jac_stoi_aq - stoichiometry for Jacobian for 
	#	aqueous-phase
	# self.njac_aq - number of Jacobian elements per 
	#	aqueous-phase reaction 
	# self.jac_den_indx_aq - index of Jacobian denominators 
	# self.jac_indx_aq - index of  Jacobian for aqueous-phase
	# self.y_arr_aq - y indices for aqueous-phase
	# self.y_rind_aq - reactant indices for aqueous-phase 
	# self.uni_y_rind_aq - y indices for reactants for aqueous-phase
	# self.y_pind_aq - y indices for products for aqueous-phase
	# self.uni_y_pind_aq - y indices for products for aqueous-phase
	# self.reac_col_aq - columns of sparse matrix for aqueous-phase
	# self.prod_col_aq - columns of sparse matrix for aqueous-phase
	# self.rstoi_flat_aq - reactant stoichiometries for Jacobian for 
	# 	aqueous-phase
	# self.pstoi_flat_aq - product stoichiometries for Jacobian for
	#	aqueous-phase
	# self.rr_arr_aq - aqueous-phase reaction rate indices
	# self.rr_arr_p_aq - aqueous-phase reaction rate indices
	# self.eqn_num - number of reactions in gas- and aqueous-phase
	# self.partit_cutoff - the product of saturation vapour pressure
	#	and activity coefficient above which gas-particle
	#	partitioning assumed negligible (Pa)
	# diff_vol - diffusion volumes of components according to 
	#		Fuller et al. (1969)
	# DStar_org - diffusion coefficient of components at initial 
	#	temperature (cm2/s)
	# corei - index of core component
	# ser_H2O - whether to serialise the gas-particle partitioning 
	#	of water
	# self.C_p2w - concentration of components on the wall due to 
	#	particle deposition to wall (# molecules/cm3)
	# the following inputs are used only for the saving module:
	# self.sch_name - path to chemical scheme
	# sav_nam - name of folder to save in
	# self.comp_namelist - chemical scheme name of components
	# self.dydt_trak - name of components to track change tendencies
	# rbou00 - original particle size bin bounds
	# ub_rad_amp - amplificatin factor for upper bin size bound
	# indx_plot - indices of components to plot the gas-phase 
	#	temporal profile of
	# self.comp0 - names of components to plot the gas-phase temporal 
	#	profile of
	# self.inname - path to model variables file
	# self.rel_SMILES - SMILES strings of components in chemical 
	#	scheme
	# self.Psat_Pa_rec - pure component saturation vapour pressures 
	#	(Pa) at 298.15 K
	# self.Psat_Pa - pure component saturation vapour pressures (Pa) 
	#	at starting temperature in chamber
	# self.OC - oxygen to carbon ratio of components
	# wat_hist - flag for history of particle-phase with respect to 
	#	water partitioning,
	# 	where 0 is dry (therefore on the deliquescence curve) 
	#	and 1 is wet 
	#	(therefore on the efflorescence curve)
	# self.Pybel_objects - the pybel objects for components
	# self.pcont - flag for whether particle injection continuous or 
	#	instantaneous
	# self.dil_fac - chamber dilution factor (fraction of chamber/s)
	# NOi - NO index
	# HO2i - HO2 index
	# NO3i - NO3 index
	# z_prt_coeff - fraction of total gas-particle partitioning 
	#	coefficient 
	#	below which partitioning to a particle size bin is 
	#	treated as zero,
	#	e.g. because surface area of that size bin is tiny 
	# self.con_C_indx - index of components with constant 
	# 	gas-phase concentration
	# self.seed_eq_wat - whether seed particles to be equilibrated 
	#	with water prior to ODE solver
	# self.Vwat_inc - whether suppled seed particle volume contains 
	#	equilibrated water
	# tot_in_res - record of total input of injected components 
	#	(ug/m3)
	# Compti - index for total injection record for instantaneously 
	#	injected components
	# self.cont_inf_reci - index of components with continuous 
	#	influx in record
	# self.con_infl_indx - index of components with continuous 
	#	influx in concentration array
	# tot_in_res_indx - index of components with recorded influx
	# chamSA - chamber surface area (m2)
	# chamV - chamber volume (m3)
	# self.tf_UVC - transmission factor for 254 nm wavelength light
	# tempt_cnt # count on chamber temperatures
	# ------------------------------------------------------------
	
	# start timer
	if (self.spin_up == 0): # if no spin-up
		self.st_time = time.time()
	
	step_no = 0 # track number of time steps
	sumt = 0. # track time through simulation (s)
	self.sumt = 0. # track time through simulation (s)
	# counters on updates
	light_time_cnt = 0 # light time status count
	gasinj_cnt = 0 # count on injection times of components
	if (self.pconct[0, 0] == 0. and len(self.pconct[0, :]) > 1 and 
	self.pcont[0, 0] == 0):
		# count on injection times of particles
		seedt_cnt = 1
		self.seedx_tcnt = 1 
	else:
		seedt_cnt = 0
		self.seedx_tcnt = 0
	
	# current status of lights
	self.light_stat_now = self.light_stat[light_time_cnt]
	
	# current status of whether injection of particles instantaneous 
	# or continuous, if not stated assume instantaneous
	pcontf = 0
	if (self.pconct[0, 0] == 0 and self.pcont[0, 0] == 1):
		pcontf = 1
	infx_cnt = 0 # count on constant gas-phase influx occurrences
	infx_cnt0 = 0 # remember count at start of integration step
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
	y0[:] = y[:]
	# remember initial particle number concentrations (# particles/cm3)
	N_perbin0 = np.zeros((N_perbin.shape[0], N_perbin.shape[1]))
	N_perbin0[:, :] = N_perbin[:, :]
	x0 = np.zeros((len(x)))# remember initial particle sizes (um)
	x0[:] = x[:]
	t0 = self.update_stp # remember initial integration step (s)
	# flag for changing integration time step due to changing initial values	
	ic_red = 0
	tnew = self.update_stp # the time to integrate over (s)
	# fraction of newly injected seed particles
	pconcn_frac = 0.
	# numpy array version of chemical scheme names
	self.comp_namelist_np = np.array(self.comp_namelist)
	# turn off flag for ongoing injection of particles
	self.pcont_ongoing = 0

	# find out what to do with the gas-wall partitioning 
	# coefficient,
	# note that self.kw and self.Cw are spread over wall bins in 
	# rows and components 
	# in columns (the latter spread is done in partit_var_prep.py)
	if (self.wall_on > 0):
		if (sum(sum(self.kw == -1)) > 0 ):
			self.kwf = -1 # Huang et al. 2018 treatment
		else:
			# standard PyCHAM (GMD paper) treatment with 
			# same gas-wall 
			# partitioning coefficient for all components
			self.kwf = 0
	else:
		self.kwf = 0 # filler when no wall
	
	# prepare recording matrices, including recording of initial
	# conditions, note initial change tendencies not recorded 
	# in this call but are recorded below
	[trec, yrec, Cfactor_vst, Nres_dry, Nres_wet, x2, 
	seedt_cnt, rbou_rec, Cfactor, infx_cnt, 
	temp_now, cham_env, Pnow, 
	RHn, Cinfl_now] = rec_prep.rec_prep(nrec_steps, y, y0, 
	num_sb, num_comp, N_perbin, core_diss, mfp,
	accom_coeff, y_mw, surfT, R_gas, NA,
	x, therm_sp, H2Oi, act_coeff,
	sumt, Pnow, light_time_cnt, 
	Jlen, Cfactor, 
	Vbou, tnew, 
	np_sum, update_count, injectt, gasinj_cnt, 
	inj_indx, Ct, seedt_cnt, corei, 
	lowsize, uppsize, x, std, rbou, 
	infx_cnt, MV, diff_vol, DStar_org, 
	tempt_cnt, RHt_cnt, nuci, 
	t0, pcontf, NOi, HO2i, NO3i, z_prt_coeff,
	tot_in_res, Compti, 
	tot_in_res_indx, chamSA, chamV, wat_hist, self, vol_Comp, volP)
	
	import ode_solv
	import dydt_rec
	importlib.reload(ode_solv) # import most recent version
	importlib.reload(ode_solv_wat) # import most recent version
	importlib.reload(dydt_rec) # import most recent version
	
	while (self.tot_time-sumt) > (self.tot_time/1.e10):

		# remembering variables at the start of the 
		# integration step -----------------------------------
		# remember initial concentrations (# molecules/cm3 (air))
		y0[:] = y[:]
		# remember initial particle number concentration (# particles/cm3)
		N_perbin0[:] = N_perbin[:]
		x0[:] = x[:] # remember initial particle sizes (um)
		temp_now0 = temp_now # remember temperature (K)
		# remember water history flag at start of integration step
		wat_hist0 = wat_hist
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
		
		# ------------------------------------------------------------
		
		# flag for stability in gas-particle partitioning for solver while loop
		gpp_stab = 0
		# flag for stability in gas-particle partitioning for time interval reset
		stab_red = 0
		lin_int = 0 # flag to linearly interpolate changes to chamber
		t00 = tnew # remember the initial integration step for this integration step (s)
		save_cntf = 0 # flag for updating count on number of recordings
		
		while (gpp_stab != 1): # whilst ode solver flagged as unstable

			# if integration interval decreased, reset 
			# concentrations to those at start of interval
			if (gpp_stab == -1):
				y[:] = y0[:] # (# molecules/cm3)
			
			# for change tendencies, t=0 recording done
			#  inside rec_prep
			# record any change tendencies of specified 
			# components after t=0
			if (len(self.dydt_vst) > 0 and save_cntf == 0):
				if ((sumt-(self.save_step*(save_cnt-1))
					 > -1.e-10)):
					
					if (sumt-
					(self.save_step*(save_cnt-1)) > 
					-1.e-10):
						dydt_cnt = save_cnt-1

					# before solving ODEs for 
					# chemistry, dilution, 
					# gas-particle partitioning and
					# gas-wall partitioning, 
					# estimate and record any 
					# change tendencies (# 
					# molecules/cm3/s) resulting 
					# from 
					# these processes
					if (self.testf != 5):
						self = dydt_rec.dydt_rec(y, rrc, dydt_cnt, 
							num_sb, 
							num_comp, core_diss, kelv_fac, kimt, 
							act_coeff, dydt_erh_flag, H2Oi, 
							wat_hist, self)
						
			# record output if on the first attempt at solving this time interval, 
			# note that recording here in this way means we include any
			# instantaneous changes at this time step without interpolation to
			# smaller instantaneous changes (when interpolation forced due to 
			# instability)
			# note also that recording before calling cham_up means that 
			# any instantaneous 
			# changes occurring during this upcoming time step are not 
			# recorded at the very start 
			# of the time step
			if (save_cntf == 0 and (sumt-(self.save_step*(save_cnt-1)) > 
				-1.e-10) and self.testf != 5):
				
				[trec, yrec, Cfactor_vst, save_cntf, Nres_dry, Nres_wet, 
				x2, rbou_rec, cham_env] = rec.rec(save_cnt-1, trec, 
				yrec, Cfactor_vst, y, sumt, num_sb, num_comp, N_perbin, 
				core_diss, kelv_fac, kimt, act_coeff, Cfactor, Nres_dry, 
				Nres_wet, x2, x, MV, H2Oi, Vbou, rbou, rbou_rec, 
				cham_env, temp_now, Pnow, tot_in_res, self)
				# prepare for recording next point
				save_cnt += 1
			
			
			

			# update chamber variables
			[temp_now, Pnow, light_time_cnt, tnew, ic_red, 
			update_count, Cinfl_now, seedt_cnt, Cfactor, 
			infx_cnt, 
			gasinj_cnt, DStar_org, y, tempt_cnt, RHt_cnt, 
			N_perbin, x,
			pconcn_frac,  pcontf, tot_in_res, 
			self] = cham_up.cham_up(sumt, 
			Pnow0, light_time_cnt0, 
			tnew, np_sum, update_count, 
			injectt, gasinj_cnt0, inj_indx, Ct,
			seedt_cnt0, num_comp, y0, y, N_perbin0, 
			corei, 
			lowsize, uppsize, num_sb, MV, x0, std, 
			H2Oi, rbou, 
			infx_cnt0, Cfactor, diff_vol, 
			DStar_org, tempt_cnt0, RHt_cnt0, nuci,
			y_mw, temp_now0, gpp_stab, t00, x0,
 			pcontf, Cinfl_now, surfT,
			act_coeff, tot_in_res, Compti, self, vol_Comp, 
			volP)

			# aligning time interval with pre-requisites --		
			# ensure end of time interval does not surpass 
			# recording time
			if ((sumt+tnew) > self.save_step*(save_cnt-1)):
				tnew = (self.save_step*(save_cnt-1))-sumt
				# temporarily set the update step for operator-split processes
				# to align with recording time step, this ensures that
				# recording and operator-split intervals don't fall out of sync
				self.update_stp = tnew
				update_count = 0.
				ic_red = 1

			# ensure update to operator-split processes interval not surpassed
			if (update_count+tnew > self.update_stp):
				tnew = (self.update_stp-update_count)
				ic_red = 1
			
			# ensure simulation end time not surpassed
			if (sumt+tnew > self.tot_time):
				tnew = (self.tot_time-sumt)
				ic_red = 1
			
			# ------------------------------------------------------------
			# if particles and/or wall present		
			if ((num_sb-self.wall_on) > 0 or self.wall_on > 0):
				
				# update partitioning variables
				[kimt, kelv_fac] = partit_var.kimt_calc(y, mfp, num_sb, 
				num_comp, 
				accom_coeff, y_mw, surfT, R_gas, temp_now, NA, N_perbin, 
				x.reshape(1, -1)*1.e-6, therm_sp, H2Oi, act_coeff, 1,
				Pnow, DStar_org, z_prt_coeff, chamSA, chamV, self)
				
				# update particle-phase activity coefficients, note the output,
				# note that if ODE solver unstable, then y resets to y0 via
				# the cham_up module prior to this call
				[act_coeff, wat_hist, RHn, y, 
				dydt_erh_flag] = act_coeff_update.ac_up(y, H2Oi, RH0, temp_now, 
				wat_hist0, act_coeff, num_comp, (num_sb-self.wall_on))
				
			else: # fillers
			
				kimt = np.zeros((num_sb+self.wall_on, 
					num_comp))
				kelv_fac = np.zeros((
					num_sb-self.wall_on, 1))
				dydt_erh_flag = 0
			
			# reaction rate coefficient going into this 
			# time step
			[rrc, erf, err_mess] = rrc_calc.rrc_calc(
				y[H2Oi], temp_now, y, 
				Pnow, Jlen, y[NOi], y[HO2i], y[NO3i], 
				sumt, self)

			# if error message from reaction rate 
			# calculation
			if (erf == 1): 
				yield(err_mess)
			
			# update Jacobian inputs based on 
			# particle-phase fractions of components
			[rowvalsn, colptrsn, jac_mod_len, 
			jac_part_hmf_indx, rw_indx, 
				jac_part_H2O_indx] = jac_up.jac_up(
				y[num_comp:num_comp*
				((num_sb-self.wall_on+1))], rowvals, 
				colptrs, (num_sb-self.wall_on), num_comp, 
				H2Oi, y[H2Oi], ser_H2O, self)
			
			
			# if water gas-particle partitioning serialised
			if (ser_H2O == 1 and (num_sb-self.wall_on) > 0
				 and (sum(N_perbin) > 0)): 
				
				# if on the deliquescence curve rather 										# than the 
				# efflorescence curve in terms of water 									# gas-particle partitioning
				if (wat_hist == 1):
					# flag that water gas-particle 
					# partitioning solved separately
					self.odsw_flag = 1

					# call on ode solver for water
					[y, res_t] = ode_solv_wat.ode_solv(y, 
					tnew,
					Cinfl_now, rowvalsn, colptrsn, 
					num_comp, 
					num_sb, act_coeff, core_diss, 
					kelv_fac, kimt, 
					(num_sb-self.wall_on), 
					jac_mod_len, jac_part_hmf_indx,
 					rw_indx, N_perbin, 
					jac_part_H2O_indx, H2Oi, self)

					
					
					# check on stability of water 
					# partitioning	
					if (any(y[H2Oi::num_comp] < 0.)): 

						# identify components with negative 
						# concentrations
						neg_comp_indx = y < 0.
						# transform into components in columns, 
						# locations in rows
						neg_comp_indx = neg_comp_indx.reshape(
						num_sb+1, num_comp)
						# get component indices with negative 
						# concentration
						neg_comp_indx = np.unique((
						np.where(neg_comp_indx == 1))[1])
						# get chemical scheme names of components 
						# with negative concentration
						neg_names = self.comp_namelist_np[neg_comp_indx]

						# isolate just water concentrations 
						# (molecules/cm3)
						y_H2O = y[H2Oi::num_comp]
						# sum the negative concentrations and 
						# convert to absolute value (molecules/cm3)
						neg_H2O = np.abs(sum(y_H2O[y_H2O<0.]))
						
						# allow a given fraction of water 
						# concentrations to be negative
						if (neg_H2O/sum(
						np.abs(y[H2Oi::num_comp])) > 0. ):
				
							gpp_stab = -1 # maintain unstable flag
							# tell user what's happening
							yield (str('Note: negative water concentration generated following call to ode_solv_wat module, the program assumes this is because of a change in relative humidity in chamber air, and will automatically half the integration time interval and linearly interpolate any change to chamber conditions supplied by the user.  To stop this the simulation must be cancelled using the Quit button in the PyCHAM graphical user interface.  Current update time interval is ' + str(tnew) + ' seconds'))
							
							if (tnew < 1.e-20): # if time step has decreased to unreasonably low and solver still unstable then break
								ode_brk_err_mess.ode_brk_err_mess(y0, neg_names, rrc, num_comp, 
									(num_sb-self.wall_on), act_coeff, neg_comp_indx, 
									N_perbin, core_diss, kelv_fac, kimt, 0, H2Oi, y, self)

								yield (str('Error: negative concentrations generated following call to ode_solv_wat module, the program has assumed this is because of a change in chamber condition (e.g. injection of components), and has automatically halved the integration time interval and linearly interpolated any change to chamber conditions supplied by the user.  However, the integration time interval has now decreased to ' + str(tnew) + ' seconds, which is assumed too small to be useful, so the program has been stopped.  The components with negative concentrations are : ' + str(neg_names) + '.  The problem could be too stiff for the solver and the relevant fluxes (change tendencies) have been output to the file ODE_solver_break_relevant_fluxes.txt for your analysis of problem stiffness.  You could identify the maximum and minimum fluxes to gain indication of the components and/or processes making the problem stiff.  Therefafter you could modify the relevant model variables (supplied by the user) and the chemical scheme (supplied by the user).' ))
							# half the update and 
							# integration time step (s) 
							# if necessary
							tnew = tnew/2.
							stab_red = 1 # remember that time step temporarily reduced due to instability
							continue

						else: # if acceptable
							gpp_stab = 1 # change to stable flag
						
					else: # if solution stable, change stability flag to represent this
						gpp_stab = 1 # change to stable flag
				else:
					# water gas-particle partitioning not solved separately
					self.odsw_flag = 0
				# zero partitioning of water to particles 
				# for integration without 
				# water gas-particle partitioning
				# if particles present
				if (num_sb > self.wall_on):
					kimt[0:num_sb-self.wall_on, 
						H2Oi] = 0.
			
			else:
				# water gas-particle partitioning not 
				# solved separately
				self.odsw_flag = 0

			

			# model component concentration changes to 
			# get new concentrations molecules/cm3 (air))
			[y, res_t] = ode_solv.ode_solv(y, tnew, rrc,
				Cinfl_now, rowvalsn, colptrsn, num_comp, 
				num_sb, act_coeff,
				core_diss, kelv_fac, kimt, 
				(num_sb-self.wall_on),
				jac_mod_len, jac_part_hmf_indx, rw_indx, 
				N_perbin, 
				jac_part_H2O_indx, H2Oi, self)
			
			
		
			# if any components set to have constant 
			# gas-phase concentration
			if (any(self.con_C_indx)): # then keep constant
				y[self.con_C_indx] = y0[self.con_C_indx] # (# molecules/cm3)
			
			# if negative, suggests ODE solver instability, 
			# but could also be numerical 
			# limits, especially if concentrations are
			# relatively 
			# close to zero, so allow 
			# some leeway
			if (any(y/np.sum(np.abs(y))<-1.e-30)):
			
				# identify components with negative concentrations
				neg_comp_indx = y < 0.
				# transform into components in columns, locations in rows
				neg_comp_indx = neg_comp_indx.reshape(num_sb+1, num_comp)
				# get component indices with negative concentration
				neg_comp_indx = np.unique((np.where(neg_comp_indx == 1))[1])
				# get chemical scheme names of components with 
				# negative concentration
				neg_names = self.comp_namelist_np[neg_comp_indx]				

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
				
				# if time step has decreased to unreasonably 
				# low and solver still unstable then break
				if (tnew < 1.e-20):
					# estimate gas-phase reaction fluxes 
					# for all reactions and partitioning 
					# fluxes for troublesome components
					ode_brk_err_mess.ode_brk_err_mess(y0, neg_names, rrc, 
						num_comp, (num_sb-self.wall_on), act_coeff, 
						neg_comp_indx, N_perbin, core_diss, kelv_fac, 
						kimt, 1, H2Oi, y, self)

					yield (str('Error: negative concentrations generated following call to ode_solv module, the program has assumed this is because of a change in chamber condition (e.g. injection of components), and has automatically halved the integration time interval and linearly interpolated any change to chamber conditions supplied by the user.  However, the integration time interval has now decreased to ' + str(tnew) + ' seconds, which is assumed too small to be useful, so the program has been stopped.  The components with negative concentrations are : ' + str(neg_names) + '.  The problem could be too stiff for the solver and the relevant fluxes (change tendencies) have been output to the file ODE_solver_break_relevant_fluxes.txt for your analysis of problem stiffness.  You could identify the maximum and minimum fluxes to gain indication of the components and/or processes making the problem stiff.  Therefafter you could modify the relevant model variables (supplied by the user) and the chemical scheme (supplied by the user).' ))
						
				# half the update and integration time step (s) if necessary	
				tnew = tnew/2.
				# remember that time step temporarily 
				# reduced due to instability
				stab_red = 1 
				
			else: # if solution stable, change stability flag to represent this
				
				# account for any partial addition of 
				# newly injected seed particles
				self.pconc[:, seedt_cnt] -= self.pconc[:, seedt_cnt]*pconcn_frac
				# reset fraction of newly injected seed particles
				pconcn_frac = 0.
				gpp_stab = 1 # change to stable flag
		
		
		# end of integration stability condition section ------
		step_no += 1 # track number of steps
		sumt += tnew # total time through simulation (s)
		self.sumt += tnew

		# if any gas-phase components constrained 
		# to observations
		if (any(self.obs_comp_i)):
				# get observed concentrations now
				# loop through components
				for ci in range(len(self.obs_comp_i)):
					y[self.obs_comp_i[ci]] = np.interp(
					sumt, self.obs[:, 0], self.obs[:, ci+1])
				
		# dilute particle number following an integration 
		# time step, e.g. for flow-reactor -----
		# note that concentrations of components inside 
		# particles (and in the gas-phase) will have been
		# reduced due to dilution inside the ODE solver
		# note self.pp_dil set in def_mod_var and obs_file_open
		if (self.dil_fac_now > 0 and self.pp_dil == 1):
			N_perbin -= N_perbin*(self.dil_fac_now*tnew)
		
		# if particle size bins present, rebin
		if ((num_sb-self.wall_on) > 0):

			# update particle sizes
			# if multiple particle size bins present 
			# containing particles
			if (((num_sb-self.wall_on) > 1) and 
			(any(N_perbin > 1.e-10))):

				if (siz_str == 0): # moving centre
					(N_perbin, Varr, y, x, redt, t, bc_red) = mov_cen.mov_cen_main(N_perbin, 
					Vbou, num_sb, num_comp, y_mw, x, Vol0, tnew, 
					y0, MV, ic_red, y, res_t, self)
					
				if (siz_str == 1): # full-moving
					(Varr, x, y[num_comp:(num_comp*(num_sb-self.wall_on+1))], 
					N_perbin, Vbou, rbou) = fullmov.fullmov((num_sb-self.wall_on), N_perbin,
 					num_comp, y[num_comp:(num_comp)*(num_sb-self.wall_on+1)], MV*1.e12, 
					Vol0, Vbou, rbou)
			
			# time since operator-split processes 
			# last called (s)
			update_count += tnew
			
			# if time met to implement operator-split 
			# processes
			if (update_count >= (
			self.update_stp*9.999999e-1)):

				if (any(N_perbin > 1.e-10)):
				
					# particle-phase 
					# concentration(s) 
					# (# molecules/cm3)
					Cp = np.transpose(y[num_comp:(num_comp)*(
					num_sb-self.wall_on+1)].reshape(
					num_sb-self.wall_on, num_comp))
					
					# coagulation
					[N_perbin, y[num_comp:(num_comp)*(num_sb-self.wall_on+1)], x, Gi, eta_ai, 
						Varr, Vbou, rbou] = coag.coag(self.RH[RHt_cnt], temp_now, x*1.e-6, 
						(Varr*1.0e-18).reshape(1, -1), 
						y_mw.reshape(-1, 1), x*1.e-6, 
						Cp, (N_perbin).reshape(1, -1), update_count, 
						(Vbou*1.0e-18).reshape(1, -1), rbou,
						num_comp, 0, Vol0, rad0, Pnow, 0,
						Cp, (N_perbin).reshape(1, -1),
						(Varr*1.e-18).reshape(1, -1),
						coag_on, siz_str, self)
										

					# if particle loss to walls turned on, 
					# account for this now
					if ((McMurry_flag > -1) and (self.wall_on > 0)):
						[N_perbin, y[num_comp:(num_comp)*(
						num_sb-self.wall_on+1)]] = wallloss.wallloss(
							N_perbin.reshape(-1, 1), 
							y[num_comp:(num_comp)*(
							num_sb-self.wall_on+1)], Gi, eta_ai,
 							x*2.e-6, y_mw, 
							Varr*1.e-18, (num_sb-self.wall_on),
							num_comp, temp_now, update_count, 
							inflectDp, pwl_xpre, pwl_xpro, inflectk,
							chamR, McMurry_flag, 
							0, p_char, e_field, 
							(num_sb-self.wall_on), self)
			

				if (self.nucv1 > 0.): # nucleation
					
					[N_perbin, y, x, Varr, np_sum, 
					rbou, Vbou] = nuc.nuc(sumt,
					np_sum, N_perbin, y, 
					y_mw.reshape(-1, 1),  
					num_comp, Varr, x, new_partr, 
					MV, siz_str, rbou, Vbou, 
					(num_sb-self.wall_on), self)
				
				# reset count that tracks when next 
				# operator-split should be 
				# called (s)
				update_count = 0.
		
		# update the percentage time through simulation 
		# in the GUI progress bar
		yield (sumt/self.tot_time*100.)
		

		# if ozone isopleth being made, then store ozone result
		if (self.testf == 5):
			
			if (sumt + t0 >= self.tot_time): # ensure simulation doesn't end
				self.tot_time += t0*2.

			# keep [VOC] constant
			y[self.VOCi] = y0[self.VOCi]
			# keep [NO]+[NO2] constant
			y_NO_reset = (y[self.NOi]/(y[self.NOi]+y[self.NO2i]))*(y0[self.NOi]+y0[self.NO2i])
			y_NO2_reset = (y[self.NO2i]/(y[self.NOi]+y[self.NO2i]))*(y0[self.NOi]+y0[self.NO2i])
			y[self.NOi] = y_NO_reset
			y[self.NO2i] = y_NO2_reset
			
			
			# change to ozone concentration between start 
			# and finish of this integration step
			O3_changen = (y[self.O3i]-self.O3equil)
			
			dfracn = np.abs(O3_changen)/y[self.O3i] # new fractional change
			
			# remember this [O3] (before integration 
			# generates the next new attempt)
			self.O3equil = y[self.O3i]

			# if O3 concentration still changing significantly (not at equilibrium)
			if ((np.abs(O3_changen)/y[self.O3i])/tnew > 1.e-4 or step_no < 2 or sumt < 1800. or dfracn > dfrac0):
				
				if (step_no > 1):

					# new ozone concentration to begin integration with
					if (np.abs(O3_changen)/y[self.O3i] <= 3.e-4 and (np.abs(dfracn-dfrac0)/dfrac0) < 0.05): # if iteration going in right direction
				
						y[self.O3i] += O3_changen*((np.abs(O3_changen)/y[self.O3i])*1.e4)

					elif (np.abs(O3_changen)/y[self.O3i] < 0.05): # if change by integration is modest
						y[self.O3i] += O3_changen/2.
					elif (np.abs(O3_changen)/y[self.O3i] >= 0.05): # if change by integration is large
						y[self.O3i] = y[self.O3i]
				else:
						y[self.O3i] += O3_changen/2.

				# remember this fraction change
				dfrac0 = dfracn
				
				
			else: # if O3 close enough to equilibrium
				# remember new ozone concentration 
				# for next set of [NOx] and 
				# [VOC] values
				self.O3equil = y[self.O3i]
				
				return() # end this call to simulation

		# record output at experiment end
		if (sumt >= (self.tot_time-self.tot_time/1.e10) and self.testf != 5):

			# update chamber variables, note this 
			# ensures that any changes made
			# coincidentally with the experiment end 
			# are captured
			[temp_now, Pnow, light_time_cnt, tnew, ic_red, 
			update_count, Cinfl_now, seedt_cnt, Cfactor, infx_cnt, 
			gasinj_cnt, DStar_org, y, tempt_cnt, RHt_cnt, N_perbin, x,
			pconcn_frac,  pcontf, tot_in_res, 
			self] = cham_up.cham_up(sumt, 
			Pnow0, light_time_cnt0, 
			tnew, np_sum, 
			update_count, 
			injectt, gasinj_cnt0, inj_indx, Ct,
			seedt_cnt0, num_comp, y0, y, N_perbin0, 
			corei, 
			lowsize, uppsize, num_sb, MV, x0, std, 
			H2Oi, rbou, 
			infx_cnt0, Cfactor, diff_vol, 
			DStar_org, tempt_cnt0, RHt_cnt0, nuci,
			y_mw, temp_now0, gpp_stab, t00, x0,  
			pcontf, Cinfl_now, surfT,
			act_coeff, tot_in_res, Compti, self, vol_Comp, 
			volP)
			
			[trec, yrec, Cfactor_vst, save_cnt, Nres_dry, 
			Nres_wet,
			x2, rbou_rec, cham_env] = rec.rec(save_cnt-1, 
			trec, yrec, Cfactor_vst, y, sumt, num_sb, 
			num_comp, N_perbin, core_diss, 
			kelv_fac, kimt, act_coeff, Cfactor, Nres_dry, 
			Nres_wet, x2, x, MV, H2Oi, Vbou, rbou, 
			rbou_rec, cham_env, temp_now, Pnow, 
			tot_in_res, self)
		
		# if time step was temporarily reduced, then reset
		if (ic_red == 1 or stab_red == 1):
			self.update_stp = t0
			tnew = self.update_stp
			ic_red = 0 # reset flag
			stab_red = 0 # reset flag
		
		# remember the gas-phase water concentration from 
		# previous integration step (# molecules/cm3)
		y_H2O0 = y[H2Oi]
			
	time_taken = time.time()-self.st_time
	
	# re-merge continuous influx of water with that of other 
	# components, this will allow further commands from the GUI
	if (self.H2Oin == 1):

		# index
		self.con_infl_indx = np.concatenate((self.con_infl_indx, np.array((H2Oi)).reshape(1)))
		# influx rate
		self.con_infl_C = np.concatenate((self.con_infl_C, self.con_infl_H2O), axis=0)
	
	# save results
	save.saving(yrec, Nres_dry, Nres_wet, trec, sav_nam, 
		num_comp, Cfactor_vst, 0, num_sb, y_mw, MV, time_taken, 
		x2, rbou_rec, rbou00, ub_rad_amp, 
		indx_plot, H2Oi, siz_str, cham_env, self)
	
	return() # end of function
