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
import mov_cen
import fullmov
import wallloss
import nuc
import coag

def ode_updater(update_stp, 
	tot_time, save_stp, y, rindx, 
	pindx, rstoi, pstoi, nreac, nprod, jac_stoi, njac, 
	jac_den_indx, jac_indx, RO2_indx, H2Oi,	temp, tempt, 
	Pnow, light_stat, light_time, daytime, lat, lon, af_path, 
	dayOfYear, photo_path, Jlen, con_infl_C, nrec_steps, 
	dydt_vst, siz_str, num_sb, num_comp, corei, core_diss, Psat, 
	mfp, therm_sp,
	accom_coeff, y_mw, surfT, R_gas, NA, y_dens, DStar_org, 
	x, Varr, act_coeff, Cw, kw, Cfactor, tf, light_ad, y_arr, 
	y_rind, 
	uni_y_rind, y_pind, uni_y_pind, reac_col, prod_col, 
	rstoi_flat, pstoi_flat, rr_arr, rr_arr_p, rowvals, 
	colptrs, wall_on, jac_wall_indx, jac_part_indx, Vbou,
	N_perbin, Vol0, rad0, np_sum, new_partr, nucv1, nucv2, 
	nucv3, nuc_comp, nuc_ad, RH, coag_on, inflectDp, pwl_xpre, 
	pwl_xpro, inflectk, chamR, Rader, p_char, e_field, 
	injectt, inj_indx, Ct, pconc, pconct, mean_rad, lowsize, 
	uppsize, std, rbou, const_infl_t):

	import ode_solv # import most updated version
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
	# H2Oi - index of water
	# temp - temperature in chamber (K)
	# tempt - times that temperatures reached (s)
	# Pnow - pressure inside chamber (Pa)
	# light_stat - lights status
	# light_time - time light status attained (s)
	# daytime - time of day experiment starts (s)
	# lat - latitude of experiment (degrees)
	# lon - longitude of experiment (degrees)
	# af_path - path to actinic flux values
	# dayOfYear - number of days since 31st December experiment held on
	# photo_path - path to file with absorption cross-sections
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
	# corei - index of seed material
	# core_diss - dissociation constant of seed
	# Psat - pure component saturation vapour pressure 
	# 	(molecules/cc (air))
	# mfp - mean free path (m)
	# accom_coeff - accommodation coefficient
	# y_mw - molecular weight (g/mol)
	# surfT - surface tension (g/s2)
	# R_gas - ideal gas constant (kg.m2.s-2.K-1.mol-1)
	# NA - Avogadro's constant (molecules/mol)
	# y_dens - component densities (kg/m3)
	# DStar_org - gas-phase diffusion coefficient (m2/s)
	# x - particle radii (um)
	# Varr - particle volume (um3)
	# therm_sp - thermal speed (m/s)
	# act_coeff - activity coefficient
	# Cw - effective absorbing mass of wall (molecules/cc (air))
	# kw - gas-wall mass transfer coefficient (/s)
	# Cfactor - conversion factor for concentrations (ppb/molecules/cc)
	# tf - transmission factor for natural sunlight
	# light_ad - marker for whether to adapt time interval for 
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
	# Cw - effective absorbing mass concentration of wall 
	#	(molecules/cc (air))
	# kw - mass transfer coefficient of wall partitioning (/s)
	# jac_wall_indx - index of inputs to Jacobian from wall 
	# 	partitioning
	# jac_part_indx - index of inputs to Jacobian from particle
	#	partitioning
	# Vbou - volume boundary of size bins (um3)
	# N_perbin - number concentration of particles per size bin 
	#	(#/cc (air))
	# Vol0 - initial single particle volumes per size bin (um3)
	# rad0 - initial radius at particle centres (um)
	# np_sum - number concentration of newly nucleated particles 
	#		(#/cc (air))
	# new_partr - radius of newly nucleated particles (cm)
	# nucv1, v2, v3 - nucleation parameters
	# nuc_comp - the nucleating component
	# nuc_ad - marker for whether to reduce time step to allow 
	#	for accurate capture of nucleation
	# RH - relative humidity
	# coag_on - whether coagulation to be modelled
	# inflectDp - particle diameter at which wall loss inflection occurs (m)
	# pwl_xpre - x value preceding inflection point
	# pwl_xpro - x value proceeding inflection point
	# inflectk - deposition rate at inflection (/s)
	# chamR - spherical-equivalent radius of chamber (m2)
	# Rader - marker for treamnt of particle deposition to walls
	# p_char - average number of charges per particle (/particle)
	# e_field - average electric field inside chamber (g.m/A.s3)
	# injectt - time of injection of components (s)
	# inj_indx - index of component(s) being injected after 
	#	experiment start
	# Ct - concentration(s) (ppb) of component(s) injected 
	#	instantaneously after experiment start
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
	# ------------------------------------------------------------
	
	step_no = 0 # track number of time steps
	sumt = 0.0 # track time through simulation (s)
	# counters on updates
	light_time_cnt = 0 # light time status count
	gasinj_cnt = 0 # count on injection times of components
	if pconct[0, 0] == 0.: 
		seedt_cnt = 1 # count on injection times of particles
	else:
		seedt_cnt = 0
	infx_cnt = 0 # count on constant gas-phase influx occurrences
	save_cnt = 1 # count on recording results
	# count on time since update to integration initial values/constants last called (s)
	update_count = 0.0
	y0 = np.zeros((len(y))) # remember initial concentrations (molecules/cc (air))
	t0 = update_stp # remember initial integration step (s)
	# flag for changing integration time step due to changing initial values	
	ic_red = 0
	tnew = update_stp # the time to integrate over (s)

	# prepare recording matrices, including recording of initial
	# conditions
	[trec, yrec, dydt_vst, Cfactor_vst, Nres_dry, Nres_wet, x2, MV, seedt_cnt, rbou_rec] = rec_prep.rec_prep(nrec_steps, 
	y, rindx, 
	rstoi, pindx, pstoi, nprod, dydt_vst, nreac, 
	num_sb, num_comp, N_perbin, core_diss, Psat, mfp,
	accom_coeff, y_mw, surfT, R_gas, temp, tempt, NA,
	y_dens*1.e-3, DStar_org, x, therm_sp, H2Oi, act_coeff,
	RO2_indx, sumt, Pnow, light_stat, light_time, 
	light_time_cnt, daytime, lat, lon, af_path, 
	dayOfYear, photo_path, Jlen, Cw, kw, Cfactor, tf, 
	light_ad, wall_on, Vbou, tnew, nuc_ad, nucv1, nucv2, nucv3, 
	np_sum, update_stp, update_count, injectt, gasinj_cnt, 
	inj_indx, Ct, pconc, pconct, seedt_cnt, mean_rad, corei, 
	lowsize, uppsize, rad0, x, std, rbou, const_infl_t, 
	infx_cnt, con_infl_C)
	
	print('Starting loop through update steps')	
	while (tot_time-sumt)>(tot_time/1.e10):
		
		y0[:] = y[:] # remember initial concentrations (molecules/cc (air))
		# update chamber variables
		[temp_now, Pnow, lightm, light_time_cnt, tnew, ic_red, update_stp, 
			update_count, Cinfl_now, seedt_cnt] = cham_up.cham_up(sumt, temp, tempt, 
			Pnow, light_stat, light_time, light_time_cnt, light_ad, 
			tnew, nuc_ad, nucv1, nucv2, nucv3, np_sum, 
			update_stp, update_count, lat, lon, dayOfYear, photo_path, 
			af_path, injectt, gasinj_cnt, inj_indx, Ct, pconc, pconct, 
			seedt_cnt, num_comp, y, N_perbin, mean_rad, corei, lowsize, 
			uppsize, num_sb, MV, rad0, x, std, y_dens, H2Oi, rbou, 
			const_infl_t, infx_cnt, con_infl_C, wall_on)
		
		# ensure end of time interval does not surpass recording time
		if ((sumt+tnew)>save_stp*save_cnt):
			tnew = (save_stp*save_cnt)-sumt
			ic_red = 1
		
		# ensure update step interval not surpassed
		if (update_count+tnew>update_stp):
			tnew = (update_stp-update_count)
			ic_red = 1
		
		if (num_sb-wall_on)>0: # if particles present
			# update partitioning variables
			[kimt, kelv_fac] = partit_var.kimt_calc(y, mfp, num_sb, num_comp, accom_coeff, y_mw,   
			surfT, R_gas, temp_now, NA, y_dens, N_perbin, DStar_org, 
			x.reshape(1, -1)*1.0e-6, Psat, therm_sp, H2Oi, act_coeff, wall_on, 1)
			
		else: # fillers
			kimt = kelv_fac = 0.
		
		# reaction rate coefficient
		
		rrc = rrc_calc.rrc_calc(RO2_indx, 
			y[H2Oi], temp_now, lightm, y, daytime+sumt, 
			lat, lon, af_path, dayOfYear, Pnow, 
			photo_path, Jlen, tf)
		
		
		# model component concentration changes to get new concentrations
		# (molecules/cc)
		[res, res_t] = ode_solv.ode_solv(y, tnew, rindx, pindx, rstoi, pstoi,
			nreac, nprod, rrc, jac_stoi, njac, jac_den_indx, jac_indx,
			Cinfl_now, y_arr, y_rind, uni_y_rind, y_pind, uni_y_pind, 
			reac_col, prod_col, rstoi_flat, 
			pstoi_flat, rr_arr, rr_arr_p, rowvals, colptrs, num_comp, 
			num_sb, wall_on, Psat, Cw, act_coeff, kw, jac_wall_indx,
			corei, core_diss, kelv_fac, kimt, (num_sb-wall_on), 
			jac_part_indx)
		
		# take last installment from res
		y = res[-1, :]

		step_no += 1 # track number of steps
		sumt += tnew # total time through simulation (s)
		update_count += tnew # time since operator-split processes last called
		
		if num_sb-wall_on>0:
			# update particle sizes
			if (num_sb-wall_on)>1 and (N_perbin>1.0e-10).sum()>0: # if particles present
				if (siz_stru == 0): # moving centre
					(N_perbin, Varr, y, x, redt, t, bc_red) = mov_cen.mov_cen_main(N_perbin, 
					Vbou, num_sb, num_comp, y_mw, x, Vol0, tnew, 
					update_stp, y0, MV, Psat[0, :], ic_red, res, res_t)
				if (si_stru == 1): # full-moving
					(Varr, x, y[num_comp:(num_comp*(num_sb-wall_on+1))], 
					N_perbin, Vbou, rbou) = fullmov.fullmov((num_sb-wall_on), N_perbin,
 					num_comp, y, MV, Vol0, Vbou, rbou)

			# if time met to implement operator-split processes
			if (update_count>=update_stp*9.999999e-1):
				if ((N_perbin>1.e-10).sum()>0):
					# particle-phase concentrations (molecules/cc (air))
					Cp = np.transpose(y[num_comp:(num_comp)*(num_sb-wall_on+1)].reshape(num_sb-wall_on, num_comp))
					# coagulation
					[N_perbin, y[num_comp:(num_comp)*(num_sb-wall_on+1)], x, Gi, eta_ai, 
						Varr, Vbou, rbou] = coag.coag(RH, temp_now, x*1.0e-6, 
						(Varr*1.0e-18).reshape(1, -1), 
						y_mw.reshape(-1, 1), x*1.0e-6, 
						Cp, (N_perbin).reshape(1, -1), update_count, 
						(Vbou*1.0e-18).reshape(1, -1), rbou,
						num_comp, 0, (np.squeeze(y_dens*1.0e-3)), Vol0, rad0, Pnow, 0,
						Cp, (N_perbin).reshape(1, -1), (Varr*1.0e-18).reshape(1, -1),
						coag_on, siz_stru)
				
					if ((Rader>-1) and (wall_on == 1)): #if particle loss to walls turned on
						# particle loss to walls
						[N_perbin, 
						y[num_comp:(num_comp)*(num_sb-wall_on+1)]] = wallloss.wallloss(
							N_perbin.reshape(-1, 1), 
							y[num_comp:(num_comp)*(num_sb-wall_on+1)], Gi, eta_ai,
 							x*2.0e-6, y_mw, 
							Varr*1.0e-18, num_sb, num_comp, temp_now, update_count, 
							inflectDp, pwl_xpre, pwl_xpro, inflectk, chamR, Rader, 
							0, p_char, e_field, (num_sb-wall_on))
				if (nucv1>0.): # nucleation
					[N_perbin, y, x[0], Varr[0], np_sum] = nuc.nuc(sumt, np_sum, 
						N_perbin, y, y_mw.reshape(-1, 1), 
						np.squeeze(y_dens*1.0e-3),  
						num_comp, Varr[0], x[0], new_partr, MV, nucv1, nucv2, 
						nucv3, nuc_comp[0])
				
				# reset count to that since original operator-split processes interval met (s)
				update_count = sumt%t0

		
		print('time through simulation (s): ', sumt)
		
		# record output
		if sumt-(save_stp*save_cnt)>-1.e-10:

			[trec, yrec, dydt_vst, Cfactor_vst, save_cnt, 
				Nres_dry, Nres_wet, x2, rbou_rec] = rec.rec(save_cnt, trec, yrec, 
				dydt_vst, Cfactor_vst, y, sumt, rindx, rstoi, rrc, pindx, pstoi, 
				nprod, nreac, num_sb, num_comp, N_perbin, core_diss, 
				Psat, kelv_fac, kimt, kw, Cw, act_coeff, Cfactor, Nres_dry, 
				Nres_wet, x2, x, MV, H2Oi, Vbou, rbou, wall_on, rbou_rec)		
		
		# if time step was temporarily reduced, then return
		if ic_red == 1:
			update_stp = t0
			tnew = update_stp
			ic_red = 0
	
	return(trec, yrec, dydt_vst, Cfactor_vst, Nres_dry, Nres_wet, x2, rbou_rec)
