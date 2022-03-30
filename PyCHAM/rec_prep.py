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
'''module to prepare the recording matrices'''
# recording matrices are prepared and initial conditions stored

import numpy as np
import partit_var
import rrc_calc
import cham_up
import scipy.constants as si
import importlib

# define function
def rec_prep(nrec_step, 
	y, y0, rindx, 
	rstoi, pindx, pstoi, nprod, nreac, 
	num_sb, num_comp, N_perbin, core_diss, Psat, mfp,
	accom_coeff, y_mw, surfT, R_gas, temp, tempt, NA, 
	y_dens, x, therm_sp, H2Oi, act_coeff, 
	sumt, Pnow, light_time_cnt, 
	Jlen, Cw, kw, Cfactor, 
	wall_on, Vbou, tnew, nuc_ad, nucv1, nucv2, nucv3, 
	np_sum, update_stp, update_count, injectt, gasinj_cnt, 
	inj_indx, Ct, pmode, pconc, pconct, seedt_cnt, mean_rad, corei, 
	seed_name, seedx, lowsize, uppsize, rad0, radn, std, rbou, 
	const_infl_t, infx_cnt, con_infl_C, MV, partit_cutoff, diff_vol, 
	DStar_org, C_p2w, RH, RHt, tempt_cnt, RHt_cnt, 
	Pybel_objects, nuci, nuc_comp, t0, pcont, pcontf, 
	NOi, HO2i, NO3i, z_prt_coeff, seed_eq_wat, Vwat_inc,
	tot_in_res, Compti, tot_time, cont_inf_reci, cont_inf_i, 
	tot_in_res_indx, chamSA, chamV, kwf, self):
	
	# inputs: --------------------------------------------------------
	# nrec_step - number of steps to record on
	# y - initial concentrations (# molecules/cm3 (air))
	# y0 - component concentrations prior to integration step (# molecules/cm3 (air))
	# rindx - indices of reactants
	# rstoi - stoichiometries of reactants
	# pindx - indices of products
	# pstoi - stoichiometries of products
	# nprod - number of products
	# nreac - number of reactants
	# num_sb - number of size bins
	# num_comp - number of components
	# N_perbin - particle concentration per size bin (#/cc (air))
	# core_diss - dissociation constant of seed
	# Psat - pure component saturation vapour pressure 
	#	(# molecules/cm3 (air))
	# mfp - mean free path (m)
	# accom_coeff - accommodation coefficient
	# y_mw - molecular weight (g/mol)
	# surfT - surface tension (g/s2)
	# R_gas - ideal gas constant (kg.m2.s-2.K-1.mol-1)
	# temp - temperature (K)
	# tempt - times temperatures achieved (s)
	# NA - Avogadro constants (molecules/mol)
	# y_dens - component densities (g/cm3)
	# x - particle radii (um)
	# therm_sp - thermal speed (m/s)
	# H2Oi - index for water
	# act_coeff - activity coefficient
	# self.RO2_indx - index of alkyl peroxy radicals
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
	# Cw - effective absorbing mass of wall (molecules/cc (air))
	# kw - gas-wall partitioning coefficient (/s)
	# Cfactor - conversion factor (ppb/molecules/cc)
	# self.tf - transmission factor for natural sunlight
	# self.light_ad - marker for whether to adapt time interval to 
	#		changing natural light intensity
	# wall_on - marker for whether wall present
	# Vbou - volume boundary of particle size bins (um3)
	# tnew - the proposed integration time interval (s)
	# nuc_ad - marker for whether to adapt time step based on 
	#	nucleation
	# nucv1 - nucleation parameter one
	# nucv2 - nucleation parameter two
	# nucv3 - nucleation parameter three
	# np_sum - cumulative number concentration of new 
	#	particles so far (#/cc (air))
	# update_stp - time interval between operator-split processes
	# update_count - count on time since last operator-split (s)
	# injectt - time of instantaneous injection of components (s)
	# gasinj_cnt - count on injection times of components
	# inj_indx - index of component(s) being injected after 
	#	experiment start
	# Ct - concentration(s) (ppb) of component(s) injected 
	# instantaneously after experiment start
	# pmode - whether number size distributions expressed as modes or explicitly
	# pconc - concentration of injected particles (#/cc (air))
	# pconct - times of particle injection (s)
	# seedt_cnt - count on injection of seed particles
	# mean_rad - mean radius for particle number size 
	#	distribution (um)
	# corei - index of core component
	# seed_name - name(s) of component(s) comprising seed particles
	# seedx - mole ratio of components comprising seed particles
	# lowsize - lower size bin boundary (um)
	# uppsize - upper size bin boundary (um)
	# rad0 - initial radius at size bin centres (um)
	# radn - current radius at size bin centres (um)	
	# std - standard deviation for injected particle number size 
	#	distributions
	# rbou - size bin radius bounds (um)
	# const_infl_t - times for constant influxes (s)
	# infx_cnt - count on constant influx occurrences
	# con_infl_C - influx rates of components with constant influx 
	#	(ppb/s)
	# MV - molar volume (cc/mol)
	# partit_cutoff - the product of saturation vapour pressure and
	#		activity coefficient above which gas-particle
	#		partitioning assumed zero (Pa)
	# diff_vol - diffusion volume of components according to Fuller et al. (1969)
	# DStar_org - gas-phase diffusion coefficient of components (cm2/s)
	# self.seedi - index of seed component(s)
	# C_p2w - concentration of components on the wall due to particle
	# deposition to wall (molecules/cc)
	# RH - relative humidities (fraction 0-1)
	# RHt - times through experiment at which relative humidities reached (s)
	# tempt_cnt - count on chamber temperatures
	# RHt_cnt - chamber relative humidity counts
	# Pybel_objects - the pybel objects for components
	# nuci - index of nucleating component
	# nuc_comp - name of nucleating component
	# t0 - initial integration step (s)
	# pcont - flag for whether particle injection instantaneous or continuous
	# pcontf - current status of particle injection (instantaneous or continuous)
	# NOi - index of NO
	# HO2i - index of HO2
	# NO3i - index of NO3
	# z_prt_coeff - fraction of total gas-particle partitioning coefficient 
	#	below which partitioning to a particle size bin is treated as zero,
	#	e.g. because surface area of that size bin is tiny
	# seed_eq_wat - whether seed particles to be equilibrated with water prior to ODE solver
	# Vwat_inc - whether suppled seed particle volume contains equilibrated water
	# tot_in_res - record of total input of injected components (ug/m3)
	# Compti - index for total injection record for instantaneously injected components
	# tot_time - total experiment time (s)
	# cont_inf_reci - index of components with continuous influx in record
	# cont_inf_i - index of components with continuous influx in concentration array
	# tot_in_res_indx - index of components with recorded influx
	# chamSA - chamber surface area (m2)
	# chamV - chamber volume (m3)
	# kwf - flag for treatment of gas-wall partitioning coefficient
	# tf_UVC - transmission factor for 254 nm wavelength light (0-1)
	# ----------------------------------------------------------------

	# note that instaneous changes occur before recording --------------------
	# update chamber variables
	
	# array to record cumulative influxes of influxing components
	tot_in_res_ft = np.zeros((nrec_step+1, len(tot_in_res)))
	# index of components being stored for influx
	tot_in_res_ft[0, :] = tot_in_res_indx
	tot_in_res_ft[1, :] = tot_in_res # influx at start

	[temp_now, Pnow, lightm, light_time_cnt, tnew, ic_red, update_stp, 
		update_count, Cinfl_now, seedt_cnt, Cfactor, infx_cnt, 
		gasinj_cnt, DStar_org, y, tempt_cnt, RHt_cnt, 
		Psat, N_perbin, x, pconcn_frac,  pcontf, tot_in_res, Cinfl_nowp_indx, 
		Cinfl_nowp] = cham_up.cham_up(sumt, 
		temp, tempt, Pnow, light_time_cnt, 0, 
		nuc_ad, nucv1, nucv2, nucv3, np_sum, 
		update_stp, update_count, 
		injectt, gasinj_cnt, inj_indx, Ct, pmode, pconc, pconct, 
		seedt_cnt, num_comp, y0, y, N_perbin, mean_rad, corei, seedx, seed_name, 
		lowsize, uppsize, num_sb, MV, rad0, radn, std, y_dens, H2Oi, rbou, 
		const_infl_t, infx_cnt, con_infl_C, wall_on, Cfactor, diff_vol, 
		DStar_org, RH, RHt, tempt_cnt, RHt_cnt, Pybel_objects, nuci, nuc_comp,
		y_mw, temp[0], Psat, 0, t0, x, pcont,  pcontf, 0., surfT, act_coeff,
		seed_eq_wat, Vwat_inc, tot_in_res, Compti, tot_time, cont_inf_reci, 
		cont_inf_i, self)
	
	# note that recording occurs after any instaneous changes--------------------
	# array to record time through simulation (s)
	trec = np.zeros((nrec_step))
	# array to record concentrations with time (# molecules/cm3)
	yrec = np.zeros((nrec_step, len(y)))
	yrec[0, :] = y # record initial concentrations (# molecules/cm3)
	# array to record conversion factor
	Cfactor_vst = np.zeros((nrec_step, 1))
	Cfactor_vst[0] = Cfactor

	# arrays to record particle radii and number size distributions	
	if ((num_sb-wall_on) > 0): # if particle size bins present
		x2 = np.zeros((nrec_step, (num_sb-wall_on)))
		Nres_dry = np.zeros((nrec_step, (num_sb-wall_on)))
		Nres_wet = np.zeros((nrec_step, (num_sb-wall_on)))
		rbou_rec = np.zeros((nrec_step, (num_sb-wall_on+1)))
		# concentration of components on the wall due to 
		# particle-wall loss, stacked by component first then by
		# size bin (molecules/cc)
		yrec_p2w = np.zeros((nrec_step, (num_sb-wall_on)*num_comp))
	else:
		x2 = 0.
		Nres_dry = 0.
		Nres_wet = 0.
		rbou_rec = 0.
		yrec_p2w = 0.
	
	if ((num_sb-wall_on) > 0 or wall_on == 1): # if particles or wall present
		
		# update partitioning variables
		[kimt, kelv_fac, kw] = partit_var.kimt_calc(y, mfp, num_sb, num_comp, accom_coeff, y_mw,   
		surfT, R_gas, temp_now, NA, y_dens*1.e3, N_perbin, 
		x.reshape(1, -1)*1.0e-6, Psat, therm_sp, H2Oi, act_coeff, wall_on, 1, partit_cutoff, 
		Pnow, DStar_org, z_prt_coeff, chamSA, chamV, kwf, self)
		
	if (num_sb-wall_on) > 0: # if particles present
		# single particle radius (um) at size bin centre 
		# including contribution of water
		x2[0, :] = x

		# single particle size bin bounds including water (um)
		rbou_rec[0, :] = rbou
		
		# estimate particle number size distributions ----------------------------------
		Vnew = np.zeros((num_sb-wall_on))
		ish = N_perbin[:, 0] > 1.e-10
		
		if ish.sum()>0:
			# rearrange particle concentrations into size bins in rows, components in columns
			Cn = y[num_comp:num_comp*(num_sb-wall_on+1)].reshape(num_sb-wall_on, num_comp)
			# new volume of single particle per size bin (um3) excluding volume of water	
			Vnew[ish] = (np.sum((Cn[ish, :]/(si.N_A*N_perbin[ish]))*MV[:, 0]*1.e12, 1)-
					((Cn[ish, H2Oi]/(si.N_A*N_perbin[ish, 0]))*MV[H2Oi, 0]*1.e12))
			# loop through size bins to find number of particles in each 
			# (# particle/cc (air))
			for Ni in range(0, (num_sb-wall_on)):
				ish = (Vnew>=Vbou[Ni])*(Vnew<Vbou[Ni+1])
				Nres_dry[0, Ni] = N_perbin[ish, 0].sum()
			Nres_wet[0, :] = N_perbin[:, 0] # record with water
		
		# end of number size distribution part ----------------------------------------
		
		# concentration of components on the wall due to particle deposition to wall (# molecules/cm3)
		yrec_p2w[0, :] = C_p2w

	else: # fillers
		kimt = kelv_fac = 0.
	
	# update reaction rate coefficients
	rrc = rrc_calc.rrc_calc(y[H2Oi], temp_now, lightm, y, Pnow, 
			Jlen, y[NOi], y[HO2i], y[NO3i], 
			0., self)

	# chamber environmental conditions ----------------------------------
	# initiate the array for recording chamber temperature (K), pressure (Pa) 
	# and relative humidity (fraction (0-1))
	cham_env = np.zeros((nrec_step, 3))
	
	cham_env[0, 0] = temp_now # temperature (K)
	cham_env[0, 1] = Pnow # pressure (Pa)
	cham_env[0, 2] = y[H2Oi]/Psat[0, H2Oi] # relative humidity (fraction (0-1))
	
	# --------------------------------------------------------------------------------

	return(trec, yrec, Cfactor_vst, Nres_dry, Nres_wet, x2, seedt_cnt, rbou_rec, Cfactor, 
		infx_cnt, yrec_p2w, temp_now, cham_env, Pnow, Psat, cham_env[0, 2], 
		Cinfl_now, tot_in_res_ft)