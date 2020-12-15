'''module to prepare the recording matrices'''
# recording matrices are prepared and initial conditions stored

import numpy as np
import partit_var
import rrc_calc
import cham_up
import scipy.constants as si

# define function
def rec_prep(nrec_step, 
	y, rindx, 
	rstoi, pindx, pstoi, nprod, dydt_vst, nreac, 
	num_sb, num_comp, N_perbin, core_diss, Psat, mfp,
	accom_coeff, y_mw, surfT, R_gas, temp, tempt, NA, 
	y_dens, x, therm_sp, H2Oi, act_coeff, 
	RO2_indx, sumt, Pnow, light_stat, light_time, 
	light_time_cnt, daytime, lat, lon, af_path, 
	dayOfYear, photo_path, Jlen, Cw, kw, Cfactor, tf, 
	light_ad, wall_on, Vbou, tnew, nuc_ad, nucv1, nucv2, nucv3, 
	np_sum, update_stp, update_count, injectt, gasinj_cnt, 
	inj_indx, Ct, pmode, pconc, pconct, seedt_cnt, mean_rad, corei, 
	seed_name, seedVr, lowsize, uppsize, rad0, radn, std, rbou, 
	const_infl_t, infx_cnt, con_infl_C, MV, partit_cutoff, diff_vol, DStar_org, seedi):
	
	# inputs: --------------------------------------------------------
	# nrec_step - number of steps to record on
	# y - initial concentrations (molecules/cc (air))
	# rindx - indices of reactants
	# rstoi - stoichiometries of reactants
	# pindx - indices of products
	# pstoi - stoichiometries of products
	# nprod - number of products
	# dydt_vst - recording for change tendencies
	# nreac - number of reactants
	# num_sb - number of size bins
	# num_comp - number of components
	# N_perbin - particle concentration per size bin (#/cc (air))
	# core_diss - dissociation constant of seed
	# Psat - pure component saturation vapour pressure 
	#	(molecules/cm3 (air))
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
	# RO2_indx - index of peroxy radicals
	# sumt - cumulative time through simulation (s)
	# Pnow - chamber pressure (Pa)
	# light_stat - light status
	# light_time - times corresponding to light status
	# light_time_cnt - count on light status array
	# daytime - start time of simulation (s)
	# lat - latitude
	# lon - longitude
	# af_path - path to actinic flux file
	# dayOfYear - day of year
	# photo_path - path to photochemistry file
	# Jlen - length of photochemistry equations
	# Cw - effective absorbing mass of wall (molecules/cc (air))
	# kw - gas-wall partitioning coefficient (/s)
	# Cfactor - conversion factor (ppb/molecules/cc)
	# tf - transmission factor for natural sunlight
	# light_ad - marker for whether to adapt time interval to 
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
	# seedVr - volume ration of component(s) comprising seed particles
	# seed_name - name(s) of component(s) comprising seed particles
	# seedVr - volume ratio of component(s) comprising seed particles
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
	# seedi - index of seed component(s)
	# ----------------------------------------------------------------

	# array to record time through simulation (s)
	trec = np.zeros((nrec_step))
	# array to record concentrations with time (molecules/cc/s)
	yrec = np.zeros((nrec_step, len(y)))
	yrec[0, :] = y # record initial concentrations (molecules/cc/s)
	# array to record conversion factor
	Cfactor_vst = np.zeros((nrec_step, 1))
	Cfactor_vst[0] = Cfactor

	# arrays to record particle radii and number size distributions	
	if ((num_sb-wall_on) > 0):
		x2 = np.zeros((nrec_step, (num_sb-wall_on)))
		Nres_dry = np.zeros((nrec_step, (num_sb-wall_on)))
		Nres_wet = np.zeros((nrec_step, (num_sb-wall_on)))
		rbou_rec = np.zeros((nrec_step, (num_sb-wall_on+1)))
	else:
		x2 = 0.
		Nres_dry = 0.
		Nres_wet = 0.
		rbou_rec = 0.

	# update chamber variables
	[temp_now, Pnow, lightm, light_time_cnt, tnew, ic_red, update_stp, 
		update_count, Cinfl_now, seedt_cnt, Cfactor, infx_cnt, 
		gasinj_cnt, DStar_org] = cham_up.cham_up(sumt, temp, tempt, 
		Pnow, light_stat, light_time, light_time_cnt, light_ad, 0, 
		nuc_ad, nucv1, nucv2, nucv3, np_sum, 
		update_stp, update_count, lat, lon, dayOfYear, photo_path, 
		af_path, injectt, gasinj_cnt, inj_indx, Ct, pmode, pconc, pconct, 
		seedt_cnt, num_comp, y, N_perbin, mean_rad, corei, seedVr, seed_name, 
		lowsize, uppsize, num_sb, MV, rad0, radn, std, y_dens, H2Oi, rbou, 
		const_infl_t, infx_cnt, con_infl_C, wall_on, Cfactor, seedi, diff_vol, DStar_org)
	
	
	if (num_sb-wall_on)>0: # if particles present
		
		# update partitioning variables
		[kimt, kelv_fac] = partit_var.kimt_calc(y, mfp, num_sb, num_comp, accom_coeff, y_mw,   
		surfT, R_gas, temp_now, NA, y_dens*1.e3, N_perbin, 
		x.reshape(1, -1)*1.0e-6, Psat, therm_sp, H2Oi, act_coeff, wall_on, 1, partit_cutoff, 
		Pnow, DStar_org)
		
		# single particle radius (um) at size bin centre 
		# including contribution of water
		x2[0, :] = x

		# single particle size bin bounds including water (um)
		rbou_rec[0, :] = rbou
		
		# estimate particle number size distributions ----------------------------------
		Vnew = np.zeros((num_sb-wall_on))
		ish = N_perbin[:, 0]>1.0e-10
		
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

	else: # fillers
		kimt = kelv_fac = 0.
	
	# update reaction rate coefficients
	rrc = rrc_calc.rrc_calc(RO2_indx, 
			y[H2Oi], temp_now, lightm, y, daytime+sumt, 
			lat, lon, af_path, dayOfYear, Pnow, 
			photo_path, Jlen, tf)

	if len(dydt_vst)>0:		
		# record any change tendencies of specified components
		import dydt_rec
		dydt_vst = dydt_rec.dydt_rec(y, rindx, rstoi, rrc, pindx, pstoi, nprod, 0, 
					dydt_vst, nreac, num_sb, num_comp, N_perbin, core_diss, Psat, kelv_fac, 
					kimt, kw, Cw, act_coeff, corei)


	return(trec, yrec, dydt_vst, Cfactor_vst, Nres_dry, Nres_wet, x2, seedt_cnt, rbou_rec, Cfactor, infx_cnt)
