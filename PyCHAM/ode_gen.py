'''module to solve the gas-phase reactions odes, and odes for partitioning of gas-particle and gas-wall'''

import numpy as np
from assimulo.problem import Explicit_Problem
from assimulo.solvers import CVode
import numba
from numba import jit, f8
import matplotlib.pyplot as plt
import ipdb
from kimt_calc import kimt_calc
from recording import recording
from mov_cen_main import mov_cen_main as movcen # moving centre method for rebinning
from coag import coag
from wallloss import wallloss
from nuc import nuc
import scipy.constants as si
from rate_valu_calc import rate_valu_calc # function to update rate coefficients
import math
from pp_dursim import pp_dursim

def ode_gen(t, y, num_speci, num_eqn, rindx, pindx, rstoi, pstoi, H2Oi, 
			TEMP, RO2_indices, num_sb, 
			Psat, mfp, accom_coeff, surfT, y_dens, 
			N_perbin, DStar_org, y_mw, x, core_diss, Varr, Vbou, RH, rad0, 
			Vol0, end_sim_time, pconc, 
			save_step, rbou, therm_sp,
			Cw, light_time, light_stat, nreac, 
			nprod, prodn, reacn, new_partr, MV, nucv1, nucv2, nucv3, inflectDp, 
			pwl_xpre, pwl_xpro, inflectk, nuc_comp, ChamR, Rader, PInit, testf, kgwt,
			dydt_vst, daytime, lat, lon, act_flux_path, DayOfYear, Ct, injectt, inj_indx,
			corei, const_compi, const_comp, const_infli, Cinfl, act_coeff, p_char, 
			e_field, const_infl_t, int_tol, photo_par_file, Jlen, dil_fac, pconct,
			lowersize, uppersize, mean_rad, std):
	
	# ----------------------------------------------------------
	# inputs
	
	# t - suggested time step length (s)
	# num_speci - number of components
	# num_eqn - number of equations
	# Psat - saturation vapour pressures (molecules/cm3 (air))
	# y_dens - components' densities (kg/m3)
	# y_mw - components' molecular weights (g/mol)
	# x - radii of particle size bins (um) (excluding walls)
	# therm_sp - thermal speed of components (m/s)
	# DStar_org - gas-phase diffusion coefficient of components (m2/s)
	# Cw - concentration of wall (molecules/cm3 (air))
	# light_time - times (s) of when lights on and lights off (corresponding to light 
	# 				status in light_stat)
	# light_stat - order of lights on (1) and lights off (0)
	# chamA - chamber area (m2)
	# nreac - number of reactants per equation
	# nprod - number of products per equation
	# pindx - indices of equation products (cols) in each equation (rows)
	# prodn - pindx no. of columns
	# reacn - rindx no. of columns
	# rindx - index of reactants per equation
	# rstoi - stoichometry of reactants per equation
	# pstoi - stoichometry of products
	# pconc - concentration of seed particles (#/cc (air)) (1)
	# new_partr - radius of two ELVOC molecules together in a newly nucleating 
	# particle (cm)
	# rad0 - original radius at size bin centres (um)
	# MV - molar volume (cc/mol) (1D array)
	# nuc_comp - index of the nucleating component
	# ChamR - spherical equivalent radius of chamber (below eq. 2 Charan (2018)) (m)
	# Rader - flag of whether to use Rader and McMurry approach (1), or manual inputs
	# 		for wall loss (0) or to ignore particle wall loss (-1)
	# PInit - pressure inside chamber (Pa)
	# testf - flag to say whether in normal mode (0) or testing mode (1)
	# kgwt - mass transfer coefficient for vapour-wall partitioning (/s)
	# dydt_vst - dictionary containing record of tendency to change due to box model 
	#			mechanisms
	# daytime - time of day experiment starts (s)
	# DayOfYear - day of the year (1-365)
	# Ct - concentrations achieved due to injections after experiment start 
	# (molecules/cc (air))
	# injectt - time of injections (s)
	# inj_indx - index of component being injected after experiment start
	# corei - index of core component
	# const_compi - index of components with constant gas-phase concentration
	# const_infli - index of components with constant influx to chamber
	# Cinfl - concentration of constant influx(es) (molecules/cc.s)
	# act_coeff - activity coefficient of components (dimensionless)
	# p_char - average number of charges per particle (/particle)
	# e_field - average electric field inside chamber (g.m/A.s3)
	# const_infl_t - times of constant influx (s)
	# int_tol - absolute (0 index) and relative (1 index) tolerances for integration
	# photo_par_file - name of file with with estimates for photolysis absorption
	# 					cross-sections and quantum yields
	# Jlen - number of photolysis reactions
	# dil_fac - dilution factor rate (/s)
	# pconct - times (s) at which seep particles injected into chamber
	# lowersize - smallest radius bound (um)
	# uppersize - greatest radius bound (um)
	# mean_rad - mean radius of particles (relevant if only one size bin or number size
	# distribution being calculated (um)
	# std - standard deviation for lognormal size distribution (dimensionless)
	# ------------------------------------------------------------------------------------

	# count on injection times of seed particles
	seedt_count = 0

	# ------------------------------------------------------------------------------------
	# testing mode
	if testf==1:
		return(0, 0, 0, 0) # return dummies
		
	if testf==2:
		# called from test_kimt_calc.py
		# recreate the solute effect used in dydt function below
		sol_eff = np.zeros((num_sb-1))
		for ibin in range(num_sb-1): # size bin loop
			Csit = y[num_speci*(ibin+1):num_speci*(ibin+2)]
			
			# sum of molecular concentrations per bin (molecules/cc (air))
			conc_sum = np.zeros((1))
			conc_sum[0] = ((Csit.sum()-Csit[corei])+Csit[corei]*core_diss)

			sol_eff[ibin] = Csit[H2Oi]/conc_sum # mole fraction of water in particle
		return(sol_eff)
		
	# ------------------------------------------------------------------------------------

	R_gas = si.R # ideal gas constant (kg.m2.s-2.K-1.mol-1)
	NA = si.Avogadro # Avogadro's number (molecules/mol)
	
	step = 0 # ode time interval step number
	t0 = t # remember original suggested time step (s)
	# final +1 for ELVOC in newly nucleating particles
	y0 = np.zeros((num_speci+num_sb*num_speci))	
	y0[:] = y[:] # initial concentrations (molecules/cc (air))
	y00 = np.zeros((num_speci+num_sb*num_speci))	
	y00[:] = y[:] # initial concentrations (molecules/cc (air))
	
	
	# initial volumes of particles in size bins at start of time steps
	if num_sb>1:
		Vstart = np.zeros((num_sb-1))
		Vstart[:] = Vol0[:]*N_perbin
	else:
		Vstart = 0.0
	sumt = 0.0 # total time integrated over (s)
	
	# record initial values
	if num_sb>0:
		# particle-phase concentrations (molecules/cc (air))
		yp = np.transpose(y[num_speci:-(num_speci)].reshape(num_sb-1, num_speci)) 
	else:
		yp = 0.0
	
	if len(light_time)>0:
		# check whether lights have changed
		timediff = (sumt)-np.array(light_time)
		timedish = (timediff == np.min(timediff[timediff>=0])) # index of reference time
		lightm = (np.array(light_stat))[timedish] # whether lights on or off now
	else: # if no input provided default to lights off
		lightm = 0
		
	
	if num_sb>1:
		# update partitioning coefficients
		[kimt, kelv_fac] = kimt_calc(y, mfp, num_sb, num_speci, accom_coeff, y_mw,   
							surfT, R_gas, TEMP, NA, y_dens, N_perbin, DStar_org, 
							x.reshape(1, -1)*1.0e-6, Psat, therm_sp, H2Oi)
	else:
		kimt = 0.0
		kelv_fac = 0.0
	
	save_count = int(1) # count on number of times saving code called
	
	# reaction rate coefficients at experiment time = 0s
	reac_coef = rate_valu_calc(RO2_indices, y[H2Oi], TEMP, lightm, y, 
								daytime+sumt, 
								lat, lon, act_flux_path, DayOfYear, PInit, 
								photo_par_file, Jlen)
	
	[t_out, y_mat, Nresult_dry, Nresult_wet, x2, dydt_vst] = recording(y, N_perbin, x, 
				save_count-1, sumt,
    			0, 0,  0, 0, 0, math.ceil(end_sim_time/save_step), 
				num_speci, num_sb, y_mw[:, 0], y_dens[:, 0]*1.0e-3, yp, Vbou, rindx, 
				rstoi, pindx, nprod, dydt_vst, RO2_indices, H2Oi, TEMP, lightm, nreac,
				pconc[:, seedt_count], core_diss, Psat, kelv_fac, kimt, kgwt, Cw, daytime+sumt, lat, lon, 
				act_flux_path, DayOfYear, act_coeff, PInit, photo_par_file, Jlen, 
				reac_coef)
	
	
	tnew = t0 # initial time step (s)
	# number concentration of nucleated particles formed (# particles/cc (air))
	new_part_sum1 = 0.0
	
	# count in number of time steps since time interval was last reduced
	tinc_count = 10
	
	# count on injection times of components
	inj_count = 0
	
	# number of components with constant gas-phase concentration and with constant influx
	num_const_compi = len(const_compi)
	const_infli_len = len(const_infli)
	# needs to be a numpy array to be used in integrator
	const_compi = np.array(const_compi)
	const_infli = np.array(const_infli)
	
	print('starting ode solver')
	
	while sumt < end_sim_time: # step through time intervals to do ode
		
		if (sumt+tnew)>end_sim_time: # ensure we finish at correct time
			tnew = end_sim_time-sumt # integration time step (s)
		t = tnew # reset integration time (s)
		
		if len(light_time)>0:
			# check whether lights have changed
			light_ind = ((sumt)>np.array(light_time)).sum()-1
			# whether lights on (1) or off (0) during this step
			lightm = (np.array(light_stat))[light_ind]
		else: # if no input provided default to lights off
			lightm = 0
		
		# --------------------------------------------------------------------------------
		# component injections check
		
		if len(injectt)>0: # if any injections occur
			# check whether latest component injection occurs
			inj_count_new = (sumt >= injectt).sum()
			# update injections if new injection time reached
			if inj_count_new>inj_count:
				for i in range(len(inj_indx)):
					# account for injection in gas-phase concentration 
					# (molecules/cc (air))
					y[int(inj_indx[i])] += Ct[i, inj_count]
				inj_count = inj_count_new # update count on injections
				
				
# 			# only need this if using inj_count_new = ((sumt+tnew)>injectt).sum() above
#			# rather than inj_count_new = ((sumt)>injectt).sum()
# 			# withdraw injections if previous injection time no longer reached, this is in
# 			# case integration step has been reduced by moving centre method
# 			if inj_count_new<inj_count:
# 				for i in range(len(inj_indx)):
# 					# account for injection in gas-phase concentration 
# 					# (molecules/cc (air))
# 					y[int(inj_indx[i])] -= Ct[i, inj_count]
# 				inj_count = inj_count_new # update count on injections
		# --------------------------------------------------------------------------------
		# constant influxes check
		
		# update index counter for constant influxes - used in integrator below
		if len(const_infl_t)>0:
			inf_ind = int(((sumt)>=const_infl_t).sum())-1
			
		# --------------------------------------------------------------------------------
		# seed particle influx check
		
		if int((sumt>=pconct).sum())-1 > seedt_count:
			
			seedt_count += 1
			
			[y[num_speci:-num_speci], N_perbin, x, 
						Varr] = pp_dursim(y[num_speci:-num_speci], N_perbin, 
									mean_rad[0, seedt_count],
									pconc[:, seedt_count], corei, lowersize, 
									uppersize, num_speci, num_sb, MV, rad0, 
									std[0, seedt_count], y_dens, H2Oi)

		# --------------------------------------------------------------------------------
		
		# update reaction rate coefficients
		reac_coef = rate_valu_calc(RO2_indices, y[H2Oi], TEMP, lightm, y, 
									daytime+sumt, 
									lat, lon, act_flux_path, DayOfYear, PInit, 
									photo_par_file, Jlen)
		
		y0[:] = y[:] # update initial concentrations (molecules/cc (air))
		# update particle volumes at start of time step (um3)
		Vstart = Varr*N_perbin
		redt = 1 # reset time reduction flag
		
		
		if num_sb>1:
			# update partitioning coefficients
			[kimt, kelv_fac] = kimt_calc(y, mfp, num_sb, num_speci, accom_coeff, y_mw,   
							surfT, R_gas, TEMP, NA, y_dens, N_perbin, DStar_org, 
							x.reshape(1, -1)*1.0e-6, Psat, therm_sp, H2Oi)
		
		
		# ensure no confusion that components are present due to low value fillers for  
		# concentrations (molecules/cc (air))
		y0[y0==1.0e-40] = 0.0
		
		# enter a while loop that continues to decrease the time step until particle
		# size bins don't change by more than one size bin (given by moving centre)
		# note, need to have rstoi and pstoi multiplication in the gas-phase reaction part
		while redt == 1:
			print('cumulative time (s)', sumt, lightm)

			# numba compiler to convert to machine code
			@jit(f8[:](f8, f8[:]), nopython=True)
			# ode solver -------------------------------------------------------------
			def dydt(t, y):
				
				# empty array to hold rate of change (molecules/cc(air).s)
				dydt = np.zeros((len(y)))
				# gas-phase rate of change ------------------------------------
				for i in range(num_eqn): # equation loop
			
					# gas-phase rate of change (molecules/cc (air).s)
					gprate = ((y[rindx[i, 0:nreac[i]]]**rstoi[i, 0:nreac[i]]).prod())*reac_coef[i] 
					dydt[rindx[i, 0:nreac[i]]] -= gprate*rstoi[i, 0:nreac[i]] # loss of reactants
					dydt[pindx[i, 0:nprod[i]]] += gprate*pstoi[i, 0:nprod[i]] # gain of products
				
				# honour the constant concentration of components with this property
				if num_const_compi>0:
					dydt[const_compi[:]] = 0.0
				# the constant gas-phase influx of components with this property
				if const_infli_len>0:
					for i in range(const_infli_len):
						dydt[const_infli[i]] += Cinfl[i, inf_ind]
					
				if num_sb>1: # as num_sb includes 1 for wall
					# gas-particle partitioning, based on eqs. 3 and 4 of Zaveri et al.
					# (2008): doi:10.1029/2007JD008782
					# and eq. 3 of Riipinen et al.
					# (2010): doi:10.1016/j.atmosenv.2009.11.022
					# -----------------------------------------------------------
					for ibin in range(num_sb-1): # size bin loop
							
						Csit = y[num_speci*(ibin+1):num_speci*(ibin+2)]
						# sum of molecular concentrations per bin (molecules/cc (air))
						conc_sum = np.zeros((1))
						
						conc_sum[0] = ((Csit.sum()-Csit[corei])+Csit[corei]*core_diss)
						
						# prevent numerical error due to division by zero
						ish = conc_sum==0.0
						conc_sum[ish] = 1.0e-40
							
						# particle surface gas-phase concentration (molecules/cc (air))
						Csit = (Csit/conc_sum)*Psat[:, 0]*kelv_fac[ibin]*act_coeff[:, 0]
						# partitioning rate (molecules/cc.s)
						dydt_all = kimt[:, ibin]*(y[0:num_speci]-Csit)
							
						# gas-phase change
						dydt[0:num_speci] -= dydt_all
						# particle-phase change
						dydt[num_speci*(ibin+1):num_speci*(ibin+2)] += dydt_all
				
				if (kgwt*Cw)>1.0e-10:
					# -----------------------------------------------------------
					# gas-wall partitioning (dydt is molecules/cc.s (air))
						
					# concentration at wall (molecules/cc (air))
					Csit = y[num_speci*num_sb:num_speci*(num_sb+1)]
					Csit = (Psat[:,0]*(Csit/Cw)*act_coeff[:, 0]) # with Raoult term
					
					dydt_all = (kgwt)*(y[0:num_speci]-Csit)
							
					# gas-phase change
					dydt[0:num_speci] -= dydt_all
					# wall concentration change
					dydt[num_speci*num_sb:num_speci*(num_sb+1)] += dydt_all
					
				return(dydt)
				
			   
			mod = Explicit_Problem(dydt, y0)
			mod_sim = CVode(mod) # define a solver instance
			# absolute tolerance, going higher than can 1.0e-3 cause issues with water 
			# vapour
			mod_sim.atol = int_tol[0]
			# relative tolerance, going higher than 1.0e-4 can cause issues with water 
			# vapour
			mod_sim.rtol = int_tol[1]
			mod_sim.discr = 'BDF' # the integration approach, default is 'Adams'
			t_array, res = mod_sim.simulate(t)
		
			y = res[-1, :]
			
			# low value filler for concentrations (molecules/cc (air)) to prevent 
			# numerical errors
			y0[y0==0.0] = 1.0e-40
			y[y==0.0] = 1.0e-40
			
			if num_sb>1 and (N_perbin>1.0e-10).sum()>0:
				
				# call on the moving centre method for rebinning particles
				(N_perbin, Varr, y[num_speci::], x, redt, t, tnew, 
				y[0:num_speci]) = movcen(N_perbin, 
				Vbou, np.transpose(y[num_speci::].reshape(num_sb, num_speci)), 
				(np.squeeze(y_dens*1.0e-3)), num_sb, num_speci, y_mw, x, Vol0, t, 
				t0, tinc_count, y0[num_speci::], MV, Psat[:, 0], y[0:num_speci], 
				y0[0:num_speci])
			else: # if bypassing moving centre
				redt = 0
				if t<t0 and tinc_count<=0:
					tnew = t*2.0
				if tnew>t0: # in case tnew exceeds user-defined maximum for time step
					tnew = t0
			
			
			# if time step needs reducing then reset gas-phase concentrations to their
			# values preceding the ode, this will already have been done inside moving
			# centre module for particle-phase
			if redt == 1:
				y[0:num_speci] = y0[0:num_speci]
			
			# start counter to determine when to next try increasing time interval
			if redt == 1:
				tinc_count = 10
			if redt == 0 and tinc_count>-1:
				tinc_count -= 1
			if tinc_count==-1:
				tinc_count = 10
		
		sumt += t # total time covered (s)
		step += 1 # ode time interval step number

		if num_sb>1:
			if (N_perbin>1.0e-10).sum()>0:
				# coagulation
				# y indices due to final element in y being number of ELVOC molecules
				# contributing to newly nucleated particles
				[N_perbin, y[num_speci:-(num_speci)], x, Gi, eta_ai, Varr] = coag(RH, 
						TEMP, x*1.0e-6, (Varr*1.0e-18).reshape(1, -1), 
						y_mw.reshape(-1, 1), x*1.0e-6, 
						np.transpose(y[num_speci::].reshape(num_sb, num_speci)), 
						(N_perbin).reshape(1, -1), t, (Vbou*1.0e-18).reshape(1, -1), 
						num_speci, 0, (np.squeeze(y_dens*1.0e-3)), rad0, PInit, 0,
						np.transpose(y[num_speci::].reshape(num_sb, num_speci)),
						(N_perbin).reshape(1, -1), (Varr*1.0e-18).reshape(1, -1))

				
				if Rader > -1:
					
					# particle loss to walls
					[N_perbin, 
					y[num_speci:-(num_speci)]] = wallloss(N_perbin.reshape(-1, 1), 
														y[num_speci:-(num_speci)], Gi, 
														eta_ai, x*2.0e-6, y_mw, 
														Varr*1.0e-18, num_sb, num_speci, 
														TEMP, t, inflectDp, pwl_xpre,
														pwl_xpro, inflectk, ChamR, Rader, 
														0, p_char, e_field)
					
			# particle nucleation
			if len(nuc_comp)>0:
				
				[N_perbin, y, x[0], Varr[0], new_part_sum1] = nuc(sumt, new_part_sum1, 
							N_perbin, y, y_mw.reshape(-1, 1), np.squeeze(y_dens*1.0e-3),  
							num_speci, x[0], new_partr, t, MV, nucv1, nucv2, 
							nucv3, nuc_comp[0])
							
			# dilution of aerosol (gases and particles), most likely due to extraction
			# from chamber
			y -= y*(dil_fac*t) # dilution of gases (molecules/cc (air))
			N_perbin -= N_perbin*(dil_fac*t) # dilution of particle phase (#/cc (air))
# 			print('y af', y)
			
		# save at every time step given by save_step (s) and at end of experiment
		if sumt>=save_step*save_count or sumt == end_sim_time:
			
			if num_sb>0:
				
				# particle-phase concentrations (molecules/cc (air))
				yp = np.transpose(y[num_speci:-(num_speci)].reshape(num_sb-1, 
										num_speci))
			else:
				yp = 0.0
			
			# record values at the end of this time step, note that these will correspond
			# to sumt, which is the cumulative time (s) at the end of this time step.  
			# This makes sense because the integration results are what is being saved 
			# and these correspond to sumt (time at the end of the time step)
			[t_out, y_mat, Nresult_dry, Nresult_wet, x2, dydt_vst] = recording(y, 
				N_perbin, x, save_count, 
				sumt, y_mat, Nresult_dry, Nresult_wet, x2, t_out, int(end_sim_time/save_step), 
				num_speci, num_sb, y_mw[:, 0], y_dens[:, 0]*1.0e-3, yp, Vbou, rindx, 
				rstoi, pindx, nprod, dydt_vst, RO2_indices, H2Oi, TEMP, lightm, nreac,
				pconc, core_diss, Psat, kelv_fac, kimt, kgwt, Cw, daytime+sumt, lat, lon, 
				act_flux_path, DayOfYear, act_coeff, PInit, photo_par_file, Jlen, reac_coef)
				
			save_count += int(1) # track number of times saved at
			
	return(t_out, y_mat, Nresult_dry, Nresult_wet, x2, dydt_vst)