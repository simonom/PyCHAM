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

def ode_gen(t, y, num_speci, num_eqn, rindx, pindx, rstoi, pstoi, H2Oi, 
			TEMP, RO2_indices, num_sb, 
			Psat, mfp, accom_coeff, surfT, y_dens, 
			N_perbin, DStar_org, y_mw, x, core_diss, Varr, Vbou, RH, rad0, 
			Vol0, end_sim_time, pconc, 
			save_step, rbou, therm_sp,
			Cw, light_time, light_stat, nreac, 
			nprod, prodn, reacn, new_partr, MV, nucv1, nucv2, nucv3, inflectDp, 
			pwl_xpre, pwl_xpro, inflectk, nuc_comp, ChamR, Rader, PInit, testf, kwgt):
	
	# ----------------------------------------------------------
	# inputs
	
	# t - suggested time step length (s)
	# num_speci - number of components
	# num_eqn - number of equations
	# Psat - saturation vapour pressures (molecules/cm3 (air))
	# y_dens - components' densities (g/cc)
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
	# MV - molar volume (cc/mol) (1D array)
	# nuc_comp - index of the nucleating component
	# ChamR - spherical equivalent radius of chamber (below eq. 2 Charan (2018)) (m)
	# Rader - flag of whether or not to use Rader and McMurry approach
	# PInit - pressure inside chamber (Pa)
	# testf - flag to say whether in normal mode (0) or testing mode (1)
	# kgwt - mass transfer coefficient for vapour-wall partitioning (cm3/molecule.s)
	# ----------------------------------------------------------
	
	if testf==1:
		return(0,0,0,0) # return dummies

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
	[t_out, y_mat, Nresult, x2] = recording(y, N_perbin, x, step, sumt,
    			0, 0, 0, 0, int(end_sim_time/save_step), 
				num_speci, num_sb, y_mw[:, 0], y_dens[:, 0]*1.0e-3, yp, Vbou)
	
	
	tnew = 0.46875 # relatively small time step (s) for first bit
	# number concentration of nucleated particles formed (# particles/cc (air))
	new_part_sum1 = 0.0
	
	save_count = int(1) # count on number of times saving code called
	# count in number of time steps since time interval was last reduced
	tinc_count = 10
	
	from rate_valu_calc import rate_valu_calc # function to update rate coefficients
	print('starting ode solver')
	
	while sumt < end_sim_time: # step through time intervals to do ode
		
		# increase time step if time step not been decreased for 10 steps
		if tinc_count == 0 and tnew<t0:
			tnew = tnew*2.0

		# check whether lights on or off
		timediff = sumt-np.array(light_time)
		timedish = (timediff == np.min(timediff[timediff>=0])) # reference time index
		lightm = (np.array(light_stat))[timedish] # whether lights on or off now
		 
		
		# update reaction rate coefficients
		reac_coef = rate_valu_calc(RO2_indices, y[H2Oi], TEMP, lightm, y)
		
		
		y0[:] = y[:] # update initial concentrations (molecules/cc (air))
		# update particle volumes at start of time step (um3)
		Vstart = Varr*N_perbin
		redt = 1 # reset time reduction flag
		t = tnew # reset integration time (s)
		
		if num_sb>1:
			# update partitioning coefficients
			[kimt, kelv_fac] = kimt_calc(y, mfp, num_sb, num_speci, accom_coeff, y_mw,   
							surfT, R_gas, TEMP, NA, y_dens, N_perbin, DStar_org, 
							x.reshape(1, -1)*1.0e-6, Psat, therm_sp, H2Oi)
		
		# ensure no confusion that components are present due to low value fillers for  
		# concentrations (molecules/cc (air))
		y0[y0==1.0e-40] = 0.0
		print()
		# enter a while loop that continues to decrease the time step until particle
		# size bins don't change by more than one size bin (given by moving centre)
		while redt == 1:
			print('cumulative time', sumt)
			
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
					dydt[rindx[i, 0:nreac[i]]] -= gprate # loss of reactants
					dydt[pindx[i, 0:nprod[i]]] += gprate # gain of products
	
				if num_sb>1:
						
					# -----------------------------------------------------------
					for ibin in range(num_sb-1): # size bin loop
							
						Csit = y[num_speci*(ibin+1):num_speci*(ibin+2)]
						# sum of molecular concentrations per bin (molecules/cc (air))
						conc_sum = np.zeros((1))
						if pconc>0.0: # if seed particles present
							conc_sum[0] = ((Csit[0:-1].sum())+Csit[-1]*core_diss)
						else:
							conc_sum[0] = Csit.sum()
						# prevent numerical error due to division by zero
						ish = conc_sum==0.0
						conc_sum[ish] = 1.0e-40
							
						# particle surface gas-phase concentration (molecules/cc (air))
						Csit = (Csit/conc_sum)*Psat[:, 0]*kelv_fac[ibin]
						# partitioning rate (molecules/cc.s)
						dydt_all = kimt[:, ibin]*(y[0:num_speci]-Csit)
							
						# gas-phase change
						dydt[0:num_speci] -= dydt_all
						# particle-phase change
						dydt[num_speci*(ibin+1):num_speci*(ibin+2)] += dydt_all
						
				# -----------------------------------------------------------
				# gas-wall partitioning eq. 14 of Zhang et al.  
				# (2015) (https://www.atmos-chem-phys.net/15/4197/2015/
				# acp-15-4197-2015.pdf) (molecules/cc.s (air))
					
				# concentration at wall (molecules/cc (air))
				Csit = y[num_speci*num_sb:num_speci*(num_sb+1)]
				Csit = (Psat[:,0]*(Csit/Cw))
					
				# eq. 14 of Zhang et al. (2015) 
				# (https://www.atmos-chem-phys.net/15/4197/2015/
				# acp-15-4197-2015.pdf) (molecules/cc.s (air))
				dydt_all = (kwgt*Cw)*(y[0:num_speci]-Csit)
						
				# gas-phase change
				dydt[0:num_speci] -= dydt_all
				# wall concentration change
				dydt[num_speci*num_sb:num_speci*(num_sb+1)] += dydt_all
					
				return(dydt)
				
			   
			mod = Explicit_Problem(dydt, y0)
			mod_sim = CVode(mod) # define a solver instance
			mod_sim.atol = 1.0e-3 # absolute tolerance
			mod_sim.rtol = 1.0e-3 # relative tolerance
			mod_sim.discr = 'BDF' # the integration approach, default is 'Adams'
			t_array, res = mod_sim.simulate(t)
		
			y = res[-1, :]
			
			# low value filler for concentrations (molecules/cc (air)) to prevent 
			# numerical errors
			y0[y0==0.0] = 1.0e-40
			y[y==0.0] = 1.0e-40
			
			
			if num_sb>1 and (N_perbin>1.0e-10).sum()>0:
				
				# call on the moving centre method for rebinning particles
				(N_perbin, Varr, y[num_speci::], x, redt, t, tnew) = movcen(N_perbin, 
				Vbou, np.transpose(y[num_speci::].reshape(num_sb, num_speci)), 
				(np.squeeze(y_dens*1.0e-3)), num_sb, num_speci, y_mw, x, Vol0, t, 
				t0, tinc_count, y0[num_speci::], MV)
			else:
				redt = 0
				
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
						num_speci, 0, (np.squeeze(y_dens*1.0e-3)), rad0, PInit)
				
				# particle loss to walls
				[N_perbin, y[num_speci:-(num_speci)]] = wallloss(N_perbin.reshape(-1, 1), 
													y[num_speci:-(num_speci)], Gi, eta_ai, 
													x*2.0e-6, y_mw, Varr*1.0e-18, 
													num_sb, num_speci, TEMP, 
													t, inflectDp, pwl_xpre,
													pwl_xpro, inflectk, ChamR, Rader)
			
			
			# particle nucleation
			if sumt<3600.0 and pconc==0.0:
				[N_perbin, y, x[0], Varr[0], new_part_sum1] = nuc(sumt, new_part_sum1, 
							N_perbin, y, y_mw.reshape(-1, 1), np.squeeze(y_dens*1.0e-3),  
							num_speci, x[0], new_partr, t, MV, nucv1, nucv2, 
							nucv3, nuc_comp)
			
		if sumt>=save_step*save_count: # save at every time step given by save_step (s)
			
			if num_sb>0:
				# particle-phase concentrations (molecules/cc (air))
					yp = np.transpose(y[num_speci:-(num_speci)].reshape(num_sb-1, 
										num_speci))
			else:
				yp = 0.0
			
			# record new values
			[t_out, y_mat, Nresult, x2] = recording(y, N_perbin, x, save_count, 
				sumt, y_mat, Nresult, x2, t_out, int(end_sim_time/save_step), 
				num_speci, num_sb, y_mw[:, 0], y_dens[:, 0]*1.0e-3, yp, Vbou)
			save_count += int(1)
		
	return(t_out, y_mat, Nresult, x2)

