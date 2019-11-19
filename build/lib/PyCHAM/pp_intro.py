'''module to set up particle phase part of box model, calls on init_water_partit to initiate water partitioning with seed particles and wall'''

import numpy as np
import Size_distributions # custom library - see source code
from init_water_partit import init_water_partit
import scipy.constants as si

def pp_intro(y, num_speci, spec_list, Pybel_objects, TEMP, H2Oi,
			mfp, accom_coeff, y_mw, surfT, 
			DStar_org, RH, num_sb, lowersize, uppersize, total_pconc, tmax, 
			nuc_comp, voli, volP, testf, std, loc, scale, therm_sp,
			Cw, y_dens, Psat, core_diss, kgwt):
	
			
	# inputs -----------------------------------
	# y_mw - molecular weight (g/mol) of components (num_speci, 1)
	# num_sb - number of size bins (excluding wall)
	# lowersize - lowest size bin radius bound (um)
	# uppersize - largest size bin radius bound (um)
	# total_pconc - starting particle concentration (# particle/cc (air))
	# tmax - maximum time step used in ode solver (s)
	# nuc_comp - index of the nucleating component (integer)
	# voli - index of components with vapour pressures given in volP (integer)
	# volP - vapour pressures of components for manual setting (Pa)
	# testf - test flag to say whether in normal mode (0) or test mode for front.py (1)
	#       or test mode for pp_intro.py
	# std - geometric standard deviation of the particle number concentration 
	# 		(dimensionless)
	# loc - shift of lognormal probability distribution function for seed particles 
	# number-size distribution (um)
	# scale - scaling factor of lognormal probability distribution function for seed 
	# particles (dimensionless)
	# Cw - concentration of wall (g/m3 (air))
	# y_dens - liquid density of components (kg/m3) (num_speci, 1)
	# Psat - saturation vapour pressure of components (molecules/cc (air))
	# core_diss - core dissociation constant
	# kgwt - mass transfer coefficient for vapour-wall partitioning (cm3/molecule.s)
	# ------------------------------------------
	if testf==1:
		return(0,0,0,0,0,0,0,0,0,0,0,0) # return dummies
	
	
	# if elements of nuc_comp are relative index (-n), change to absolute
	if nuc_comp<0:
		nuc_comp = num_speci+nuc_comp
	
	R_gas = si.R # ideal gas constant (kg.m2.s-2.K-1.mol-1)
	NA = si.Avogadro # Avogadro's number (molecules/mol)
	
	# create dummy variables if no size bins
	if num_sb==0:
		N_perbin = np.zeros((1,1))
		x = np.zeros((1,1))
		Varr = 0.0
		Vbou = 0.0
		rad0 = 0.0
		Vol0 = 0.0
		rbou = 0.0
	
	# create a number concentration for a lognormal distribution (particles/cc (air))
	# this is where gas partitioning to wall set up
	if testf==2:
		print('calling Size_distributions.lognormal')
	if num_sb>1:
		
		[N_perbin, x, rbou, rwid, Vbou, Varr] = Size_distributions.lognormal(num_sb, 
									total_pconc, std, lowersize, uppersize, loc, scale)
		if testf==2:
			print('finished with Size_distributions.lognormal')

	if num_sb == 1:
		N_perbin = np.array((total_pconc)) # particles per cc (air)
		x = np.zeros(1)
		x[0] = meansize
		# volume bounds of size bin (um3)
		Vbou = np.array(((lowersize**3.0)*(4.0/3.0)*np.pi, 
						(uppersize**3.0)*(4.0/3.0)*np.pi))
		# volume of single particle (um3)
		Varr = np.array(((meansize**3.0)*(4.0/3.0)*np.pi))
	
	if num_sb>0:
		# remember the radii (um) and volumes (um3) at size bin centre before water 
		# partitioning 
		rad0 = np.zeros((len(x)))
		rad0 = x
		Vol0 = np.zeros((len(Varr)))
		Vol0[:] = Varr[:]
	
	num_sb += 1 # add one to size bin number to account for wall

	# append particle-phase concentrations of species to y (molecules/cc (air))
	# note, there's a very small amount of each species in the particle phase to
	# prevent nan errors later on
	y = np.append(y, np.ones((num_sb*num_speci))*1.0e-40)
	
	# molar volume (multiply y_dens by 1e-3 to convert from kg/m3 to g/cc and give
	# MV in units cc/mol)
	MV = (y_mw/(y_dens*1.0e-3)).reshape(num_speci, 1)	
	
	if (total_pconc>0.0):
	
		# core concentration in size bins (molecules/cc (air)):
		# core mass concentration in each size bin (molecules/cc (air))
		y[num_speci*2-1:(num_speci*(num_sb+1)-1):num_speci] = ((y_dens[-1]*1.0e-3)*
								(Varr*1.0e-12*N_perbin)*(1.0/y_mw[-1])*NA)
	
	#if (total_pconc>0.0):
	if testf==2:
		print('calling init_water_partit.py')
	# allow water to equilibrate with particles
	[y, Varr, x] = init_water_partit(x, y, H2Oi, Psat, mfp, num_sb, num_speci, 
					accom_coeff, y_mw, surfT, R_gas, TEMP, NA, y_dens, 
					N_perbin, DStar_org, RH, core_diss, Varr, Vbou, Vol0, tmax, MV,
					therm_sp, Cw, total_pconc, kgwt)
	if testf==2:
		print('finished with init_water_partit.py')		

	# radius of two molecules together (cm)
	new_partr = ((((MV[nuc_comp,0]/si.N_A)*2.0)*3/(4.0*np.pi))**(1/3))
	
	
	return(y, N_perbin, x, Varr, Vbou, rad0, Vol0, rbou, new_partr, MV, num_sb, nuc_comp)