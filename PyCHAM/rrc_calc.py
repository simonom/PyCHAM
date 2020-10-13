'''module to link to calculation of reaction coefficients'''
# this module sets up the final details for calculating
# reaction rate coefficient, which is done via the called module Rate_coeffs 


import numpy as np
import scipy.constants as si


def rrc_calc(RO2_indices, H2O, TEMP, lightm, y, time, lat, lon, act_flux_path, 
		DayOfYear, PInit, photo_par_file, Jlen, tf):

	# ---------------------------------------------
	# inputs:
	
	# TEMP - temperature (K)
	# lightm - whether lights off (0) or on (1)
	# y - component concentrations (molecules/cc (air))
	# time - time of day (for natural light photolysis)
	# lat - latitude
	# lon - longitude
	# act_flux_path - path to file stating actinic flux
	# DayOfYear - day number of the year (1-365)
	# PInit - chamber pressure (Pa)
	# photo_par_file - name of file with with estimates for photolysis absorption
	# 					cross-sections and quantum yields
	# Jlen - number of photolysis reactions
	# tf - sunlight transmission factor
	# ---------------------------------------------
	
	# calculate total RO2 concentration
	if (RO2_indices.size == 0):
		RO2 = 0
	else:
		RO2 = np.sum(y[RO2_indices[:,1]])
		
	# calculate concentrations of third body (M), nitrogen and oxygen
	# calculate gas-phase concentrations of M, N2 and O2 (molecules/cc (air))
	# 1.0e-6 converts from molecules/m3 to molecules/cc
	# R and Avogadro's constant set the same as in atmosphereFunctions.f90 of AtChem2
	M_val = (PInit/(8.3144621*TEMP)*si.N_A)*1.e-6
	# N2 and O2 given the same multiplication as in atmosphereFunctions.f90 of AtChem2
	N2_val = M_val*0.7809
	O2_val = M_val*0.2095
	
	import rate_coeffs
	# calculate the new rate coefficient array (/s) 
	rrc = rate_coeffs.evaluate_rates(RO2, H2O, TEMP, lightm, time, lat, lon, 
						act_flux_path, DayOfYear, M_val, N2_val, 
						O2_val, photo_par_file, Jlen, tf)

		
	return(rrc)

