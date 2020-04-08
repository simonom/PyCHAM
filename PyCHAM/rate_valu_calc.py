'''module to link ode_gen with Rate_coefficients for calculation of gas-phase reaction coefficients'''

import numpy as np


def rate_valu_calc(RO2_indices, H2O, TEMP, lightm, y, time, lat, lon, act_flux_path, 
					DayOfYear, PInit, photo_par_file, Jlen):

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
	M_val = (PInit/(8.3144621*TEMP)*6.02214129e+23)*1.0e-6
	# N2 and O2 given the same multiplication as in atmosphereFunctions.f90 of AtChem2
	N2_val = M_val*0.7809
	O2_val = M_val*0.2095
	
	# note, using __import__ rather than import allows opening in run time, thereby using
	# updated module
	Rate_coeffs = __import__('Rate_coeffs') # calculate the rate coef. for each equation
	# Calculate the new rate coef. array (/s) 
	reac_coef = Rate_coeffs.evaluate_rates(RO2, H2O, TEMP, lightm, time, lat, lon, 
											act_flux_path, DayOfYear, M_val, N2_val, 
											O2_val, photo_par_file, Jlen)
	
	return reac_coef
