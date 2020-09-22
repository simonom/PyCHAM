'''module for calculating reaction rate coefficients (automatically generated)'''
# module to hold expressions for calculating rate coefficients # 
# created at 2020-09-22 17:16:23.009955

import numpy
import photolysisRates

def evaluate_rates(RO2, H2O, TEMP, lightm, time, lat, lon, act_flux_path, DayOfYear, M, N2, O2, photo_par_file, Jlen, tf):

	# inputs: ------------------------------------------------------------------
	# RO2 - names of components included in peroxy radical list
	# M - third body concentration (molecules/cc (air))
	# N2 - nitrogen concentration (molecules/cc (air))
	# O2 - oxygen concentration (molecules/cc (air))
	# H2O, TEMP: given by the user
	# lightm: given by the user and is 0 for lights off and 1 for on
	# reaction rate coefficients and their names parsed in eqn_parser.py 
	# Jlen - number of photolysis reactions
	# tf - sunlight transmission factor
	# ------------------------------------------------------------------------

	# calculate generic reaction rate coefficients given by chemical scheme

	# estimate and append photolysis rates
	J = photolysisRates.PhotolysisCalculation(time, lat, lon, TEMP, act_flux_path, DayOfYear, photo_par_file, Jlen, tf)

	if lightm == 0:
		J = [0]*len(J)
	rate_values = numpy.zeros((11))
	# reac_coef has been formatted so that python can recognize it
	rate_values[0] = J[1]
	rate_values[1] = J[2]
	rate_values[2] = 6.3e-16*numpy.exp(-580/TEMP)*0.57
	rate_values[3] = 6.3e-16*numpy.exp(-580/TEMP)*0.37
	rate_values[4] = 1.2e-12*numpy.exp(490/TEMP)*0.65
	rate_values[5] = 1.2e-12*numpy.exp(490/TEMP)*0.35
	rate_values[6] = 1.2e-11*numpy.exp(440/TEMP)*0.482
	rate_values[7] = 1.2e-11*numpy.exp(440/TEMP)*0.293
	rate_values[8] = 1.2e-11*numpy.exp(440/TEMP)*0.065
	rate_values[9] = 1.2e-11*numpy.exp(440/TEMP)*0.08
	rate_values[10] = 6.3e-16*numpy.exp(-580/TEMP)*0.06
	
	return rate_values
