'''automatically produces a module for calculating reacion rate coefficients'''
# function to generate a module for calculation of reaction rate coefficients

import datetime

def write_rate_file(reac_coef_g, reac_coef_aq, rrc, rrc_name, testf): # define function

	# inputs: ----------------------------------------------------------------------------
	# reac_coef_g - gas-phase reaction rate coefficient expression from the equation file
	# reac_coef_aq - aqueous-phase reaction rate coefficient expression from the equation file
	# rrc - expression for generic reaction rate coefficients
	# rrc_name - name given to generic reaction rate coefficients	
	# testf - flag for mode: 0 in gas-phase equation mode, 2 for test mode, 3 for
	#			aqueous-phase equation mode
	# ------------------------------------------------------------------------------------

	# open/create relevant file to write module to
	if (testf == 0):
		f = open('PyCHAM/rate_coeffs.py', mode='w')
	if (testf == 3):
		f = open('PyCHAM/rate_coeffs_aq.py', mode='w')
	if (testf == 2):
		f = open('rate_coeffs.py', mode='w')
		
	f.write('\'\'\'module for calculating reaction rate coefficients (automatically generated)\'\'\'\n')
	f.write('# module to hold expressions for calculating rate coefficients # \n') # python will convert \n to os.linesep
	f.write('# created at %s\n' %(datetime.datetime.now()))
	f.write('\n')
	f.write('import numpy\n')
	f.write('import photolysisRates\n')
	f.write('\n')

	# following part is the function (there should be an indent at the start of each line)
	# suggest using one tab
	f.write('def evaluate_rates(RO2, H2O, TEMP, lightm, time, lat, lon, act_flux_path, DayOfYear, M, N2, O2, photo_par_file, Jlen, tf):\n')
	f.write('\n')
	f.write('	# inputs: ------------------------------------------------------------------\n')
	f.write('	# RO2 - names of components included in peroxy radical list\n')	
	f.write('	# M - third body concentration (molecules/cc (air))\n')
	f.write('	# N2 - nitrogen concentration (molecules/cc (air))\n')
	f.write('	# O2 - oxygen concentration (molecules/cc (air))\n')
	f.write('	# H2O, TEMP: given by the user\n')
	f.write('	# lightm: given by the user and is 0 for lights off and 1 for on\n')
	f.write('	# reaction rate coefficients and their names parsed in eqn_parser.py \n')
	f.write('	# Jlen - number of photolysis reactions\n')
	f.write('	# tf - sunlight transmission factor\n')
	f.write('	# ------------------------------------------------------------------------\n')
	f.write('\n')
	f.write('	# calculate generic reaction rate coefficients given by chemical scheme\n')
	# code to calculate rate coefficients given by chemical scheme file
	for line in rrc:
		f.write('	%s \n' %line)
	f.write('\n')
	f.write('	# estimate and append photolysis rates\n')
	f.write('	J = photolysisRates.PhotolysisCalculation(time, lat, lon, TEMP, act_flux_path, DayOfYear, photo_par_file, Jlen, tf)\n')
	f.write('\n')
	f.write('	if lightm == 0:\n')
	f.write('		J = [0]*len(J)\n')

	# calculate the rate coefficient for each equation
	f.write('	rate_values = numpy.zeros((%i))\n' %(len(reac_coef_g)+len(reac_coef_aq)))
	# BE NOTIFIED!!!: before writing the script, 'reac_coef' must be converted to 
	# python-compatible format
	f.write('	\n')	
	f.write('	# reac_coef has been formatted so that python can recognize it\n')
	f.write('	# gas-phase reactions\n')
	for eqn_key in range (len(reac_coef_g)):
		f.write('	rate_values[%s] = %s\n' %(eqn_key, reac_coef_g[eqn_key]))
	f.write('	\n')
	f.write('	# aqueous-phase reactions\n')
	for eqn_key_aq in range (1, len(reac_coef_aq)+1):
		f.write('	rate_values[%s] = %s\n' %(eqn_key+eqn_key_aq, reac_coef_aq[eqn_key_aq-1]))
	f.write('	\n')
	f.write('	return(rate_values)\n')
	f.close()

	return()
