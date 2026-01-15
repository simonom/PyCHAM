##########################################################################################
#                                                                                        #
#    Copyright (C) 2018-2026 Simon O'Meara : simon.omeara@manchester.ac.uk               #
#                                                                                        #
#    All Rights Reserved.                                                                #
#    This file is part of PyCHAM                                                         #
#                                                                                        #
#    PyCHAM is free software: you can redistribute it and/or modify it under             #
#    the terms of the GNU General Public License as published by the Free Software       #
#    Foundation, either version 3 of the License, or (at your option) any later          #
#    version.                                                                            #
#                                                                                        #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT               #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              #
#    details.                                                                            #
#                                                                                        #
#    You should have received a copy of the GNU General Public License along with        #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                #
#                                                                                        #
##########################################################################################
'''automatically produces a module for calculating reacion rate coefficients'''
# function to generate a module for calculation of reaction rate coefficients

import datetime

def write_rate_file(testf, self): # define function
	
	# inputs: ----------------------------------------------------------------------------
	# self.reac_coef_g - gas-phase reaction rate coefficient expression 
	#	from the equation file
	# self.reac_coef_aq - aqueous-phase reaction rate coefficient expression from the 
	# equation file
	# self.reac_coef_su - surface (e.g. wall) reaction rate coefficient expression from 
	# the equation file
	# self.rrc - expression for generic reaction rate coefficients
	# self.rrc_name - name given to generic reaction rate coefficients	
	# testf - flag for mode: 0 in gas-phase equation mode, 2 for test mode, 3 for
	#			aqueous-phase equation mode
	# ------------------------------------------------------------------------------------

	# open/create relevant file to write module to
	if (testf == 0):
		f = open(str(self.PyCHAM_path + '/PyCHAM/rate_coeffs.py'), mode='w')
	if (testf == 3):
		f = open(str(self.PyCHAM_path + '/PyCHAM/rate_coeffs_aq.py'), mode='w')
	if (testf == 2):
		f = open('rate_coeffs.py', mode='w')
	f.write('##########################################################################################\n')
	f.write('#                                                                                        #\n')
	f.write('#    Copyright (C) 2018-2026 Simon O\'Meara : simon.omeara@manchester.ac.uk               #\n')
	f.write('#                                                                                        #\n')
	f.write('#    All Rights Reserved.                                                                #\n')
	f.write('#    This file is part of PyCHAM                                                         #\n')
	f.write('#                                                                                        #\n')
	f.write('#    PyCHAM is free software: you can redistribute it and/or modify it under             #\n')
	f.write('#    the terms of the GNU General Public License as published by the Free Software       #\n')
	f.write('#    Foundation, either version 3 of the License, or (at your option) any later          #\n')
	f.write('#    version.                                                                            #\n')
	f.write('#                                                                                        #\n')
	f.write('#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT               #\n')
	f.write('#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       #\n')
	f.write('#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              #\n')
	f.write('#    details.                                                                            #\n')
	f.write('#                                                                                        #\n')
	f.write('#    You should have received a copy of the GNU General Public License along with        #\n')
	f.write('#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                #\n')
	f.write('#                                                                                        #\n')
	f.write('##########################################################################################\n')
	f.write('\'\'\'module for calculating reaction rate coefficients (automatically generated)\'\'\'\n')
	f.write('# module to hold expressions for calculating rate coefficients # \n') # python will convert \n to os.linesep
	f.write('# created at %s\n' %(datetime.datetime.now()))
	f.write('\n')
	f.write('import numpy\n')
	f.write('import photolysisRates\n')
	f.write('\n')

	# following part is the function (there should be an indent at the start of each line)
	# suggest using one tab
	f.write('def evaluate_rates(RO2, H2O, RH, TEMP, time, M, N2, O2, Jlen, NO, HO2, NO3, sumt, self):\n')
	f.write('\n')
	f.write('	# inputs: ------------------------------------------------------------------\n')
	f.write('	# RO2 - total concentration of alkyl peroxy radicals (# molecules/cm3) \n')	
	f.write('	# M - third body concentration (# molecules/cm3 (air))\n')
	f.write('	# N2 - nitrogen concentration (# molecules/cm3 (air))\n')
	f.write('	# O2 - oxygen concentration (# molecules/cm3 (air))\n')
	f.write('	# H2O, RH, TEMP: given by the user\n')
	f.write('	# self.light_stat_now: given by the user and is 0 for lights off and >1 for on\n')
	f.write('	# reaction rate coefficients and their names parsed in eqn_parser.py \n')
	f.write('	# Jlen - number of photolysis reactions\n')
	f.write('	# self.tf - sunlight transmission factor\n')
	f.write('	# NO - NO concentration (# molecules/cm3 (air))\n')
	f.write('	# HO2 - HO2 concentration (# molecules/cm3 (air))\n')
	f.write('	# NO3 - NO3 concentration (# molecules/cm3 (air))\n')
	f.write('	# self.tf_UVC - transmission factor for 254 nm wavelength light (0-1) \n')
	f.write('	# ------------------------------------------------------------------------\n')
	f.write('\n')
	
	f.write('	erf = 0; err_mess = \'\' # begin assuming no errors\n')

	# determine whether to bypass calculations because no chemical reactions were found
	if (len(self.reac_coef_g)+len(self.reac_coef_aq)+len(self.reac_coef_su) > 0):
		
		f.write('\n')
		f.write('	# calculate any generic reaction rate \n')
		f.write('	# coefficients given by chemical scheme \n')
		f.write('\n')
		if self.rrc:
			f.write('	try:\n')
			f.write('		gprn=0\n')
			# code to calculate any generic rate coefficients 
			# given by chemical scheme file
			for line in self.rrc:
				f.write('		# keep count on reaction number \n')
				f.write('		gprn += 1 \n')
				f.write('		%s \n' %line)
			f.write('\n')
			f.write('	except:\n')
			f.write('		erf = 1 # flag error\n')
			f.write('		err_mess = str(\'Error: generic reaction rates failed to be calculated inside rate_coeffs.py at number \' + str(gprn) + \', please check chemical scheme and associated chemical scheme markers, which are stated in the model variables input file\') # error message\n')
			f.write('		return([], erf, err_mess)\n')
	
		f.write('	# estimate and append photolysis rates\n')
		f.write('	J = photolysisRates.PhotolysisCalculation(TEMP, Jlen, sumt, self)\n')
		f.write('\n')
		f.write('	if (self.light_stat_now == 0):\n')
		f.write('		J = [0]*len(J)\n')

		# calculate the rate coefficient for each equation
		f.write('	rate_values = numpy.zeros((%i))\n' %(len(self.reac_coef_g)+len(self.reac_coef_aq)+len(self.reac_coef_su)))
		
		f.write('	\n')
		f.write('	# if reactions have been found in the chemical scheme\n')
		f.write('	# gas-phase reactions\n')
		f.write('	gprn = 0 # keep count on reaction number\n')
		# in case pause needed
		#f.write('	import ipdb; ipdb.set_trace()\n')
		f.write('	try:\n') # in case there are any issues with calculating a rate coefficient
		for eqn_key in range (len(self.reac_coef_g)):
			f.write('		gprn += 1 # keep count on reaction number\n')
			f.write('		# remember equation in case needed for error reporting\n')
			f.write('		rc_eq_now = \'%s\' \n' %(self.reac_coef_g[eqn_key]))
			f.write('		rate_values[%s] = %s\n' %(eqn_key, self.reac_coef_g[eqn_key]))
		f.write('	except:\n') # in case there are any issues with calculating a rate coefficient
		f.write('		erf = 1 # flag error\n')
		f.write('		err_mess = (str(\'Error: Could not calculate \'+ \n')
		f.write('		\'rate coefficient for equation number \' \n')
		f.write('		+ str(gprn) + \' \' + rc_eq_now + \n')
		f.write('		\' (message from rate coeffs.py)\'))\n')
		
		f.write('	\n')
		f.write('	# aqueous-phase reactions\n')
		for eqn_key_aq in range (len(self.reac_coef_aq)):
			f.write('	rate_values[%s] = %s\n' %(len(self.reac_coef_g)+eqn_key_aq, self.reac_coef_aq[eqn_key_aq]))
		f.write('	\n')
		f.write('	# surface (e.g. wall) reactions\n')
		for eqn_key_su in range (len(self.reac_coef_su)):
			f.write('	rate_values[%s] = %s\n' %(len(self.reac_coef_g)+len(self.reac_coef_aq)+eqn_key_su, self.reac_coef_su[eqn_key_su]))
		f.write('	\n')
	
	else: # if no reactions found
		f.write('	rate_values = numpy.zeros((%i))\n' %(len(self.reac_coef_g)+len(self.reac_coef_aq)+len(self.reac_coef_su)))		
	f.write('	return(rate_values, erf, err_mess)\n')
	f.close()

	return()
