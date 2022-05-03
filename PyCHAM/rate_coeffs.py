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
'''module for calculating reaction rate coefficients (automatically generated)'''
# module to hold expressions for calculating rate coefficients # 
# created at 2022-05-03 15:06:40.250481

import numpy
import photolysisRates

def evaluate_rates(RO2, H2O, TEMP, lightm, time, M, N2, O2, Jlen, NO, HO2, NO3, sumt, self):

	# inputs: ------------------------------------------------------------------
	# RO2 - total concentration of alkyl peroxy radicals (molecules/cm3) 
	# M - third body concentration (molecules/cm3 (air))
	# N2 - nitrogen concentration (molecules/cm3 (air))
	# O2 - oxygen concentration (molecules/cm3 (air))
	# H2O, TEMP: given by the user
	# lightm: given by the user and is 0 for lights off and 1 for on
	# reaction rate coefficients and their names parsed in eqn_parser.py 
	# Jlen - number of photolysis reactions
	# self.tf - sunlight transmission factor
	# NO - NO concentration (# molecules/cm3 (air))
	# HO2 - HO2 concentration (# molecules/cm3 (air))
	# NO3 - NO3 concentration (# molecules/cm3 (air))
	# self.tf_UVC - transmission factor for 254 nm wavelength light (0-1) 
	# ------------------------------------------------------------------------

	erf = 0; err_mess = '' # begin assuming no errors
	# calculate any generic reaction rate coefficients given by chemical scheme

	# estimate and append photolysis rates
	J = photolysisRates.PhotolysisCalculation(TEMP, Jlen, sumt, self)

	if (lightm == 0):
		J = [0]*len(J)
	rate_values = numpy.zeros((11))
	
	# reac_coef has been formatted so that python can recognize it
	# gas-phase reactions
	gprn = 0 # keep count on reaction number
	try:
		gprn += 1 # keep count on reaction number
		rate_values[0] = 6.3e-16*numpy.exp(-580/TEMP)*0.57
		gprn += 1 # keep count on reaction number
		rate_values[1] = 6.3e-16*numpy.exp(-580/TEMP)*0.37
		gprn += 1 # keep count on reaction number
		rate_values[2] = 1.2e-12*numpy.exp(490/TEMP)*0.65
		gprn += 1 # keep count on reaction number
		rate_values[3] = 1.2e-12*numpy.exp(490/TEMP)*0.35
		gprn += 1 # keep count on reaction number
		rate_values[4] = 1.2e-11*numpy.exp(440/TEMP)*0.482
		gprn += 1 # keep count on reaction number
		rate_values[5] = 1.2e-11*numpy.exp(440/TEMP)*0.293
		gprn += 1 # keep count on reaction number
		rate_values[6] = 1.2e-11*numpy.exp(440/TEMP)*0.065
		gprn += 1 # keep count on reaction number
		rate_values[7] = 1.2e-11*numpy.exp(440/TEMP)*0.08
		gprn += 1 # keep count on reaction number
		rate_values[8] = 6.3e-16*numpy.exp(-580/TEMP)*0.06
		gprn += 1 # keep count on reaction number
		rate_values[9] = 6.3e-16*numpy.exp(-580/TEMP)*0.06
		gprn += 1 # keep count on reaction number
		rate_values[10] = 6.3e-16*numpy.exp(-580/TEMP)*0.06
	except:
		erf = 1 # flag error
		err_mess = str('Error: estimating reaction rate for reaction number ' + str(gprn) + ' failed, please check chemical scheme (including whether definitions for generic rate coefficients have been included), and associated chemical scheme markers, which are stated in the model variables input file') # error message
	
	# aqueous-phase reactions
	
	return(rate_values, erf, err_mess)
