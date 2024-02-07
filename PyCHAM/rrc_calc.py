##########################################################################################
#                                                                                        #
#    Copyright (C) 2018-2024 Simon O'Meara : simon.omeara@manchester.ac.uk               #
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
'''module to link to calculation of reaction coefficients'''
# this module sets up the final details for calculating
# reaction rate coefficient, which is done via the called module Rate_coeffs 


import numpy as np
import scipy.constants as si
try:
	import rate_coeffs
except:
	import os
	if os.path.exists('rate_coeffs'): # remove any bad functions
		os.remove(rate_coeffs)
import importlib


def rrc_calc(H2O, TEMP, y, PInit, Jlen, NO, HO2, NO3, sumt, self):

	import rate_coeffs # in case failure to import previous version using import command above

	# ---------------------------------------------
	# inputs:
	# self.RO2_indices - indices of RO2 components
	# H2O - concentrations of water (# molecules/cm3)
	# TEMP - temperature (K)
	# self.light_status_now - whether lights off (0) or on (1)
	# y - component concentrations (# molecules/cm3 (air))
	# self.lat - latitude
	# self.lon - longitude
	# self.af_path - path to file stating actinic flux
	# self.DayOfYear - day number of the year (1-365)
	# PInit - chamber pressure (Pa)
	# self.photo_file - name of file with with estimates for photolysis absorption
	# 					cross-sections and quantum yields
	# Jlen - number of photolysis reactions
	# self.tf - sunlight transmission factor
	# NO - concentration of NO (# molecules/cm3)
	# HO2 - concentration of HO2 (# molecules/cm3)
	# NO3 - concentration of NO3 (# molecules/cm3)
	# self.tf_UVC - transmission factor for 254 nm wavelength light (0-1)
	# ---------------------------------------------
	
	time = self.daytime+sumt # time of day (for natural light photolysis) (s)
	
	# start by assuming no error message
	erf = 0
	err_mess = ''
	
	# calculate total RO2 concentration
	if (self.RO2_indices.size == 0):
		RO2 = 0
	else:
		RO2 = np.sum(y[self.RO2_indices[:, 1]])
		
	# calculate concentrations of third body (M), nitrogen and oxygen
	# calculate gas-phase concentrations of M, N2 and O2 (# molecules/cm3 (air))
	# 1.0e-6 is the 1 cm3 volume expressed as m3
	# R and Avogadro's constant set the same as in atmosphereFunctions.f90 of AtChem2
	M_val = ((PInit*1.e-6)/(8.3144621*TEMP))*si.N_A
	
	# N2 and O2 given the same multiplication as in atmosphereFunctions.f90 of AtChem2
	N2_val = M_val*0.7809
	O2_val = M_val*0.2095
	
	importlib.reload(rate_coeffs) # ensure latest version uploaded
	# calculate the new rate coefficient array (/s) 
	[rrc, erf, err_mess] = rate_coeffs.evaluate_rates(RO2, H2O, TEMP, time, 
					M_val, N2_val, O2_val, Jlen, NO, HO2, NO3, sumt, self)
	
	return(rrc, erf, err_mess)
