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
'''module to calculate gas-phase concentration of water'''
# based on input relative humidity and temperature the
# concentration of water in the gas-phase is estimated

import numpy as np

def water_calc(TEMP, RH, NA):

	# saturation vapour pressure (Pa) (TEMP in K.)
	Psat_water = (np.exp((-0.58002206E4/TEMP)+0.13914993E1-(0.48640239E-1*TEMP) 
		+ (0.41764768E-4*(TEMP**2.e0))-(0.14452093E-7*(TEMP**3.e0))+ 
		(0.65459673e1*np.log(TEMP))))

	# convert saturation vapour pressures from Pa to molecules/cc (air) using ideal
	# gas law, R has units cc.Pa/K.mol
	H2Oc = (RH*Psat_water)*(NA/(8.3144598e6*TEMP))
	
	# convert Psat_water to log10(atm) from Pa
	Psat_water = np.log10(Psat_water/101325.0)
	H2O_mw = 18.0 # state molecular weight of water (g/mol)
	
	return(H2Oc, Psat_water, H2O_mw)
