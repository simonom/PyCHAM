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
'''estimating the natural light intensity for photochemistry'''
# equivalent to the method for AtChem2 for equations related to 
# solar zenith angle; note follows the equations of 
# Chapter 1 ("The Atmosphere and UV-B Radiation at 
# Ground Level" by S. Madronich) of the textbook
# Environmental UV Photobiology (1993), Young, Antony R., editor.

import numpy as np


def zenith(self):
    
	# inputs ------------------------------------------------------
	# self - reference to program
	# -------------------------------------------------------------

	time = self.sumt+self.daytime # time (s) through day

	pi = 4.*np.arctan(1.) # 1pi=180deg=3.14. For future calculation
	radian = 1.8e2/pi # 1 unit rad; use this to convert [\deg] to [\rad]

	# latitude conversion from degrees to radian
	lat_rad = self.lat/radian # in [\rad]
	
	# the day angle (radians), from calcTheta inside solarFunctions.f90 file of AtChem2,
	# assuming not a leap year (correct for 2010)
	theta = 2.*pi*self.dayOfYear/365.
	
	# inputs for solar declination angle (see reference above)
	b0 =  0.006918
	b1 = -0.399912
	b2 =  0.070257
	b3 = -0.006758
	b4 =  0.000907
	b5 = -0.002697
	b6 =  0.001480

	dec = b0 + b1*np.cos(theta) + b2*np.sin(theta) + b3*np.cos(2.*theta) + b4*np.sin(2.*theta) + b5*np.cos(3.*theta) + b6*np.sin(3.*theta)
	
	# initially used (prior to 08/07/2022) value for solar declination angle
	#dec = 0.41 # solar declination angle (rad), consistent with AtChem2 default
	
	# equation of time accounts for the discrepancy between the apparent and the mean
	# solar time at a given location
	c0 = 0.000075
	c1 = 0.001868
	c2 = -0.032077
	c3 = -0.014615
	c4 = -0.040849
	eqtime = c0 + c1*np.cos(theta)+c2*np.sin(theta)+c3*np.cos(2.0*theta)+c4*np.sin(2.0*theta)

	# get the fraction hour by dividing by seconds per hour, then remove any hours
	# from previous days by using the remainder function and 24 hours
	currentFracHour = np.remainder(time/3600., 24.)

	# local hour angle (lha): representing the time of the day - taken from 
	# solarFunctions.f90 of AtChem2
	lha = pi*((currentFracHour/12.)-(1.+self.lon/180.))+eqtime

	sinld = np.sin(lat_rad)*np.sin(dec)
	cosld = np.cos(lat_rad)*np.cos(dec)
   
	cosx = (np.cos(lha)*cosld)+sinld
	secx = 1.E+0/(cosx+1.E-30)

	# as in solarFunctions.f90 of AtChem2, set negative cosx to 0
	if (cosx < 0.):
		cosx = 0.
		secx = 100.
	
	return (secx, cosx)