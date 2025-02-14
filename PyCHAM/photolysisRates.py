########################################################################
#								       #
# Copyright (C) 2018-2025					       #
# Simon O'Meara : simon.omeara@manchester.ac.uk			       #
#								       #
# All Rights Reserved.                                                 #
# This file is part of PyCHAM                                          #
#                                                                      #
# PyCHAM is free software: you can redistribute it and/or modify it    #
# under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or    #
# (at  your option) any later version.                                 #
#                                                                      #
# PyCHAM is distributed in the hope that it will be useful, but        #
# WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# General Public License for more details.                             #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with PyCHAM.  If not, see <http://www.gnu.org/licenses/>.      #
#                                                                      #
########################################################################
'''module for calculating photolysis rates'''
# can use ambient sunlight expression or artificial chamber lights for 
# calculating photolysis rates

import scipy
import os
import numpy as np
import shutil
from lamp_photo import lamp_photo
import zenith

def PhotolysisCalculation(TEMP, Jlen, sumt, self):
	
	# inputs:------------------------------------------------------
	# self.daytime - time of day experiment starts (for natural light photolysis) (s)
	# self.lat - latitude
	# self.lon - longitude
	# TEMP - temperature inside chamber (K)
	# self.af_path - name of path to file containing known actinic flux 
	#	(only used if lights on inside chamber)
	# self.DayOfYear - number of days through the calendar year
	# self.photo_path - name of file containing estimates for wavelength-dependent
	# 	absorption cross-sections and quantum yields
	# Jlen - number of photolysis reactions
	# self.tf - the transmission factor (for light intensity)
	# sumt - total time through experiment (s)
	# self.tf_UVC - transmission factor for 254 nm wavelength light (0-1)
	# -----------------------------------------------------------
	
	self.sumt = sumt
	
	# if dark
	if (sum(self.light_stat == 3) == 0):
		J = np.zeros((Jlen)) # prepare output
	else: # if a stored array of photolysis rates
		J = self.stored_J[:, sum(self.light_time<=sumt)-1]
		# if a wavelength-independent transmission factor
		if (self.tf_range == 0):
			J = J*self.tf
		return(J)  

	cwd = os.getcwd() # address of current working directory
	
	# if using MCM chemical scheme and natural light
	if (self.photo_path == str(cwd + '/PyCHAM/photofiles/MCMv3.2') 
		and self.af_path == 'no'):
	
		# get solar zenith angle following the equations of 
		# Chapter 1 ("The Atmosphere and UV-B Radiation at 
		# Ground Level" by S. Madronich) of the textbook
		# Environmental UV Photobiology (1993)
		
		# if secx and cosx specified by user in the model variables 
		# file (to keep solar intensity constant), then use these
		if hasattr(self, 'secx') and hasattr(self, 'cosx'):
			secx = self.secx
			cosx = self.cosx
		
		else: # otherwise estimate current secx and cosx
			(secx, cosx) = zenith.zenith(self)
		
		# The Hayman (1997) parameterisation for MCM reactions as described in
		# Saunders et al. (2003): https://doi.org/10.5194/acp-3-161-2003
		#J          L           M          N
		J[1] = 6.073E-05*cosx**(1.743)*np.exp(-1.0*0.474*secx)
		J[2] = 4.775E-04*cosx**(0.298)*np.exp(-1.0*0.080*secx)
		J[3] = 1.041E-05*cosx**(0.723)*np.exp(-1.0*0.279*secx)
		J[4] = 1.165E-02*cosx**(0.244)*np.exp(-1.0*0.267*secx)
		J[5] = 2.485E-02*cosx**(0.168)*np.exp(-1.0*0.108*secx)
		J[6] = 1.747E-01*cosx**(0.155)*np.exp(-1.0*0.125*secx)
		J[7] = 2.644E-03*cosx**(0.261)*np.exp(-1.0*0.288*secx)
		J[8] = 9.312E-07*cosx**(1.230)*np.exp(-1.0*0.307*secx)
		J[11] = 4.642E-05*cosx**(0.762)*np.exp(-1.0*0.353*secx)
		J[12] = 6.853E-05*cosx**(0.477)*np.exp(-1.0*0.323*secx)
		J[13] = 7.344E-06*cosx**(1.202)*np.exp(-1.0*0.417*secx)
		J[14] = 2.879E-05*cosx**(1.067)*np.exp(-1.0*0.358*secx)
		J[15] = 2.792E-05*cosx**(0.805)*np.exp(-1.0*0.338*secx)
		J[16] = 1.675E-05*cosx**(0.805)*np.exp(-1.0*0.338*secx)
		J[17] = 7.914E-05*cosx**(0.764)*np.exp(-1.0*0.364*secx)
		J[18] = 1.140E-05*cosx**(0.396)*np.exp(-1.0*0.298*secx)
		J[19] = 1.140E-05*cosx**(0.396)*np.exp(-1.0*0.298*secx)
		J[21] = 7.992E-07*cosx**(1.578)*np.exp(-1.0*0.271*secx)
		J[22] = 5.804E-06*cosx**(1.092)*np.exp(-1.0*0.377*secx)
		J[23] = 1.836E-05*cosx**(0.395)*np.exp(-1.0*0.296*secx)
		J[24] = 1.836E-05*cosx**(0.395)*np.exp(-1.0*0.296*secx)
		J[31] = 6.845E-05*cosx**(0.130)*np.exp(-1.0*0.201*secx)
		J[32] = 1.032E-05*cosx**(0.130)*np.exp(-1.0*0.201*secx)
		J[33] = 3.802E-05*cosx**(0.644)*np.exp(-1.0*0.312*secx)
		J[34] = 1.537E-04*cosx**(0.170)*np.exp(-1.0*0.208*secx)
		J[35] = 3.326E-04*cosx**(0.148)*np.exp(-1.0*0.215*secx)
		J[41] = 7.649E-06*cosx**(0.682)*np.exp(-1.0*0.279*secx)
		J[51] = 1.588E-06*cosx**(1.154)*np.exp(-1.0*0.318*secx)
		J[52] = 1.907E-06*cosx**(1.244)*np.exp(-1.0*0.335*secx)
		J[53] = 2.485E-06*cosx**(1.196)*np.exp(-1.0*0.328*secx)
		J[54] = 4.095E-06*cosx**(1.111)*np.exp(-1.0*0.316*secx)
		J[55] = 1.135E-05*cosx**(0.974)*np.exp(-1.0*0.309*secx)
		J[56] = 7.549E-06*cosx**(1.015)*np.exp(-1.0*0.324*secx)
		J[57] = 3.363E-06*cosx**(1.296)*np.exp(-1.0*0.322*secx)
		J[61] = 7.537E-04*cosx**(0.499)*np.exp(-1.0*0.266*secx)
		
		# if tranmission factor of light
		# that is not wavelength dependent
		if (self.tf_range == 0):
			J = J*self.tf
		# if tranmission factor of light
		# is wavelength dependent, then assume the effect
		# on J as calculated by Hayman (1997) is that of
		# solar transmittance through a window (clear 
		# glass in Figure 1b of 
		# doi.org/10.1007/s44223-022-00009-6). Simon 
		# O'Meara (11/06/2024) estimated this effect 
		# by getting the ratio
		# of J w/ and w/o transmittance impediment where
		# w/ impedimenent is using the window 
		# transmittance cited above and w/o is 
		# transmission factor=1. The J was calculated using 
		# the MAC_Actinic_Flux_Spectrum_wUVCx3p5.csv file
		# from PyCHAM/photofiles, which is estimated in
		# lamp_photo.py. Note that J[57] and J[61]
		# are not estimated in lamp_photo, but nor are they 
		# present in the full MCM chemical scheme, so their
		# transmittance factor is assumed to be 0.
		# The resulting ratio for every J 
		# is: 
		if (self.tf_range == 2):
			Hayman_wd_tf_clear_glass = [0, 0.173121721, 0.707269159, 
			0.826452633, 0.870432426, 0.838243408, 0.892998655, 
			0.858393338, 0.620252851, 0, 0, 0.442839373, 
			0.598900199, 0.258290129, 0.232623684, 0.446158007, 
			0.446158007, 0.416693792, 0.795799329, 0.795799329, 
			0.795799329, 0.170721083, 0.287727411, 0.559692167, 
			0.559692166, 0, 0, 0, 0, 0, 0, 0.413105058, 0.544067465, 
			0.773449079, 0.74310927, 0.894902116, 0, 0, 0, 0, 0, 
			0.813654941, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.808833912, 
			0.733110394, 0.801542645, 0.741646608, 0.829273355, 
			0.821086824, 0, 0, 0, 0, 0]
			J = J*self.tf*Hayman_wd_tf_clear_glass

	# if a file path for user-supplied absorption cross-sections
	# and quantum yields are supplied and actinic flux is
	# based on the Madronich equations for actinic flux from
	# solar irradiation
	
	# if a file path for the actinic flux is supplied,
	# (e.g. if chamber lamps are turned on)
	# note this option is for user-supplied actinic flux and either a
	# user-supplied absorption cross-section and quantum yield 
	# file or MCM recommended absorption cross-section and 
	# quantum yields in combination with MCM photolysis reactions 
	# (http://mcm.leeds.ac.uk/MCMv3.3.1/parameters/photolysis.htt)
	if (self.af_path != 'no'):
		# call on module to process photolysis files and estimate J values
		J = lamp_photo(J, TEMP, self)
	
	# in case a print out of photolysis rate to command line needed
	#Jcn = 0
	#for Jn in J:
	#	print(Jcn, Jn)
	#	Jcn += 1
	#import ipdb; ipdb.set_trace()
	
	return(J)
