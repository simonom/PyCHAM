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
'''interpreting and implementing particle number size distribution inputs'''
# a module that takes the variables relevant to particle number size 
# distributions and outputs the particle number concentration per size bin

import numpy as np
import size_distr # for allocating particle number size distributions

def part_nsd(lowersize, num_asb, uppersize, mean_radn, stdn, pmode, pconcn, space_mode, testf):

	# inputs: ----------------------------
	# lowersize - the smallest bound on the particle size range (um)
	# num_asb - number of size bins (excluding wall)
	# uppersize - the largest bound on the particle size range (um)
	# mean_radn - mean radius of modes (um)
	# stdn - standard deviation of modes (must be greater than 1)
	# pmode - whether distribution described as modes or explicitly
	# pconcn - particle number concentration (# particles/cm3)
	# space_mode - how to space out particles (logarithmically or 
	#		linearly)
	# testf - flag for whether in testing mode
	# ------------------------------------

	# if lower bound of particle sizes set to 0, this will cause an error 
	# when taking log10, so change to very small value (um)
	if (lowersize == 0.):
		lowersize = 9.e-4
	
	# if mean radius not stated explicitly calculate from size ranges (um)
	if (num_asb > 0):
		
		if (any(mrn == -1.e6 for mrn in mean_radn)):
			if (lowersize > 0.):
				mean_radn[mean_radn == -1.e6] = [10**((np.log10(lowersize)+np.log10(uppersize))/2.)]
			if (lowersize == 0.):
				mean_radn[mean_radn == -1.e6] = [10**((np.log10(uppersize))/2.)]

	# if multiple size bins, this call will assume a lognormal distribution if initial 
	# particle concentration is described by mode, or will assign particles to size bins if
	# initial particle concentration per size bin provided
	if (num_asb > 1):
			
		# set scale and standard deviation input for lognormal probability distribution 
		# function, following guidance here: 
		# http://all-geo.org/volcan01010/2013/09/how-to-use-lognormal-distributions-in-python/
		# if mean_radn and stdn are not already arrays then transform to a list 
		if type(mean_radn) != np.ndarray:
			scale = [np.exp(np.log(mean_radn))]
		else:
			scale = np.exp(np.log(mean_radn))
		
		if type(stdn) != np.ndarray:
			stdn = [np.log(stdn)]
		else:
			stdn = np.log(stdn)
			
		loc = 0. # no shift
		
		[N_perbin, x, rbou, Vbou, Varr, upper_bin_rad_amp] = size_distr.lognormal(num_asb, 
			pmode, pconcn, stdn, lowersize, uppersize, loc, scale, space_mode)
		
		if (testf == 2):
			print('finished with size_distr.lognormal')
		
	if (num_asb == 1):

		N_perbin = np.array((sum(pconcn))).reshape(-1, 1) # (# particles/cm3 (air) at experiment start)
		x = np.zeros(1) # radii at size bin centre
		# mean radius of this one size bin (um)
		try: # if mean_radn an array
			if (len(mean_radn[0]) == 1): # if a scalar
				meansize = mean_radn[0][0]
			else : # if an array
				meansize = sum(mean_radn)/len(mean_radn)
		except: # if mean_radn a scalar
			meansize = mean_radn
			
		x[0] = meansize

		# extend uppersize to reduce chance of particles growing beyond this
		upper_bin_rad_amp = 1.e6
		uppersize = uppersize*upper_bin_rad_amp
		# volume bounds of size bin (um3)
		Vbou = np.array(((lowersize**3.)*(4./3.)*np.pi, 
						(uppersize**3.)*(4./3.)*np.pi))
		# volume of single particle (um3)
		Varr = np.zeros((1, 1))
		Varr[0] = (4./3.)*np.pi*(x[0]**3.)
		# radius bounds of size bin (um)
		rbou = ((Vbou*3.0)/(4.0*np.pi))**(1.0/3.0)	

	return(N_perbin, x, rbou, Vbou, Varr, upper_bin_rad_amp)