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
'''plot particle number size distributions to be used in simulation'''
# this module receives the seed particle number size distribution inputs 
# from the GUI and demonstrates these distributions graphically.  E.g.
# when a user wants to check their supplied values

import numpy as np
import matplotlib.pyplot as plt
import part_nsd # determining particle number size distributions

def plotter_nsd(lowersize, num_asb, uppersize, mean_rad, std, pmode, pconc, 
		space_mode, testf, pconct): # define module
	
	# inputs: ----------------------------
	# lowersize - the smallest bound on the particle size range (um)
	# num_asb - number of size bins (excluding wall)
	# uppersize - the largest bound on the particle size range (um)
	# mean_rad - mean radius of modes (um)
	# std - standard deviation of modes (must be greater than 1)
	# pmode - whether distribution described as modes or explicitly
	# pconc - particle number concentration (# particles/cm3)
	# space_mode - how to space out particles (logarithmically or 
	#		linearly)
	# testf - flag for whether in testing mode
	# pconct - the timings of particle injection
	# ------------------------------------
	
	# prepare figure
	plt.ion() # show results to screen and turn on interactive mode
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))

	# loop through number of supplied number size distributions
	for ti in range(len(pconct[0])):
		
		# get inputs at this time
		mean_radn = mean_rad[:, ti]
		stdn = std[:, ti]
		pconcn = pconc[:, ti]
	
		# get particle number concentration per size bin (# particles/cm3)
		# and radius at particle size bin centres (um)
		[N_perbin, x, rbou, Vbou, Varr, 
			upper_bin_rad_amp] = part_nsd.part_nsd(lowersize, 
			num_asb, uppersize, mean_radn, stdn, pmode, pconcn, space_mode, testf)

		# prepare dN/dlogDp
		# don't use the first boundary as it could be zero, which will error when log10 taken
		log10D = np.log10(rbou[1::]*2.)

		if (num_asb > 1) :
			# note, can't append zero to start of log10D to cover first size bin as the log10 of the
			# non-zero boundaries give negative results due to the value being below 1, so instead
			# assume same log10 distance as the next pair
			log10D = np.append((log10D[0]-(log10D[1]-log10D[0])), log10D)
			# radius distance covered by each size bin (log10(um))
			dlog10D = (log10D[1::]-log10D[0:-1])
		if (num_asb == 1): # single particle size bin
			# assume lower radius bound is ten times smaller than upper
			dlog10D = (log10D[0]-np.log10((rbou[1]/10.)*2.))	
		
		dNdlog10D = np.zeros((num_asb))
		dNdlog10D[:] = N_perbin[:, 0]/dlog10D

		# plot this number size distribution against diameter
		ax0.loglog(x*2., dNdlog10D, label = str('Time = ' + str(pconct[0][ti]) + ' s'))
		
	ax0.legend() # show legend
	ax0.set_xlabel(str('Particle Diameter (' + u'\u03BC' + 'm)'), size = 14)
	ax0.set_ylabel(str('dN (' + u'\u0023' + ' particles cm' + u'\u207B' + u'\u00B3' + ')/dlog' + u'\u2081' + u'\u2080' + '(Dp)'), size=14)
	#   ), size=14)

	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')	


	return() 