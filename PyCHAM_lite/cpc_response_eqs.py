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
'''solving the weighting of particles of different ages over response time of instrument'''
# module to estimate the weighting of particles of different ages during the response time of instrument to represent the mixing of particles of different ages due to differential flow prior to counter 
# File Created at 2021-08-16 12:38:02.703258

import numpy as np

# function for weighting
def cpc_response(delays, wfuncs):
	
	# inputs: -----------------
	# ---------------------------
	
	# remember all times (s) 
	t_all = np.zeros((100)) 
	sht = 1.0 # shortest delay
	resp_timepeak = 1.0 # delay at peak weighting
	# time range for increasing weighting (s)
	if ((resp_timepeak-sht) > 0):
		t = np.arange(sht, (resp_timepeak), (resp_timepeak-sht)/50.)
	if (resp_timepeak == sht):
		t = np.ones((50))*resp_timepeak
	wpre = 1.*t 
	t_all[0:50] = t[:] 
	
	lot = 1.0 # longest delay (s)
	# time range for decreasing weighting (s)
	if ((lot-resp_timepeak) > 0):
		t = np.arange((resp_timepeak+(lot-resp_timepeak)/50.), (lot), (lot-resp_timepeak)/51.)
	if (resp_timepeak == lot):
		t = np.ones((50))*resp_timepeak
	wpro = 1.*t 
	t_all[50::] = t[:] 
	
	# join weighting
	w = np.append(wpre, wpro)
	
	# integrate weight curve
	area = np.trapz(w, t_all)
	# normalise so that integral is one
	if (area > 0):
		w = w/area
	
	return(w, t_all)