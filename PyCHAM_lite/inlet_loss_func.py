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
'''solving the size-dependent loss rate of particles in instrument inlet'''
# module to estimate loss rate (fraction /s) of particles during passage through an instrument inlet 
# File Created at 2021-08-16 12:38:02.687302

import numpy as np

# function for loss rate
def inlet_loss_func(Dp_all):
	
	# inputs: -------------------
	# Dp_all - diameter (um) of size bins (columns) per time step (rows)
	# -----------------------------
	
	sd_lrate = np.zeros((Dp_all.shape[0], Dp_all.shape[1]))
	try: # in case loss function string is acceptable
		# estimate loss rate (fraction/s)
		for it in range(Dp_all.shape[0]):
			Dp = Dp_all[it, :]
			sd_lrate[it, :] = 0.
	except: # in case of issue
		sd_lrate = 'Error, function of loss rate of particles during passage through inlet failed, please revise'
	
	return(sd_lrate)