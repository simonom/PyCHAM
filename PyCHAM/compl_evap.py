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
'''module to allow complete evaporation of the smallest size bin when conditions set in Vchange_check.py met'''
# rather than reduce the time step to unpractical length, allow the smallest size bin
# to completely evaporate using this module

def compl_evap(ytest, n0, Vnew, Vol0, nc, sbi):

	# inputs: ----------------------------------------
	# ytest - concentration (molecules/cc (air)) of components
	# n0 - number concentration of particles per size bin (# particles/cc (air))
	# Vnew - new estimates of single particle volume (um3)
	# Vol0 - default volume per size bin (um3), volume at centre of size bin bounds
	# nc - number of components
	# sbi - index of size bins in question
	# ------------------------------------------------

	# allow complete evaporation of the smallest size bin
	# components move to gas-phase (molecules/cc (air))
	ytest[0:nc] += ytest[(sbi+1)*nc:(sbi+2)*nc]
	# effectively zero concentration in particle-phase (molecules/cc (air))
	ytest[(sbi+1)*nc:(sbi+2)*nc] = 1.0e-40
	# remove particle concentration from size bin (# particles/cc (air))
	n0[sbi] = 1.0e-40
	# with no particles in smalles size bin now, default to volume at mid-point between 
	# size bin bounds (um3)
	Vnew[0] = Vol0[0]
	
	return(ytest, n0, Vnew)
