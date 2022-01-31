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
'''radius response to humidity'''
# a module to alter the sizes of particles due to conditioning outside the chamber,
# for example, for an inlet to an instrument with a different relative humidity to the 
# chamber

import numpy as np
import scipy.constants as si

# define function
def rad_resp_hum(yrec, x, dryf, H2Oi, num_comp, num_asb, Nwet, y_MV):

	# inputs: -------------------------------------------
	# yrec - concentration (molecules/cm3) of components per particle size bin (rows) with time (columns)
	# x - original wet radius of particles (um)
	# dryf - relative humidity to equilibrate to (due to conditioning inside instrument)
	# H2Oi - index of water
	# num_comp - number of components
	# num_asb - number of particle size bins (excluding wall)
	# Nwet - number concentration of particles (# particles/cm3)
	# y_MV - molar volumes of components (cm3/mol)
	# -----------------------------------------------------
	
	# holder for particle-phase concentrations
	yn = np.zeros((yrec.shape[0], yrec.shape[1]))
	yn[:, :] = yrec[:, :]
	# zero particle-phase water (molecules/cm3)
	yn[:, H2Oi::num_comp] = 0.
	
	# tile molar volumes (cm3/mol) over times
	y_MV = np.array(y_MV)
	y_MV = np.tile(y_MV.reshape(1, -1), [yn.shape[0], 1])
	# tile molar volumes over size bins
	y_MV = np.tile(y_MV, [1, num_asb])
	# convert from cm3/mol to um3/mol
	y_MV = y_MV*1e12
	
	# repeat particle number concentration over components (# particles/cm3)
	Nwet_rep = np.repeat(Nwet, num_comp, axis = 1)
	
	
	# factor to convert component particle-phase concentrations to volume (um3/(molecules/cm3))
	vol_fac = np.zeros((Nwet_rep.shape[0], Nwet_rep.shape[1]))
	nz_indx = Nwet_rep != 0. # non-zero index
	vol_fac[nz_indx] = (y_MV[nz_indx]/(Nwet_rep[nz_indx]*si.N_A))*(3./(4.*np.pi))
	
	# empty array to hold radius of particles equilibrated 
	# with instrument relative humidity (um)
	xn = np.zeros((x.shape[0], x.shape[1]))
	
	for isb in range(num_asb): # loop through size bins
		# sum all non-water components in this size bin for all time steps (molecules/cm3)
		ysum = yn[:, (isb*num_comp):((isb+1)*num_comp)].sum(axis=1)
		# equilibrate water in instrument (molecules/cm3)
		yn[:, (isb*num_comp+H2Oi)] = ysum*(dryf/(1.-dryf))
		# new volume of single particles (um3)
		voln = ((yn[:, (isb*num_comp):((isb+1)*num_comp)]*vol_fac[:, (isb*num_comp):((isb+1)*num_comp)]).sum(axis=1))
		# new radius of particles (um)
		xn[:, isb] = voln**(1./3.)
		
		# fill radii of size bins without particles with default values (um)
		xn[:, isb][Nwet[:, isb] < 1.e-3] = x[:, isb][Nwet[:, isb] < 1.e-3] 
	
	return(xn, yn)