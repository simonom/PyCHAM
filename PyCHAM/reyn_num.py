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
'''module to calculate Reynold number'''

import numpy as np
import scipy.constants as si

def Reyn_num(sbr, eta_a, rho_a, kin_visc, sbrho, Kni):

	# --------------------------------------------------------
	# inputs:
	# sbr - particle radius (m)
	# eta_a - dynamic viscosity of air (4.54) (g/m.s)	
	# rho_a - density of air (g/m3)
	# kin_visc - kinematic viscosity of air (m2/s)
	# sbrho - average density of particles (g/m3)
	# Kni - Knudsen number (dimensionless)

	# --------------------------------------------------------
	# outputs:
	# Re - Reynold number (dimensionless)
	# Vf - terminal fall speed (m/s)
	# --------------------------------------------------------


	# Cunningham slip-flow correction (15.30) with constant taken
	# from text below
	G = 1.0+Kni*(1.249+0.42*(np.exp(-0.87/Kni)))
	
	# terminal fall speed (m/s) (20.4)
	Vf = ((2.0*sbr**2.0*(sbrho-rho_a)*si.g)/(9.0*eta_a))*G
	# particle Reynolds number (15.26) (dimensionless)
	Re = (2.0*sbr*Vf)/kin_visc
	
	return Re, Vf

