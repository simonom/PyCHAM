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
''' module to estimate variables related to the fluidity of particles in air'''

# eq. nos. refer to Jacobson (2005)

import numpy as np 
import scipy.constants as si
import pdb

def reg_determ(RH, T, sbr, self):

	# ------------------------------------------------
	# inputs:
	# RH - relative humidity at this time (fraction) 
	# T - temperature (K)
	# sbr - size bin radii (m)
	# self.Pressn - air pressure now (Pa)
	
	# ------------------------------------------------
	# outputs:
	# Kni - Knudsen number of particles
	# rho_a - mass density of moist air
	# eta_a - dynamic viscosity (4.54) (g/m.s)
	# kin_visc - kinematic viscosity of air (4.55) (m2/s)	
	# ------------------------------------------------

	# saturation vapour pressure of water vapour (Pa (kg/m.s2) 
	# (2.61))
	Pvs = (6.112e2*np.exp(6816*(1.0/273.15-1.0/T)+
		(5.1309*np.log(273.15/T))))
	
	# vapour pressure (partial pressure) of water vapour 
	# (Pa (kg/m.s2) (2.66)
	Pv = (RH*self.Pressn*Pvs)/(RH*Pvs+self.Pressn-Pvs)
	# pressure of dry air (Pa) (2.22)
	Pd = self.Pressn-Pv
	
	# ratio of molecular weights of water and dry air (p. 30 for Mw air
	# & CRC for Mw water vapour) (2.31)
	epsilon = (18.015/28.966)
	# mass mixing ratio of water vapour (2.31)
	omega = epsilon*(Pv/Pd)
	# dry-air gas constant (m2/s2.K)  
	# (2.24 (R* value from p. 29 (converted to units 
	# m2/s.K}))
	R_dash = (8314.51/28.966)
	# gas constant for moist air (2.37) (m3Pa/g.K)	
	Rm = R_dash*((1.0+omega/epsilon)/(1.0+omega))	
	# mass density of moist air (g/m3) (2.36 Jacobson 2005)
	# (multiply by 1.0e3) to convert from kg/m3 to g/m3
	rho_a = self.Pressn/(Rm*T)*1.0e3
	
	
	# dynamic viscosity of air (g/m.s) (4.54) 
	# (p. 30 for Mw air and p. 102 for av. diameter of air molecule (m))
	eta_a = (5.0/(16.0*si.N_A*(3.673e-10**2.0)))*(((28.966*8314.51*T)
		/(np.pi))**0.5)
	# kinematic viscosity of air (m2/s) (4.55)
	kin_visc = eta_a/rho_a

	# thermal speed of air molecule (m/s) (2.3) (15.25) 
	# (4.8096e-26 is average mass of one air molecule (p. 18))
	nu_a = ((8.0*si.k*T)/(np.pi*4.8096e-26))**(0.5)	
	# mean free path of an air molecule (m) (15.24)
	lamb = (2.0*kin_visc)/nu_a

	# Knudsen number for particles in air (15.23) - determines flow regime
	Kni = lamb/sbr	
	
	return(Kni, eta_a, rho_a, kin_visc)

