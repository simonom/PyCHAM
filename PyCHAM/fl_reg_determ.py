''' module to estimate variables related to the fluidity of particles in air'''

# eq. nos. refer to Jacobson (2005)

import numpy as np 
import scipy.constants as si
import pdb

def reg_determ(RH, T, sbr, Pa):

	# ------------------------------------------------
	# inputs:
	# RH - relative humidity at this time (fraction) 
	# T - temperature (K)
	# sbr - size bin radii (m)
	# Pa - pressure inside chamber (Pa)
	
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
	Pv = (RH*Pa*Pvs)/(RH*Pvs+Pa-Pvs)
	# pressure of dry air (Pa) (2.22)
	Pd = Pa-Pv
	
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
	rho_a = Pa/(Rm*T)*1.0e3
	
	
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
		
	return Kni, eta_a, rho_a, kin_visc

