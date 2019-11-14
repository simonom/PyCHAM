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
