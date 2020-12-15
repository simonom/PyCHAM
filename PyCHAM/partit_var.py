'''module to estimate the particle and wall partitioning coefficient'''
# the kimt_calc module is called at the start of the model loop time interval to update
# the mass transfer coefficient of gases to particles

import numpy as np
from part_prop import part_prop
import scipy.constants as si

def kimt_calc(y, mfp, num_sb, num_comp, accom_coeff, y_mw, surfT, R_gas, TEMP, NA, 
		y_dens, N_perbin, radius, Psat, therm_sp,
		H2Oi, act_coeff, wall_on, caller, partit_cutoff, Press, DStar_org):
	
	# inputs:---------------------------------------------------------------------------
	
	# y - concentration of components' molecules (molecules/cc (air))
	# mfp - mean free path of gas molecules (m) (num_comp, 1)
	# accom_coeff - accommodation coefficients of components in each size bin
	# y_mw - molecular weight of components (g/mol)
	# radius - particle radius (m)
	# surfT - surface tension (g/s2==mN/m==dyn/cm)
	# R_gas - the universal gas constant (cc.Pa/K.mol)
	# TEMP - current temperature in chamber (K)
	# NA - Avogadro's constant (molecules/mol) 
	# y_dens - density of components (kg/m3)
	# N_perbin - number of particles per size bin (excluding wall)
	# Psat - liquid-phase saturation vapour pressures of components (molecules/cc (air))	
	# therm_sp - thermal speed of components (m/s) (num_comp, 1)
	# H2Oi - water index (integer)
	# act_coeff - activity coefficient of components (dimensionless)
	# wall_on - marker for whether wall present
	# caller - marker for the calling function
	# partit_cutoff - the product of Psat and act_coeff above which gas-particle 
	# 		partitioning assumed zero (Pa)
	# Press - air pressure (Pa)
	# DStar_org - gas-phase diffusion coefficient of components (cm2/s)
	# ------------------------------------------------------------------------------------
	
	if num_sb == 0: # fillers
		kimt = 0.
		kelv = 0.

		return(kimt, kelv)

	if (num_sb > 0) and (wall_on > 0): # if wall present
		y_part = y[num_comp:-(num_comp)]
	if (num_sb > 0) and (wall_on == 0): # if wall absent
		y_part = y[num_comp::]
	
	# density (g/cm3) and average molecular weight (g/mol) of particles (excluding wall)
	[tot_rho, ish, avMW] = part_prop(y_part, num_comp, (num_sb-wall_on), NA, y_mw, y_dens, 
					N_perbin)
	
	# Knudsen number (dimensionless)
	Kn = np.repeat(mfp, (num_sb-wall_on), 1)/np.repeat(radius, num_comp, 0)
	
	# update accommodation coefficients if necessary
	# note, using __import__ rather than import allows opening in run time, thereby using
	# updated module
	accom_coeff_calc = __import__('accom_coeff_calc')
	accom_coeff_now = accom_coeff_calc.accom_coeff_func(accom_coeff, radius)

	# Non-continuum regime correction 
	# calculate a correction factor according to the continuum versus non-continuum 
	# regimes
	# expression taken from Jacobson et al (2000), page 457, or Jacobson (2005), page 530 
	# (eq. 16.19).
	# They reference:
	# Fuchs and Sutugin 1971
	# Pruppacher and Klett 1997
	Inverse_Kn = Kn**-1.
	correct_1 = (1.33+0.71*Inverse_Kn)/(1.+Inverse_Kn)
	correct_2 = (4.*(1.-accom_coeff_now))/(3.*accom_coeff_now)
	correct_3 = 1.e0+(correct_1+correct_2)*Kn
	correction = correct_3**-1.
	
	# kelvin factor for each size bin (excluding wall), eq. 16.33 Jacobson et al. (2005)
	# note that avMW has units g/mol, surfT (g/s2==mN/m==dyn/cm), R_gas is multiplied by 
	# 1e7 for units g cm2/s2.mol.K, 
	# TEMP is K, radius is multiplied by 1e2 to give cm and tot_rho is g/cm3
	kelv = np.zeros((num_sb-wall_on, 1))
	kelv[ish, 0] = np.exp((2.e0*avMW[ish]*surfT)/(R_gas*1.e7*TEMP*radius[0, ish]*1.e2*tot_rho[ish]))


	# ------------------------------------------------
	# gas phase diffusion coefficient*Fuch-Sutugin correction (cm2/s)
	# eq. 5 Zaveri et al. (2008)
	kimt = (DStar_org)*correction
	# final partitioning coefficient (converting radius from m to cm)
	# eq. 16.2 of Jacobson (2005) and eq. 5 Zaveri et al. (2008)
	# components in rows and size bins in columns (/s)
	kimt = (4.*np.pi*(radius*1.e2)*N_perbin.reshape(1, -1))*kimt
	# transpose kimt ready for multiplication inside ode solver, so size bins
	# in rows and components in columns
	kimt = np.transpose(kimt)

	# zero partitioning to particles for any components with low condensability
	if (partit_cutoff): # if a value provided (default is empty list)

		# convert partit_cutoff from Pa to molecules/cc (air), note README states
		# that just one value accepted for partit_cutoff input
		partit_cutoff_Pa = partit_cutoff[0]*(NA/((R_gas*1.e6)*TEMP))
		highVPi = (Psat*act_coeff) > partit_cutoff_Pa
		highVPi[:, H2Oi] = 0 # mask water to allow its partitioning
		kimt[highVPi] = 0.

	# zero partitioning coefficient for size bins where no particles - enables significant
	# computation time acceleration and is physically realistic
	ish = N_perbin[:, 0]<=1.e-10
	kimt[ish, :] = 0.
	
	return(kimt, kelv)
