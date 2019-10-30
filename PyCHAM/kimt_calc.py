'''module to estimate the particle and wall partitioning coefficient'''

import numpy as np
from part_prop import part_prop
import ipdb

def kimt_calc(y, mfp, num_sb, num_speci, accom_coeff, y_mw, surfT, R_gas, TEMP, NA, 
				y_dens, N_perbin, DStar_org, radius, Psat, therm_sp,
				H2Oi):

	# ------------------------------------------------------------------------------------
	# inputs:
	
	# y - concentration of components' molecules (molecules/cc (air))
	# N_perbin - number of particles in a size bin (excluding wall)
	# mfp - mean free path of gas molecules (m) (num_speci, 1)
	# DStar_org - gas molecule diffusion coefficient (m2/s) (num_speci, 1)
	# radius - particle radius (m)
	# Psat - vapour pressures (molecules/cc (air))
	# surfT - surface tension (g/s2==mN/m==dyn/cm)
	# therm_sp - thermal speed of components (m/s) (num_speci, 1)
	# nuccor - Fuchs-Sutugin correction for newly nucleating particles
	# H2Oi - water index (integer)
	# ------------------------------------------------------------------------------------
	
	# density (g/cm3) and average molecular weight (g/mol) of particles (excluding wall)
	[tot_rho, ish, avMW] = part_prop(y[num_speci:-(num_speci)], num_speci, num_sb-1, NA, 
										y_mw, y_dens, N_perbin)
	
	# Knudsen number (dimensionless)
	Kn = np.repeat(mfp, num_sb-1, 1)/np.repeat(radius, num_speci, 0)

	# Non-continuum regime correction 
    # calculate a correction factor according to the continuum versus non-continuum 
    # regimes
    # expression taken from Jacobson et al (2000), page 457, or Jacobson (2005), page 530 
    # (eq. 16.19).
    # They reference:
    # Fuchs and Sutugin 1971
    # Pruppacher and Klett 1997
	Inverse_Kn = np.power(Kn, -1.0E0)
	correct_1 = (1.33E0+0.71*Inverse_Kn)/(1.0+Inverse_Kn)
	correct_2 = (4.0E0*(1.0E0-accom_coeff))/(3.0E0*accom_coeff)
	# repeat over particle bins
	correct_2 = np.repeat(correct_2.reshape(num_speci, 1), num_sb-1, 1)
	correct_3 = 1.0E0+(correct_1+correct_2)*Kn
	correction = np.power(correct_3, -1.0E0)

	# kelvin factor for each size bin (excluding wall), eq. 16.33 Jacobson et al. (2005)
	# note that avMW has units g/mol, surfT (g/s2==mN/m==dyn/cm), R_gas is multiplied by 
	# 1e7 for units
	# g cm2/s2.mol.K, TEMP is K, radius is multiplied by 1e2 to give cm and tot_rho is
	# g/cm3
	kelv_fac = np.zeros((num_sb-1))
	kelv_fac[ish] = np.exp((2.0E0*avMW[ish]*surfT)/(R_gas*1.0e7*TEMP*(radius[0, ish]*
					1.0e2)*tot_rho[ish]))

	
	# gas phase diffusion coefficient*Fuch-Sutugin correction (cm2/s)
	# eq. 5 Zaveri et al. (2008), scale by 1e4 to convert from m2/s to cm2/s
	kimt = (DStar_org*1e4)*correction
	# final partitioning coefficient (converting m to cm)
	# eq. 16.2 of Jacobson (2005) and eq. 5 Zaveri et al. (2008)
	# species in rows and size bins in columns (cm3/s)
	kimt = (4.0E0*np.pi*(radius*1.0e2)*N_perbin.reshape(1, -1))*kimt
	
	# zero partitioning for any components with vapour pressures greater than
	# 10^11 molecules/cc (air) - accelerates computation time with negligible change to 
	# output
	highVPi = np.squeeze(Psat>1.0e11)
	highVPi[H2Oi] = 0 # mask water
	kimt[highVPi, :] = 0.0
	
	# zero partitioning coefficient for size bins where no particles - enables significant
	# computation time acceleration and is physically realistic
	ish = N_perbin<=1.0e-10
	kimt[:, ish] = 0.0
	
	return kimt, kelv_fac