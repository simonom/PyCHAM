'''module to set up particle phase part of box model, calls on init_water_partit to initiate water partitioning with seed particles and wall'''

import numpy as np
from init_water_partit import init_water_partit
import scipy.constants as si
from scipy import stats # import the scipy.stats module

def pp_dursim(y, N_perbin, mean_rad, pconc, corei, lowersize, uppersize, num_speci, 
				num_sb, MV, rad0, std, y_dens, H2Oi):
	
			
	# inputs -----------------------------------
	# y - concentrations of components in particle phase (molecules/cc (air))
	# N_perbin - number concentration of particles (# /cc (air))
	# mean_rad - mean radius of seed particles at this time (um)
	# pconc - number concentration of seed particles (#/cc (air))
	# corei - index of seed material
	# lowersize - smallest radius bound (um)
	# uppersize - greatest radius bound (um)
	# num_speci - number of components
	# num_sb - number of size bins
	# MV - molar volume (cc/mol)
	# rad0 - original radius at particle centres (um)
	# std - standard deviation for lognormal size distribution calculation (dimensionless)
	# y_dens  - density of components (kg/m3)
	# H2Oi - index of water
	# ------------------------------------------
	
	
	# if mean radius not stated explicitly calculate from size ranges (um)
	if mean_rad == -1.0e6 and num_sb>0:
		if lowersize>0.0:
			mean_rad = 10**((np.log10(lowersize)+np.log10(uppersize))/2.0)
		if lowersize == 0.0:
			mean_rad = 10**((np.log10(uppersize))/2.0)
	
	R_gas = si.R # ideal gas constant (kg.m2.s-2.K-1.mol-1)
	NA = si.Avogadro # Avogadro's number (molecules/mol)
	print(mean_rad, rad0)
	if num_sb == 2:
		N_perbin += np.array((pconc)) # (# particles/cc (air))
	if num_sb > 2 and len(pconc)==num_sb-1:
		N_perbin += np.array((pconc)) # (# particles/cc (air))
	if num_sb > 2 and len(pconc)==1:
		
		# set scale and standard deviation input for lognormal probability distribution 
		# function, following guidance here: 
		# http://all-geo.org/volcan01010/2013/09/how-to-use-lognormal-distributions-in-python/
		scale = np.exp(np.log(mean_rad))
		std = np.log(std)
		loc = 0.0 # no shift
	
		# number fraction-size distribution - enforce high resolution to ensure size
		# distribution of seed particles fully captured
		if lowersize>0.0:
			hires = 10**(np.linspace(np.log10(lowersize), np.log10(uppersize), (num_sb-1)*1.0e2))
		else:
			hires = 10**(np.linspace(-4, np.log10(uppersize), (num_sb-1)*1.0e2))
		pdf_output = stats.lognorm.pdf(hires, std, loc, scale)
		pdf_out = np.interp(rad0, hires, pdf_output)	
		# number concentration of seed in all size bins (# particle/cc (air))
		pconc_new = (pdf_out/sum(pdf_out))*pconc
		mean_rad = rad0
		N_perbin += pconc_new # (# particles/cc (air))
	
	# molecular concentration of seed required to comprise these additional seed particle
	# (molecules/cc (air))
	print(pconc)
	print(pconc_new)
	print((pdf_out/sum(pdf_out)))
	
	# volume concentration of seed particles (cm3/cc (air)), note mean_rad scaled up to 
	# cm from um
	Vperbin = ((pconc_new*(4.0/3.0)*np.pi*(mean_rad*1.0e-4)**3.0))
	# corresponding molecular concentration of seed material (molecules/cc (air))
	y[corei:(num_speci*(num_sb-1)+corei):num_speci] += Vperbin/(MV[corei]/NA)
	
	# loop through size bins to estimate new total volume concentrations (um3/cc (air))
	Vtot = np.zeros((num_sb-1))
	Varr = np.zeros((num_sb-1)) # volume concentration of single particles (um3/cc (air))
	mass_conc = 0.0 # mass concentration of particles (g/cm3)
	for i in range(num_sb-1):
		Vtot[i] = np.sum(y[num_speci*i:num_speci*(i+1)]/NA*MV*1.0E12)
		# new volume concentration of single particles (um3/cc (air))
		if N_perbin[i]>0.0:
			Varr[i] = Vtot[i]/N_perbin[i]
		else:
			Varr[i] = (4.0/3.0*np.pi)*rad0[i]**3.0
		# multiply y_dens by 1e-3 to get g/cm3 from kg/m3
		mass_conc += np.sum((y_dens[corei, 0]*1.0e-3)*Vperbin[i])
		
	mass_conc = mass_conc*1.0e12 # convert from g/cc (air) to ug/m3 (air)
	if mass_conc < 1.0e-10:
		mass_conc = 0.0
	print(str('Total dry (no water) mass concentration of particles after injection of seed particles is ' + str(mass_conc) + ' ug/m3 (air)'))

	# new radius of single particles (um)
	x = ((3.0/(4.0*np.pi))*Varr)**(1.0/3.0)
	
	return(y, N_perbin, x, Varr)