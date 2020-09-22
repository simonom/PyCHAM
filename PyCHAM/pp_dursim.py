'''module to set up particle phase part of box model, calls on init_water_partit to initiate water partitioning with seed particles and wall'''

import numpy as np
from init_water_partit import init_water_partit
import scipy.constants as si
from scipy import stats # import the scipy.stats module

def pp_dursim(y, N_perbin, mean_rad, pconc, corei, lowersize, uppersize, num_speci, 
				num_sb, MV, rad0, radn, std, y_dens, H2Oi, rbou):
	
			
	# inputs -----------------------------------
	# y - concentrations of components in particle phase (molecules/cc (air))
	# N_perbin - number concentration of particles (#/cc (air))
	# mean_rad - mean radius of seed particles at this time (um)
	# pconc - number concentration of seed particles (#/cc (air))
	# corei - index of seed material
	# lowersize - smallest radius bound (um)
	# uppersize - greatest radius bound (um)
	# num_speci - number of components
	# num_sb - number of size bins (exlcuding wall)
	# MV - molar volume (cc/mol)
	# rad0 - original radius at size bin centres (um)
	# radn - current radius at size bin centres (um)
	# std - standard deviation for lognormal size distribution calculation (dimensionless)
	# y_dens  - density of components (kg/m3)
	# H2Oi - index of water
	# rbou - radius bounds per size bin (um)
	# ------------------------------------------
	
	
	# if mean radius not stated explicitly calculate from size ranges (um)
	if mean_rad == -1.0e6 and num_sb>0:
		if lowersize>0.0:
			mean_rad = 10**((np.log10(lowersize)+np.log10(uppersize))/2.0)
		if lowersize == 0.0:
			mean_rad = 10**((np.log10(uppersize))/2.0)
	
	R_gas = si.R # ideal gas constant (kg.m2.s-2.K-1.mol-1)
	NA = si.Avogadro # Avogadro's number (molecules/mol)
	
	if num_sb == 1:
		N_perbin += np.array((pconc)) # (# particles/cc (air))
		pconc_new = pconc

	# number concentration stated per size bin in multi size bin simulation
	if num_sb > 1 and (len(pconc) == num_sb): 
		N_perbin += np.array((pconc)) # (# particles/cc (air))
		pconc_new = pconc

	# total number concentration stated in multi size bin simulation
	if num_sb > 1 and len(pconc)==1:
		
		# set scale and standard deviation input for lognormal probability distribution 
		# function, following guidance here: 
		# http://all-geo.org/volcan01010/2013/09/how-to-use-lognormal-distributions-in-python/
		scale = np.exp(np.log(mean_rad))
		std = np.log(std)
		loc = 0.0 # no shift
	
		# number fraction-size distribution - enforce high resolution to ensure size
		# distribution of seed particles fully captured
		hires = 10**(np.linspace(np.log10((rad0[0]-(rbou[1]-rbou[0])/2.1)), np.log10(uppersize), (num_sb)*1.0e2))
		pdf_output = stats.lognorm.pdf(hires, std, loc, scale)
		pdf_out = np.interp(radn, hires, pdf_output)	
		# number concentration of seed in all size bins (# particle/cc (air))
		pconc_new = (pdf_out/sum(pdf_out))*pconc
		N_perbin[:, 0] += pconc_new # (# particles/cc (air))
	
	# molecular concentration of seed required to comprise these additional seed particle
	# (molecules/cc (air)):
	
	# volume concentration of seed particles (cm3/cc (air)), note radn scaled up to 
	# cm from um
	Vperbin = ((pconc_new*(4.0/3.0)*np.pi*(radn*1.0e-4)**3.0))
	# corresponding molecular concentration of seed material (molecules/cc (air))
	y[corei:(num_speci*(num_sb)+corei):num_speci] += Vperbin/(MV[corei]/NA)
	
	# loop through size bins to estimate new total volume concentrations (um3/cc (air))
	Vtot = np.zeros((num_sb))
	Varr = np.zeros((num_sb)) # volume concentration of single particles (um3/cc (air))
	mass_conc = 0.0 # mass concentration of particles (g/cm3)
	for i in range(num_sb):
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
	print(str('Total dry (no water) mass concentration of newly injected seed particles: ' + str(mass_conc) + ' ug/m3 (air)'))

	# new radius of single particles (um)
	radn = ((3.0/(4.0*np.pi))*Varr)**(1.0/3.0)
	
	return(y, N_perbin, radn, Varr)
