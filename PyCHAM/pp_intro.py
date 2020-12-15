'''module to set up particle phase part of box model, calls on init_water_partit to initiate water partitioning with seed particles and wall'''
# using the user-defined or default values, the initial number size distribution is determined here

import numpy as np
import size_distr # custom library - see source code
from init_water_partit import init_water_partit
import scipy.constants as si

def pp_intro(y, num_comp, Pybel_objects, TEMP, H2Oi,
		mfp, accom_coeff, y_mw, surfT, 
		RH, siz_str, num_asb, lowersize, uppersize, pmode, pconc, 
		pconct, nuc_comp, testf, std, mean_rad, therm_sp,
		y_dens, Psat, core_diss, kgwt, space_mode, seedVr, 
		spec_namelist, act_coeff, wall_on, partit_cutoff, Press,
		seedi):
	
	# inputs -----------------------------------
	# TEMP - temperature (K) in chamber at start of experiment
	# y_mw - molecular weight (g/mol) of components (num_comp, 1)
	# siz_str - the size structure
	# num_asb - number of size bins (excluding wall)
	# lowersize - lowest size bin radius bound (um)
	# uppersize - largest size bin radius bound (um)
	# pmode - whether particle number concentrations given as modes or explicitly
	# pconc - starting particle concentration (# particle/cc (air)) - if scalar then
	# gets split between size bins in Size_distributions call, or if an array, elements 
	# are allocated to corresponding size bins
	# pconct - time(s) through experiment that particle number concentrations 
	#	correspond to (s)
	# nuc_comp - name of the nucleating component
	# testf - test flag to say whether in normal mode (0) or test mode for front.py (1)
	#       or test mode for pp_intro.py
	# std - geometric standard deviation of the particle number concentration 
	# 		(dimensionless)
	# mean_rad - either the mean radius (um) of particles in lognormal number-size 
	#			distribution (in which case pconc should be scalar), or mean radius of
	#			particles where just one size bin present (in which case pconc is also
	#			scalar)
	# y_dens - liquid density of components (kg/m3) (num_comp, 1)
	# Psat - saturation vapour pressure of components (molecules/cc (air))
	# core_diss - core dissociation constant
	# kgwt - mass transfer coefficient for vapour-wall partitioning (/s)
	# space_mode - string specifying whether to space size bins logarithmically or 
	# linearly
	# seedVr - volume ratio of component(s) comprising seed particles
	# spec_namelist - names of components noted in chemical scheme file
	# act_coeff - activity coefficient of components
	# wall_on - whether or not to consider wall
	# partit_cutoff - product of vapour pressure and activity coefficient
	#		at which gas-particle partitioning assumed zero (Pa)
	# Press - pressure inside chamber
	# seedi - index of seed components
	# ------------------------------------------
	
	if testf==1: # in test mode
		return(0,0,0,0,0,0,0,0,0,0,0,0) # return dummies
	
	# isolate the starting number size distribution information
	i = (pconct[0, :] == 0) # index of initial information
	
	if (sum(i) == 0): # if no initial information provide fillers
		pconcn = np.zeros((1))
		stdn = [1.e20]
		mean_radn = [-1.e6]
	else: # obtain initial information
		pconcn = pconc[:, i]

		try: # if mean_rad already an array
			if (std[:, i].ndim > 1):
				stdn = std[:, i]
			if (std[:, i].ndim == 1):
				stdn = std[:, i]
			if (std[:, i].ndim == 0):
				stdn = [std[:, i]]
		except:
			stdn = [std[:, i]]
		
		try: # if mean_rad already an array
			if (mean_rad[:, i].ndim > 1):
				mean_radn = mean_rad[:, i]
			if (mean_rad[:, i].ndim == 1):
				mean_radn = mean_rad[:, i]
			if (mean_rad[:, i].ndim == 0):
				mean_radn = [mean_rad[:, i]]
		except:
			mean_radn = [mean_rad[:, i]]
	
	# if mean radius not stated explicitly calculate from size ranges (um)
	if (num_asb > 0):
		if (any(mrn == -1.e6 for mrn in mean_radn)):
			if (lowersize > 0.):
				mean_radn[mean_radn == -1.e6] = [10**((np.log10(lowersize)+np.log10(uppersize))/2.0)]
			if (lowersize == 0.):
				mean_radn[mean_radn == -1.e6] = [10**((np.log10(uppersize))/2.0)]
	
	# index of nucleating component
	if len(nuc_comp)>0:
		nuc_compi = spec_namelist.index(nuc_comp[0])
		nuc_comp = np.empty(1, dtype=int)
		nuc_comp[0] = nuc_compi

	R_gas = si.R # ideal gas constant (kg.m2.s-2.K-1.mol-1)
	NA = si.Avogadro # Avogadro's number (molecules/mol)
	
	
	if num_asb == 0: # create dummy variables if no size bins
		N_perbin = np.zeros((1,1))
		x = np.zeros((1,1))
		Varr = np.zeros((1,1))
		Vbou = np.zeros((1,1))
		rad0 = np.zeros((1,1))
		Vol0 = np.zeros((1,1))
		rbou = np.zeros((1,1))
		upper_bin_rad_amp = 1.0e6
	
	# create a number concentration for a lognormal distribution (particles/cc (air))
	# this is where gas partitioning to wall set up
	if (testf == 2):
		print('calling Size_distributions.lognormal')
		
	# if multiple size bins, this call will assume a lognormal distribution if initial 
	# particle concentration is described by mode, or will assign particles to size bins if
	# initial particle concentration per size bin provided
	if (num_asb > 1):
			
		# set scale and standard deviation input for lognormal probability distribution 
		# function, following guidance here: 
		# http://all-geo.org/volcan01010/2013/09/how-to-use-lognormal-distributions-in-python/
		# if mean_radn and stdn are not already arrays then transform to a list 
		if type(mean_radn) != np.ndarray:
			scale = [np.exp(np.log(mean_radn))]
		else:
			scale = np.exp(np.log(mean_radn))
		
		if type(stdn) != np.ndarray:
			stdn = [np.log(stdn)]
		else:
			stdn = np.log(stdn)
			
		loc = 0. # no shift
		
		[N_perbin, x, rbou, Vbou, Varr, upper_bin_rad_amp] = size_distr.lognormal(num_asb, 
			pmode, pconcn, stdn, lowersize, uppersize, loc, scale, space_mode)
		

		if (testf == 2):
			print('finished with Size_distributions.lognormal')
			
	if (num_asb == 1):

		N_perbin = np.array((sum(pconcn))).reshape(-1, 1) # (# particles/cc (air))
		x = np.zeros(1) # radii at size bin centre
		# mean radius of this one size bin (um)
		if (len(mean_radn[0]) == 1): # if a scalar
			meansize = mean_radn[0][0]
		else : # if an array
			meansize = sum(mean_radn)/len(mean_radn)
			
		x[0] = meansize

		# extend uppersize to reduce chance of particles growing beyond this
		upper_bin_rad_amp = 1.0e6
		uppersize = uppersize*upper_bin_rad_amp
		# volume bounds of size bin (um3)
		Vbou = np.array(((lowersize**3.0)*(4.0/3.0)*np.pi, 
						(uppersize**3.0)*(4.0/3.0)*np.pi))
		# volume of single particle (um3)
		Varr = np.zeros((1, 1))
		Varr[0] = (4./3.)*np.pi*(x[0]**3.0)
		# radius bounds of size bin (um)
		rbou = ((Vbou*3.0)/(4.0*np.pi))**(1.0/3.0)

	# set first volume and radius bound to zero, thereby allowing shrinkage to zero in the 
	# smallest bin
	# remember initial first radius bound for saving
	rbou00 = rbou[0]
	Vbou[0] = 0.
	rbou[0] = 0. # this reversed in saving.py back to rbou00
	
	if (num_asb > 0):
		# remember the radii (um) and volumes (um3) at size bin centre before water 
		# partitioning
		rad0 = np.zeros((len(x)))
		rad0[:] = x[:]
		Vol0 = np.zeros((len(Varr)))
		Vol0[:] = Varr[:]
	
	if (wall_on > 0):
		num_asb += 1 # add one to size bin number to account for wall

	# append particle and wall concentrations of components to y (molecules/cc (air))
	y = np.append(y, np.zeros((num_asb*num_comp)))
	
	# molar volume (multiply y_dens by 1e-3 to convert from kg/m3 to g/cc and give
	# MV in units cc/mol)
	MV = (y_mw/(y_dens*1.0e-3)).reshape(num_comp, 1)
	if (sum(pconcn) > 0.0): # account for concentration of components comprising seed
		for ci in range(len(seedi)): # loop through indices of seed components
			# concentration in all size bins (molecules/cc (air)):
			y[num_comp+seedi[ci]:(num_comp*(num_asb)+seedi[ci]):num_comp] = (NA/MV[seedi[ci]])*(Varr*1.e-12*(seedVr[ci]/sum(seedVr)))*N_perbin[:, 0]
	
	# print mass concentration of particles (scale y_dens by 1e-3 to convert from kg/m3
	# to g/cm3)
	if (num_asb > 0): # with particles
		mass_conc = 0. # start cumulation
		for i in range(num_asb-1): # as size bin now account for wall too
			mass_conc += sum((y_dens[:, 0]*1.0e-3)*((y[num_comp*(i+1):num_comp*(i+2)]/si.N_A)*MV[:,0]))
			mass_conc -= (y_dens[int(H2Oi), 0]*1.0e-3)*((y[num_comp*(i+1)+int(H2Oi)]/si.N_A)*MV[int(H2Oi), 0])
		mass_conc = mass_conc*1.e12 # convert from g/cc (air) to ug/m3 (air)
		
		print(str('Total dry (no water) mass concentration of particles at start of simulation is ' + str(mass_conc) + ' ug/m3 (air)'))
	else: # no particles
		print('No particle size bins detected, simulation will not include particles')

	# start counter on number concentration of newly nucleated particles (#/cc (air))
	np_sum = 0.
	
	return(y, N_perbin, x, Varr, Vbou, rad0, Vol0, rbou, MV, num_asb, nuc_comp, rbou00, 
			upper_bin_rad_amp, np_sum)
