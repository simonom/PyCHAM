'''module to estimate component volatilities and liquid densities'''
# This module is responsible for
# setting key properties of components, including liquid-phase saturation vapour pressures
# and liquid-phase densities.  It does this using either UManSysProp (default), or with
# user settings

import numpy as np
import sys
import os
from git import Repo
import shutil
import scipy.constants as si
import errno
import stat

def prop_calc(comp_list, Pybel_objects, TEMP, H2Oi, num_comp, Psat_water, vol_Comp, 
				volP, testf, corei, pconc, umansysprop_update, core_dens, spec_namelist,
				ode_gen_flag, nuci, nuc_comp, num_asb, dens_comp, dens, seed_name):

	# inputs: ------------------------------------------------------------
	# comp_list - array of SMILE strings for components 
	# (omitting water and core, if present)
	# Pybel_objects - list of Pybel objects representing the species in comp_list
	# (omitting water and core, if present)
	# TEMP - temperature (K) in chamber at time function called
	# vol_Comp - names of components (corresponding to those in chemical scheme file)
	# 			that have vapour pressures manually set in volP
	# testf - flag for whether in normal mode (0) or testing mode (1)
	# corei - index of seed particle component
	# pconc - initial number concentration of particles (#/cc (air))
	# umansysprop_update - marker for cloning UManSysProp so that latest version used
	# core_dens - density of core material (g/cc (liquid/solid density))
	# spec_namelist - list of components' names in chemical equation file
	# ode_gen_flag - whether or not called from middle or ode_gen
	# nuci - index of nucleating component
	# nuc_comp - name of nucleating component
	# num_asb - number of actual size bins (excluding wall)
	# dens_comp - chemical scheme names of components with manually assigned 
	# 	densities
	# dens - manually assigned densities (g/cc)
	# seed_name - chemical scheme name(s) of component(s) comprising seed particles
	# ------------------------------------------------------------
	
	
	
	if testf==1:
		return(0,0,0) # return dummies
		
	cwd = os.getcwd() # address of current working directory
	if umansysprop_update == 1:
		print('Cloning latest version of UManSysProp via prop_calc module')
		# download latest version of umansysprop
		
		# check if there is an existing umansysprop folder
		if os.path.isdir(cwd + '/umansysprop'): 
			def handleRemoveReadonly(func, path, exc):
				excvalue = exc[1]
				if not os.access(path, os.W_OK):
					# Is the error an access error ?
					os.chmod(path, stat.S_IWUSR)
					func(path)
				else:
					raise
			# remove existing folder, onerror will change permission of directory if 
			# needed
			shutil.rmtree(cwd + '/umansysprop', ignore_errors=False, onerror=handleRemoveReadonly)
		
		git_url = 'https://github.com/loftytopping/UManSysProp_public.git'
		Repo.clone_from(git_url, (cwd + '/umansysprop'))
	
	# point to umansysprop folder
	sys.path.insert(1, (cwd + '/umansysprop')) # address for updated version
	
	from umansysprop import boiling_points
	from umansysprop import vapour_pressures
	from umansysprop import liquid_densities

	NA = si.Avogadro # Avogadro's number (molecules/mol)
	y_dens = np.zeros((num_comp, 1)) # components' liquid density (kg/m3)
	# vapour pressures of components, ensures any seed component called 
	# core has zero vapour pressure
	Psat = np.zeros((1, num_comp))

	
	if (ode_gen_flag == 0): # estimate densities if called from middle
		
		for i in range (num_comp): # loop through components
			
			# density estimation ---------------------------------------------------------
			if (i == H2Oi): # liquid-phase density of water
				y_dens[i] = 1.*1.E3 # (kg/m3 (particle))
				continue
			# core component properties
			if (i == corei[0]): # density of core
				y_dens[i] = core_dens*1.e3 # core density (kg/m3 (particle))
				continue
			if comp_list[i] == '[HH]': # omit H2 as unliked by liquid density code
				# liquid density code does not like H2, so manually input kg/m3
				y_dens[i] = 1.0e3
			else:
				# density (convert from g/cc to kg/m3)
				y_dens[i] = liquid_densities.girolami(Pybel_objects[i])*1.0E3
			# ----------------------------------------------------------------------------
		
	# account for any manually assigned component densities (kg/m3)
	if (len(dens_comp) > 0  and ode_gen_flag == 0):
		for i in range (len(dens_comp)):
			# index of component in list of components
			dens_indx = spec_namelist.index(dens_comp[i])
			y_dens[dens_indx] = dens[i]
	
	# estimate vapour pressures (log10(atm))
	for i in range (num_comp):
		
		if (i == corei[0]): # if this component is 'core'
			# core component not included in Pybel_objects, continuing
			# here means its vapour pressure is 0 Pa, which is fine
			continue
		
		# water vapour pressure already given by Psat_water (log10(atm))
		if i == H2Oi:
			Psat[0, i] = Psat_water
			continue # water not included in Pybel_objects
		
		# vapour pressure (log10 atm) (# eq. 6 of Nannoolal et al. (2008), with dB of 
		# that equation given by eq. 7 of same reference)
		Psat[0, i] = ((vapour_pressures.nannoolal(Pybel_objects[i], TEMP, 
						boiling_points.nannoolal(Pybel_objects[i]))))
	
	ish = (Psat == 0.)
	Psat = (np.power(10.0, Psat)*101325.0) # convert to Pa from atm
	# retain low volatility where wanted following unit conversion
	Psat[ish] = 0.
	
	# manually assigned vapour pressures (Pa)
	if (len(vol_Comp) > 0 and ode_gen_flag == 0):
		for i in range (len(vol_Comp)):
			# index of component in list of components
			vol_indx = spec_namelist.index(vol_Comp[i])
			Psat[0, vol_indx] = volP[i]

	# ensure if nucleating component is core that it is involatile
	if (nuc_comp == 'core'):
		Psat[0, nuci] = 0.
	
	Psat_Pa = np.zeros((1, num_comp)) # for storing vapour pressures in Pa (Pa)
	Psat_Pa[0, :] = Psat[0, :]
    
	# convert saturation vapour pressures from Pa to molecules/cc (air) using ideal
	# gas law, R has units cc.Pa/K.mol
	Psat = Psat*(NA/((si.R*1.e6)*TEMP))
	# now, in preparation for ode solver, repeat over number of size bins
	if (num_asb > 0):
		Psat = np.repeat(Psat, num_asb, axis=0)
	
	return(Psat, y_dens, Psat_Pa)
