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
from water_calc import water_calc

def prop_calc(rel_SMILES, Pybel_objects, TEMP, H2Oi, num_comp, Psat_water, vol_Comp, 
			volP, testf, corei, pconc, umansysprop_update, core_dens, spec_namelist,
			ode_gen_flag, nuci, nuc_comp, num_asb, dens_comp, dens, seed_name,
			y_mw, self):

	# inputs: ------------------------------------------------------------
	# rel_SMILES - array of SMILE strings for components 
	# (omitting water and core, if present)
	# Pybel_objects - list of Pybel objects representing the components in rel_SMILES
	# (omitting water and core, if present)
	# TEMP - temperature (K) in chamber at time function called
	# vol_Comp - names of components (corresponding to those in chemical scheme file)
	# 			that have vapour pressures manually set in volP
	# testf - flag for whether in normal mode (0) or testing mode (1)
	# corei - index of seed particle component
	# pconc - initial number concentration of particles (# particles/cm3 (air))
	# umansysprop_update - marker for cloning UManSysProp so that latest version used
	# core_dens - density of core material (g/cm3 (liquid/solid density))
	# spec_namelist - list of component names in chemical equation file
	# ode_gen_flag - whether or not called from middle or ode_gen
	# nuci - index of nucleating component
	# nuc_comp - name of nucleating component
	# num_asb - number of actual size bins (excluding wall)
	# dens_comp - chemical scheme names of components with manually assigned 
	# 	densities
	# dens - manually assigned densities (g/cm3)
	# seed_name - chemical scheme name(s) of component(s) comprising seed particles
	# y_mw - molar mass of components (g/mol)
	# self - reference to PyCHAM
	# ------------------------------------------------------------
	
	
	if (testf == 1):
		return(0, 0, 0) # return dummies
		
	cwd = os.getcwd() # address of current working directory
	
	# if update required, note this update flag is set when model variables are checked
	if (umansysprop_update == 1):
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
	# oxygen:carbon ratio of components
	OC = np.zeros((1, num_comp))
	# hydrogen:carbon ratio of components
	self.HC = np.zeros((1, num_comp))
	# nominal molar mass of components
	self.nom_mass = np.zeros((1, num_comp))

	
	if (ode_gen_flag == 0): # estimate densities if called from middle
		
		for i in range (num_comp): # loop through components
			
			# density estimation ---------------------------------------------------------
			if (i == H2Oi): # liquid-phase density of water
				y_dens[i] = 1.*1.e3 # (kg/m3 (particle))
				continue
			# core component properties
			if (i == corei[0]): # density of core
				y_dens[i] = core_dens*1.e3 # core density (kg/m3 (particle))
				continue
			if rel_SMILES[i] == '[HH]': # omit H2 as unliked by liquid density code
				# liquid density code does not like H2, so manually input kg/m3
				y_dens[i] = 1.e3
			else:
				# density (convert from g/cm3 to kg/m3)
				y_dens[i] = liquid_densities.girolami(Pybel_objects[i])*1.e3
			# ----------------------------------------------------------------------------
		
	# account for any manually assigned component densities (kg/m3)
	if (len(dens_comp) > 0  and ode_gen_flag == 0):
		for i in range (len(dens_comp)):
			# index of component in list of components
			dens_indx = spec_namelist.index(dens_comp[i])
			y_dens[dens_indx] = dens[i]

	# for records (e.g. plotting volatility basis set), 
	# estimate and list the pure component saturation vapour 
	# pressures (Pa) at standard temperature (298.15 K), though note that
	# any manually assigned vapour pressures overwrite these later in
	# this module
	Psat_Pa_rec = np.zeros((num_comp))
	
	# estimate vapour pressures (log10(atm)) and O:C ratio
	# note when the O:C ratio and vapour pressure at 298.15 K are
	# combined, one can produce the two-dimensional volatility
	# basis set, as shown in Fig. 1 of https://doi.org/10.5194/acp-20-1183-2020
	for i in range (num_comp):
		
		if (i == corei[0]): # if this component is 'core'
			# core component not included in Pybel_objects
			# assign an assumed O:C ratio of 0.
			OC[0, i] = 0.
			self.HC[0, i] = 0.
			self.nom_mass[0, i] = 132.
			# continuing
			# here means its vapour pressure is 0 Pa, which is fine; if a 
			# different vapour pressure is specified it is accounted for below
			continue
		
		# water vapour pressure already given by Psat_water (log10(atm))
		# and water not included in Pybel_objects
		if (i == H2Oi):
			Psat[0, i] = Psat_water
			if (TEMP == 298.15):
				Psat_Pa_rec[i] = Psat[0, i]
			else:
				[_, Psat_Pa_rec[i], _] = water_calc(298.15, 0.5, si.N_A)
			OC[0, i] = 0.
			self.HC[0, i] = 0.
			self.nom_mass[0, i] = 2.*1.+1.*16.
			continue
		
		if (spec_namelist[i] == 'O3'):
			# vapour pressure of ozone from https://doi.org/10.1063/1.1700683
			Psat[0, i] =  np.log10((8.25313-(814.941587/TEMP)-0.001966943*TEMP)*1.31579e-3)
			if (TEMP == 298.15):
				Psat_Pa_rec[i] = Psat[0, i]
			else:
				Psat_Pa_rec[i] =  np.log10((8.25313-(814.941587/298.15)-0.001966943*298.15)*1.31579e-3)
			OC[0, i] = 0.
			self.HC[0, i] = 0.
			self.nom_mass[0, i] = 0.*1.+3.*16.
			continue

		# possibly use different method for vapour pressure (log10(atm)) of HOMs
		#if 'api_' in spec_namelist[i] or 'API_' in spec_namelist[i]:
			#Psat[0, i] = -0.12*rel_SMILES[i].count('O') + rel_SMILES[i].count('C')*-0.22-2.5
			
			#	Psat[0, i] = ((vapour_pressures.myrdal_and_yalkowsky(Pybel_objects[i], TEMP, boiling_points.nannoolal(Pybel_objects[i]))))
			
			#if (TEMP == 298.15):
				#Psat_Pa_rec[i] = Psat_Pa[0, i]
			#else:
				#Psat_Pa_rec[i] = -0.1*rel_SMILES[i].count('O') + rel_SMILES[i].count('C')*-0.2+0.5
				#Psat_Pa_rec[i] = ((vapour_pressures.myrdal_and_yalkowsky(Pybel_objects[i], 298.15, boiling_points.nannoolal(Pybel_objects[i]))))
		
		else: # for non-HOM components
			# vapour pressure (log10(atm)) (eq. 6 of Nannoolal et al. (2008), with dB of 
			# that equation given by eq. 7 of same reference)
			Psatnow = ((vapour_pressures.nannoolal(Pybel_objects[i], TEMP, 
						boiling_points.nannoolal(Pybel_objects[i]))))

			# in case you want to ensure small molecules don't contribute to particle mass
			#if rel_SMILES[i].count('C')<=5:
			#	 Psatnow += 10 # ensure no condensation of small molecules
			
			try: # in case array
				Psat[0, i] = Psatnow[0]
			except: # in case float
				Psat[0, i] = Psatnow
			
			if (TEMP == 298.15):
				try: # in case array
					Psat_Pa_rec[i] = Psatnow[0] # note transfer to Pa is below
				except: # in case float
					Psat[0, i] = Psatnow # note transfer to Pa is below
			else: 
				Psatnow= ((vapour_pressures.nannoolal(Pybel_objects[i], 298.15, 
						boiling_points.nannoolal(Pybel_objects[i]))))
			
				try: # in case array
					Psat_Pa_rec[i]  = Psatnow[0]
				except: # in case float
					Psat_Pa_rec[i] = Psatnow

		# if component is chlorine, then H:C is 0 and can continue
		if (rel_SMILES[i] == 'ClCl'):
			self.HC[0, i] = 0.
			self.nom_mass[0, i] = 70.
			continue

		# if hydrogen is present in this molecule
		if ('H' in Pybel_objects[i].formula):

			

			Hindx_start = Pybel_objects[i].formula.index('H')+1
			Hindx_end = Hindx_start
			for Hnum_test in Pybel_objects[i].formula[Hindx_start::]:
				try:
					float(Hnum_test) # only continue if this character is a number
					Hindx_end += 1
					if (Hindx_end == len(Pybel_objects[i].formula)):
						Hcount = float(Pybel_objects[i].formula[Hindx_start:Hindx_end])
				except:
					if (Hindx_end != Hindx_start):
						Hcount = float(Pybel_objects[i].formula[Hindx_start:Hindx_end])
					else:
						Hcount = 1. # if no number then Hydrogen must be alone
					break

		else: # if no hydrocarbons
			Hcount = 0.
			self.HC[0, i] = 0.

		self.nom_mass[0, i] = Hcount*1.+rel_SMILES[i].count('O')*16.+rel_SMILES[i].count('C')*12.+rel_SMILES[i].count('N')*14.+rel_SMILES[i].count('S')*32.
		
		#if (rel_SMILES[i] == '[N+](=O)(O)[O-]'):
		#	print(Hcount, self.nom_mass[0, i])
		#	import ipdb; ipdb.set_trace()

		# O:C ratio determined from SMILES string
		if (rel_SMILES[i].count('C') > 0):
			OC[0, i] = (rel_SMILES[i].count('O'))/(rel_SMILES[i].count('C'))
			self.HC[0, i] = Hcount/(rel_SMILES[i].count('C'))
		else: # if no carbons in this component
			OC[0, i] = 0.
	
	ish = (Psat == 0.)
	
	Psat = (10.**Psat)*101325. # convert to Pa from atm
	Psat_Pa_rec = (10.**Psat_Pa_rec)*101325 # convert to Pa from atm
	# retain low volatility where wanted following unit conversion
	Psat[ish] = 0.
	
	# list to remember which components have vapour pressures specified
	vi_rec = []
	
	# manually assigned vapour pressures (Pa)
	if (len(vol_Comp) > 0 and ode_gen_flag == 0):
		for i in range (len(vol_Comp)):
			# index of component in list of components
			vol_indx = spec_namelist.index(vol_Comp[i])
			Psat[0, vol_indx] = volP[i]
			Psat_Pa_rec[vol_indx] = volP[i]
			vi_rec.append(vol_indx)

	# ensure if nucleating component is core that it is involatile
	if (nuc_comp == 'core'):
		Psat[0, nuci] = 0.
	
	Psat_Pa = np.zeros((1, num_comp)) # for storing vapour pressures in Pa (Pa)
	Psat_Pa[0, :] = Psat[0, :]
    
	# convert saturation vapour pressures from Pa to # molecules/cm3 (air) using ideal
	# gas law, R has units cc.Pa/K.mol
	Psat = Psat*(NA/((si.R*1.e6)*TEMP))
	# now, in preparation for ode solver, repeat over number of size bins
	if (num_asb > 0):
		Psat = np.repeat(Psat, num_asb, axis=0)
	
	return(Psat, y_dens, Psat_Pa, Psat_Pa_rec, OC, self)