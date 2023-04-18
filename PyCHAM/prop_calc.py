##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2023 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
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

def prop_calc(rel_SMILES, Pybel_objects, H2Oi, num_comp, Psat_water, vol_Comp, 
			volP, testf, corei, pconc, umansysprop_update, core_dens,
			ode_gen_flag, nuci, nuc_comp, dens_comp, dens, seed_name,
			y_mw, tempt_cnt, self):

	# inputs: ------------------------------------------------------------
	# rel_SMILES - array of SMILE strings for components 
	# (omitting water and core, if present)
	# Pybel_objects - list of Pybel objects representing the components in rel_SMILES
	# (omitting water and core, if present)
	# self.TEMP - temperature (K) in chamber at all times
	# vol_Comp - names of components (corresponding to those in chemical scheme file)
	# 			that have vapour pressures manually set in volP
	# testf - flag for whether in normal mode (0) or testing mode (1)
	# corei - index of seed particle component
	# pconc - initial number concentration of particles (# particles/cm3 (air))
	# umansysprop_update - marker for cloning UManSysProp so that latest version used
	# core_dens - density of core material (g/cm3 (liquid/solid density))
	# self.comp_namelist - list of component names from the chemical equation file
	# ode_gen_flag - whether or not called from middle or ode_gen
	# nuci - index of nucleating component
	# nuc_comp - name of nucleating component
	# self.num_asb - number of actual size bins (excluding wall)
	# dens_comp - chemical scheme names of components with manually assigned 
	# 	densities
	# dens - manually assigned densities (g/cm3)
	# seed_name - chemical scheme name(s) of component(s) comprising seed particles
	# y_mw - molar mass of components (g/mol)
	# self - reference to PyCHAM
	# tempt_cnt - count on temperatures
	# ------------------------------------------------------------
	
	if (testf == 1):
		return(0, 0, 0) # return dummies
	
	# default values for error message and error flag
	err_mess = ''
	erf = 0
		
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
	self.Psat = np.zeros((1, num_comp))
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
			dens_indx = self.comp_namelist.index(dens_comp[i])
			y_dens[dens_indx] = dens[i]

	# for records (e.g. plotting volatility basis set), 
	# estimate and list the pure component saturation vapour 
	# pressures (Pa) at standard temperature (298.15 K), though note that
	# any manually assigned vapour pressures overwrite these later in
	# this module
	self.Psat_Pa_rec = np.zeros((num_comp))
	
	# prepare carbon and oxygen atom count
	self.Cnum = np.zeros((num_comp, 1))
	self.Onum = np.zeros((num_comp, 1))
	
	# prepare for gathering the indices of any HOM-RO2-MCM-RO2 accretion products
	self.aoRO2_indx = []
	
	# estimate vapour pressures (log10(atm)) and O:C ratio
	# note when the O:C ratio and vapour pressure at 298.15 K are
	# combined, one can produce the two-dimensional volatility
	# basis set, as shown in Fig. 1 of https://doi.org/10.5194/acp-20-1183-2020
	for i in range (num_comp):
	
		# note, rec_now_flag is only changed below if alternative vapour pressure estimation 
		# method uncommented for HOMs
		rec_now_flag = 0

		if ('_ao' in self.comp_namelist[i]): # if part of the extension
			# if it's the peroxy radical formed from autoxidation
			if (self.comp_namelist[i][-4::] == '_ao1'):
				self.aoRO2_indx.append(i)
				
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
			self.Psat[0, i] = Psat_water
			if (self.TEMP[tempt_cnt] == 298.15):
				self.Psat_Pa_rec[i] = self.Psat[0, i]
			else:
				[_, self.Psat_Pa_rec[i], _] = water_calc(298.15, 0.5, si.N_A)
			OC[0, i] = 0.
			self.HC[0, i] = 0.
			self.nom_mass[0, i] = 2.*1.+1.*16.
			continue
		
		if (self.comp_namelist[i] == 'O3'):
			# vapour pressure of ozone from https://doi.org/10.1063/1.1700683
			self.Psat[0, i] =  np.log10((8.25313-(814.941587/self.TEMP[tempt_cnt])-0.001966943*self.TEMP[tempt_cnt])*1.31579e-3)
			if (self.TEMP[tempt_cnt] == 298.15):
				self.Psat_Pa_rec[i] = self.Psat[0, i]
			else:
				self.Psat_Pa_rec[i] =  np.log10((8.25313-(814.941587/298.15)-0.001966943*298.15)*1.31579e-3)
			OC[0, i] = 0.
			self.HC[0, i] = 0.
			self.nom_mass[0, i] = 0.*1.+3.*16.
			continue

		# possibly use different method for vapour pressure (log10(atm)) of HOMs
		
		if ('_ao' in self.comp_namelist[i] or 'C10H' in self.comp_namelist[i]  or 'C18H' in self.comp_namelist[i] or 'C19H' in self.comp_namelist[i] or 'C20H' in self.comp_namelist[i]):
			Psatnow = -2.63+-0.50*rel_SMILES[i].count('O') + (rel_SMILES[i].count('C')-5)*-0.80
			rec_now_flag = 1 # record this vapour pressure
				
		# vapour pressure (log10(atm)) (eq. 6 of Nannoolal et al. (2008), with dB of 
		# that equation given by eq. 7 of same reference)
		else: 
			Psatnow = ((vapour_pressures.nannoolal(Pybel_objects[i], self.TEMP[tempt_cnt], 
					boiling_points.nannoolal(Pybel_objects[i]))))

		# in case you want to ensure small molecules don't contribute to particle mass
		#if rel_SMILES[i].count('C')<=5:
		#	 Psatnow += 10 # ensure no condensation of small molecules


		
		try: # in case array
			self.Psat[0, i] = Psatnow[0]
		except: # in case float
			self.Psat[0, i] = Psatnow
	
		if (self.TEMP[tempt_cnt] == 298.15 or rec_now_flag == 1):
			try: # in case array
				self.Psat_Pa_rec[i] = Psatnow[0] # note transfer to Pa is below
			except: # in case float
				self.Psat_Pa_rec[i] = Psatnow # note transfer to Pa is below
		else: 
			Psatnow = ((vapour_pressures.nannoolal(Pybel_objects[i], 298.15, 
					boiling_points.nannoolal(Pybel_objects[i]))))
			
			try: # in case array
				self.Psat_Pa_rec[i]  = Psatnow[0]
			except: # in case float
				self.Psat_Pa_rec[i] = Psatnow
		
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

		# carbon and oxygen numbers in this component
		self.Cnum[i, 0] = rel_SMILES[i].count('C')
		self.Onum[i, 0] = rel_SMILES[i].count('O')

		# O:C ratio determined from SMILES string
		if (rel_SMILES[i].count('C') > 0):
			OC[0, i] = self.Onum[i, 0]/self.Cnum[i, 0] 
			self.HC[0, i] = Hcount/self.Cnum[i, 0] 
		else: # if no carbons in this component
			OC[0, i] = 0.
	
	ish = (self.Psat == 0.) # non-volatiles
	
	self.Psat = (10.**self.Psat)*101325. # convert to Pa from atm
	self.Psat_Pa_rec = (10.**self.Psat_Pa_rec)*101325. # convert to Pa from atm
	# retain low volatility where wanted following unit conversion
	self.Psat[ish] = 0.
	
	# in preparation for ode solver, tile over size and wall bins if present
	if (self.num_asb+self.wall_on > 0):
		self.Psat = np.tile(self.Psat, (self.num_asb+self.wall_on, 1))
	else:
		self.Psat = np.tile(self.Psat, (1, 1))
	
	# list to remember which components have vapour pressures specified
	vi_rec = []
	
	# list to remember which walls affected by wall-specific vapour pressures
	self.P_wfunc_wi = []
	# list to remember which components affected by wall-specific vapour pressures
	self.P_wfunc_ci = []
	# list to remember the user-defined vapour pressure
	self.P_wfunc = []
	
	# manually assigned vapour pressures (Pa)
	if (len(vol_Comp) > 0 and ode_gen_flag == 0):
		for i in range (len(vol_Comp)):

			if '_wall' in vol_Comp[i]: # this is specific to a wall
				# get wall number
				# get location of wall number
				wn = vol_Comp[i].rfind('l')
				wn = int(float(vol_Comp[i][wn+1::]))

				try: # first see if an individual component has been named
					vol_indx = [self.comp_namelist.index(vol_Comp[i][0:-6])]
				except: # could be a group of components
					group_name = vol_Comp[i][0:-6]

					# check if an inequality present
					if '<' or '>' or '==' in group_name:
						# get locations of underscores
						us_indx = group_name.rfind('_') 
						# get inequality
						inequal =  group_name[us_indx+1::]

						# get index of components in this group
						if '<' in inequal:
							if '=' in inequal: # less than or equal to
								# get index of all components in this group
								vol_indx = self.Psat[0, :] <= float(inequal[2::])
							else: # just less than
								# get index of all components in this group
								vol_indx = self.Psat[0, :] < float(inequal[1::])
						if '>' in inequal:
							if '=' in inequal: # greater than or equal to
								# get index of all components in this group
								vol_indx = self.Psat[0, :] >= float(inequal[2::])
							else: # just greater than
								# get index of all components in this group
								vol_indx = self.Psat[0, :] > float(inequal[1::])
						if '==' in inequal:
							# get index of all components in this group
							vol_indx = self.Psat[0, :] == float(inequal[2::])

					if 'RO2' in group_name: # if RO2, as categorised by the chemical scheme
						vol_indx = self.RO2_indices[:, 1]
				
				# assign user-defined vapour pressure for this wall (Pa)
				self.Psat[(self.num_asb-1)+wn, vol_indx] = volP[i]

				# remember which components affected
				self.P_wfunc_ci.append(vol_indx)

				# remember which walls affected
				self.P_wfunc_wi.append([wn*i for i in [1]*len(vol_indx)])
				# remember the user-defined vapour pressure
				self.P_wfunc.append(volP[i])
				
			else: # not specific to a wall

				# index of component in list of components
				try: # first see if an individual component has been named
					vol_indx = self.comp_namelist.index(vol_Comp[i])
				except: # could be a group of components
					group_name = vol_Comp[i][0:-6]
					
					# check if an inequality present
					if '<' or '>' or '==' in group_name:
						# get locations of underscores
						us_indx = group_name.rfind('_') 
						# get inequality
						inequal =  group_name[us_indx+1::]

						# get index of components in this group
						if '<' in inequal:
							if '=' in inequal: # less than or equal to
								# get index of all components in this group
								vol_indx = self.Psat[0, :] <= float(inequal[2::])
							else: # just less than
								# get index of all components in this group
								vol_indx = self.Psat[0, :] < float(inequal[1::])
						if '>' in inequal:
							if '=' in inequal: # greater than or equal to
								# get index of all components in this group
								vol_indx = self.Psat[0, :] >= float(inequal[2::])
							else: # just greater than
								# get index of all components in this group
								vol_indx = self.Psat[0, :] > float(inequal[1::])
						if '==' in inequal:
							# get index of all components in this group
							vol_indx = self.Psat[0, :] == float(inequal[2::])

					if 'RO2' in group_name: # if RO2, as categorised by the chemical scheme
						vol_indx = self.RO2_indices

				# assign user-defined vapour pressure (Pa)
				self.Psat[:, vol_indx] = volP[i]
				self.Psat_Pa_rec[vol_indx] = volP[i]
				vi_rec.append(vol_indx)

	# ensure if nucleating component is core that it is involatile
	if (nuc_comp == 'core'):
		self.Psat[0, nuci] = 0.
	
	self.Psat_Pa = np.zeros((1, num_comp)) # for storing vapour pressures in Pa (Pa)
	self.Psat_Pa[0, :] = self.Psat[0, :]
    
	# convert saturation vapour pressures from Pa to # molecules/cm3 (air) using ideal
	# gas law, R has units cm3.Pa/K.mol
	self.Psat = self.Psat*(NA/((si.R*1.e6)*self.TEMP[tempt_cnt]))

	# remember Psat (# molecules/cm3) in case it is altered 
	# by user-defined inputs in partit_var.py
	self.Psat_num_rec = np.zeros((self.Psat.shape))
	self.Psat_num_rec[:, :] = self.Psat[:, :]
	
	# if vapour pressure plot requested then make this now --------------------------------------------------------------------
	if (self.testf == 3.2): 
		
		import matplotlib.pyplot as plt
		
		plt.ion() # show results to screen and turn on interactive mode
		
		# prepare plot
		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
		
		# get indices of vapour pressure in descending order
		des_ind = np.flip(np.argsort(self.Psat_Pa_rec, axis = 0))
		
		array_names = np.squeeze(np.array(self.comp_namelist))

		# plot vapour pressures against component names in descending order of volatility
		ax0.semilogy(np.arange(len(self.Psat_Pa_rec)), self.Psat_Pa_rec[des_ind], '+')
		
		ax0.set_ylabel(r'Pure Component Saturation Vapour Pressure at 298.15 K (Pa)', fontsize = 14)
		ax0.set_xlabel(r'Component name', fontsize = 14)
		# set location of x ticks
		ax0.set_xticks(np.arange(len(self.comp_namelist)))
		ax0.set_xticklabels(array_names[des_ind], rotation = 45)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		ax0.set_title(str('Pure Component Saturation Vapour Pressures at 298.15 K of All Components'), fontsize = 14)
		
		# tell middle that plot made and need to stop running code now
		err_mess = 'Stop'
		erf = 1
	# end of plotting section ----------------------------------------------------------------------------------------------------------
	
	return(y_dens, OC, self, err_mess, erf)