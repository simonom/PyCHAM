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

# called/returned from/to the front.py and ode_gen.py modules, 
# this module is responsible for
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

def volat_calc(comp_list, Pybel_objects, TEMP, H2Oi, num_comp, Psat_water, vol_Comp, 
				volP, testf, corei, seed_name, pconc, umansysprop_update, core_dens, comp_namelist,
				ode_gen_flag, nuci, nuc_comp, self):

	# inputs: ------------------------------------------------------------
	# comp_list - array of SMILE strings for components 
	# (omitting water and core, if present)
	# Pybel_objects - list of Pybel objects representing the species in comp_list
	# (omitting water and core, if present)
	# TEMP - temperature (K) in chamber at time function called
	# Psat_water - pure component saturation vapour pressure of water (log10(atm)) 
	# vol_Comp - names of components (corresponding to those in chemical scheme file)
	# 			that have vapour pressures manually set in volP
	# testf - flag for whether in normal mode (0) or testing mode (1)
	# corei - index of seed particle component
	# seed_name - name(s) of components(s) comprising seed particles
	# pconc - initial number concentration of particles (#/cc (air))
	# umansysprop_update - marker for cloning UManSysProp so that latest version used
	# core_dens - density of core material (g/cc (liquid/solid density))
	# comp_namelist - list of components' names in chemical equation file
	# ode_gen_flag - whether or not called from front or ode_gen
	# nuci - index of nucleating component
	# nuc_comp - name of nucleating component
	# self - reference to PyCHAM
	# ------------------------------------------------------------
	
	
	
	if (testf == 1):
		return(0,0,0) # return dummies
		
	cwd = os.getcwd() # address of current working directory
	if umansysprop_update == 1:
		print('Cloning latest version of UManSysProp in volat_calc module')
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

	NA = si.Avogadro # Avogadro's number (# molecules/mol)
	y_dens = np.zeros((num_comp, 1)) # components' liquid density (kg/m3)
	
	if ode_gen_flag == 0: # estimate densities
		
		for i in range (num_comp):
			
			# density estimation ---------------------------------------------------------
			if i == H2Oi:
				y_dens[i] = 1.0*1.e3 # (kg/m3 (particle))
				continue
			# core properties
			if (i == corei[0]):
				y_dens[i] = core_dens*1.e3 # core density (kg/m3 (particle))
				continue
			# nucleating component density, if component is core (kg/m3 (particle))
			if i == nuci and nuc_comp[0] == 'core': 
				y_dens[i] = 1.*1.e3
				continue
			if comp_list[i] == '[HH]': # omit H2 as unliked by liquid density code
				# liquid density code does not like H2, so manually input kg/m3
				y_dens[i] = 1.e3
			else:
				# density (convert from g/cc to kg/m3)
				y_dens[i] = liquid_densities.girolami(Pybel_objects[i])*1.e3
			# ----------------------------------------------------------------------------
	
	# note that self.Psat already tiled over size bins and wall bins

	# estimate vapour pressures (log10(atm))
	for i in range (num_comp):
		
		if (i == corei[0]): # if this core component
			continue # core component not included in Pybel_objects
		if (i == nuci and nuc_comp[0] == 'core'):
			continue # core component not included in Pybel_objects
		
		# water vapour pressure already given by Psat_water (log10(atm))
		if (i == H2Oi):
			self.Psat[:, i] = Psat_water
			continue # water not included in Pybel_objects
		
		# vapour pressure (log10 atm) (# eq. 6 of Nannoolal et al. (2008), with dB of 
		# that equation given by eq. 7 of same reference)
		self.Psat[:, i] = ((vapour_pressures.nannoolal(Pybel_objects[i], TEMP, 
						boiling_points.nannoolal(Pybel_objects[i]))))
	
	ish = (self.Psat == 0.) # non-volatiles
	
	self.Psat = (np.power(10.0, self.Psat)*101325.) # convert to Pa from atm

	# retain low volatility where wanted
	self.Psat[ish] = 0.
	
	# note that self.Psat alreday tiled over size bins and wall bins

	# list to remember which components have vapour pressures specified
	vi_rec = []
	
	# list to remember which walls affected by wall-specific vapour pressures
	self.P_wfunc_wi = []
	# list to remember which components affected by wall-specific vapour pressures
	self.P_wfunc_ci = []
	# list to remember the user-defined vapour pressure
	self.P_wfunc = []
	
	# manually assigned vapour pressures (Pa)
	if (len(vol_Comp) > 0):
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
	
	# --------------------------------

	# ensure if nucleating component is core that it is involatile
	if (nuc_comp[0] == 'core'):
		self.Psat[:, nuci] = 0.
	    
	# convert saturation vapour pressures from Pa to # molecules/cm3 (air) using ideal
	# gas law, R has units cm3.Pa/K.mol
	self.Psat = self.Psat*(NA/((si.R*1.e6)*TEMP))
	
	# remember Psat (# molecules/cm3) in case it is altered 
	# by user-defined inputs in partit_var.py
	self.Psat_num_rec[:, :] = self.Psat[:, :]

	return(self, y_dens)