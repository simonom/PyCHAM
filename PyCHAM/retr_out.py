########################################################################
#                                                                      #
# Copyright (C) 2018-2026                                              #
# Simon O'Meara : simon.omeara@manchester.ac.uk                        #
#                                                                      #
# All Rights Reserved.                                                 #
# This file is part of PyCHAM                                          #
#                                                                      #
# PyCHAM is free software: you can redistribute it and/or modify it    #
# under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or    #
# (at  your option) any later version.                                 #
#                                                                      #
# PyCHAM is distributed in the hope that it will be useful, but        #
# WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# General Public License for more details.                             #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with PyCHAM.  If not, see <http://www.gnu.org/licenses/>.      #
#                                                                      #
########################################################################
'''code to open saved files and return useful variables'''
# e.g. called by PyCHAM plotting codes to obtain the required model 
# outputs for evaluation

import numpy as np
import os
import ast
import pickle
import scipy.constants as si
import matplotlib.pyplot as plt

# define function, note that retr_out is called by the click202 function in gui.py 
# and it can return progress updates on loading results
def retr_out(self):
	
	# inputs: -------------------------------
	# self.dir_path - path of directory requested by the calling
	# code to be looked at
	# ---------------------------------------
	
	# name of file where experiment constants saved
	fname = str(self.dir_path + '/model_and_component_constants')
	
	try: # try opening file
		const_in = open(fname)
	except: # if still can't find a valid file return error message
		err_mess = str('Error - no such file ' + fname + 
			', please check it still exists')
		self.l203a.setText(err_mess)
		# set border around error message
		if (self.bd_pl == 1):
			self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
			self.bd_pl = 2
		else:
			self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
			self.bd_pl = 1
		return()

	# output format changed in version 4 of PyCHAM to accelerate 
	# loading of results, therefore check which version saved 
	# the results in question and treat accordingly
	try:
		# if v4 or later this will work
		load_path = str(self.dir_path + '/nom_mass.npy') # path
		nom_mass = np.load(load_path, allow_pickle=True)
		v4_flag = 1
	except:
		# if earlier than v4 then flag
		v4_flag = 0
	
	# prepare to hold constants
	const = {} # create empty dictionary to hold constants
	# empty dictionary to contain indices of certain groups of components
	group_indx = {}
	group_indx['RO2i'] = [] # filler in case of no RO2i
	group_indx['RO2pooli'] = [] # filler in case of no RO2pooli
	group_indx['ROi'] = [] # filler in case of no ROi
	group_indx['HOMRO2'] = [] # filler in case of no HOMRO2
	group_indx['HOMs'] = [] # filler	
	group_indx['OOH'] = [] # filler
	group_indx['HOM_OOH'] = [] # filler
	group_indx['OH'] = [] # filler
	group_indx['HOM_OH'] = [] # filler	
	group_indx['carbonyl'] = [] # filler	
	group_indx['HOM_carbonyl'] = [] # filler
	group_indx['NO3'] = [] # filler
	group_indx['HOM_NO3'] = [] # filler
	group_indx['PRAMpr_indx'] = [] # filler
	group_indx['PRAMcsmon_indx'] = [] # filler
	group_indx['PRAMcsacc_indx'] = [] # filler	
	
	for line in const_in.readlines():
		
		dlist = [] # empty list to hold values

		# get chemical scheme names
		if (str(line.split(',')[0]) == 'chem_scheme_names'):
			# find index of first [ and index of last ]
			icnt = 0 # count on characters
			for i in line:
				if i == '[':
					st_indx = icnt
					break
				icnt += 1 # count on characters

			for cnt in range(10):
				if line[-cnt] == ']':
					fi_indx = -cnt+1
					break

			comp_names = ast.literal_eval(line[st_indx:fi_indx])
			yield (37.)

		if (str(line.split(',')[0]) == str('molecular_weights_g/'+
			'mol_corresponding_to_component_names')):
			# find index of first [ and index of last ]
			icnt = 0 # count on characters
			for i in line:
				if i == '[':
					st_indx = icnt
					break
				icnt += 1 # count on characters

			for cnt in range(10):
				if line[-cnt] == ']':
					fi_indx = -cnt+1
					break

			y_MM = ast.literal_eval(line[st_indx:fi_indx])
			yield (7.)
		if (str(line.split(',')[0]) == 'molar_volumes_cm3/mol'):
			# find index of first [ and index of last ]
			icnt = 0 # count on characters
			for i in line:
				if i == '[':
					st_indx = icnt
					break
				icnt += 1 # count on characters

			for cnt in range(10):
				if line[-cnt] == ']':
					fi_indx = -cnt+1
					break

			MV = ast.literal_eval(line[st_indx:fi_indx])
			yield (17.)
		if (str(line.split(',')[0]) == 'nominal_molar_mass_g/mol'):
			# find index of first [ and index of last ]
			icnt = 0 # count on characters
			for i in line:
				if i == '[':
					st_indx = icnt
					break
				icnt += 1 # count on characters

			for cnt in range(10):
				if line[-cnt] == ']':
					fi_indx = -cnt+1
					break

			# nominal molar masses (g/mol)
			nom_mass = ast.literal_eval(line[st_indx:fi_indx])
			yield (12.)

		if (str(line.split(',')[0]) == 'pure_component_saturation_vapour_pressures_at_298.15K_Pa'):
			# find index of first [ and index of last ]
			icnt = 0 # count on characters
			for i in line:
				if i == '[':
					st_indx = icnt
					break
				icnt += 1 # count on characters

			for cnt in range(10):
				if line[-cnt] == ']':
					fi_indx = -cnt+1
					break

			PsatPa = ast.literal_eval(line[st_indx:fi_indx])
			yield (55.)
		if (str(line.split(',')[0]) == 'oxygen_to_carbon_ratios_of_components'):
			
			# find index of first [ and index of last ]
			icnt = 0 # count on characters
			for i in line:
				if line[icnt:icnt+2] == '[[':
					st_indx = icnt+1
					break
				icnt += 1 # count on characters

			for cnt in range(10):
				if line[-cnt:-cnt+2] == ']]':
					fi_indx = -cnt+1
					break
			
			OC = ast.literal_eval(line[st_indx:fi_indx])
			yield (60.)

		if (str(line.split(',')[0]) == 'hydrogen_to_carbon_ratios_of_components'):
			# find index of first [ and index of last ]
			icnt = 0 # count on characters
			for i in line:
				if line[icnt:icnt+2] == '[[':
					st_indx = icnt+1
					break
				icnt += 1 # count on characters
		
			for cnt in range(10):
				if line[-cnt:-cnt+2] == ']]':
					fi_indx = -cnt+1
					break
		
			HC = ast.literal_eval(line[st_indx:fi_indx])
			yield (65.)	

		if (str(line.split(',')[0]) == 'SMILES'):
			# find index of first [ and index of last ]
			icnt = 0 # count on characters
			for i in line:
				if line[icnt] == '[':
					st_indx = icnt
					break
				icnt += 1 # count on characters

			for cnt in range(10):
				if line[-cnt:] == ']':
					fi_indx = -cnt+1
					break
			rel_SMILES = ast.literal_eval(line[st_indx:fi_indx])
			yield (42.)

		# get indices of organic peroxy radicals
		if (str(line.split(',')[0]) == 'organic_peroxy_radical_index'):		
			# find index of first [ and index of last ]
			icnt = 0 # count on characters
			for i in line:
				if i == '[':
					st_indx = icnt+1
					break
				icnt += 1 # count on characters
			for cnt in range(10):
				if line[-cnt] == ']':
					fi_indx = -cnt
					break

			if (st_indx == len(line)+fi_indx): # if empty list
				continue
			else: # if list has contents
				group_indx['RO2i'] = list(np.array((
				line[st_indx:fi_indx].strip(' ').split(','))).astype('int'))			
			yield (22.)

		if (str(line.split(',')[0]) == 'organic_alkoxy_radical_index'):
			# find index of first [ and index of last ]
			icnt = 0 # count on characters
			for i in line:
				if i == '[':
					st_indx = icnt+1
					break
				icnt += 1 # count on characters
			for cnt in range(10):
				if line[-cnt] == ']':
					fi_indx = -cnt
					break

			if (st_indx == len(line)+fi_indx): # if empty list
				continue
			else: # if list has contents
				group_indx['ROi'] = list(np.array((line[st_indx:fi_indx].strip(' ').split(','))).astype('int'))			
			yield (27.)

		if (str(line.split(',')[0]) == 'organic_HOM_peroxy_radical_index'):
			
			# find index of first [ and index of last ]
			icnt = 0 # count on characters
			for i in line:
				if i == '[':
					st_indx = icnt+1
					break
				icnt += 1 # count on characters
			for cnt in range(10):
				if line[-cnt] == ']':
					fi_indx = -cnt
					break

			if (st_indx == len(line)+fi_indx): # if empty list
				continue
			else: # if list has contents
				group_indx['HOMRO2'] = list(np.array((line[
				st_indx:fi_indx].strip(' ').split(','))).astype('int'))
			yield (32.)

		if (str(line.split(',')[0]) == 'organic_HOMs_index'):
			
			# find index of first [ and index of last ]
			icnt = 0 # count on characters
			for i in line:
				if i == '[':
					st_indx = icnt+1
					break
				icnt += 1 # count on characters
			for cnt in range(10):
				if line[-cnt] == ']':
					fi_indx = -cnt
					break

			if (st_indx == len(line)+fi_indx): # if empty list
				continue
			else: # if list has contents
				group_indx['HOMs'] = list(np.array((line[st_indx:fi_indx].strip(' ').split(','))).astype('int'))	
			yield (32.)


		if (str(line.split(',')[0]) == 
			'factor_for_multiplying_ppb_to_get_molec/cm3_with_time'):

			# find index of first [ and index of last ]
			icnt = 0 # count on characters
			for i in line:
				if i == '[':
					st_indx = icnt
					break
				icnt += 1 # count on characters
			for cnt in range(10):
				if line[-cnt] == ']':
					fi_indx = -cnt+1
					break

			# conversion factor to change gas-phase 
			# concentrations from ppb into # molecules/cm3 
			# (air)
			Cfactor = ast.literal_eval(line[st_indx:fi_indx])
			
			yield (47.)

		for i in line.split(',')[1::]:
			
			if str(line.split(',')[0]) == 'number_of_size_bins':
				dlist.append(int(i))
				yield (1.)
			if str(line.split(',')[0]) == 'output_by_sim_sch_ext':
				dlist.append(str(i))
				yield (69.)
			if str(line.split(',')[0]) == 'output_by_sim_mv_ext':
				dlist.append(str(i))
				yield (70.)
			if str(line.split(',')[0]) == 'number_of_components':
				dlist.append(int(i))
				yield (2.)
			if str(line.split(',')[0]) == 'wall_on_flag_0forNO_>0forYES':
				dlist.append(int(i))
				yield (49.)
			if (str(line.split(',')[0]) == 'index_of_water'):
				i = i.strip('\n')
				i = i.strip('[[')
				i = i.strip(']]')
				i = i.strip('[')
				i = i.strip(']')
				i = i.strip(' ')
				dlist.append(int(i))
				yield (66.)
			if (str(line.split(',')[0]) == 'index_of_seed_components'):
				i = i.strip('\n')
				i = i.strip('[[')
				i = i.strip(']]')
				i = i.strip('[')
				i = i.strip(']')
				i = i.strip(' ')
				dlist.append(int(i))
				yield (67.)
			if (str(line.split(',')[0]) == 'space_mode'):
				i = i.strip('\n')
				i = i.strip('[')
				i = i.strip(']')
				i = i.strip(' ')
				i = i.strip('\'')
				dlist.append(str(i))
				yield (50.)
			if str(line.split(',')[0]) == 'simulation_computer_time(s)':
				i = i.strip('\n')
				i = i.strip('[')
				i = i.strip(']')
				i = i.strip(' ')
				dlist.append(float(i))
				yield (48.)
			if (str(line.split(',')[0]) == str('size_structure_0_for_'+
				'moving_centre_1_for_full_moving')):
				i = i.strip('\n')
				i = i.strip('[[')
				i = i.strip(']]')
				i = i.strip('[')
				i = i.strip(']')
				i = i.strip(' ')
				dlist.append(int(i))
				yield (68.)

		const[str(line.split(',')[0])] = dlist
	const_in.close()
	
	# extract required data from dictionary, note this prepared above
	num_sb = int((const['number_of_size_bins'])[0]) # number of size bins
	num_comp = int((const['number_of_components'])[0]) # number of components
	
	wall_on = const['wall_on_flag_0forNO_>0forYES'][0]
	
	space_mode = const['space_mode'][0]
	
	H2Oi = int((const['index_of_water'])[0]) # index of water
	
	seedi = const['index_of_seed_components'] # index of seed components
	siz_str = const['size_structure_0_for_moving_centre_1_for_full_moving']	
	
	if (v4_flag == 1): # if results in question saved in version 4 or later

		# load component molar masses (g/mol)
		load_path = str(self.dir_path + '/y_mw.npy') # path
		y_MM = np.load(load_path, allow_pickle=True)

		# cm^3/mol
		load_path = str(self.dir_path + '/MV.npy') # path
		MV = np.load(load_path, allow_pickle=True)
		
		load_path = str(self.dir_path + '/comp_namelist.npy') # path
		comp_names = (np.load(load_path, allow_pickle=True)).tolist()
			
		load_path = str(self.dir_path + '/rel_SMILES.npy') # path
		rel_SMILES = (np.load(load_path, allow_pickle=True)).tolist()
		
		# path
		load_path = str(self.dir_path + 
		'/pure_component_saturation_vapour_pressures_at_298p15K_Pa.npy')
		PsatPa = (np.load(load_path, allow_pickle=True)).tolist()

		# path
		load_path = str(self.dir_path + 
		'/pure_component_saturation_vp_at_startT_molec_percm3.npy')
		Psatmolecpercm3_0 = (np.load(load_path, allow_pickle=True)).tolist()
		
		load_path = str(self.dir_path + 
			'/oxygen_to_carbon_ratios_of_components.npy') # path
		OC = (np.load(load_path, allow_pickle=True)).tolist()
		
		load_path = str(self.dir_path + 
			'/hydrogen_to_carbon_ratios_of_components.npy') # path
		HC = (np.load(load_path, allow_pickle=True)).tolist()
		
		
		try: # this output added after others in version 4 of PyCHAM

			# path
			# note this definition of organic peroxy radicals is based
			# on SMILES strings, so not limited to the list of RO2 in the RO2 pool
			load_path = str(self.dir_path + '/organic_peroxy_radical_index.npy')
			group_indx['RO2i'] = (np.load(load_path, allow_pickle=True)).tolist()

			# path
			# note this definition of organic peroxy radicals is based
			# on the RO2 pool, so is limited to the list of RO2 in the RO2 pool
			try:
				load_path = str(self.dir_path + 
					'/organic_peroxy_radical_pool_index.npy')
				group_indx['RO2pooli'] = (np.load(load_path, 
					allow_pickle=True)).tolist()
			except:
				group_indx['RO2pooli'] = []

			# path
			load_path = str(self.dir_path + '/organic_alkoxy_radical_index.npy')
			group_indx['ROi'] = (np.load(load_path, allow_pickle=True)).tolist()
		
			# path
			load_path = str(self.dir_path + '/organic_HOM_peroxy_radical_index.npy')
			group_indx['HOMRO2'] = (np.load(load_path, allow_pickle=True)).tolist()

			load_path = str(self.dir_path + '/organic_HOMs_index.npy') # path
			group_indx['HOMs'] = (np.load(load_path, allow_pickle=True)).tolist()

			load_path = str(self.dir_path + '/organic_ROOR_index.npy') # path
			group_indx['ROOR'] = (np.load(load_path, allow_pickle=True)).tolist()

			load_path = str(self.dir_path + '/OOH_index.npy') # path
			group_indx['OOH'] = (np.load(load_path, allow_pickle=True)).tolist()
		
			load_path = str(self.dir_path + '/HOM_OOH_index.npy') # path
			group_indx['HOM_OOH'] = (np.load(load_path, allow_pickle=True)).tolist()	
		
			load_path = str(self.dir_path + '/OH_index.npy') # path
			group_indx['OH'] = (np.load(load_path, allow_pickle=True)).tolist()	
		
			load_path = str(self.dir_path + '/HOM_OH_index.npy') # path
			group_indx['HOM_OH'] = (np.load(load_path, allow_pickle=True)).tolist()	
		
			load_path = str(self.dir_path + '/carbonyl_index.npy') # path
			group_indx['carbonyl'] = (np.load(load_path, 
				allow_pickle=True)).tolist()	
		
			load_path = str(self.dir_path + '/HOM_carbonyl_index.npy') # path
			group_indx['HOM_carbonyl'] = (np.load(load_path, 
				allow_pickle=True)).tolist()	
		
			load_path = str(self.dir_path + '/NO3_index.npy') # path
			group_indx['NO3'] = (np.load(load_path, allow_pickle=True)).tolist()
		
			load_path = str(self.dir_path + '/HOM_NO3_index.npy') # path
			group_indx['HOM_NO3'] = (np.load(load_path, allow_pickle=True)).tolist()

			load_path = str(self.dir_path + '/PRAMpr_indx.npy') # path
			group_indx['PRAMpr_indx'] = (np.load(load_path, 
				allow_pickle=True)).tolist()	
		
			load_path = str(self.dir_path + '/PRAMcsmon_indx.npy') # path
			group_indx['PRAMcsmon_indx'] = (np.load(load_path, 
				allow_pickle=True)).tolist()
		
			load_path = str(self.dir_path + '/PRAMcsacc_indx.npy') # path
			group_indx['PRAMcsacc_indx'] = (np.load(load_path,
				allow_pickle=True)).tolist()


		except:
			group_indx['RO2i'] = []
			group_indx['RO2pooli'] = []
			group_indx['ROi'] = []
			group_indx['HOMRO2'] = []
			group_indx['HOMs'] = []
			group_indx['ROOR'] = []
			group_indx['OOH']  = []
			group_indx['HOM_OOH'] = []
			group_indx['OH'] = []
			group_indx['HOM_OH'] = []
			group_indx['carbonyl'] = []
			group_indx['HOM_carbonyl'] = []
			group_indx['NO3'] = []
			group_indx['HOM_NO3'] = []
			group_indx['PRAMpr_indx'] = []
			group_indx['PRAMcsmon_indx'] = []
			group_indx['PRAMcsacc_indx'] = []
	
	try:
		comp_time = (const["simulation_computer_time(s)"])[0]
	except:
		comp_time = 0.

	try:
		output_by_sim_sch_ext = (const["output_by_sim_sch_ext"])[0][0:-1]
	except:
		output_by_sim_sch_ext = 0.
	try:
		output_by_sim_mv_ext = (const["output_by_sim_mv_ext"])[0][0:-1]
	except:
		output_by_sim_mv_ext = 0.
	
	yield (75.)
	# withdraw index and names of components to plot the 
	# gas-phase concentration temporal profile of
	fname = str(self.dir_path + 
		'/components_with_initial_gas_phase_concentrations_specified')
	# check file size (bytes) to see if file contains more than just the header
	if (os.stat(fname).st_size > 123):
		indx_plot = np.loadtxt(fname, delimiter=',', skiprows=1, dtype='str')
		# chemical scheme names of components to plot
		try: # in case components with initial concentrations provided
			comp0 = indx_plot[1].tolist()
			# indices of components
			indx_plot = indx_plot[0].tolist()
			indx_plot = [int(i) for i in indx_plot]
		except: # in case components with initial concentrations not provided
			comp0 = []
			indx_plot = []
		
	else:
		comp0 = []
		indx_plot = []

	# withdraw times (s)
	fname = str(self.dir_path + '/time')
	t_array = np.loadtxt(fname, delimiter=',', skiprows=1)
	timehr = t_array/3600.0 # convert from s to hr
	
	
	try: # this output added on 23/02/2021
		# withdraw environmental conditions (temperature in first column (K)),
		# pressure in second column (Pa), relative humidity in third
		# column (0-1 fraction), transmission factor of light in fourth
		# column (fraction 0-1), times in rows
		fname = str(self.dir_path + '/chamber_environmental_conditions')
		cham_env = np.loadtxt(fname, delimiter=',', skiprows=1)
	except:
		cham_env = []
		
	# withdraw generation number of components, note this output added on 28/04/2022
	try:
		fname = str(self.dir_path + '/component_generation')
		if (os.path.getsize(fname) > 483): # only open if not empty
			gen_num = np.loadtxt(fname, delimiter=',', skiprows=1)
		else:
			gen_num = []
	except:
		gen_num = []
	
	yield (80.)

	# withdraw concentrations (ppb in gas, 
	# # molecules/cm^3 in particle and wall)
	fname = str(self.dir_path + 
	'/concentrations_all_components_all_times_gas_particle_wall')
	y = np.loadtxt(fname, delimiter=',', skiprows=1)
	
	# following will only load for certain simulation setups 
	# (mostly whether particles included)
	try:
		# withdraw the wall concentration of components due to 
		# particle deposition to wall
		fname = str(self.dir_path + 
			'/concentrations_all_components_all_times_on_wall_due_to_particle' +
			'_deposition_to_wall')
		yrec_p2w = np.loadtxt(fname, delimiter = ',', skiprows = 2)
	except:
		yrec_p2w = []
	
	try:
		# withdraw number-size distributions (# particles/cm3 (air))
		fname = str(self.dir_path + '/particle_number_concentration_dry')
		N = np.loadtxt(fname, delimiter=',', skiprows=1)
		if ((num_sb-wall_on) == 1): # if just one size bin, ensure two dimensions
			N = N.reshape(-1, 1)
	except:
		N = np.zeros((0, 0))
	
	yield (85.)

	try:
		# withdraw number-size distributions (# particles/cm3 (air))
		fname = str(self.dir_path + '/particle_number_concentration_wet')
		Nwet = np.loadtxt(fname, delimiter=',', skiprows=1)
		if ((num_sb-wall_on)<= 1): # if just one size bin, ensure two dimensions
			Nwet = Nwet.reshape(-1, 1)
	except:
		Nwet = np.zeros((0, 0))
	
	try:
		# particle radii (um)
		fname = str(self.dir_path + '/size_bin_radius')
		# skiprows=1 omits header
		x = np.loadtxt(fname, delimiter=',', skiprows=1)
	except:
		x = []

	try:
		# particle size bin bounds (radii) (um^3)
		fname = str(self.dir_path + '/size_bin_bounds')
		rbou_rec = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header
	except:
		rbou_rec =  np.zeros((0, 0))

	yield (90.)

	try: # in case this output is saved for a given simulation
		# withdraw consumptions (ug/m3)
		fname = str(self.dir_path + '/total_concentration_of_injected_components')
		if (os.path.getsize(fname) > 483): # only open if not empty
			tot_in_res = np.loadtxt(fname, delimiter=',', skiprows=1) # ug/m3
		else:
			tot_in_res = []
	except: # in case not saved, e.g. for older outputs
		tot_in_res = []
	
	yield (95.)
	
	# convert vapour pressure at starting simulation temperature 
	# from molecules/cm3 to Pa using ideal gas law
	PsatPa0 = Psatmolecpercm3_0/(si.Avogadro/((si.R*1.e6)*cham_env[0, 0]))
	
	# create a class to hold outputs
	class ro_outputs:
		
		sp = output_by_sim_sch_ext # chemical scheme path
		vp = output_by_sim_mv_ext # model variables path
		gi = group_indx # indices of groups of components
		# for each component, the generation number
		gen_numbers = gen_num 
		# hydrogen:carbon ratios for each component, this 
		# output added on 31/05/2022
		HyC = HC
		nominal_mass = nom_mass 
		nsb = num_sb
		nc = num_comp
		cfac = Cfactor
		yrec = y
		Nrec_dry = N
		rad = rbou_rec
		cen_size = x
		thr = timehr
		rSMILES = rel_SMILES
		comp_MW = y_MM
		Nrec_wet = Nwet
		names_of_comp = comp_names
		comp_MV = MV
		proc_time = comp_time
		wf = wall_on
		spacing = space_mode
		plot_indx = indx_plot
		init_comp = comp0
		part_to_wall = yrec_p2w
		vpPa = PsatPa
		vpPa0 = PsatPa0
		O_to_C = OC
		H2O_ind = H2Oi
		seed_ind = seedi
		siz_struc = siz_str
		env_cond = cham_env
		total_influx = tot_in_res
		
	self.ro_obj = ro_outputs() # create object to hold outputs
	
	return(self)

def retr_out_noncsv(output_by_sim, self): # similar to above function but for when non-csv files need interrogating

	import netCDF4 as nc # e.g. for PyCHAM output saved to netcdf file, or for EASY outputs
	# inputs: -------------------------------
	# output_by_sim - name of folders requested by the calling code to be looked at
	# ---------------------------------------
	# if a .dat file used (e.g. FACSIMILE output)
	
	if (output_by_sim[-4::] == '.dat'):
	
		datafile = open(output_by_sim)
		# flag stating whether column titles yet reached
		col_title = 0
		data_cnt = -1 # count on lines of data
		# prepare to create dictionary containing component names, 
		# times and concentrations
		data_dic = {}
		for line in datafile.readlines():
			dlist = [] # prepare to convert to python list
			# identify whether tabs used to separate columns (if not then 
			# spaces are)
			if '\t' in line:
				sep = '\t'
			if ' ' in line:
				sep = ' '
			if (col_title == 2): # ready to read in concentrations and times
				data_cnt += 1 # count on lines of data
			# loop through sections of line separated by a space and/or tab
			for i in line.split(sep):
				if (i == ''): # ignore spaces
					continue
				if (i.strip() == 'PRINT'): # identifier for header
					break # skip header
				if (i.strip()[0:4] == 'TIME'): # column titles
					col_title = 1 # flag that column titles being read
				if (col_title == 1):
					if str(i)[-1::] == '\n': # in case new line symbol included
						i = str(i)[0:-1]
					dlist.append(str(i)) # list column headers
				if (col_title == 2): # ready to read in concentrations and times
					dlist.append(float(i.strip()))
			if (col_title == 2):
				data_dic[str('data'+str(data_cnt))] = dlist
			if (col_title == 1):
				data_dic['col_title'] = dlist
				col_title = 2 # ready to read in concentrations and times
		# extract times (s), component names and concentrations with time (molecules/cm^3)
		# from dictionary
		comp_names = [i for i in data_dic['col_title'][1::]]
		Crec = np.zeros((data_cnt+1, len(comp_names))) # empty array for concentrations with time (molecules/cm^3)
		time_s = np.zeros((data_cnt+1, 1)) # empty array for times (s)
		for key in data_dic: # loop through dictionary keys
			if (key[0:4] == 'data'): # if this a useful entry
				rn = int(key[4::])# get row number
				time_s[rn] = data_dic[key][0] # times (s)
				Crec[rn, :] = data_dic[key][1::] # concentrations with time (molecules/cm^3)

	# if a .nc file used (e.g. EASY output)
	
	if (output_by_sim[-3::] == '.nc'):

		yield (0.)
		
		ds = nc.Dataset(output_by_sim) # open file
		
		timehr = ds['time'][:]/3600. # get time (seconds to hours)

		yield (10.)

		# get time saved at
		self.tsaved = ds.history
		# get model version used
		self.modv = ds.source

		yield (20.)
		
		# check whether SOA (anhydrous particle) mass concentration saved
		try:
			self.SOAmass = ds.variables[
			str('mass_concentration_of_secondary_particulate' +
				'_organic_matter_dry_aerosol_particles_in_air')]
		except:
			self.SOAmass = 'not saved'

		yield (30.)

		# number of particle size bins
		num_sb = ds['number_of_particle_size_bins'][:]

		# get total number of components
		num_comp = np.array((ds['number_of_components']))[0]

		# get wall on flag
		wall_on = ds['wall_on_flag'][:]

		# get chemical scheme names of components
		comp_names = ds['component_chemical_scheme_name'][:].tolist()

		# get SMILES strings of components in chemical scheme
		rel_SMILES = ds['SMILES_strings'][:]

		# get molar masses of components in chemical scheme
		y_MM = ds['molar_masses'][:]

		# get molar volumes of components in chemical scheme (cm^3/mol)
		y_MV = ds['molar_volumes'][:]

		# pure-component saturation vapour pressure of components at 298.15 K (Pa)
		PsatPa = ds['pure_component_saturation_vapour_pressure_of_components_at_298.15_K']

		# oxygen to carbon ratios of components 
		OC = ds['oxygen_to_carbon_ratios']

		# index of water component
		H2Oi = ds['water_component_index'][:]

		# index of seed component
		seedi = ds['seed_component_indices'][:]

		# index of organic peroxy radicals
		ro2i = ds['organic_peroxy_radical_component_indices'][:]

		# empty dictionary to contain indices of certain groups of components
		group_indx = {}
		try:
			group_indx['RO2i'] = ro2i
		except:
			group_indx['RO2i'] = [] # filler in case of no RO2i

		# chemical scheme names of components with initial concentration specified
		comp0 = ds['names_of_components_with_initial_gas_phase_concentrations_specified'][:]

		# indices of components with initial concentration specified
		indx_plot = ds['indices_of_components_with_initial_gas_phase_concentrations_specified'][:]
		
		yield (40.)
		
		
		# get the factor for converting molecules/cm^3 to ppb
		Cfactor = ds['factor_for_multiplying_ppb_to_get_molec/cm3_with_time'][:]

		yield (50.)
		
		# prepare to hold concentrations of components (molecules/cm^3)
		y = np.zeros((len(timehr), len(comp_names)))

		# loop through these components to get their gas-phase concentrations (molecules/cm^3)
		for ci in range(len(comp_names)):
			
			# set variable name
			var_name = str('number_concentration_of_' +
				'gas_phase_' + comp_names[ci] + '_molecules_in_air')

			# get values in molecules/cm^3
			try:
				y[:, ci] = ds[str('/concentrations_g/' + var_name)][:]
				# convert value to ppb to be consistent with gas-phase concetrations
				# returned by default saving type
				y[:, ci] = y[:, ci]/Cfactor
			except: # in case this component not saved, then continue to next component
				continue
		
		yield (60.)

		# hydrous and anhydrous particle number concentrations (# particles/cm^3 (air))
		N = ds['particle_number_concentration_anhydrous'][:] # anhydrous
		Nwet = ds['particle_number_concentration_anhydrous'][:] # hydrous

		# radius bounds (um)
		rbou_rec = ds['particle_radius_bounds'][:]

		# hydrous particle radius (um)
		x_rec = ds['hydrous_particle_radius'][:]

		# mode of spacing between particle size bins
		space_mode = ds['particle_size_bin_spacing_mode'][:]

		yield (70.)

		# open records needed to return rates of production and loss of
		# individual components

		# rates of reactions (molecules/cm^3/s), time in rows
		# and reactions in columns
		reactionRates = ds['rates_of_reaction_of_all_reactions'][:, :]

		# indices of components acting as reactants (columns) per reaction (rows)
		rindx_g = ds['reactant_indices'][:]

		# reactant stoichiometries per reaction (row) per reactant (column)
		rstoi_g = ds['reactant_stoichiometries'][:]

		# number of reactants per reaction
		nreac_g = ds['number_of_reactants_per_reaction'][:]

		# indices of components acting as products (columns) per reaction (rows)
		pindx_g = ds['product_indices'][:]

		# product stoichiometries per reaction (row) per reactant (column)
		pstoi_g = ds['product_stoichiometries'][:]

		# number of products per reaction
		nprod_g = ds['number_of_products_per_reaction'][:]

		# strings of gas-phase chemical reaction equation per reaction
		eq_str = ds['equations_per_reaction'][:]

		yield (80.)
		
		# create a class to hold outputs
		class ro_outputs:
	
			#sp = output_by_sim_sch_ext # chemical scheme path
			#vp = output_by_sim_mv_ext # model variables path
			gi = group_indx # indices of groups of components
			# for each component, the generation number
			#gen_numbers = gen_num 
			# hydrogen:carbon ratios for each component, this 
			# output added on 31/05/2022
			#HyC = HC
			#nominal_mass = nom_mass 
			nsb = num_sb
			nc = num_comp
			cfac = Cfactor
			yrec = y
			Nrec_dry = N
			rad = rbou_rec
			cen_size = x_rec
			thr = timehr
			rSMILES = rel_SMILES
			comp_MW = y_MM
			Nrec_wet = Nwet
			names_of_comp = comp_names
			comp_MV = y_MV
			#proc_time = comp_time
			wf = wall_on
			spacing = space_mode
			plot_indx = indx_plot
			init_comp = comp0
			#part_to_wall = yrec_p2w
			vpPa = PsatPa
			#vpPa0 = PsatPa0
			O_to_C = OC
			H2O_ind = H2Oi
			seed_ind = seedi
			#siz_struc = siz_str
			#env_cond = cham_env
			#total_influx = tot_in_res
			reactionRates_ret = np.array((reactionRates))
			rindx_g_ret = np.array((rindx_g))
			rstoi_g_ret = np.array((rstoi_g))
			nreac_g_ret = np.array((nreac_g))
			pindx_g_ret = np.array((pindx_g))
			pstoi_g_ret = np.array((pstoi_g))
			nprod_g_ret = np.array((nprod_g))
			eq_str_ret = eq_str
		
		self.ro_obj = ro_outputs() # create object to hold outputs
	
		yield (100.)
	
	return(self) # end of retr_out_noncsv function