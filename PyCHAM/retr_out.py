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
'''code to open saved files and return useful variables'''
# e.g. called by PyCHAM plotting codes to obtain the required model outputs for evaluation

import numpy as np
import os
import ast

# define function, note that retr_out is called by the click202 function in gui.py 
# and it can return progress updates on loading results
def retr_out(self):
	
	# inputs: -------------------------------
	# self.dir_path - path of directory requested by the calling code to be looked at
	# ---------------------------------------

	# name of file where experiment constants saved
	fname = str(self.dir_path + '/model_and_component_constants')
	
	try: # try opening file
		const_in = open(fname)
	except:
		err_mess = str('Error - no such file ' + fname + ', please check it still exists')
		self.l203a.setText(err_mess)
		# set border around error message
		if (self.bd_pl == 1):
			self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
			self.bd_pl = 2
		else:
			self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
			self.bd_pl = 1
		return()
		
	const = {} # create empty dictionary to hold constants
	# empty dictionary to contain indices of certain groups of components
	group_indx = {}
	group_indx['RO2i'] = [] # filler in case of no RO2i
	group_indx['ROi'] = [] # filler in case of no ROi
	group_indx['HOMRO2'] = [] # filler in case of no HOMRO2

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

		if (str(line.split(',')[0]) == 'molecular_weights_g/mol_corresponding_to_component_names'):
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

			y_MW = ast.literal_eval(line[st_indx:fi_indx])
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
				group_indx['RO2i'] = list(np.array((line[st_indx:fi_indx].strip(' ').split(','))).astype('int'))			
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
				group_indx['HOMRO2'] = list(np.array((line[st_indx:fi_indx].strip(' ').split(','))).astype('int'))
			yield (32.)

		if str(line.split(',')[0]) == 'factor_for_multiplying_ppb_to_get_molec/cm3_with_time':

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

			# conversion factor to change gas-phase concentrations from # molecules/cm3 
			# (air) into ppb
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
			if (str(line.split(',')[0]) == 'size_structure_0_for_moving_centre_1_for_full_moving'):
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
	# withdraw index and names of components to plot the gas-phase concentration temporal profile of
	fname = str(self.dir_path + '/components_with_initial_gas_phase_concentrations_specified')
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
		# withdraw chamber environmental conditions (s)
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

	# withdraw concentrations (ppb in gas, # molecules/cm3 in particle and wall)
	fname = str(self.dir_path + '/concentrations_all_components_all_times_gas_particle_wall')
	y = np.loadtxt(fname, delimiter=',', skiprows=1)
	
	# following will only load for certain simulation setups (mostly whether particles included)
	
	try:
		# withdraw the wall concentration of components due to particle deposition to wall
		fname = str(self.dir_path + '/concentrations_all_components_all_times_on_wall_due_to_particle_deposition_to_wall')
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
		# particle sizes (um)
		fname = str(self.dir_path + '/size_bin_radius')
		x = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header
	except:
		x = []

	try:
		# particle size bin bounds (radii) (um3)
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

	# create a class to hold outputs
	class ro_outputs:
		sp = output_by_sim_sch_ext # chemical scheme path
		vp = output_by_sim_mv_ext # model variables path
		gi = group_indx # indices of groups of components
		gen_numbers = gen_num # for each component, the generation number
		HyC = HC # hydrogen:carbon ratios for each component, this output added on 31/05/2022
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
		comp_MW = y_MW
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
		O_to_C = OC
		H2O_ind = H2Oi
		seed_ind = seedi
		siz_struc = siz_str
		env_cond = cham_env
		total_influx = tot_in_res

	self.ro_obj = ro_outputs() # create object to hold outputs

	yield('Results Loaded')
	
	return(self)

def retr_out_noncsv(output_by_sim, comp_of_int): # similar to above function but for when non-csv files need interrogating
	
	import netCDF4 as nc # for EASY outputs
	
	# inputs: -------------------------------
	# output_by_sim - name of folders requested by the calling code to be looked at
	# comp_of_int - components of interest
	# ---------------------------------------
	
	# if a .dat file used (e.g. FACSIMILE output)
	if (output_by_sim[-4::] == '.dat'):
	
		datafile = open(output_by_sim)
	
		# flag stating whether column titles yet reached
		col_title = 0
		data_cnt = -1 # count on lines of data
		
		# prepare to create dictionary containing component names, times and concentrations
		data_dic = {}
		for line in datafile.readlines():
			
			dlist = [] # prepare to convert to python list
			
			# identify whether tabs used to separate columns (if not then spaces are)
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
			
	
		# extract times (s), component names and concentrations with time (molecules/cm3)
		# from dictionary
		comp_names = [i for i in data_dic['col_title'][1::]]
	
		Crec = np.zeros((data_cnt+1, len(comp_names))) # empty array for concentrations with time (molecules/cm3)
		time_s = np.zeros((data_cnt+1, 1)) # empty array for times (s)
	
		for key in data_dic: # loop through dictionary keys
			if (key[0:4] == 'data'): # if this a useful entry
				rn = int(key[4::])# get row number
				time_s[rn] = data_dic[key][0] # times (s)
				Crec[rn, :] = data_dic[key][1::] # concentrations with time (molecules/cm3)
	
		
		return(time_s, comp_names, Crec, [], [])
	
	# if a .nc file used (e.g. EASY output)
	if (output_by_sim[-3::] == '.nc'):
		
		ds = nc.Dataset(output_by_sim) # open file
		
		time_s = ds['time'][:] # get time (seconds)
		
		# in case want to see what variables are present
		#print(ds.variables); import ipdb; ipdb.set_trace()
		# empty array ready for component concentrations (molecules/cm3)
		Crec = np.zeros((len(time_s), len(comp_of_int)))
		
		# retrieve concentrations (molecules/cm3) of components of interest
		c_cnt = 0 # count on components
		for comp_name in comp_of_int:
			#if comp_name == 'OCRESOL':
			#	Crec[:, c_cnt] = ds[str(comp_name[1::]+'_0_0')][:]
			#else:
			Crec[:, c_cnt] = ds[str(comp_name+'_0')][:]
			c_cnt += 1 # count on components
			
		return(time_s, comp_of_int, Crec, [], [])
		
	else: # if file type unrecognised return fillers
		return([], [], [], [], [])