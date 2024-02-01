'''code to open saved PyCHAM files and return useful variables'''
# called by PyCHAM plotting codes to obtain the required model outputs for evaluation

import numpy as np

def retr_out(output_by_sim):
	
	# inputs: -------------------------------
	# output_by_sim - name of folder requested by the calling code to be looked at
	# ---------------------------------------
	
	# name of file where experiment constants saved
	fname = str(output_by_sim + '/model_and_component_constants')
	
	const_in = open(fname)

	const = {} # prepare to create dictionary
	for line in const_in.readlines():
		# convert to python list
		dlist = []
		for i in line.split(',')[1::]:
			
			if str(line.split(',')[0]) == 'number_of_size_bins':
				dlist.append(int(i))
			if str(line.split(',')[0]) == 'number_of_components' or str(line.split(',')[0]) == 'wall_on_flag_0forNO_1forYES':
				dlist.append(int(i))
			if str(line.split(',')[0]) == 'molecular_weights_g/mol_corresponding_to_component_names' or  str(line.split(',')[0]) == 'molecular_volumes_cm3/mol' or  str(line.split(',')[0]) == 'molar_volumes_cm3/mol':
				i = i.strip('\n')
				i = i.strip('[')
				i = i.strip(']')
				i = i.strip(' ')
				dlist.append(float(i))
			if str(line.split(',')[0]) == 'component_names' or str(line.split(',')[0]) == 'space_mode':
				i = i.strip('\n')
				i = i.strip('[')
				i = i.strip(']')
				i = i.strip(' ')
				i = i.strip('\'')
				dlist.append(str(i))
			if str(line.split(',')[0]) == 'factor_for_multiplying_ppb_to_get_molec/cm3_with_time':
				i = i.strip('\n')
				i = i.strip(' ')				
				i = i.strip('[[')
				i = i.strip('[')
				i = i.strip(']')
				dlist.append(float(i))
			if str(line.split(',')[0]) == 'simulation_computer_time(s)':
				i = i.strip('\n')
				i = i.strip('[')
				i = i.strip(']')
				i = i.strip(' ')
				dlist.append(float(i))
			
		const[str(line.split(',')[0])] = dlist
	const_in.close()
	num_sb = int((const['number_of_size_bins'])[0]) # number of size bins
	num_speci = int((const['number_of_components'])[0]) # number of species
	# conversion factor to change gas-phase concentrations from molecules/cc 
	# (air) into ppb 
	Cfactor = const['factor_for_multiplying_ppb_to_get_molec/cm3_with_time']
	PyCHAM_names = const['component_names']
	y_MW = const['molecular_weights_g/mol_corresponding_to_component_names']
	spec_namelist = const['component_names']
	try:
		wall_on = const['wall_on_flag_0forNO_1forYES'][0]
	except:
		wall_on = 1
	try:
		space_mode = const['space_mode'][0]
	except:
		space_mode = 0
	try:
		MV = const["molecular_volumes_cm3/mol"]
	except:
		MV = const["molar_volumes_cm3/mol"]

	try:
		speed = (const["simulation_computer_time(s)"])[0]
	except:
		speed = 10.
		
	# withdraw times (s)
	fname = str(output_by_sim+'/time')
	t_array = np.loadtxt(fname, delimiter=',', skiprows=1)
	timehr = t_array/3600.0 # convert from s to hr
	
	# withdraw concentrations
	fname = str(output_by_sim+'/concentrations_all_components_all_times_gas_particle_wall')
	y = np.loadtxt(fname, delimiter=',', skiprows=1)
	
	try:
		# withdraw number-size distributions (# particles/cc (air))
		fname = str(output_by_sim+'/particle_number_concentration_dry')
		N = np.loadtxt(fname, delimiter=',', skiprows=1)
	except:
		N = []
	
	try:
		# withdraw number-size distributions (# particles/cc (air))
		fname = str(output_by_sim+'/particle_number_concentration_wet')
		Nwet = np.loadtxt(fname, delimiter=',', skiprows=1)
	except:
		Nwet = []
	
	try:
		# particle sizes (um)
		fname = str(output_by_sim+'/size_bin_radius')
		x = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header
	except:
		x = []

	try:
		# particle size bin bounds (radii) (um3)
		fname = str(output_by_sim+'/size_bin_bounds')
		rbou_rec = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header
	except:
		rbou_rec = []
	
	return(num_sb, num_speci, Cfactor, y, N, rbou_rec, x, timehr, PyCHAM_names, y_MW, 
		Nwet, spec_namelist, MV, speed, wall_on, space_mode)