'''code to print out the PyCHAM results that are equivalent to the Figure 4 Baker et al. 2024 results: doi.org/10.5194/acp-24-4789-2024'''


# import dependencies start ------------------------
import numpy as np
import ast
# import dependencies end --------------------------

# define function
def BakerFig4_reproduction():


	# name of file where experiment constants saved
	fname = str('model_and_component_constants')
	const_in = open(fname)
	const = {} # create empty dictionary to hold constants

	for line in const_in.readlines():
		
		dlist = [] # empty list to hold values

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


		for i in line.split(',')[1::]:
			if str(line.split(',')[0]) == 'number_of_components':
				dlist.append(int(i))
		const[str(line.split(',')[0])] = dlist
	const_in.close()
	# extract required data from dictionary, note this prepared above
	num_comp = int((const['number_of_components'])[0]) # number of components

	# open concentrations of all components (columns)
	# at all times (rows)
	fname = str('concentrations_all_components_all_times_gas_particle_wall')
	y = np.loadtxt(fname, delimiter=',', skiprows=1)

	# open times through simulation (s)
	fname = str('time')
	t_array = np.loadtxt(fname, delimiter=',', skiprows=1)
	
	# get indices of Baker product classes
	load_path = str('HOMMonBaker_indx.npy') # path
	HOMMonBaker_indx = (np.load(load_path, allow_pickle=True)).tolist()

	load_path = str('HOMFragBaker_indx.npy') # path
	HOMFragBaker_indx = (np.load(load_path, allow_pickle=True)).tolist()

	load_path = str('HOMRO2Baker_indx.npy') # path
	HOMRO2Baker_indx = (np.load(load_path, allow_pickle=True)).tolist()

	load_path = str('ROORBaker_indx.npy') # path
	ROORBaker_indx = (np.load(load_path, allow_pickle=True)).tolist()

	# get index of times needed (near end of low 
	# HO2/RO2 period (11 hours) and 
	# near end of high HO2/RO2 period (23 hours))
	t_indx_loHO2 = (t_array == 11.*3600.)
	t_indx_hiHO2 = (t_array == 23.*3600.)

	# get gas-phase concentrations (ppb)
	# of all species at this times
	y_loHO2 = np.squeeze(y[t_indx_loHO2, 0:num_comp])
	y_hiHO2 = np.squeeze(y[t_indx_hiHO2, 0:num_comp])

	# sum abundances of of product classes (ppb)
	loHO2_HOMMonBaker = sum(y_loHO2[HOMMonBaker_indx])
	loHO2_HOMFragBaker = sum(y_loHO2[HOMFragBaker_indx])
	loHO2_HOMRO2Baker = sum(y_loHO2[HOMRO2Baker_indx])
	loHO2_ROORBaker = sum(y_loHO2[ROORBaker_indx])
	hiHO2_HOMMonBaker = sum(y_hiHO2[HOMMonBaker_indx])
	hiHO2_HOMFragBaker = sum(y_hiHO2[HOMFragBaker_indx])
	hiHO2_HOMRO2Baker = sum(y_hiHO2[HOMRO2Baker_indx])
	hiHO2_ROORBaker = sum(y_hiHO2[ROORBaker_indx])

	# print results to command line
	print('PyCHAM estimates for Baker product classes at low HO2:RO2: ')
	print('All products: ', (loHO2_HOMFragBaker+loHO2_HOMMonBaker+loHO2_ROORBaker+
	loHO2_HOMRO2Baker))
	print('Fragments (C5-C9): ', loHO2_HOMFragBaker)
	print('Monomers: ', loHO2_HOMMonBaker)
	print('Accretion products: ', loHO2_ROORBaker)
	print('C10-HOM-RO2: ', loHO2_HOMRO2Baker)
	print(' ')
	print('PyCHAM estimates for Baker product classes at high HO2:RO2: ')
	print('All products: ', (hiHO2_HOMFragBaker+hiHO2_HOMMonBaker+hiHO2_ROORBaker+
	hiHO2_HOMRO2Baker))
	print('Fragments (C5-C9): ', hiHO2_HOMFragBaker)
	print('Monomers: ', hiHO2_HOMMonBaker)
	print('Accretion products: ', hiHO2_ROORBaker)
	print('C10-HOM-RO2: ', hiHO2_HOMRO2Baker)
	


	return() # end function
BakerFig4_reproduction() # call function