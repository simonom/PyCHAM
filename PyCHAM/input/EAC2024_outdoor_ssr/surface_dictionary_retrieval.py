'''code to retrieve surface deposition velocities from the surface_dictionary module of INCHEM-Py'''
def surface_dictionary_retrieval(): # define function
	# import dependencies
	import surface_dictionary
	import numpy as np # for containing and saving data

	# prepare a list to hold component names
	comp_names = []
	# prepare a list to hold deposition velocities
	dv = []

	# get data from surface_dictionary.py
	surface_dict = surface_dictionary.surface_deposition(1.,'TRUE','TRUE')

	for key, value in surface_dict.items(): # loop through keys and their value from dcitionary
		comp_names.append(key[0:key.index("_")])
		dv.append(value)

	# now convert to numpy matrix
	mat_form = np.concatenate((np.array((comp_names)).reshape(-1, 1), np.array((dv)).reshape(-1, 1)), axis = 1)

	# now save as csv
	np.savetxt("chems_and_dep_vel.csv", mat_form, fmt = '%s', delimiter=",")

	return

surface_dictionary_retrieval() # call function