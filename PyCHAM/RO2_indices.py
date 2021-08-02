'''generates an array of component indices for the components that constitue a particular component type'''
# for peroxy radicals makes a two column array, with the first column giving the index of components
# included in the peroxy radical list based on component name and the
# second column containing the index based on component SMILE

import numpy as np

def RO2_indices(comp_namelist, RO2_names):
    
	# store the names of RO2 species which are present in the equation file
	# get a list of INDICES of RO2 that present in the equation file 
	# (or total species dict)
	# empty list for RO2_indices
	RO2_indices0 = []
	RO2_indices = []
    
	for name in RO2_names:
        
		if (name in comp_namelist):
			# get the RO2 index
			index0 = RO2_names.index(name)
			RO2_indices0.append(index0)
			# get the comp_namelist index for this RO2 species
			index1 = comp_namelist.index(name)
			RO2_indices.append(index1)
    
	# ensure elements in RO2_indices are integer (iterable)
	RO2_indices0 = (np.asarray(RO2_indices0, dtype=int)).reshape(-1, 1)
	RO2_indices = (np.asarray(RO2_indices, dtype=int)).reshape(-1, 1)
	RO2_indices = np.hstack((RO2_indices0, RO2_indices))
    
	return (RO2_indices)