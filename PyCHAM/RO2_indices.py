'''generates an array of component indices for the components that constitue a particular component type'''
# for peroxy radicals makes a two column array, with the first column giving the index of components
# included in the peroxy radical list based on component name and the
# second column containing the index based on component SMILE

import numpy as np

def RO2_indices(comp_namelist, RO2_names):

	# inputs: -----------------------------------------------------------------------------------------
	# comp_namelist - all chemical scheme names
	# RO2_names - all RO2 (non-HOM) names
	# ---------------------------------------------------------------------------------------------------
    
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
	
def HOMRO2_indices(comp_namelist):

	# inputs: -----------------------------------------------------------------------------------------
	# comp_namelist - all chemical scheme names
	# ---------------------------------------------------------------------------------------------------
    
	# store the names of HOMRO2 species which are present in the equation file
	# get a list of INDICES of RO2 that present in the equation file 
	# (or total species dict)
	# empty list for HOMRO2 indices
	HOMRO2_indices = []
    
	cin = 0 # count on components
	for name in comp_namelist: # loop through names of all components
        
		if ('API_' in name) or ('api_' in name):
			if ('RO2 in name'):
				# store the HOMRO2 index
				HOMRO2_indices.append(cin)
				
		cin += 1 # count on components
	# ensure elements in HOMRO2_indices are integer (iterable)
	HOMRO2_indices = (np.asarray(HOMRO2_indices, dtype=int)).reshape(-1, 1)
    
	return (HOMRO2_indices)