'''generates an array of component indices for the components contained in the peroxy radical list'''
# makes a two column array, with the first giving the index of components
# included in the peroxy radical list based on component name and the
# second column containing the index based on component SMILE

import numpy as np

def write_RO2_indices(smiles_array, RO2_names):
    
    # store the names of RO2 species which are present in the equation file
    # get a list of INDICES of RO2 that present in the equation file 
    # (or total species dict)
    # empty list for RO2_indices
    RO2_indices0 = []
    RO2_indices = []
    
    for name in RO2_names:
        
        if name in smiles_array:
            # get the RO2 index
            index0 = RO2_names.index(name)
            RO2_indices0.append(index0)
            # get the smiles_array index for this RO2 species
            index1 = smiles_array.index(name)
            RO2_indices.append(index1)
    
    # Ensure elements in RO2_indices are int (iterable)
    RO2_indices0 = (np.asarray(RO2_indices0, dtype=int)).reshape(-1, 1)
    RO2_indices = (np.asarray(RO2_indices, dtype=int)).reshape(-1, 1)
    RO2_indices = np.hstack((RO2_indices0, RO2_indices))
    
    return (RO2_indices)
