'''estimate of the diffusion volume of components'''
# following the method of Fuller et al. (1966) (doi.org/10.1021/ie50677a007) and Fuller et al. (1969) (doi.org/10.1021/j100845a020), the calculated, with the volume increments given in Table 4.1 
# of the Taylor (1993) textbook Multicomponent Mass Transfer, ISBN: 0-471-57417-1
# diffusion volumes will be used to estimate gas-phase diffusion coefficients
# documents for pybel: http://openbabel.org/docs/dev/UseTheLibrary/Python_PybelAPI.html#pybel.Molecule.sssr

import numpy as np
import pybel
import math

# define function
def diff_vol_est(Pybel_object):

	# inputs: -------------------------------------------------------------------------
	# Pybel_object - Pybel objects of components
	# -----------------------------------------------------------------------------------
	
	# simple molecules given in Table 4.1 of Taylor (1993)
	# SMILE strings for simple molecules
	simp_molec = ['[He]', '[Ne]', '[Ar]', '[Kr]', '[Xe]', '[H][H]', r'N#N', 'O=O', 'CO', 'O=C=O', 'NO=N', '[H]N([H])[H]', 'O', 'S(F)(F)(F)(F)(F)F', 'ClCl', 'BrBr', 'O=S=O']
	
	# diffusion volumes for simple molecules given in Table 4.1 of  the Taylor (1993) textbook 
	# Multicomponent Mass Transfer, ISBN: 0-471-57417-1
	simp_molec_vol = [2.67, 5.98, 16.2, 24.5, 32.7, 6.12, 18.5, 16.3, 18., 26.7, 35.9, 20.7, 13.1, 71.3, 38.4, 69., 41.8]
	
	Pybel_object_ref = [] # empty list
	# convert to SMARTS (pybel_objects)
	# generate pybel objects from SMILES
	for i in simp_molec:	
		Pybel_object_ref.append(pybel.readstring('smi', i))
	
	diff_vol = np.zeros((len(Pybel_object))) # empty results array
	
	for compi in range(len(Pybel_object)): # component loop
	
		simp_cnt = 0 # count
	
		# check if component is a simple molecule
		for i_smrt in Pybel_object_ref: # loop through simple molecules
			if (Pybel_object[compi].formula == i_smrt.formula):
				diff_vol[compi] = simp_molec_vol[simp_cnt]
				break
			simp_cnt += 1
			
		# if a simple molecule move onto next component
		if (simp_cnt < len(Pybel_object_ref)):
			continue
			
		# SMARTS for the atoms given in Table 4.1 of  the Taylor (1993) textbook 
		# Multicomponent Mass Transfer, ISBN: 0-471-57417-1
		smarts = []
		smarts.append(pybel.Smarts('[#6]'))
		smarts.append(pybel.Smarts('[H]'))
		smarts.append(pybel.Smarts('[#8]'))
		smarts.append(pybel.Smarts('[#7]'))
		smarts.append(pybel.Smarts('[#16]'))
		smarts.append(pybel.Smarts('[#9]'))
		smarts.append(pybel.Smarts('[#17]'))
		smarts.append(pybel.Smarts('[#35]'))
		smarts.append(pybel.Smarts('[#53]'))
	
		# count atoms present, final element will be used for number of rings
		atm_cnt = np.zeros((len(smarts)+1))
	
		# atomic and molecular diffusion volume increments (Table 4.1 Taylor 1993)
		vol_inc = np.array((15.9, 2.31, 6.11, 4.54, 22.9, 14.7, 21., 21.9, 29.8, -18.3))
	
		i_cnt = 0
		for i_smrt in smarts: # loop through reference groups
			atm_cnt[i_cnt] = len(i_smrt.findall(Pybel_object[compi]))
			i_cnt += 1
			
		# smarts for rings
		smarts = []
		smarts.append(pybel.Smarts('[r3]'))
		smarts.append(pybel.Smarts('[r4]'))
		smarts.append(pybel.Smarts('[r5]'))
		smarts.append(pybel.Smarts('[r6]'))
		smarts.append(pybel.Smarts('[r7]'))
		smarts.append(pybel.Smarts('[r8]'))
		smarts.append(pybel.Smarts('[r9]'))
		smarts.append(pybel.Smarts('[r10]'))
		
		# check for rings
		num_atm = 3 # count on number of atoms in rings
		for i_smrt in smarts: # loop through reference groups
			atm_cnt[i_cnt]  += math.ceil(len(i_smrt.findall(Pybel_object[compi]))/num_atm)
			num_atm += 1 # count on number of atoms in rings

		# diffusion volume
		diff_vol[compi] = np.sum((atm_cnt*vol_inc))

	return(diff_vol)