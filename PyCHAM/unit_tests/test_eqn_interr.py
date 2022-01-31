'''module to test that equation interrogation functions correctly'''
# module that tests the functionality of interrogation of chemical reaction equations, 
# assumes calling from the PyCHAM home directory -
# note not all function outputs tested below - needs updating

import sys
import os
import numpy as np
# ensure modules can be seen 
# (assumes calling function is in the home folder)
sys.path.append(str(os.getcwd() + '/PyCHAM'))

import eqn_interr

def test_eqn_interr(): # define function

	# prepare inputs to testing module -----------------------
	
	# prepare correct answers to compare against --------------
	corr_rindx_g = np.array(([[0, 1], [3, 0]]))
	corr_rstoi_g = np.array(([[1., 1.], [1., 1.]]))
	corr_pindx_g = np.array(([[2, 0], [4, 5]]))
	corr_pstoi_g = np.array(([[1., 0.], [1., 1.]]))
	corr_y_arr_g = np.array(([0, 1, 2]))
	corr_y_arr_aq = np.array(([0, 1, 2]))
	
	# prepare function inputs --------------------------
	
	# number of gas-phase (first index) and particle-phase reactions (second index)
	num_eqn = [2, 2]
	# gas-phase equations in list of strings
	eqn_list = ['% 8.05D-16*EXP(-640/TEMP)*0.58 : APINENE + O3 = MGLYOX ;', '% 8.05D-16*EXP(-640/TEMP)*0.58 : ELVOC_O3 = ELVOC_OH + OH;']
	# particle-phase equations in list of strings
	aqeqn_list = ['$% 1.0D-14 : MGLYOX + AMM_SUL = C3H7N2O ;', '$% 1.0D-14 : C3H7N2O + C3H7N2O = C6H8N2 ;']
	# chemical scheme markers list
	chem_scheme_markers = ['%', 'RO2', '+', '', '', ';', '+', ';', '$', '%', ':', ';']
	# chemical scheme name list
	comp_name = ['APINENE', 'O3', 'MGLYOX', 'ELVOC_O3', 'ELVOC_OH', 'OH', 'AMM_SUL', 'C3H7N2O', 'C6H8N2']
	# SMILES strings
	comp_smil = ['CC1=CCC2CC1C2(C)C', '[O-][O+]=O', 'O=CC(=O)C', 'OCC1=COC(OO)CO2C(OO)C1C2(C(O))COOCC1=CCC2CC1C2(C)C', 'OCC1=COC(OO)CO2C(OO)C1C2(C(O))COOCC1=CC(O)C2C(O)C1C2(C)CO', '[OH]', 'O=S(=O)(O)O.N.N', 'CC(=O)CN', 'CC1=CN=C(C)C=N1']
	num_sb = 20 # number of size bins
	wall_on = 1# whether wall on
	
	
	# call on function to get outputs
	[rindx_g, rstoi_g, pindx_g, pstoi_g, reac_coef_g, 
			nreac_g, nprod_g, jac_stoi_g, 
			jac_den_indx_g, njac_g, jac_indx_g, 				
			y_arr_g, y_rind_g, uni_y_rind_g, y_pind_g, 
			uni_y_pind_g, reac_col_g, prod_col_g, rstoi_flat_g, pstoi_flat_g, 
			rr_arr_g, rr_arr_p_g, rindx_aq, rstoi_aq, pindx_aq, pstoi_aq, reac_coef_aq, 
			nreac_aq, nprod_aq, jac_stoi_aq, 
			jac_den_indx_aq, njac_aq, jac_indx_aq, 				
			y_arr_aq, y_rind_aq, uni_y_rind_aq, y_pind_aq, 
			uni_y_pind_aq, reac_col_aq, prod_col_aq, rstoi_flat_aq, pstoi_flat_aq, 
			rr_arr_aq, rr_arr_p_aq, comp_namelist, comp_list, Pybel_objects, 
			comp_num] = eqn_interr.eqn_interr(num_eqn, eqn_list, aqeqn_list, chem_scheme_markers, comp_name, 
		comp_smil, num_sb, wall_on)
	
	# run checks
	if ((rindx_g - corr_rindx_g).any()):
		print('rindx_g wrong')
	if ((rstoi_g - corr_rstoi_g).any()):
		print('rstoi_g wrong')
	if ((pindx_g - corr_pindx_g).any()):
		print('pindx_g wrong')
	if ((pstoi_g - corr_pstoi_g).any()):
		print('pstoi_g wrong')
	if ((y_arr_g - corr_y_arr_g).any()):
		print('y_arr_g wrong')
	if ((y_arr_aq - corr_y_arr_aq).any()):
		print('y_arr_aq wrong')
	
	

test_eqn_interr() # call on test