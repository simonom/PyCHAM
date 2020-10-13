'''preparing the matrices for aqueous-phase reactions'''
# mostly just adjusting indices to account for number of size bins

import numpy as np

def aq_mat_prep(rindx, rstoi, pindx, pstoi, reac_coef, 
	nprod, jac_stoi, njac,
	jac_den_indx, jac_indx, 				
	y_arr, y_rind, uni_y_rind, y_pind, 
	uni_y_pind, reac_col, prod_col, rstoi_flat, pstoi_flat, 
	rr_arr, rr_arr_p, num_sb, wall_on, num_eqn, comp_num):

	# inputs: --------------------------------------------------------------
	# rindx - index for reactants
	# rstoi - stoichiometries of reactants
	# pindx - index for products
	# pstoi - stoichiometries of products
	# nprod - number of products per equation
	# njac - number of Jacobian elements affected per equation
	# jac_stoi - stoichiometries for jacobian
	# jac_den_indx_aq - indices for Jacobian denominators per equation
	# y_arr - indices of array for holding reactant concentrations
	# y_rind - indices of concentration array for reactants
	# uni_y_rind - indices of reactants in concentration array
	# y_pind - indices of concentration array for products
	# uni_y_pind - indices of products in concentration array
	# reac_col - columns for sparse matrix 
	# prod_col - columns for sparse matrix
	# rstoi_flat - flattened reactant stoichiometries
	# pstoi_flat - falttened product stoichiometries
	# rr_arr - reaction rate array indices for reactants
	# rr_arr_p - reaction rate indices for products
	# num_sb - number of size bins
	# wall_on - marker for whether wall on or off
	# num_eqn - number of aqueous reactions
	# comp_num - number of components
	# ----------------------------------------------------------------------

	num_asb = (num_sb-wall_on) # number of particle size bins
	# shapes of arrays
	rindxs = rindx.shape
	pindxs = pindx.shape
	y_arrl = len(y_arr)
	y_rindl = len(y_rind)
	y_pindl = len(y_pind)
	rr_arrl = len(rr_arr)
	rr_arr_pl = len(rr_arr_p)
	reac_coll = len(reac_col)
	prod_coll = len(prod_col)
	uni_y_rindl = len(uni_y_rind)
	uni_y_pindl = len(uni_y_pind)
	
	# tile where possible (ahead of vectorised operation)
	rstoi = np.tile(rstoi, (num_asb, 1))
	pstoi = np.tile(pstoi, (num_asb, 1))
	rstoi_flat = np.tile(rstoi_flat, (1, num_asb))
	pstoi_flat = np.tile(pstoi_flat, (1, num_asb))
	
	for sbi in range(num_sb-wall_on): # size bin loop
		
		if (sbi == 0): # first size bin - account for gas-phase indices prior
			rindx += (sbi+1)*(comp_num+2)
			pindx += (sbi+1)*(comp_num+2)
			y_rind = y_rind+(sbi+1)*(comp_num+2)
			y_pind = y_pind+(sbi+1)*(comp_num+2)
			uni_y_rind = uni_y_rind+(comp_num+2)
			uni_y_pind = uni_y_pind+(comp_num+2)

		if (sbi > 0): # larger size bin
			rindx = np.append(rindx, rindx[0:rindxs[0], 0:rindxs[1]]+(sbi)*(comp_num+2), axis=0)
			pindx = np.append(pindx, pindx[0:pindxs[0], 0:pindxs[1]]+(sbi)*(comp_num+2), axis=0)
			y_arr = np.append(y_arr, y_arr[0:y_arrl]+(sbi)*(max(y_arr[0:y_arrl])+1), axis=0)
			y_rind = np.append(y_rind, y_rind[0:y_rindl]+(sbi)*(comp_num+2))
			y_pind = np.append(y_pind, y_pind[0:y_pindl]+(sbi)*(comp_num+2))
			rr_arr = np.append(rr_arr, rr_arr[0:rr_arrl]+(sbi)*(max(rr_arr[0:rr_arrl])+1), axis=0)
			rr_arr_p = np.append(rr_arr_p, rr_arr_p[0:rr_arr_pl]+(sbi)*(max(rr_arr_p[0:rr_arr_pl])+1), axis=0)
			reac_col = np.append(reac_col, reac_col[reac_coll-1]+reac_col[-reac_coll+1::])
			prod_col = np.append(prod_col, prod_col[prod_coll-1]+prod_col[-prod_coll+1::])
			uni_y_rind = np.append(uni_y_rind, uni_y_rind[0:uni_y_rindl]+(sbi)*(comp_num+2), axis=0)
			uni_y_pind = np.append(uni_y_pind, uni_y_pind[0:uni_y_pindl]+(sbi)*(comp_num+2), axis=0)
			jac_stoi = np.append(jac_stoi, np.zeros((jac_stoi.shape[0], int(max(njac)))), axis=1)
			jac_den_indx = np.append(jac_den_indx, np.zeros((jac_den_indx.shape[0], int(max(njac)))), axis=1)
			for eqi in range(njac.shape[0]): # equation loop
				jac_stoi[eqi, int(sbi*njac[eqi]):int((sbi+1)*njac[eqi])] = jac_stoi[eqi, 0:int(njac[eqi])]
				jac_den_indx[eqi, int(sbi*njac[eqi]):int((sbi+1)*njac[eqi])] = jac_den_indx[eqi, 0:int(njac[eqi])]+(sbi)*(comp_num+2)
	
	# account for size bins
	njac = njac*(num_sb-wall_on)
	
	# ensure integer type
	njac = njac.astype(int)
	jac_den_indx = jac_den_indx.astype(int)

	return(rindx, rstoi, pindx, pstoi, reac_coef, 
		nprod, jac_stoi, njac,
		jac_den_indx, jac_indx, 				
		y_arr, y_rind, uni_y_rind, y_pind, 
		uni_y_pind, reac_col, prod_col, rstoi_flat, pstoi_flat, 
		rr_arr, rr_arr_p)
