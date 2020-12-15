'''parses the input files to automatically create the solver file'''
# input files are interpreted and used to create the necessary
# arrays and python files to solve problem

import numpy as np
import sch_interr
import xml_interr
import eqn_interr
import photo_num
import write_RO2_indices
import write_dydt_rec
import write_ode_solv
import write_rate_file
import jac_setup
import aq_mat_prep

# define function to extract the chemical mechanism
def extr_mech(sch_name, chem_sch_mrk, xml_name, photo_path, 
		con_infl_nam, int_tol, wall_on, num_sb, const_comp):

	# inputs: ----------------------------------------------------
	# sch_name - file name of chemical scheme
	# chem_sch_mrk - markers to identify different sections of 
	# 			the chemical scheme
	# xml_name - name of xml file
	# photo_path - path to file containing absorption 
	# 	cross-sections and quantum yields
	# con_infl_nam - chemical scheme names of components with 
	# 		constant influx
	# int_tol - integration tolerances
	# wall_on - marker for whether to include wall partitioning
	# num_sb - number of size bins (including any wall)
	# const_comp - chemical scheme name of components with 
	#	constant concentration
	# ------------------------------------------------------------

	print('Starting equation parsing')
	f_open_eqn = open(sch_name, mode='r') # open the chemical scheme file
	
	# read the file and store everything into a list
	total_list_eqn = f_open_eqn.readlines()
	f_open_eqn.close() # close file
	
	# interrogate scheme to list equations
	[eqn_list, aqeqn_list, eqn_num, rrc, rrc_name, 
		RO2_names] = sch_interr.sch_interr(total_list_eqn, chem_sch_mrk)
	
	# interrogate xml to list component names and SMILES
	[comp_smil, comp_name] = xml_interr.xml_interr(xml_name)

	# get equation information for chemical reactions
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
		comp_num] = eqn_interr.eqn_interr(eqn_num, 
		eqn_list, aqeqn_list, chem_sch_mrk, comp_name, comp_smil, num_sb, wall_on)
	
	print(('Number of gas-phase chemical reactions identified: ' + str(eqn_num[0])))
	print(('Number of unique components identified in chemical scheme file: ' + str(comp_num)))

	[rowvals, colptrs, jac_indx_g, jac_indx_aq, jac_part_indx, jac_wall_indx] = jac_setup.jac_setup(jac_den_indx_g, njac_g, comp_num, num_sb, eqn_num, nreac_g, nprod_g, rindx_g, pindx_g, jac_indx_g, wall_on, nreac_aq, nprod_aq, rindx_aq, pindx_aq, jac_indx_aq, (num_sb-wall_on))
	
	# prepare aqueous-phase reaction matrices for applying to reaction rate calculation
	if (eqn_num[1] > 0): # if aqueous-phase reactions present
		[rindx_aq, rstoi_aq, pindx_aq, pstoi_aq, reac_coef_aq, 
			nprod_aq, jac_stoi_aq, njac_aq,
			jac_den_indx_aq, jac_indx_aq, 				
			y_arr_aq, y_rind_aq, uni_y_rind_aq, y_pind_aq, 
			uni_y_pind_aq, reac_col_aq, prod_col_aq, rstoi_flat_aq, pstoi_flat_aq, 
			rr_arr_aq, rr_arr_p_aq] = aq_mat_prep.aq_mat_prep(rindx_aq, rstoi_aq, 
			pindx_aq, pstoi_aq, reac_coef_aq, 
			nprod_aq, jac_stoi_aq, njac_aq, 
			jac_den_indx_aq, jac_indx_aq, 				
			y_arr_aq, y_rind_aq, uni_y_rind_aq, y_pind_aq, 
			uni_y_pind_aq, reac_col_aq, prod_col_aq, rstoi_flat_aq, pstoi_flat_aq, 
			rr_arr_aq, rr_arr_p_aq, num_sb, wall_on, eqn_num[1], comp_num) 
	
	
	print('Generating modules dependent on chemical scheme')
	# get index of components with constant influx/concentration -----------
	# empty array for storing index of components with constant influx
	con_infl_indx = np.zeros((len(con_infl_nam)))
	con_C_indx = np.zeros((len(const_comp)))
	for i in range (len(con_infl_nam)):
		# index of where constant components occur in list of components
		con_infl_indx[i] = comp_namelist.index(con_infl_nam[i])

	for i in range (len(const_comp)):
		# index of where constant concnetration components occur in list 
		# of components
		con_C_indx[i] = comp_namelist.index(const_comp[i])
	# ---------------------------------------------------------------------

	# call function to generate ordinary differential equation (ODE)
	# solver module, add two to comp_num to account for water and core component
	write_ode_solv.ode_gen(con_infl_indx, int_tol, rowvals, wall_on, comp_num+2, 
			(num_sb-wall_on), con_C_indx, 0, eqn_num)
	
	# get index of components in the peroxy radical list
	RO2_indx = write_RO2_indices.write_RO2_indices(comp_namelist, RO2_names)

	# call function to generate reaction rate calculation module
	write_rate_file.write_rate_file(reac_coef_g, reac_coef_aq, rrc, rrc_name, 0)

	# call function to generate module that tracks change tendencies
	# of certain components
	write_dydt_rec.write_dydt_rec()

	# get number of photolysis equations
	Jlen = photo_num.photo_num(photo_path)

	return(rindx_g, pindx_g, rstoi_g, pstoi_g, nreac_g, nprod_g, jac_stoi_g, 
		njac_g, jac_den_indx_g, jac_indx_g, y_arr_g, y_rind_g,
		uni_y_rind_g, y_pind_g, uni_y_pind_g, reac_col_g, prod_col_g, 
		rstoi_flat_g, pstoi_flat_g, rr_arr_g, rr_arr_p_g, rowvals, colptrs, 
		jac_wall_indx, jac_part_indx, comp_num, RO2_indx, comp_list, 
		Pybel_objects, eqn_num, comp_namelist, Jlen, 
		rindx_aq, rstoi_aq, pindx_aq, pstoi_aq, reac_coef_aq, 
		nreac_aq, nprod_aq, jac_stoi_aq, 
		jac_den_indx_aq, njac_aq, jac_indx_aq, 				
		y_arr_aq, y_rind_aq, uni_y_rind_aq, y_pind_aq, 
		uni_y_pind_aq, reac_col_aq, prod_col_aq, rstoi_flat_aq, pstoi_flat_aq, 
		rr_arr_aq, rr_arr_p_aq, comp_name, comp_smil)
