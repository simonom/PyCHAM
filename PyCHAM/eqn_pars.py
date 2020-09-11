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
	# num_sb - number of size bins
	# const_comp - chemical scheme name of components with 
	#	constant concentration
	# ------------------------------------------------------------

	print('Starting equation parsing')
	f_open_eqn = open(sch_name, mode='r') # open the chemical scheme file
	
	# read the file and store everything into a list
	total_list_eqn = f_open_eqn.readlines()
	f_open_eqn.close() # close file
	
	# interrogate scheme to list equations
	[eqn_list, eqn_num, rrc, rrc_name, RO2_names] = sch_interr.sch_interr(total_list_eqn, chem_sch_mrk)

	# interrogate xml to list component names and SMILES
	[comp_smil, comp_name] = xml_interr.xml_interr(xml_name)
	
	# get equation information for gas-phase reactions
	[rindx, rstoi, pindx, pstoi, reac_coef, comp_namelist, comp_list, 
		Pybel_objects, nreac, nprod, comp_num, jac_stoi, 
		jac_den_indx, njac, jac_indx, rowvals, colptrs, y_arr, 
		y_rind, uni_y_rind, y_pind, uni_y_pind, reac_col, prod_col, 
		rstoi_flat, pstoi_flat, rr_arr, rr_arr_p, jac_wall_indx, 
		jac_part_indx] = eqn_interr.eqn_interr(eqn_num[0], 
		eqn_list, chem_sch_mrk, comp_name, comp_smil, 0, num_sb, wall_on)
	
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
	# solver module, add two to comp_num to account for water and seed material
	write_ode_solv.ode_gen(con_infl_indx, int_tol, rowvals, wall_on, comp_num+2, 
			(num_sb-wall_on), con_C_indx)
	
	# get index of components in the peroxy radical list
	RO2_indx = write_RO2_indices.write_RO2_indices(comp_namelist, RO2_names)

	# call function to generate reaction rate calculation module
	write_rate_file.write_rate_file(reac_coef, rrc, rrc_name, 0)

	# call function to generate module that tracks change tendencies
	# of certain components
	write_dydt_rec.write_dydt_rec()

	# get number of photolysis equations
	Jlen = photo_num.photo_num(photo_path) 


	return(rindx, pindx, rstoi, pstoi, nreac, nprod, comp_num, jac_stoi, 
		njac, jac_den_indx, jac_indx, RO2_indx, comp_list, 
		Pybel_objects, eqn_num, comp_namelist, Jlen, y_arr, y_rind,
		uni_y_rind, y_pind, uni_y_pind, reac_col, prod_col, 
		rstoi_flat, pstoi_flat, rr_arr, rr_arr_p, rowvals, colptrs, 
		jac_wall_indx, jac_part_indx)
