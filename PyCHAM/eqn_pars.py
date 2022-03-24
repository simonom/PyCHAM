##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2022 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
#                                                                                       											 #
#    All Rights Reserved.                                                                									 #
#    This file is part of PyCHAM                                                         									 #
#                                                                                        											 #
#    PyCHAM is free software: you can redistribute it and/or modify it under              						 #
#    the terms of the GNU General Public License as published by the Free Software       					 #
#    Foundation, either version 3 of the License, or (at your option) any later          						 #
#    version.                                                                            										 #
#                                                                                        											 #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT                						 #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       			 #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              				 #
#    details.                                                                            										 #
#                                                                                        											 #
#    You should have received a copy of the GNU General Public License along with        					 #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                 							 #
#                                                                                        											 #
##########################################################################################
'''parses the input files to automatically create the solver file'''
# input files are interpreted and used to create the necessary
# arrays and python files to solve problem

import numpy as np
import sch_interr
import xml_interr
import eqn_interr
import photo_num
import RO2_indices
import write_dydt_rec
import write_ode_solv
import write_rate_file
import write_hyst_eq
import jac_setup
import aq_mat_prep

# define function to extract the chemical mechanism
def extr_mech(chem_sch_mrk, 
		con_infl_nam, int_tol, wall_on, num_sb,
		drh_str, erh_str, dil_fac, sav_nam, pcont, self):

	# inputs: ----------------------------------------------------
	# self.sch_name - file name of chemical scheme
	# chem_sch_mrk - markers to identify different sections of 
	# 	the chemical scheme
	# self.xml_name - name of xml file
	# self.photo_path - path to file containing absorption 
	# 	cross-sections and quantum yields
	# con_infl_nam - chemical scheme names of components with 
	# 		constant influx
	# int_tol - integration tolerances
	# wall_on - marker for whether to include wall partitioning
	# num_sb - number of size bins (including any wall)
	# self.const_comp - chemical scheme name of components with 
	#	constant concentration
	# drh_str - string from user inputs describing 
	#	deliquescence RH (fraction 0-1) as function of temperature (K)
	# erh_str - string from user inputs describing 
	#	efflorescence RH (fraction 0-1) as function of temperature (K)
	# dil_fac - fraction of chamber air extracted/s
	# sav_nam - name of folder to save results to
	# pcont - flag for whether seed particle injection is 
	#	instantaneous (0) or continuous (1)
	# self - reference to PyCHAM program
	# ------------------------------------------------------------
	
	# starting error flag and message (assumes no errors)
	erf = 0
	err_mess = ''

	f_open_eqn = open(self.sch_name, mode='r') # open the chemical scheme file
	# read the file and store everything into a list
	total_list_eqn = f_open_eqn.readlines()
	f_open_eqn.close() # close file
	
	# interrogate scheme to list equations
	[eqn_list, aqeqn_list, eqn_num, rrc, rrc_name, 
		RO2_names] = sch_interr.sch_interr(total_list_eqn, chem_sch_mrk)
	
	# interrogate xml to list all component names and SMILES
	[err_mess_new, comp_smil, comp_name] = xml_interr.xml_interr(self.xml_name)

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
		comp_num, RO_indx] = eqn_interr.eqn_interr(eqn_num, 
		eqn_list, aqeqn_list, chem_sch_mrk, comp_name, comp_smil, num_sb, wall_on)
		
	[rowvals, colptrs, jac_indx_g, jac_indx_aq, jac_part_indx, jac_wall_indx, jac_extr_indx] = jac_setup.jac_setup(jac_den_indx_g, njac_g, comp_num, num_sb, eqn_num, nreac_g, nprod_g, rindx_g, pindx_g, jac_indx_g, wall_on, nreac_aq, nprod_aq, rindx_aq, pindx_aq, jac_indx_aq, (num_sb-wall_on), dil_fac)
	
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
	
	# get index of components with constant influx/concentration -----------
	# empty array for storing index of components with constant influx
	con_infl_indx = np.zeros((len(con_infl_nam)))
	self.con_C_indx = np.zeros((len(self.const_comp))).astype('int')
	for i in range (len(con_infl_nam)):
		# water not included explicitly in chemical schemes but accounted for later in init_conc
		if (con_infl_nam[i] == 'H2O'):
			con_infl_indx[i] = comp_num
			continue
		try:
			# index of where components with constant influx occur in list of components
			con_infl_indx[i] = comp_namelist.index(con_infl_nam[i])
		except:
			erf = 1 # raise error
			err_mess = str('Error: constant influx component with name ' +str(con_infl_nam[i]) + ' has not been identified in the chemical scheme, please check it is present and the chemical scheme markers are correct')

	for i in range (len(self.const_comp)):
		try:
			# index of where constant concentration components occur in list 
			# of components
			self.con_C_indx[i] = comp_namelist.index(self.const_comp[i])
		except:
			erf = 1 # raise error
			err_mess = str('Error: constant concentration component with name ' + str(const_comp[i]) + ' has not been identified in the chemical scheme, please check it is present and the chemical scheme markers are correct')
	
	if hasattr(self, 'obs_file') and self.obs_file != []: # if observation file provided for constraint

		import openpyxl
		import os
		self.obs_file = str(os.getcwd() + self.obs_file)
		wb = openpyxl.load_workbook(filename = self.obs_file)
		sheet = wb['PyCHAMobs']
		# time (seconds) is in first column, components in later columns
		ic = 0 # count on row iteration
		for i in sheet.iter_rows(values_only=True): # loop through rows
			if (ic == 0):
				names_xlsx = i[1::] # get chemical scheme names of components
				self.obs_comp = [] # prepare for names
				for oc in names_xlsx:
					if oc != None:
						self.obs_comp.append(oc)

				# get number of components to consider
				nc_obs = len(self.obs_comp)
				# prepare to store observations, not forgetting a column for time
				self.obs = np.zeros((1, nc_obs+1))
				
			else:
				if ic > 1: # add row onto data array
					self.obs = np.concatenate((self.obs, (np.zeros((1, nc_obs+1)))), axis=0)
				self.obs[ic-1, :] = i[0:nc_obs+1]
				
			ic += 1 # count on row iteration
		
		wb.close() # close excel file
		
		# get indices of components with concentrations to fit obervations
		self.obs_comp_i = np.zeros((len(self.obs_comp))).astype('int')
		for i in range(len(self.obs_comp)):
			self.obs_comp_i[i] = comp_namelist.index(self.obs_comp[i])

	# ---------------------------------------------------------------------

	# call function to generate ordinary differential equation (ODE)
	# solver module, add two to comp_num to account for water and core component
	write_ode_solv.ode_gen(con_infl_indx, int_tol, rowvals, wall_on, comp_num+2, 
			(num_sb-wall_on), 0, eqn_num, dil_fac, sav_nam, pcont)

	# call function to generate reaction rate calculation module
	write_rate_file.write_rate_file(reac_coef_g, reac_coef_aq, rrc, rrc_name, 0)

	# call function to generate module that tracks change tendencies
	# of certain components
	write_dydt_rec.write_dydt_rec()
	
	# write the module for estimating deliquescence and efflorescence 
	# relative humidities as a function of temperature
	write_hyst_eq.write_hyst_eq(drh_str, erh_str)
	
	# get index of components in the peroxy radical list
	RO2_indx = RO2_indices.RO2_indices(comp_namelist, RO2_names)
	
	# get index of HOM-RO2 radicals
	HOMRO2_indx = RO2_indices.HOMRO2_indices(comp_namelist)
	
	# get number of photolysis equations
	Jlen = photo_num.photo_num(self.photo_path)

	return(rindx_g, pindx_g, rstoi_g, pstoi_g, nreac_g, nprod_g, jac_stoi_g, 
		njac_g, jac_den_indx_g, jac_indx_g, y_arr_g, y_rind_g,
		uni_y_rind_g, y_pind_g, uni_y_pind_g, reac_col_g, prod_col_g, 
		rstoi_flat_g, pstoi_flat_g, rr_arr_g, rr_arr_p_g, rowvals, colptrs, 
		jac_wall_indx, jac_part_indx, jac_extr_indx, comp_num, RO2_indx, RO_indx,
		HOMRO2_indx, comp_list, 
		Pybel_objects, eqn_num, comp_namelist, Jlen, 
		rindx_aq, rstoi_aq, pindx_aq, pstoi_aq, reac_coef_aq, 
		nreac_aq, nprod_aq, jac_stoi_aq, 
		jac_den_indx_aq, njac_aq, jac_indx_aq, 				
		y_arr_aq, y_rind_aq, uni_y_rind_aq, y_pind_aq, 
		uni_y_pind_aq, reac_col_aq, prod_col_aq, rstoi_flat_aq, pstoi_flat_aq, 
		rr_arr_aq, rr_arr_p_aq, comp_name, comp_smil, erf, err_mess, self)
