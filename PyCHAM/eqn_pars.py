##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2023 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
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
def extr_mech(int_tol, num_sb,
		drh_str, erh_str, sav_nam, pcont, self):
	
	# inputs: ----------------------------------------------------
	# self.sch_name - file name of chemical scheme
	# self.chem_sch_mrk - markers to identify different sections of 
	# 	the chemical scheme
	# self.xml_name - name of xml file
	# self.photo_path - path to file containing absorption 
	# 	cross-sections and quantum yields
	# self.con_infl_nam - chemical scheme names of components with 
	# 		constant influx
	# int_tol - integration tolerances
	# self.wall_on - marker for whether to include wall partitioning
	# num_sb - number of size bins (including any wall)
	# self.const_comp - chemical scheme name of components with 
	#	constant concentration
	# drh_str - string from user inputs describing 
	#	deliquescence RH (fraction 0-1) as function of temperature (K)
	# erh_str - string from user inputs describing 
	#	efflorescence RH (fraction 0-1) as function of temperature (K)
	# self.dil_fac - fraction of chamber air extracted/s
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
	[rrc, rrc_name, RO2_names, self] = sch_interr.sch_interr(total_list_eqn, self)
	
	# interrogate xml to list all component names and SMILES
	[err_mess_new, comp_smil, comp_name] = xml_interr.xml_interr(self.xml_name)

	# get equation information for chemical reactions
	[comp_list, Pybel_objects, comp_num, self] = eqn_interr.eqn_interr(comp_name, comp_smil, num_sb, self)                                                        	
	# prepare aqueous-phase and surface (e.g. wall) reaction matrices for applying to reaction rate calculation
	if (self.eqn_num[1] > 0 or self.eqn_num[2] > 0): # if aqueous-phase or surface (e.g. wall) reactions present
		[] = aq_mat_prep.aq_mat_prep(num_sb, comp_num, self)                                                  
	
	# if particle-phase equations are provided by particles not turned on then raise an error
	if (self.eqn_num[1] > 0 and num_sb-self.wall_on == 0):
		erf = 1 # raise error
		err_mess = str('Error: ' +str(self.eqn_num[1]) + ' particle-phase reactions were registered (from the chemical scheme input file), but no particle size bins have been invoked (number_size_bins variable in the model variables input file). Please ensure consistency. (message generated by eqn_pars.py module)')
	
	# get index of components with constant influx/concentration -----------
	# empty array for storing index of components with constant influx
	self.con_infl_indx = np.zeros((len(self.con_infl_nam)))
	self.con_C_indx = np.zeros((len(self.const_comp))).astype('int')
	delete_row_list = [] # prepare for removing rows of unrecognised components

	icon = 0 # count on constant influxes

	for i in range (len(self.con_infl_nam)):
		
		# water not included explicitly in chemical schemes but accounted for later in init_conc
		if (self.con_infl_nam[icon] == 'H2O'):
			self.con_infl_indx[icon] = int(comp_num)
			icon += 1 # count on constant influxes
			continue

		# if we want to remove constant influxes 
		# not present in the chemical scheme
		if (self.remove_influx_not_in_scheme == 1):

			try:
				# index of where components with constant influx occur in list of components
				self.con_infl_indx[icon] = self.comp_namelist.index(self.con_infl_nam[icon])
			except:
				# remove names of unrecognised components
				self.con_infl_nam = np.delete(self.con_infl_nam, (icon), axis=0)
				# remove emissions of unrecognised components
				self.con_infl_C = np.delete(self.con_infl_C, (icon), axis = 0)
				# remove empty indices of unrecognised components
				self.con_infl_indx = np.delete(self.con_infl_indx, (icon), axis = 0)
				icon -= 1 # count on constant influxes
		else:

			try:
				# index of where components with constant influx occur in list of components
				self.con_infl_indx[i] = self.comp_namelist.index(self.con_infl_nam[i])
			except:
				erf = 1 # raise error
				err_mess = str('Error: constant influx component with name ' +str(self.con_infl_nam[i]) + ' has not been identified in the chemical scheme, please check it is present and the chemical scheme markers are correct')
	
		icon += 1 # count on constant influxes

	if len(self.con_infl_nam) > 0:

		# get ascending indices of components with continuous influx
		si = self.con_infl_indx.argsort()
		# order constant influx indices, influx rates and names ascending by component index
		self.con_infl_indx = (self.con_infl_indx[si]).astype('int')
		self.con_infl_C = (self.con_infl_C[si])
		self.con_infl_nam = (self.con_infl_nam[si])

	for i in range (len(self.const_comp)):
		try:
			# index of where constant concentration components occur in list 
			# of components
			self.con_C_indx[i] = self.comp_namelist.index(self.const_comp[i])
		except:
			# if water then we know it will be the next component to be appended to the component list
			if (self.const_comp[i] == 'H2O'):
				self.con_C_indx[i] = len(self.comp_namelist)
			else: # if not water
				erf = 1 # raise error
				err_mess = str('Error: constant concentration component with name ' + str(self.const_comp[i]) + ' has not been identified in the chemical scheme, please check it is present and the chemical scheme markers are correct')

	
	if hasattr(self, 'obs_file') and self.obs_file != []: # if observation file provided for constraint

		import openpyxl
		import os
		#os.getcwd() + 
		self.obs_file = str(self.obs_file)
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
		
		# get indices of components with concentrations to fit observations
		self.obs_comp_i = np.zeros((len(self.obs_comp))).astype('int')
		for i in range(len(self.obs_comp)):
			self.obs_comp_i[i] = self.comp_namelist.index(self.obs_comp[i])

	# ---------------------------------------------------------------------

	if (comp_num in self.con_infl_indx):

		self.H2Oin = 1 # flag for water influx
		self.con_infl_H2O = self.con_infl_C[self.con_infl_indx==comp_num, :]
		
		# index of water in continuous influx array
		wat_indx = self.con_infl_indx==comp_num

		# do not allow continuous influx of water in the standard ode 
		# solver, instead deal with it inside the water ode-solver
		self.con_infl_C = np.delete(self.con_infl_C, wat_indx, axis=0)
		self.con_infl_indx = np.delete(self.con_infl_indx, wat_indx, axis=0)
	else:
		self.H2Oin = 0 # flag for no water influx

	
	# ensure integer
	self.con_infl_indx = self.con_infl_indx.astype('int')
	
	[rowvals, colptrs, jac_part_indx, jac_wall_indx, jac_extr_indx, self] = jac_setup.jac_setup(comp_num, 
		num_sb, (num_sb-self.wall_on), self)

	# call function to generate ordinary differential equation (ODE)
	# solver module, add two to comp_num to account for water and core component
	write_ode_solv.ode_gen(int_tol, rowvals, comp_num+2, 
			(num_sb-self.wall_on), 0, sav_nam, pcont, self)

	# call function to generate reaction rate calculation module
	write_rate_file.write_rate_file(rrc, rrc_name, 0, self)

	# call function to generate module that tracks change tendencies
	# of certain components
	write_dydt_rec.write_dydt_rec()
	
	# write the module for estimating deliquescence and efflorescence 
	# relative humidities as a function of temperature
	write_hyst_eq.write_hyst_eq(drh_str, erh_str)
	
	# get index of components in the peroxy radical list
	self = RO2_indices.RO2_indices(RO2_names, self)
	
	# get index of HOM-RO2 radicals
	self = RO2_indices.HOMRO2_indices(self)
	
	# get number of photolysis equations
	Jlen = photo_num.photo_num(self.photo_path)

	# in case equation parsing to be skipped in further simulations
	self.rowvals = rowvals; self.colptrs = colptrs; self.jac_wall_indx = jac_wall_indx 
	self.jac_part_indx = jac_part_indx; self.jac_extr_indx = jac_extr_indx; self.comp_num = comp_num
	self.rel_SMILES = comp_list; self.Pybel_objects = Pybel_objects; self.Jlen = Jlen 
	self.comp_xmlname = comp_name; self.comp_smil = comp_smil

	return(rowvals, colptrs, 
		jac_wall_indx, jac_part_indx, jac_extr_indx, comp_num, comp_list, 
		Pybel_objects, Jlen, comp_name, comp_smil, erf, err_mess, self)
