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
'''code for alterantive to eqn_pars when eqn_pars being skipped'''

def eqn_pars_skipper(self): # define function

	import numpy as np # matrix functions

	# imports: ---------------------------------
	# self - reference to PyCHAM class
	# ------------------------------------------

	if ('H2O' in self.con_infl_nam): # check if water in constant influx components 

		# index of water in continuous influx array
		wat_indx = (np.array(self.con_infl_nam) == 'H2O').reshape(-1)

		# do not allow continuous influx of water in the standard ode 
		# solver, instead deal with it inside the water ode-solver
		self.con_infl_C = np.delete(self.con_infl_C, wat_indx, axis=0)

		if self.comp_num in self.con_infl_indx:
			wat_indx = (np.array(self.con_infl_indx == self.comp_num))
			self.con_infl_indx = np.delete(self.con_infl_indx, wat_indx, axis=0)

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

	# use eqn_pars output from previous simulation
	rowvals = self.rowvals; colptrs = self.colptrs; jac_wall_indx = self.jac_wall_indx 
	jac_part_indx = self.jac_part_indx; jac_extr_indx = self.jac_extr_indx
	comp_num = self.comp_num; rel_SMILES = self.rel_SMILES; Pybel_objects = self.Pybel_objects
	Jlen = self.Jlen; comp_xmlname = self.comp_xmlname; comp_smil = self.comp_smil; 
	erf = 0.; err_mess = ''

	return(rowvals, colptrs, jac_wall_indx, 
		jac_part_indx, jac_extr_indx, comp_num, rel_SMILES, 
		Pybel_objects, Jlen, comp_xmlname, comp_smil, erf, err_mess)