##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2024 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
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
'''preparing the matrices for aqueous-phase reactions'''
# mostly just adjusting indices to account for number of size bins

import numpy as np

def aq_mat_prep(num_sb, comp_num, self):

	# inputs: --------------------------------------------------------------
	# self.rindx_aq - index for reactants
	# self.rstoi_aq - stoichiometries of reactants
	# self.pindx_aq - index for products
	# self.pstoi_aq - stoichiometries of products
	# self.nprod_aq - number of products per equation
	# self.njac_aq - number of Jacobian elements affected per equation
	# self.jac_stoi_aq - stoichiometries for jacobian
	# self.jac_den_indx_aq - indices for Jacobian denominators per equation
	# self.y_arr_aq - indices of array for holding reactant concentrations
	# self.y_rind_aq - indices of concentration array for reactants
	# self.uni_y_rind_aq - indices of reactants in concentration array
	# self.y_pind_aq - indices of concentration array for products
	# self.uni_y_pind_aq - indices of products in concentration array
	# self.reac_col_aq - columns for sparse matrix 
	# self.prod_col_aq - columns for sparse matrix
	# self.rstoi_flat_aq - flattened reactant stoichiometries
	# self.pstoi_flat_aq - falttened product stoichiometries
	# self.rr_arr_aq - reaction rate array indices for reactants
	# self.rr_arr_p_aq - reaction rate indices for products
	# num_sb - number of size bins
	# self.wall_on - marker for whether wall on or off
	# self.eqn_num - number of aqueous-phase reactions and surface (e.g. wall) reactions
	# comp_num - number of components
	# self - reference to PyCHAM
	# ----------------------------------------------------------------------

	num_asb = (num_sb-self.wall_on) # number of particle size bins

	if (self.eqn_num[1] > 0): # if particle-phase reaction present
		# aqueous (particle) phase
		
		# shapes of arrays
		rindxs = self.rindx_aq.shape
		pindxs = self.pindx_aq.shape
		y_arrl = len(self.y_arr_aq)
		y_rindl = len(self.y_rind_aq)
		y_pindl = len(self.y_pind_aq)
		rr_arrl = len(self.rr_arr_aq)
		rr_arr_pl = len(self.rr_arr_p_aq)
		reac_coll = len(self.reac_col_aq)
		prod_coll = len(self.prod_col_aq)
		uni_y_rindl = len(self.uni_y_rind_aq)
		uni_y_pindl = len(self.uni_y_pind_aq)
	
		# tile where possible (ahead of vectorised operation)
		self.rstoi_aq = np.tile(self.rstoi_aq, (num_asb, 1))
		self.pstoi_aq = np.tile(self.pstoi_aq, (num_asb, 1))
		self.rstoi_flat_aq = np.tile(self.rstoi_flat_aq, (1, num_asb))
		self.pstoi_flat_aq = np.tile(self.pstoi_flat_aq, (1, num_asb))
	
		# count number of empty elements in rindx (due to not all reactions containing
		# the maximum number of unique reactants present in a reaction)
		empty_num = sum(sum(self.rindx_aq == -2))
	
		# now change fillers to zero, and note that although this suggests the first component
		# it makes no difference as the corresponding stoichiometry is 0
		self.rindx_aq[self.rindx_aq == -2] = 0
	
		for sbi in range(num_asb): # size bin loop
		
			if (sbi == 0): # first size bin - account for gas-phase indices prior
				self.rindx_aq += (comp_num+2)
				self.pindx_aq += (comp_num+2)
				self.y_rind_aq = self.y_rind_aq+(comp_num+2)
				self.y_pind_aq = self.y_pind_aq+(comp_num+2)
				self.uni_y_rind_aq = self.uni_y_rind_aq+(comp_num+2)
				self.uni_y_pind_aq = self.uni_y_pind_aq+(comp_num+2)

			# cumulative number of empty elements in rindx
			en_cum = empty_num*sbi

			if (sbi > 0): # larger size bin
				self.rindx_aq = np.append(self.rindx_aq, self.rindx_aq[0:rindxs[0], 0:rindxs[1]]+(sbi)*(comp_num+2), axis=0)
				self.pindx_aq = np.append(self.pindx_aq, self.pindx_aq[0:pindxs[0], 0:pindxs[1]]+(sbi)*(comp_num+2), axis=0)
				self.y_arr_aq = np.append(self.y_arr_aq, self.y_arr_aq[0:y_arrl]+(sbi)*(max(y_arr_aq[0:y_arrl])+1)+en_cum, axis=0)
				self.y_rind_aq = np.append(self.y_rind_aq, self.y_rind_aq[0:y_rindl]+(sbi)*(comp_num+2))
				self.y_pind_aq = np.append(self.y_pind_aq, self.y_pind_aq[0:y_pindl]+(sbi)*(comp_num+2))
				self.rr_arr_aq = np.append(self.rr_arr_aq, self.rr_arr_aq[0:rr_arrl]+(sbi)*(max(rr_arr_aq[0:rr_arrl])+1), axis=0)
				self.rr_arr_p_aq = np.append(self.rr_arr_p_aq, self.rr_arr_p_aq[0:rr_arr_pl]+(sbi)*(max(rr_arr_p_aq[0:rr_arr_pl])+1), axis=0)
				self.reac_col_aq = np.append(self.reac_col_aq, self.reac_col_aq[reac_coll-1]+reac_col[-reac_coll+1::])
				self.prod_col_aq = np.append(self.prod_col_aq, self.prod_col_aq[prod_coll-1]+prod_col[-prod_coll+1::])
				self.uni_y_rind_aq = np.append(self.uni_y_rind_aq, self.uni_y_rind_aq[0:uni_y_rindl]+(sbi)*(comp_num+2), axis=0)
				self.uni_y_pind_aq = np.append(self.uni_y_pind_aq, self.uni_y_pind_aq[0:uni_y_pindl]+(sbi)*(comp_num+2), axis=0)
				self.jac_stoi_aq = np.append(self.jac_stoi_aq, np.zeros((self.jac_stoi_aq.shape[0], int(max(self.njac_aq)))), axis=1)
				self.jac_den_indx_aq = np.append(self.jac_den_indx_aq, np.zeros((self.jac_den_indx_aq.shape[0], int(max(self.njac_aq)))), axis=1)
				for eqi in range(self.njac_aq.shape[0]): # equation loop
					self.jac_stoi_aq[eqi, int(sbi*self.njac_Aq[eqi]):int((sbi+1)*self.njac_aq[eqi])] = self.jac_stoi_aq[eqi, 0:int(self.njac_aq[eqi])]
					self.jac_den_indx_aq[eqi, int(sbi*self.njac_aq[eqi]):int((sbi+1)*self.njac_aq[eqi])] = self.jac_den_indx_aq[eqi, 0:int(self.njac_aq[eqi])]+(sbi)*(comp_num+2)
	
		# account for size bins
		self.njac_aq = self.njac_aq*(num_asb)
	
		# ensure integer type
		self.njac_aq = self.njac_aq.astype(int)
		self.jac_den_indx_aq = self.jac_den_indx_aq.astype(int)

	if (self.eqn_num[2] > 0): # if surface (e.g. wall) reaction present
		# now do for surface (e.g. wall) interactions
		ns = (self.wall_on) # number of surface (e.g. wall) size bins
	
		# shapes of arrays
		rindxs = self.rindx_su.shape
		pindxs = self.pindx_su.shape
		y_arrl = len(self.y_arr_su)
		y_rindl = len(self.y_rind_su)
		y_pindl = len(self.y_pind_su)
		rr_arrl = len(self.rr_arr_su)
		rr_arr_pl = len(self.rr_arr_p_su)
		reac_coll = len(self.reac_col_su)
		prod_coll = len(self.prod_col_su)
		uni_y_rindl = len(self.uni_y_rind_su)
		uni_y_pindl = len(self.uni_y_pind_su)
	
		# tile where possible (ahead of vectorised operation)
		self.rstoi_su = np.tile(self.rstoi_su, (ns, 1))
		self.pstoi_su = np.tile(self.pstoi_su, (ns, 1))
		self.rstoi_flat_su = np.tile(self.rstoi_flat_su, (1, ns))
		self.pstoi_flat_su = np.tile(self.pstoi_flat_su, (1, ns))
	
		# count number of empty elements in rindx (due to not all reactions containing
		# the maximum number of unique reactants present in a reaction)
		empty_num = sum(sum(self.rindx_su == -2))
	
		# now change fillers to zero, and note that although this suggests the first component
		# it makes no difference as the corresponding stoichiometry is 0
		self.rindx_su[self.rindx_su == -2] = 0
	
		for sbi in range(ns): # surface number loop
		
			if (sbi == 0): # first surface - account for gas-phase and particle-phase indices prior
				self.rindx_su += (comp_num+2+num_asb*(comp_num+2))
				self.pindx_su += (comp_num+2+num_asb*(comp_num+2))
				self.y_rind_su = self.y_rind_su+(comp_num+2+num_asb*(comp_num+2))
				self.y_pind_su = self.y_pind_su+(comp_num+2+num_asb*(comp_num+2))
				self.uni_y_rind_su = self.uni_y_rind_su+(comp_num+2+num_asb*(comp_num+2))
				self.uni_y_pind_su = self.uni_y_pind_su+(comp_num+2+num_asb*(comp_num+2))

			# cumulative number of empty elements in rindx
			en_cum = empty_num*sbi

			if (sbi > 0): # larger surfaces
				self.rindx_su = np.append(self.rindx_su, self.rindx_su[0:rindxs[0], 0:rindxs[1]]+(sbi)*(comp_num+2), axis=0)
				self.pindx_su = np.append(self.pindx_su, self.pindx_su[0:pindxs[0], 0:pindxs[1]]+(sbi)*(comp_num+2), axis=0)
				self.y_arr_su = np.append(self.y_arr_su, self.y_arr_su[0:y_arrl]+(sbi)*(max(y_arr[0:y_arrl])+1)+en_cum, axis=0)
				self.y_rind_su = np.append(self.y_rind_su, self.y_rind_su[0:y_rindl]+(sbi)*(comp_num+2))
				self.y_pind_su = np.append(self.y_pind_su, self.y_pind_su[0:y_pindl]+(sbi)*(comp_num+2))
				self.rr_arr_su = np.append(self.rr_arr_su, self.rr_arr_su[0:rr_arrl]+(sbi)*(max(rr_arr[0:rr_arrl])+1), axis=0)
				self.rr_arr_p_su = np.append(self.rr_arr_p_su, self.rr_arr_p_su[0:rr_arr_pl]+(sbi)*(max(rr_arr_p[0:rr_arr_pl])+1), axis=0)
				self.reac_col_su = np.append(self.reac_col_su, self.reac_col_su[reac_coll-1]+reac_col[-reac_coll+1::])
				self.prod_col_su = np.append(self.prod_col_su, self.prod_col_su[prod_coll-1]+prod_col[-prod_coll+1::])
				self.uni_y_rind_su = np.append(self.uni_y_rind_su, self.uni_y_rind_su[0:uni_y_rindl]+(sbi)*(comp_num+2), axis=0)
				self.uni_y_pind_su = np.append(self.uni_y_pind_su, self.uni_y_pind_su[0:uni_y_pindl]+(sbi)*(comp_num+2), axis=0)
				self.jac_stoi_su = np.append(self.jac_stoi_su, np.zeros((self.jac_stoi_su.shape[0], int(max(self.njac_su)))), axis=1)
				self.jac_den_indx_su = np.append(self.jac_den_indx_su, np.zeros((self.jac_den_indx_su.shape[0], int(max(self.njac_su)))), axis=1)
				for eqi in range(self.njac_su.shape[0]): # equation loop
					self.jac_stoi_su[eqi, int(sbi*self.njac_su[eqi]):int((sbi+1)*self.njac_su[eqi])] = self.jac_stoi_su[eqi, 0:int(self.njac_su[eqi])]
					self.jac_den_indx_su[eqi, int(sbi*self.njac_su[eqi]):int((sbi+1)*self.njac_su[eqi])] = self.jac_den_indx_su[eqi, 0:int(self.njac_su[eqi])]+(sbi)*(comp_num+2)
	
		# account for number of surfaces
		self.njac_su = self.njac_su*(ns)
	
		# ensure integer type
		self.njac_su = self.njac_su.astype(int)
		self.jac_den_indx_su = self.jac_den_indx_su.astype(int)

	return()
