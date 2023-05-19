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
'''preparing the inputs for the ode solver Jacobian'''
# called once the gas-phase and particle-phase equations have been interrogated

# required modules
import numpy as np

def jac_setup(comp_num, num_sb, num_asb, self):

	# inputs ---------------------------------------------
	# comp_num - number of components
	# num_sb - number of size bins (including any wall)
	# self.eqn_num - number of chemical reactions
	# gas-phase reactions
	# self.nreac_g - number of reactants per reaction
	# self.nprod_g - number of products per reaction
	# self.rindx_g - index of reactants per equation
	# self.pindx_g - index of products per equation
	# self.jac_indx_g - index of Jacobian affected per equation
	# self.wall_on - flag for whether wall on or off
	# aqueous-phase reactions:
	# self.nreac_aq - number of reactants per equation
	# self.nprod_aq - number of products per equation
	# self.rindx_aq - index of reactants per equation
	# self.pindx_aq - index of prodcuts per equation
	# self.jac_indx_aq - index of jacobian affected per equation
	# surface (e.g. wall) reactions:
	# self.nreac_su - number of reactants per equation
	# self.nprod_su - number of products per equation
	# self.rindx_su - index of reactants per equation
	# self.pindx_su - index of prodcuts per equation
	# self.jac_indx_su - index of jacobian affected per equation
	# other inputs:
	# num_asb - number of actual particle size bins (excluding wall)
	# self.dil_fac - fraction of chamber air removed/s
	# self - reference to PyCHAM
	# ----------------------------------------------------

	# rows in Jacobian affected by a process		 	
	rowvals = np.empty((0))
	# indices of rowvals representing each component being 
	# differentiated by in Jacobian. Note, needs to start
	# as zero to represent first index of rowval
	colptrs = np.zeros((1))
	
	# track columns of Jacobian affected
	col_tr = 0

	# index for Jacobian - note can only be done after all equations checked as only then
	# is the total number of unique components known
	for eqni in range(self.eqn_num[0]): # equation loop
	
		# total number of components in this equation
		tot_comp = self.nreac_g[eqni]+self.nprod_g[eqni]
		# combined index of reactants and products in this equation
		totindx = np.append(self.rindx_g[eqni, 0:self.nreac_g[eqni]], self.pindx_g[eqni, 0:self.nprod_g[eqni]])
		
		# reactant loop (equivalent to columns of 2D Jacobian)
		for i in range(self.nreac_g[eqni]):
			# index of flattened Jacobian being affected in 
			# this equation, note this doesn't need to 
			# account for water and core component as below (when unique elements found) all
			# full Jacobian matrix elements that aren't indexed are
			# omitted to create the sparse Jacobian matrix (add one to num_sb to account for gas-phase)
			self.jac_indx_g[eqni, i*tot_comp:(i+1)*tot_comp] = self.rindx_g[eqni, i]*((comp_num)*(num_sb+1))+(totindx)
			# check if rowvals already sufficiently long to contain this reactant  
			if self.rindx_g[eqni, i]<=(len(colptrs)-2):

				# start index of rowvals where it is represented
				r_st = int(colptrs[self.rindx_g[eqni, i]])
				# end index of rowvals where it is represented
				r_en = int(colptrs[self.rindx_g[eqni, i]+1])
				# check if the components affected have been 
				# represented in previous equations
				for i2 in totindx:
					# if this reactant has affected this component in 
					# a previous reaction then jac_indx_g (set above) will
					# sum its effect and no need to account for
					# the Jacobian index with rowvals and 
					# colptrs				
					if (any(rowvals[r_st:r_en] == i2)):
						continue
					else: # account for in Jacobian
						new_el = np.array((float(i2))).reshape(1)
						# index of where new index enters rowvals (to ensure
						# indices per Jacobian column are present in
						# ascending order)
						r_in = r_st+int(sum(rowvals[r_st:r_en] < new_el))
						rowvals = np.concatenate([rowvals[0:r_in], new_el, rowvals[r_in::]])
						colptrs[self.rindx_g[eqni, i]+1::] += 1
						r_en += 1
						
						
			else: # if not accounted for in previous reaction
				# rows in Jacobian affected by this reactant 
				# in this equation	 	
				rowvals = np.append(rowvals, np.unique(totindx))
				while (self.rindx_g[eqni, i]>col_tr):
					colptrs = np.append(colptrs, colptrs[-1])
					col_tr += 1
				# indices of rowvals
				colptrs = np.append(colptrs, colptrs[-1]+len(np.unique(totindx)))
				# track number of columns of Jacobian affected				
				col_tr += 1
	
	# because above self.jac_indx_g contains the indices of data assuming
	# a full, rather than sparse matrix, now correct to affect
	# the sparse matrix by reducing indices where gaps occur in
	# the complete data matrix
	# get unique values in self.jac_indx_g
	uni = np.unique(self.jac_indx_g)
	# now reduce these indices based on the number of indices
	# omitted
	diff = 0 # prepare for number of empty points
	for i in range(1, len(uni)):
		diff += (uni[i]-uni[i-1])-1
		self.jac_indx_g[self.jac_indx_g==uni[i]] -= diff
	
	# colptrs responds to appendages to rowvals above and therefore may not
	# include representatives for all components across all size bins and 
	# wall, therefore extend to
	# include all components, including water and seed material.
	# Remember that any final element of colptrs should represent the final
	# rowval index, so colptrs needs to have a length one greater than the
	# number of components+number of components*number of size bins (number of size bins includes wall)
	col_shrt =  ((comp_num+2)+(comp_num+2)*num_sb+1)-len(colptrs)
	# rowvals indices for final columns of Jacobian		
	colptrs = np.append(colptrs, np.array((colptrs[-1])).repeat(col_shrt))
	
	if (self.eqn_num[1] > 0): # if aqueous-phase reactions present

		# aqueous-phase reaction part -----------------------------------------------------
		# container for Jacobian index affected by aqueous-phase reactions	
		self.jac_indx_aq = np.zeros((self.eqn_num[1], max(self.nreac_aq*(self.nreac_aq+self.nprod_aq))))
		# index of unique rows affected per reactant, if zero acts as a filler this omitted below
		col_tracker = np.unique(self.rindx_aq)
		# track rows affected per column
		row_tracker = np.ones((max(self.nreac_aq+self.nprod_aq), len(col_tracker)))*-1.
		# note that size bins are dealt with below loops
		for eqni in range(self.eqn_num[1]): # loop through reactions
			# total number of components affected per reactant
			tot_affcomp = self.nreac_aq[eqni]+self.nprod_aq[eqni]
			# combined index of reactants and products in this equation
			totindx = np.append(self.rindx_aq[eqni, 0:self.nreac_aq[eqni]], self.pindx_aq[eqni, 0:self.nprod_aq[eqni]])
			for ir in range(self.nreac_aq[eqni]):# loop through reactants (would be columns in 2D Jacobian)
				# number of full Jacobian elements passed before reaching this reactant, 
				# accounting for water and core
				st_indx = (comp_num+2)*(num_sb+1)*self.rindx_aq[eqni, ir]
				# flat Jacobian index affected by this equation, this then reduced to sparse 
				# matrix index below
				self.jac_indx_aq[eqni, ir*tot_affcomp:(ir+1)*tot_affcomp] = st_indx+totindx
			
				# the column number affected
				coli = np.where(col_tracker == self.rindx_aq[eqni, ir])[0]
				# unique rows affected per reactant
				unirow = np.unique(np.append(row_tracker[:, coli][row_tracker[:, coli]!=-1], totindx))
				row_tracker[0:len(unirow), coli] = unirow.reshape(-1, 1)
		
		# account for the the gas-phase Jacobian indices
		row_tracker[row_tracker != -1] += comp_num+2
		# flatten ready for appendage
		row_trackerf = row_tracker.flatten(order='F')
		# append to rowvals
		rowvals = np.append(rowvals, row_trackerf[row_trackerf != -1])
		# sum number of unique rows affected per component
		col_num = (row_tracker != -1).sum(axis = 0)
		# account for new rows in first size bin, add 1 to index because this represents the
		# number of elements per component
		colptrs[col_tracker+(comp_num+2)+1] += np.cumsum(col_num)
		# ensure latter inidices are consistent 
		colptrs[max(col_tracker+(comp_num+2)+1)::] = colptrs[max(col_tracker+(comp_num+2)+1)] 
		
		# because above jac_indx_aq contains the indices of data assuming
		# a full, rather than sparse matrix, now correct to affect
		# the sparse matrix by reducing indices where gaps occur in
		# the complete data matrix
		# get unique values, mask out elements equal to zero as these are 
		# unfilled from the original container
		uni = np.unique(self.jac_indx_aq[self.jac_indx_aq != 0.])
		# now reduce these unique indices based on the number of indices
		# omitted
		diff = 0 # prepare for number of empty points
		for i in range(1, len(uni)):
			diff += (uni[i]-uni[i-1])-1
			self.jac_indx_aq[self.jac_indx_aq==uni[i]] -= diff
		
		# now account for omission of full Jacobian elements between the gas-phase 
		# reaction part and the aqueous-phase reaction part
		self.jac_indx_aq[self.jac_indx_aq != 0.] -= np.amin(self.jac_indx_aq[self.jac_indx_aq != 0.])-np.amax(self.jac_indx_g)-1
		# range of aqueous-phase sparse jacobian indices
		ran_aq = np.amax(self.jac_indx_aq[self.jac_indx_aq != 0.])-np.amin(self.jac_indx_aq[self.jac_indx_aq != 0.])+1
		# shape of index matrix for just one size bin
		jsh = self.jac_indx_aq.shape

		# now repeat over particle size bins
		for sbi in range(1, num_asb):
			self.jac_indx_aq = np.concatenate((self.jac_indx_aq, self.jac_indx_aq[0:jsh[0], 0:jsh[1]]+self.ran_aq*(sbi)), axis = 0)
			
			# account for the the gas-phase Jacobian indices
			row_tracker[row_tracker != -1] += (comp_num+2)
			# flatten ready for appendage
			row_trackerf = row_tracker.flatten(order='F')
			# append to rowvals
			rowvals = np.append(rowvals, row_trackerf[row_trackerf != -1])
			# account for new rows in latter size bin, add 1 to index because this represents the
			# number of elements per component
			colptrs[col_tracker+(comp_num+2)*(sbi+1)+1] += np.cumsum(col_num)
			# ensure latter inidices are consistent 
			colptrs[max(col_tracker+(comp_num+2)*(sbi+1)+1)::] = colptrs[max(col_tracker+(comp_num+2)*(sbi+1)+1)]
	
	if (self.eqn_num[1] == 0):
		self.jac_indx_aq = np.zeros((1)) # filler


	if (self.eqn_num[2] > 0): # if surface (e.g. wall) reactions present

		# surface (e.g. wall) reaction part -----------------------------------------------------
		# container for Jacobian index affected by surface (e.g. wall) reactions	
		jac_indx_su = np.zeros((self.eqn_num[2], max(self.nreac_su*(self.nreac_su+self.nprod_su))))
		# index of unique rows affected per reactant, if zero acts as a filler this omitted below
		col_tracker = np.unique(self.rindx_su)
		# track rows affected per column
		row_tracker = np.ones((max(self.nreac_su+self.nprod_su), len(col_tracker)))*-1.
		# note that size bins are dealt with below loops
		for eqni in range(self.eqn_num[2]): # loop through reactions
			# total number of components affected per reactant
			tot_affcomp = self.nreac_su[eqni]+self.nprod_su[eqni]
			# combined index of reactants and products in this equation
			totindx = np.append(self.rindx_su[eqni, 0:self.nreac_su[eqni]], self.pindx_su[eqni, 0:self.nprod_su[eqni]])

			for ir in range(self.nreac_su[eqni]): # loop through reactants (would be columns in 2D Jacobian)
				# number of full Jacobian (i.e. not sparse Jacobian matrix) elements passed 
				# before reaching this reactant, 
				# accounting for water and core
				st_indx = (comp_num+2)*(num_sb+1)*self.rindx_su[eqni, ir]
				# flat Jacobian index affected by this equation, this then reduced to sparse 
				# matrix index below
				self.jac_indx_su[eqni, ir*tot_affcomp:(ir+1)*tot_affcomp] = st_indx+totindx
				# the column number affected
				coli = np.where(col_tracker == self.rindx_su[eqni, ir])[0]
				# unique rows affected per reactant
				unirow = np.unique(np.append(row_tracker[:, coli][row_tracker[:, coli]!=-1], totindx))
				row_tracker[0:len(unirow), coli] = unirow.reshape(-1, 1)

		# account for the the gas-phase and particle Jacobian indices
		row_tracker[row_tracker != -1] += (comp_num+2) + (comp_num+2)*num_asb
		# flatten ready for appendage
		row_trackerf = row_tracker.flatten(order='F')
		# append to rowvals
		rowvals = np.append(rowvals, row_trackerf[row_trackerf != -1])
		# sum number of unique rows affected per component
		col_num = (row_tracker != -1).sum(axis = 0)
		# account for new rows in first size bin, add 1 to index because this represents the
		# number of elements per component
		colptrs[col_tracker+((comp_num+2) + (comp_num+2)*num_asb)+1] += np.cumsum(col_num)
		# ensure latter inidices are consistent 
		colptrs[max(col_tracker+((comp_num+2) + (comp_num+2)*num_asb)+1)::] = colptrs[max(col_tracker+((comp_num+2) + (comp_num+2)*num_asb)+1)] 

		# because above jac_indx_su contains the indices of data assuming
		# a full, rather than sparse matrix, now correct to affect
		# the sparse matrix by reducing indices where gaps occur in
		# the complete data matrix
		# get unique values, mask out elements equal to zero as these are 
		# unfilled from the original container
		uni = np.unique(self.jac_indx_su[self.jac_indx_su != 0.])
		# now reduce these unique indices based on the number of indices
		# omitted
		diff = 0 # prepare for number of empty points
		for i in range(1, len(uni)):
			diff += (uni[i]-uni[i-1])-1
			self.jac_indx_su[self.jac_indx_su==uni[i]] -= diff
		
		# now account for omission of full Jacobian elements between the 
		# aqueous-phase reaction part (if aqueous-phase reactions present) or 
		# the gas-phase reaction part (if aqueous-phase reactions not present) 
		# and the surface (e.g. wall) reaction part
		self.jac_indx_su[self.jac_indx_su != 0.] -= np.amin(self.jac_indx_su[self.jac_indx_su != 0.])-max(np.amax(self.jac_indx_g), np.amax(self.jac_indx_aq))-1

		# range of surface (e.g. wall) sparse jacobian indices
		ran_su = np.amax(self.jac_indx_su[self.jac_indx_su != 0.])-np.amin(self.jac_indx_su[self.jac_indx_su != 0.])+1
		# shape of index matrix for just one size bin
		jsh = self.jac_indx_su.shape

		# now repeat over surface bins
		for sbi in range(1, self.wall_on):
			self.jac_indx_su = np.concatenate((self.jac_indx_su, self.jac_indx_su[0:jsh[0], 0:jsh[1]]+ran_su*(sbi)), axis = 0)
			
			# account for the the gas-phase and particle-phase Jacobian indices
			row_tracker[row_tracker != -1] += (comp_num+2)+((comp_num+2)*num_asb)
			# flatten ready for appendage
			row_trackerf = row_tracker.flatten(order='F')
			# append to rowvals
			rowvals = np.append(rowvals, row_trackerf[row_trackerf != -1])
			# account for new rows in latter size bin, add 1 to index because this represents the
			# number of elements per component
			colptrs[col_tracker+(comp_num+2)*(num_asb+sbi+1)+1] += np.cumsum(col_num)
			# ensure latter indices are consistent 
			colptrs[max(col_tracker+(comp_num+2)*(sbi+1)+1)::] = colptrs[max(col_tracker+(comp_num+2)*(sbi+1)+1)]
	
	if (self.eqn_num[2] == 0):
		self.jac_indx_su = np.zeros((1)) # filler
	
	# particle partitioning influence on Jacobian part -------------------------------------------------
	# loop through jacobian index to check whether the centre 
	# diagonal for gas effect on gas components is due to be 
	# affected by gas-phase reactions or whether particle 
	# partitioning is the first effect.  Similarly accommodate 
	# for particle effect for: gas on particle effect, particle on gas effect
	# and particle on particle effect

	# Jacobian index for particle effects
	jac_part_indx = np.zeros((comp_num+2)*(num_asb+1)+((comp_num+2)*(num_asb*2)))
	
	part_cnt = 0 # count on jac_part_indx inputs
	
	if (num_asb > 0): # if particle size bins are present

		# loop through components in the gas-phase (add two 
		# to account for water and seed material)
		for compi in range(comp_num+2):
			# gas effect on gas part ---------------------------------------------
			# relevant starting and finishing index in rowvals
			st_indx = int(colptrs[compi])
			en_indx = int(colptrs[compi+1])
			if ((st_indx < en_indx) == True): # rows already attributed to this column
				# if rows are already attributed, check 
				# whether the diagonal is already affected
				if (any(rowvals[st_indx:en_indx] == compi)):
					# get index
					exist_indx = st_indx+(np.where(rowvals[st_indx:en_indx] == compi)[0][0])
					# if diagonal already affected then just copy relevant
					# data index to the indexing for wall
					jac_part_indx[part_cnt] = exist_indx
					
				else: # if diagonal not already affected, then include
					# get index
					new_indx = st_indx+(sum(rowvals[st_indx:en_indx]<compi))+1
					# modify indices for sparse Jacobian matrix
					jac_part_indx[part_cnt] = new_indx # gas-particle partitioning
					self.jac_indx_g[self.jac_indx_g>=new_indx] += 1 # gas-phase reactions
					self.jac_indx_aq[self.jac_indx_aq>=new_indx] += 1 # aqueous-phase reactions
					self.jac_indx_su[self.jac_indx_su>=new_indx] += 1 # surface (e.g. wall) reactions
					new_el = np.array((compi)).reshape(1) 
					rowvals = np.concatenate([rowvals[0:new_indx], 
						new_el, rowvals[new_indx::]])
					colptrs[compi+1::] += 1
					# increase en_indx to account for new row when
					# accommodating gas effect on wall below
					en_indx += 1
					
			else: # no rows yet attributed to this column, so need to include

				# modify indices for sparse Jacobian matrix
				jac_part_indx[part_cnt] = st_indx # gas-particle partitioning
				self.jac_indx_g[self.jac_indx_g>=st_indx] += 1 # gas-phase reactions
				self.jac_indx_aq[self.jac_indx_aq>=st_indx] += 1 # aqueous-phase reactions
				self.jac_indx_su[self.jac_indx_su>=st_indx] += 1 # surface (e.g. wall) reactions
				new_el = np.array((compi)).reshape(1)
				rowvals = np.concatenate([rowvals[0:st_indx], new_el, rowvals[st_indx::]])
				colptrs[compi+1::] += 1
				# increase en_indx to account for new row when accommodating
				# gas effect on particle below
				en_indx += 1
				
			part_cnt += 1 # keep count on particle index
			
			# gas effect on particle part -------------------------------
			# already know that the gas-phase reactions will not
			# have affected this part of the Jacobian, so just need to
			# modify sparse Jacobian inputs accordingly
			jac_part_indx[part_cnt:part_cnt+num_asb] = range(en_indx, en_indx+num_asb)
			self.jac_indx_g[self.jac_indx_g>=en_indx] += num_asb # gas-phase reaction
			self.jac_indx_aq[self.jac_indx_aq>=en_indx] += num_asb # aqueous-phase reactions
			self.jac_indx_su[self.jac_indx_su>=en_indx] += num_asb # surface (e.g. wall) reactions
			new_el = (np.array(range(comp_num+2+compi, (comp_num+2)*(num_asb+1)+compi, (comp_num+2)))).reshape(-1)
			rowvals = np.concatenate([rowvals[0:en_indx], new_el, rowvals[en_indx::]])
			colptrs[compi+1::] += num_asb
			part_cnt += num_asb # keep count on particle index		
		
		
		# loop through components to include their particle effect on
		# gas phase and on particle components for the Jacobian, note if no wall this 
		# will include the final row in the final column of the Jacobian
		for sbi in range(num_asb): # size bin loop
			for compi in range(comp_num+2): # component loop
				# particle component effect on gas-phase
				
				# the flattened Jacobian index relating to this
				# relevant start index for colptrs
				stc_indx = int((comp_num+2)*(sbi+1)+(compi))
				st_indx = int(colptrs[stc_indx]) # start index for rowvals		
				jac_part_indx[part_cnt] = st_indx
				new_el = np.array((compi)).reshape(-1)
				rowvals = np.concatenate([rowvals[0:st_indx], new_el, rowvals[st_indx::]])
				colptrs[(comp_num+2)*(sbi+1)+(compi+1)::] += 1
				self.jac_indx_aq[self.jac_indx_aq>=st_indx] += 1 # aqueous-phase reactions
				self.jac_indx_su[self.jac_indx_su>=st_indx] += 1 # surface (e.g. wall) reactions
				part_cnt += 1
				
				# particle component effect on particle
				new_el = np.array(((comp_num+2)*(sbi+1)+compi)).reshape(-1)
				# account for any aqueous-phase or surface (e.g. wall) reactions
				st_indx = (st_indx) + sum(rowvals[int(colptrs[stc_indx]):int(colptrs[stc_indx+1])] < new_el)
				jac_part_indx[part_cnt] = st_indx
				# check if diagonal already in use by aqueous-phase reactions
				if any(rowvals[int(colptrs[stc_indx]):int(colptrs[stc_indx+1])] == new_el):
					part_cnt += 1
					continue # continue to next component if it is
				else:
					rowvals = np.concatenate([rowvals[0:st_indx], new_el, rowvals[st_indx::]])
					colptrs[(comp_num+2)*(sbi+1)+compi+1::] += 1
					self.jac_indx_aq[self.jac_indx_aq>=st_indx] += 1 # aqueous-phase reactions
					self.jac_indx_su[self.jac_indx_su>=st_indx] += 1 # surface (e.g. wall) reactions
					part_cnt += 1
	
	# wall partitioning influence on Jacobian part ---------------------------------------------------

	# loop through jacobian index to check whether the centre 
	# diagonal for gas effect on gas components is due to be 
	# affected by gas-phase reactions or whether wall 
	# partitioning is the first effect.  Similarly accommodate 
	# for wall effect for: gas-on-wall effect, wall-on-gas effect
	# and wall-on-wall affect
	
	# empty results array, 1st term is gas-on-gas, 2nd is gas-on-wall, 3rd is wall-on-gas, 4th is wall-on-wall
	jac_wall_indx = np.zeros(((comp_num+2)+(self.wall_on*(comp_num+2))+(self.wall_on*(comp_num+2))+(self.wall_on*(comp_num+2))))
	wall_cnt = 0 # count on jac_wall_indx inputs
	
	if (self.wall_on > 0): # if surface present
		# loop through components in the gas-phase (add two 
		# to account for water and core component)
		for compi in range(comp_num+2):
		
			# gas effect on gas part --------------------------------------------
			# relevant starting and finishing index in rowvals
			st_indx = int(colptrs[compi])
			en_indx = int(colptrs[compi+1])
			
			# check if any rows already attributed to this column
			if ((st_indx < en_indx) == True):
				# if rows are already attributed, check 
				# whether the diagonal is already affected
				if (sum(rowvals[st_indx:en_indx]==compi))>0:
					# get index
					exist_indx = st_indx+(np.where(rowvals[st_indx:en_indx]==compi)[0][0])
					# if diagonal already affected then just copy relevant
					# data index to the indexing for wall
					jac_wall_indx[wall_cnt] = exist_indx
					
				else: # if diagonal not already affected, then add
					# get index
					new_indx = st_indx+(sum(rowvals[st_indx:en_indx]<compi))+1
					# modify indices for sparse Jacobian matrix
					jac_wall_indx[wall_cnt] = new_indx
					self.jac_indx_g[self.jac_indx_g >= new_indx] += 1
					self.jac_indx_aq[self.jac_indx_aq >= new_indx] += 1
					self.jac_indx_su[self.jac_indx_su >= new_indx] += 1
					new_el = np.array((compi)).reshape(1)
					rowvals = np.concatenate([rowvals[0:new_indx], new_el, rowvals[new_indx::]])
					colptrs[compi+1::] += 1
					# increase en_indx to account for new row when
					# accommodating gas effect on wall below
					en_indx += 1
					
			else: # no rows yet attributed to this column, so need to include

				# modify indices for sparse Jacobian matrix
				jac_wall_indx[wall_cnt] = st_indx
				self.jac_indx_g[self.jac_indx_g >= st_indx] += 1
				self.jac_indx_aq[self.jac_indx_aq >= st_indx] += 1
				self.jac_indx_su[self.jac_indx_su >= st_indx] += 1
				new_el = np.array((compi)).reshape(1)
				rowvals = np.concatenate([rowvals[0:st_indx], new_el, rowvals[st_indx::]])
				colptrs[compi+1::] += 1
				# increase en_indx to account for new row when accommodating
				# gas effect on wall below
				en_indx += 1
				
			wall_cnt += 1 # keep count on wall index
			
			# gas effect on wall part -----------------------------------------
			# already know that the gas-phase photochemistry will not
			# have affected this part of the Jacobian, so just need to
			# modify sparse Jacobian inputs accordingly
			jac_wall_indx[wall_cnt:wall_cnt+self.wall_on] = range(en_indx, en_indx+self.wall_on)
			# increase gas, aqueous and particle (if present) effect indices up by number of surface bins to 
			# account for addition to Jacobian
			self.jac_indx_g[self.jac_indx_g >= en_indx] += self.wall_on
			self.jac_indx_aq[self.jac_indx_aq >= en_indx] += self.wall_on
			self.jac_indx_su[self.jac_indx_su >= en_indx] += self.wall_on
			
			if (num_asb > 0):
				jac_part_indx[jac_part_indx >= en_indx] += self.wall_on
				
			new_el = (np.array(range((comp_num+2)*(num_asb+1)+compi, (comp_num+2)*(num_asb+1)+(comp_num+2)*self.wall_on+compi, comp_num+2))).reshape(-1)
			rowvals = np.concatenate([rowvals[0:en_indx], new_el, rowvals[en_indx::]])
			colptrs[compi+1::] += self.wall_on
			wall_cnt += self.wall_on # keep count on wall index

		# wall effect on gas phase and on wall components for the Jacobian, note this will include the
		# final row in the final column of the Jacobian
		for wbi in range(self.wall_on): # wall bin loop
			for compi in range(comp_num+2): # component loop
				# wall component effect on itself in gas-phase is a location in the 
				# Jacobian that cannot be taken by surface reactions
				
				# the flattened Jacobian index relating to this
				# relevant start index for colptrs
				stc_indx = int((comp_num+2)*(wbi+num_asb+1)+(compi))
				st_indx = colptrs[stc_indx]

				jac_wall_indx[wall_cnt] = st_indx # include in wall index for Jacobian
				# include in rowvals
				rowvals = np.concatenate([rowvals[0:int(colptrs[stc_indx])], [compi], rowvals[int(colptrs[stc_indx])::]])
				# include in colptrs
				colptrs[stc_indx+1::] += 1
				self.jac_indx_su[self.jac_indx_su >= st_indx] += 1 # adjust surface (e.g. wall) reactions
				wall_cnt += 1 # move up index for wall index for Jacobian
			
				# for wall component effect on itself in wall, this is a space in the 
				# Jacobian that surface reactions might have taken, so check whether 
				# any surface (e.g. wall reactions) affect this component - note the +1 added to
				# colptrs[stc_indx] on the LHS of the equality to account for the +1 addition to colptrs
				# just above in this loop
				if ((colptrs[stc_indx]+1) == (colptrs[stc_indx+1])): # surface reaction is not present
					
					jac_wall_indx[wall_cnt] = st_indx+1 # include in wall index for Jacobian
					# include in rowvals, note that using stc_indx as the value means that 
					# we use the index for component on wall effect on itself on wall - see 
					# the stc_indx definition further up in this loop
					rowvals = np.concatenate([rowvals[0:int(colptrs[stc_indx]+1)], [stc_indx], rowvals[int(colptrs[stc_indx]+1)::]])
					
					# include in colptrs
					colptrs[stc_indx+1::] += 1
					self.jac_indx_su[self.jac_indx_su >= st_indx+1] += 1 # surface (e.g. wall) reactions
				
				else: # surface reaction is present
					# get the rows affected by surface reaction in this column, note the 
					# [1::] index removes the row affected by wall effect on gas,
					# which is added above (inside this loop)
					rows_reac = rowvals[int(colptrs[stc_indx]):int(colptrs[stc_indx+1])][1::]	
					
					# sum Jacobian indices between the wall on gas effect and this wall on wall effect
					wall_indx_add = sum(rows_reac < comp_num+2)
					# get the wall index
					jac_wall_indx[wall_cnt] = st_indx+1+wall_indx_add # include in wall index for Jacobian

					if (sum(rows_reac == stc_indx) == 0.): # if reaction does not cover wall effect on wall
						
						# include in rowvals
						rowvals = np.concatenate([rowvals[0:int(colptrs[stc_indx]+1)], [stc_indx], rowvals[int(colptrs[stc_indx]+1)::]])
						# include in colptrs
						colptrs[stc_indx+1::] += 1
						self.jac_indx_su[self.jac_indx_su > jac_wall_indx[wall_cnt]] += 1 # adjust surface (e.g. wall) reactions
										

				wall_cnt += 1 # move up index for wall index for Jacobian
				
	# end of wall influence on Jacobian part ---------------------------------------------------

	# index of the Jacobian affected by air extraction
	jac_extr_indx = np.zeros(((comp_num+2)*(num_sb+1)))
	
	# extraction effect on Jacobian part -----------------------
	if (len(self.dil_fac) > 0): # if chamber air continuously being extracted, e.g. in flow-reactor
	
		extr_cnt = 0 # count on jac_extr_indx inputs
		
		# loop through all gas- and particle-phase components
		for compi in range((comp_num+2)*(num_asb)):
		
			# gas effect on gas part --------------------------------------------
			# relevant starting and finishing index in rowvals
			st_indx = int(colptrs[compi])
			en_indx = int(colptrs[compi+1])
			
			# check whether the diagonal of Jacobian for this component is already affected
			# relevant starting and finishing index in rowvals
			st_indx = int(colptrs[compi])
			en_indx = int(colptrs[compi+1])
			
			# check if any rows already attributed to this column
			if ((st_indx < en_indx) == True):
			
				# if rows are already attributed, check 
				# whether the diagonal is already affected
				if ((sum(rowvals[st_indx:en_indx]==compi)) > 0):
				
					# get index
					exist_indx = st_indx+(np.where(rowvals[st_indx:en_indx]==compi)[0][0])
					# if diagonal already affected then just copy relevant
					# data index to the indexing for wall
					jac_extr_indx[extr_cnt] = exist_indx
					
				else: # if diagonal not already affected, then add
				
					# get index
					new_indx = st_indx+(sum(rowvals[st_indx:en_indx]<compi))+1
					# modify indices for sparse Jacobian matrix
					jac_extr_indx[extr_cnt] = new_indx
					self.jac_indx_g[self.jac_indx_g >= new_indx] += 1
					self.jac_indx_aq[self.jac_indx_aq >= new_indx] += 1
					self.jac_indx_su[self.jac_indx_su >= new_indx] += 1
					if (num_asb > 0):
						jac_part_indx[jac_part_indx >= new_indx] += 1
					jac_wall_indx[jac_wall_indx >= new_indx] += 1
					
					new_el = np.array((compi)).reshape(1)
					rowvals = np.concatenate([rowvals[0:new_indx], 
						new_el, rowvals[new_indx::]])
					colptrs[compi+1::] += 1
					
			else: # no rows yet attributed to this column, so need to include

				# modify indices for sparse Jacobian matrix
				jac_extr_indx[extr_cnt] = st_indx
				# adjust other indices
				self.jac_indx_g[self.jac_indx_g >= st_indx] += 1
				self.jac_indx_aq[self.jac_indx_aq >= st_indx] += 1
				self.jac_indx_su[self.jac_indx_su >= st_indx] += 1
				if (num_asb > 0):
					jac_part_indx[jac_part_indx >= st_indx] += 1
				jac_wall_indx[jac_wall_indx >= st_indx] += 1
				
				new_el = np.array((compi)).reshape(1)
				rowvals = np.concatenate([rowvals[0:st_indx], new_el, rowvals[st_indx::]])
				colptrs[compi+1::] += 1
				
			extr_cnt += 1 # keep count on wall index
	
	# end of extraction effect on Jacobian part -------------


	# index of the Jacobian affected by continuous influx
	self.jac_cont_infl_indx = np.zeros((comp_num+2))
	
	# continuous influx effect on Jacobian part -----------------------
	if (len(self.con_infl_indx) > 0): # if component continuously injected
	
		extr_cnt = 0 # count on jac_cont_infl_indx inputs
		
		# loop through all gas-phase components
		for compi in range(comp_num+2):
		
			# gas effect on gas part --------------------------------------------
			# relevant starting and finishing index in rowvals
			st_indx = int(colptrs[compi])
			en_indx = int(colptrs[compi+1])
			
			# check whether the diagonal of Jacobian for this component is already affected
			# relevant starting and finishing index in rowvals
			st_indx = int(colptrs[compi])
			en_indx = int(colptrs[compi+1])
			
			# check if any rows already attributed to this column
			if ((st_indx < en_indx) == True):
			
				# if rows are already attributed, check 
				# whether the diagonal is already affected
				if ((sum(rowvals[st_indx:en_indx]==compi)) > 0):
				
					# get index
					exist_indx = st_indx+(np.where(rowvals[st_indx:en_indx]==compi)[0][0])
					# if diagonal already affected then just copy relevant
					# data index to the indexing for wall
					self.jac_cont_infl_indx[extr_cnt] = exist_indx
					
				else: # if diagonal not already affected, then add
				
					# get index
					new_indx = st_indx+(sum(rowvals[st_indx:en_indx]<compi))+1
					# modify indices for sparse Jacobian matrix
					self.jac_cont_infl_indx[extr_cnt] = new_indx
					self.jac_indx_g[self.jac_indx_g >= new_indx] += 1
					self.jac_indx_aq[self.jac_indx_aq >= new_indx] += 1
					self.jac_indx_su[self.jac_indx_su >= new_indx] += 1
					if (num_asb > 0):
						jac_part_indx[jac_part_indx >= new_indx] += 1
					jac_wall_indx[jac_wall_indx >= new_indx] += 1
					
					new_el = np.array((compi)).reshape(1)
					rowvals = np.concatenate([rowvals[0:new_indx], 
						new_el, rowvals[new_indx::]])
					colptrs[compi+1::] += 1
					
			else: # no rows yet attributed to this column, so need to include

				# modify indices for sparse Jacobian matrix
				self.jac_cont_infl_indx[extr_cnt] = st_indx
				# adjust other indices
				self.jac_indx_g[self.jac_indx_g >= st_indx] += 1
				self.jac_indx_aq[self.jac_indx_aq >= st_indx] += 1
				self.jac_indx_su[self.jac_indx_su >= st_indx] += 1
				if (num_asb > 0):
					jac_part_indx[jac_part_indx >= st_indx] += 1
				jac_wall_indx[jac_wall_indx >= st_indx] += 1
				
				new_el = np.array((compi)).reshape(1)
				rowvals = np.concatenate([rowvals[0:st_indx], new_el, rowvals[st_indx::]])
				colptrs[compi+1::] += 1
				
			extr_cnt += 1 # keep count on wall index

		self.jac_cont_infl_indx = self.jac_cont_infl_indx.astype('int') # ensure integer type

	# end of continuous influx effect on Jacobian part -------------
	
	if ((num_sb == 0) and len(rowvals) >= 1): # if no particle size bins and no wall

		# if the Jacobian matrix has an empty final row, then
		# an error will be displayed during ODE solver call, so 
		# input a filler on this row;
		# add two to the number of unique components counted here to
		# account for water and seed material
		if (max(rowvals) <= (comp_num+2)-1):
			# index for final row of Jacobian
			rowvals = np.append(rowvals, (comp_num+2)-1)
			# number of columns currently short of final column
			col_shrt =  comp_num+2-len(colptrs)
			# rowval index for final column of Jacobian		
			colptrs[-1] += 1

	# ensure integer type for arrays acting as indices
	self.jac_indx_g = self.jac_indx_g.astype(int)
	self.jac_indx_aq = self.jac_indx_aq.astype(int)
	self.jac_indx_su = self.jac_indx_su.astype(int)
	jac_part_indx = jac_part_indx.astype(int)
	jac_wall_indx = jac_wall_indx.astype(int)	
	jac_extr_indx = jac_extr_indx.astype(int)
	rowvals = rowvals.astype(int)	
	colptrs = colptrs.astype(int)
	
	return(rowvals, colptrs, jac_part_indx, jac_wall_indx, jac_extr_indx, self)
