'''preparing the inputs for the ode solver Jacobian'''
# called once the gas-phase and particle-phase equations have been interrogated

# required modules
import numpy as np

def jac_setup(jac_den_indx, njac, comp_num, num_sb, num_eqn, nreac_g, nprod_g, rindx_g, pindx_g, jac_indx_g, wall_on, nreac_aq, nprod_aq, rindx_aq, pindx_aq, jac_indx_aq, num_asb):

	# inputs ---------------------------------------------
	# jac_den_indx - index of denominators for jacobian
	# njac - number of jacobian elements affected per equation
	# comp_num - number of components
	# num_sb - number of size bins (including any wall)
	# num_eqn - number of chemical reactions
	# nreac_g - number of reactants per reaction
	# nprod_g - number of products per reaction
	# rindx_g - index of reactants per equation
	# pindx_g - index of products per equation
	# jac_indx_g - index of Jacobian affected per equation
	# wall_on - flag for whether wall on or off
	# nreac_aq - number of reactants per equation
	# nprod_aq - number of products per equation
	# rindx_aq - index of reactants per equation
	# pindx_aq - index of prodcuts per equation
	# jac_indx_aq - index of jacobian affected per equation
	# num_asb - number of actual particle size bins (excluding wall) 
	# ----------------------------------------------------

	print('Preparing Jacobian inputs')

	# rows in Jacobian affected by this equation		 	
	rowvals = np.empty((0))
	# indices of rowvals representing each component being 
	# differentiated by in Jacobian, note, needs to start
	# as zero to represent first index of rowval
	colptrs = np.zeros((1))

	# ensure integer type
	jac_den_indx = jac_den_indx.astype(int)
	njac = njac.astype(int)
	
	# track columns of Jacobian affected
	col_tr = 0
	# length of Jacobian; add 1 to account for the gas phase
	len_jac = comp_num*(num_sb+1)

	# index for Jacobian - note can only be done after all equations checked as only then
	# is the total number of unique components known
	for eqni in range(num_eqn[0]): # equation loop
	
		# total number of components in this equation
		tot_comp = nreac_g[eqni]+nprod_g[eqni]
		# combined index of reactants and products in this equation
		totindx = np.append(rindx_g[eqni, 0:nreac_g[eqni]], pindx_g[eqni, 0:nprod_g[eqni]])
		
		# reactant loop (equivalent to columns of 2D Jacobian)
		for i in range(nreac_g[eqni]):
			# index of flattened Jacobian being affected in 
			# this equation, note this doesn't need to 
			# account for water and core component as below (when unique elements found) all
			# full Jacobian matrix elements that aren't indexed are
			# omitted to create the sparse Jacobian matrix  
			jac_indx_g[eqni, i*tot_comp:(i+1)*tot_comp] = rindx_g[eqni, i]*len_jac+(totindx)
			# check if rowvals already sufficiently long to contain this reactant  
			if rindx_g[eqni, i]<=(len(colptrs)-2):

				# start index of rowvals where it is represented
				r_st = int(colptrs[rindx_g[eqni, i]])
				# end index of rowvals where it is represented
				r_en = int(colptrs[rindx_g[eqni, i]+1])
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
						r_in = r_st+int(sum(rowvals[r_st:r_en]<new_el))
						rowvals = np.concatenate([rowvals[0:r_in], new_el, rowvals[r_in::]])
						colptrs[rindx_g[eqni, i]+1::] += 1
						r_en += 1
						
						
			else: # if not accounted for in previous reaction
				# rows in Jacobian affected by this reactant 
				# in this equation	 	
				rowvals = np.append(rowvals, np.unique(totindx))
				while (rindx_g[eqni, i]>col_tr):
					colptrs = np.append(colptrs, colptrs[-1])
					col_tr += 1
				# indices of rowvals
				colptrs = np.append(colptrs, colptrs[-1]+len(np.unique(totindx)))
				# track number of columns of Jacobian affected				
				col_tr += 1

	# because above jac_indx_g contains the indices of data assuming
	# a full, rather than sparse matrix, now correct to affect
	# the sparse matrix by reducing indices where gaps occur in
	# the complete data matrix
	# get unique values in jac_indx_g
	uni = np.unique(jac_indx_g)
	# now reduce these indices based on the number of indices
	# omitted
	diff = 0 # prepare for number of empty points
	for i in range(1, len(uni)):
		diff += (uni[i]-uni[i-1])-1
		jac_indx_g[jac_indx_g==uni[i]] -= diff
	
	# colptrs responds to appendages to rowvals above and therefore may not
	# include representatives for all components across all size bins and 
	# wall, therefore extend to
	# include all components, including water and seed material.
	# Remember that any final element of colptrs should represent the final
	# rowval index, so colptrs needs to have a length one greater than the
	# number of components+number of components*number of size bins
	col_shrt =  ((comp_num+2)+(comp_num+2)*num_sb+1)-len(colptrs)
	# rowvals indices for final columns of Jacobian		
	colptrs = np.append(colptrs, np.array((colptrs[-1])).repeat(col_shrt))
	
	if (num_eqn[1] > 0): # if aqueous-phase reactions present

		# aqueous-phase reaction part -----------------------------------------------------
		# container for Jacobian index affected by aqueous-phase reactions	
		jac_indx_aq = np.zeros((num_eqn[1], max(nreac_aq*(nreac_aq+nprod_aq))))
		# index of unique rows affected per reactant, if zero acts as a filler this omitted below
		col_tracker = np.unique(rindx_aq)
		# track rows affected per column
		row_tracker = np.ones((max(nreac_aq+nprod_aq), len(col_tracker)))*-1.
		# note that size bins are dealt with below loops
		for eqni in range(num_eqn[1]): # loop through reactions
			# total number of components affected per reactant
			tot_affcomp = nreac_aq[eqni]+nprod_aq[eqni]
			# combined index of reactants and products in this equation
			totindx = np.append(rindx_aq[eqni, 0:nreac_aq[eqni]], pindx_aq[eqni, 0:nprod_aq[eqni]])
			for ir in range(nreac_aq[eqni]):# loop through reactants (would be columns in 2D Jacobian)
				# number of full Jacobian elements passed before reaching this reactant, 
				# accounting for water and core
				st_indx = (comp_num+2)*(num_sb+1)*rindx_aq[eqni, ir]
				# flat Jacobian index affected by this equation, this then reduced to sparse 
				# matrix index below
				jac_indx_aq[eqni, ir*tot_affcomp:(ir+1)*tot_affcomp] = st_indx+totindx
			
				# the column number affected
				coli = np.where(col_tracker == rindx_aq[eqni, ir])[0]
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
		uni = np.unique(jac_indx_aq[jac_indx_aq != 0.])
		# now reduce these unique indices based on the number of indices
		# omitted
		diff = 0 # prepare for number of empty points
		for i in range(1, len(uni)):
			diff += (uni[i]-uni[i-1])-1
			jac_indx_aq[jac_indx_aq==uni[i]] -= diff
		
		# now account for omission of full Jacobian elements between the gas-phase 
		# reaction part and the aqueous-phase reaction part
		jac_indx_aq[jac_indx_aq != 0.] -= np.amin(jac_indx_aq[jac_indx_aq != 0.])-np.amax(jac_indx_g)-1
		# range of aqueous-phase sparse jacobian indices
		ran_aq = np.amax(jac_indx_aq[jac_indx_aq != 0.])-np.amin(jac_indx_aq[jac_indx_aq != 0.])+1
		# shape of index matrix for just one size bin
		jsh = jac_indx_aq.shape
		# now repeat over particle size bins
		for sbi in range(1, num_asb):
			jac_indx_aq = np.concatenate((jac_indx_aq, jac_indx_aq[0:jsh[0], 0:jsh[1]]+ran_aq*(sbi)), axis = 0)
			
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
	
	if (num_eqn[1] == 0):
		jac_indx_aq = np.zeros((1)) # filler

	# particle influence on Jacobian part -------------------------------------------------
	# loop through jacobian index to check whether the centre 
	# diagonal for gas effect on gas components is due to be 
	# affected by gas-phase reactions or whether particle 
	# partitioning is the first effect.  Similarly accommodate 
	# for particle effect for: gas on particle effect, particle on gas effect
	# and particle on particle effect

	# Jacobian index for particle effects
	jac_part_indx = np.zeros((comp_num+2)*(num_asb+1)+((comp_num+2)*(num_asb*2)))
	
	part_cnt = 0 # count on jac_part_indx inputs
	
	if (num_asb>0): # if particle size bins are present

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
					jac_indx_g[jac_indx_g>=new_indx] += 1 # gas-phase reactions
					jac_indx_aq[jac_indx_aq>=new_indx] += 1 # aqueous-phase reactions
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
				jac_indx_g[jac_indx_g>=st_indx] += 1 # gas-phase reactions
				jac_indx_aq[jac_indx_aq>=st_indx] += 1 # aqueous-phase reactions
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
			jac_indx_g[jac_indx_g>=en_indx] += num_asb # gas-phase reaction
			jac_indx_aq[jac_indx_aq>=en_indx] += num_asb # aqueous-phase reactions
			new_el = (np.array(range(comp_num+2+compi, 
				(comp_num+2)*(num_asb+1)+compi, (comp_num+2)))).reshape(-1)
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
				jac_indx_aq[jac_indx_aq>=st_indx] += 1 # aqueous-phase reactions
				part_cnt += 1
				
				# particle component effect on particle
				new_el = np.array(((comp_num+2)*(sbi+1)+compi)).reshape(-1)
				# account for any aqueous-phase reactions
				st_indx = (st_indx) + sum(rowvals[int(colptrs[stc_indx]):int(colptrs[stc_indx+1])] < new_el)
				jac_part_indx[part_cnt] = st_indx
				# check if diagonal already in use by aqueous-phase reactions
				if any(rowvals[int(colptrs[stc_indx]):int(colptrs[stc_indx+1])] == new_el):
					part_cnt += 1
					continue # continue to next component if it is
				else:
					rowvals = np.concatenate([rowvals[0:st_indx], new_el, rowvals[st_indx::]])
					colptrs[(comp_num+2)*(sbi+1)+compi+1::] += 1
					jac_indx_aq[jac_indx_aq>=st_indx] += 1 # aqueous-phase reactions
					part_cnt += 1
	
	# wall influence on Jacobian part ---------------------------------------------------

	# loop through jacobian index to check whether the centre 
	# diagonal for gas effect on gas components is due to be 
	# affected by gas-phase reactions or whether wall 
	# partitioning is the first effect.  Similarly accommodate 
	# for wall effect for: gas on wall effect, wall on gas effect
	# and wall on wall affect
	jac_wall_indx = np.zeros(((comp_num+2)*4))
	wall_cnt = 0 # count on jac_wall_indx inputs

	if (wall_on > 0):
		# loop through components in the gas-phase (add two 
		# to account for water and core component)
		for compi in range(comp_num+2):
			# gas effect on gas part --------------------------------------------
			# relevant starting and finishing index in rowvals
			st_indx = int(colptrs[compi])
			en_indx = int(colptrs[compi+1])
			
			# check if any rows already attributed to this column
			if ((st_indx<en_indx) == True):
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
					jac_indx_g[jac_indx_g>=new_indx] += 1
					jac_indx_aq[jac_indx_aq>=new_indx] += 1
					new_el = np.array((compi)).reshape(1)
					rowvals = np.concatenate([rowvals[0:new_indx], 
						new_el, rowvals[new_indx::]])
					colptrs[compi+1::] += 1
					# increase en_indx to account for new row when
					# accommodating gas effect on wall below
					en_indx += 1
					
			else: # no rows yet attributed to this column, so need to include

				# modify indices for sparse Jacobian matrix
				jac_wall_indx[wall_cnt] = st_indx
				jac_indx_g[jac_indx_g>=st_indx] += 1
				jac_indx_aq[jac_indx_aq>=st_indx] += 1
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
			jac_wall_indx[wall_cnt] = en_indx
			# bump gas effect and particle effect indices up by one to 
			# account for addition to Jacobian
			jac_indx_g[jac_indx_g>=en_indx] += 1
			jac_indx_aq[jac_indx_aq>=en_indx] += 1
			if (num_asb > 0):
				jac_part_indx[jac_part_indx>=en_indx] += 1
			new_el = np.array(((comp_num+2)*(num_asb+1)+compi)).reshape(1)
			rowvals = np.concatenate([rowvals[0:en_indx], new_el, rowvals[en_indx::]])
			colptrs[compi+1::] += 1
			wall_cnt += 1 # keep count on wall index
		
		# wall effect on gas phase and on wall components for the Jacobian, note this will include the
		# final row in the final column of the Jacobian
		for compi in range(comp_num+2):
			# wall component effect on itself in gas-phase
			jac_wall_indx[wall_cnt] = len(rowvals)
			new_el = np.array((compi)).reshape(1)
			rowvals = np.append(rowvals, new_el)
			colptrs[(comp_num+2)*(num_asb+1)+compi+1::] += 1
			wall_cnt += 1
			# wall component effect on itself on wall
			jac_wall_indx[wall_cnt] = len(rowvals)
			new_el = np.array(((comp_num+2)*(num_asb+1)+compi)).reshape(1)
			rowvals = np.append(rowvals, new_el)
			colptrs[(comp_num+2)*(num_asb+1)+compi+1::] += 1
			wall_cnt += 1
	# end of wall influence on Jacobian part ---------------------------------------------------
		 		
	if (num_sb == 0): # if no particle size bins and no wall

		# if the Jacobian matrix has an empty final row, then
		# an error will be displayed during ODE solver call, so 
		# input a filler on this row;
		# add two to the number of unique components counted here to
		# acount for water and seed material
		if max(rowvals)<=(comp_num+2)-1:
			# index for final row of Jacobian
			rowvals = np.append(rowvals, (comp_num+2)-1)
			# number of columns currently short of final column
			col_shrt =  comp_num+2-len(colptrs)
			# rowval index for final column of Jacobian		
			colptrs[-1] += 1

	# ensure integer type for arrays acting as indices
	jac_indx_g = jac_indx_g.astype(int)
	jac_indx_aq = jac_indx_aq.astype(int)
	jac_part_indx = jac_part_indx.astype(int)
	jac_wall_indx = jac_wall_indx.astype(int)	
	rowvals = rowvals.astype(int)	
	colptrs = colptrs.astype(int)



	return(rowvals, colptrs, jac_indx_g, jac_indx_aq, jac_part_indx, jac_wall_indx)
