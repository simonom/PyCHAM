'''module to interrogate equations to withdraw essential information for solution'''
# code to extract the equation information for chemical reactions required for 
# their solution in PyCHAM 

import numpy as np
import re
import formatting
import pybel
import sys

def eqn_interr(num_eqn, eqn_list, chem_scheme_markers, spec_name, spec_smil, phase, num_sb, wall_on):
				
	# inputs: ----------------------------------------------------------------------------
	# num_eqn - number of equations (scalar)
	# eqn_list - equations in strings
	# chem_scheme_markers - markers for separating sections of the chemical scheme
	# spec_name - name string of components in xml file (not SMILES)
	# spec_smil - SMILES from xml file
	# phase - marker for the phase being considered: 0 for gas, 1 for particulates
	# num_sb - number of size bins
	# wall_on - marker for whether to include wall partitioning
	# ------------------------------------------------------------------------------------
	
	# preparatory part ----------------------------------------------------
	comp_num = 0 # count the number of unique components
	# matrix to record indices of reactants (cols) in each equation (rows)
	rindx = np.zeros((num_eqn, 1)).astype(int)
	# matrix to arrange concentrations when reaction rate coefficient calculated
	y_arr = (np.ones((num_eqn, 1)).astype(int))*-9999
	# array to arrange reaction rates so they align with reactant stoichiometries
	rr_arr = np.empty((0))
	# same but for products
	rr_arr_p = np.empty((0))
	# index array for extracting required reactant concentrations for the
	# reaction rate coefficient calculation
	y_rind = np.empty((0))
	# index array for identifying products when assigning gains from reactions
	y_pind = np.empty((0))	
	# matrix to record indices of products (cols) in each equation (rows)
	pindx = np.zeros((num_eqn, 1)).astype(int)
	# matrix to record stoichiometries of reactants (cols) in each equation (rows)
	rstoi = np.ones((num_eqn, 1))
	jac_stoi = np.zeros((num_eqn, 1))
	# 1D array to record stoichiometries of reactants per equarion
	rstoi_flat = np.empty((0))
	# 1D array to record stoichiometries of products per equarion
	pstoi_flat = np.empty((0))
	# matrix to record stoichiometries of products (cols) in each equation (rows)
	pstoi = np.zeros((num_eqn, 1))
	# arrays to store number of reactants and products in gas-phase equations
	nreac = np.empty(num_eqn, dtype=np.int8)
	nprod = np.empty(num_eqn, dtype=np.int8)
	# list for equation reaction rate coefficients
	reac_coef = []
	# list for components' SMILE strings
	spec_list = []
	# list of Pybel objects
	Pybel_objects = []
	# a new list for the name strings of species presented in the scheme (not SMILES)
	spec_namelist = []
	# matrix containing index of components who are denominators in the
	# calculation of equation derivatives in the Jacobian
	jac_den_indx = np.zeros((num_eqn, 1))
	# total number of Jacobian elements per equation
	njac = np.zeros((num_eqn, 1))	
	# indices of Jacobian to affect per equation (rows)
	jac_indx = np.zeros((num_eqn, 1))
	# ---------------------------------------------------------------------

	max_no_reac = 0.0 # log maximum number of reactants in a reaction
	max_no_prod = 0.0 # log maximum number of products in a reaction

	# loop through equations line by line and extract the required information
	for eqn_step in range(num_eqn):
		
		line = eqn_list[eqn_step] # extract this line
		
		# work out whether equation or reaction rate coefficient part comes first
		eqn_start = str('.*\\' +  chem_scheme_markers[10])
		rrc_start = str('.*\\' +  chem_scheme_markers[9])
		# get index of these markers, note span is the property of the match object that
		# gives the location of the marker
		eqn_start_indx = (re.match(eqn_start, line)).span()[1]
		rrc_start_indx = (re.match(rrc_start, line)).span()[1]
		
		if eqn_start_indx>rrc_start_indx:
			eqn_sec = 1 # equation is second part
		else:
			eqn_sec = 0 # equation is first part
		
		# split the line into 2 parts: equation and rate coefficient
		# . means match with anything except a new line character., when followed by a * 
		# means match zero or more times (so now we match with all characters in the line
		# except for new line characters, so final part is stating the character(s) we 
		# are specifically looking for, \\ ensures the marker is recognised
		if eqn_sec == 1:
			eqn_markers = str('\\' +  chem_scheme_markers[10]+ '.*\\' +  chem_scheme_markers[11])
		else: # end of equation part is start of reaction rate coefficient part
			eqn_markers = str('\\' +  chem_scheme_markers[10]+ '.*\\' +  chem_scheme_markers[9])

		# extract the equation as a string ([0] extracts the equation section and 
		# [1:-1] removes the bounding markers)
		eqn = re.findall(eqn_markers, line)[0][1:-1].strip()
		
		eqn_split = eqn.split()
		eqmark_pos = eqn_split.index('=')
		# with stoich number; rule out the photon
		reactants = [i for i in eqn_split[:eqmark_pos] if i != '+' and i != 'hv']
		products = [t for t in eqn_split[eqmark_pos+1:] if t != '+'] # with stoich number
		
		# record maximum number of reactants across all equations
		max_no_reac = np.maximum(len(reactants), max_no_reac)
		# record maximum number of products across all equations
		max_no_prod = np.maximum(len(products), max_no_prod)

		# append columns if needed
		while max_no_reac > np.minimum(rindx.shape[1], rstoi.shape[1]): 
			rindx = np.append(rindx, (np.zeros((num_eqn, 1))).astype(int), axis=1)
			rstoi = np.append(rstoi, (np.ones((num_eqn, 1))), axis=1)
			y_arr = np.append(y_arr, (np.ones((num_eqn, 1))*-9999).astype(int), axis=1)
			y_arr_fixer = ((np.arange(0, num_eqn, dtype = 'int')).reshape(-1, 1)).repeat(max_no_reac, axis=1)
			y_arr[y_arr!=-9999] = y_arr[y_arr!=-9999]+y_arr_fixer[y_arr!=-9999] 
		while max_no_prod > np.minimum(pindx.shape[1], pstoi.shape[1]): 
			pindx = np.append(pindx, (np.zeros((num_eqn, 1))).astype(int), axis=1)
			pstoi = np.append(pstoi, (np.zeros((num_eqn, 1))), axis=1)
		while ((len(reactants)**2.0+len(reactants)*len(products))>jac_indx.shape[1]):
			jac_indx = np.append(jac_indx, (np.zeros((num_eqn, 1))), axis=1)
			jac_den_indx = np.append(jac_den_indx, (np.zeros((num_eqn, 1))), axis=1)		
			jac_stoi = np.append(jac_stoi, (np.zeros((num_eqn, 1))), axis=1)

		# .* means occurs anywhere in line and, first \ means second \ can be interpreted 
		# and second \ ensures recognition of marker
		rate_coeff_start_mark = str('\\' +  chem_scheme_markers[9])
		# . means match with anything except a new line character, when followed by a * 
		# means match zero or more times (so now we match with all characters in the line
		# except for new line characters, \\ ensures the marker
		# is recognised
		if eqn_sec == 1: # end of reaction rate coefficient part is start of equation part
			rate_coeff_end_mark = str('.*\\' +  chem_scheme_markers[10])
		else: # end of reaction rate coefficient part is end of line
			rate_coeff_end_mark = str('.*\\' +  chem_scheme_markers[11])
		
		# rate coefficient starts and end punctuation
		rate_regex = str(rate_coeff_start_mark + rate_coeff_end_mark)
		# rate coefficient expression in a string
		rate_ex = re.findall(rate_regex, line)[0][1:-1].strip()

		# convert fortran-type scientific notation to python type
		rate_ex = formatting.SN_conversion(rate_ex)
		# convert the rate coefficient expressions into Python readable commands
		rate_ex = formatting.convert_rate_mcm(rate_ex)
		if (rate_ex.find('EXP') != -1):
			print('Error in reaction rate coefficient expression: ', rate_ex)
			sys.exit()
		
		# store the reaction rate coefficient for this equation 
		# (/s once any inputs applied)
		reac_coef.append(rate_ex)
		
		# extract the stoichiometric number of the component in current equation
		reactant_step = 0
		product_step = 0
		stoich_regex = r"^\d*\.\d*|^\d*"
		numr = len(reactants) # number of reactants in this equation
		
		
		# left hand side of equations (losses)
		for reactant in reactants:
				
			if (re.findall(stoich_regex, reactant)[0] != ''):
				stoich_num = float(re.findall(stoich_regex, reactant)[0])
				# name with no stoich number
				name_only = re.sub(stoich_regex, '', reactant)
			elif (re.findall(stoich_regex, reactant)[0] == ''):
				stoich_num = 1.0
				name_only = reactant
			
			# store stoichiometry
			rstoi[eqn_step, reactant_step] = stoich_num
			jac_stoi[eqn_step, reactant_step] = -1*stoich_num

			if name_only not in spec_namelist: # if new component encountered
				spec_namelist.append(name_only) # add to chemical scheme name list
			
				# convert MCM chemical names to SMILES
				if name_only in spec_name:
					# index where xml file name matches reaction component name
					name_indx = spec_name.index(name_only)
					name_SMILE = spec_smil[name_indx] # SMILES of component
				else:
					sys.exit(str('Error: inside eqn_parser, chemical scheme name '+str(name_only)+' not found in xml file'))
			
				spec_list.append(name_SMILE) # list SMILE names
				name_indx = comp_num # allocate index to this species
				# Generate pybel
				Pybel_object = pybel.readstring('smi', name_SMILE)
				# append to Pybel object list
				Pybel_objects.append(Pybel_object)
				
				comp_num += 1 # number of unique species
				

			else: # if it's a species already encountered it will be in spec_list
				# existing index
				name_indx = spec_namelist.index(name_only)
			
			# store reactant index
			# check if index already present - i.e. component appears more than once
			if sum(rindx[eqn_step, 0:reactant_step]==int(name_indx))>0:
				# get existing index of this component
				exist_indx = (np.where(rindx[eqn_step, 0:reactant_step]==(int(name_indx))))[0]
				# add to existing stoichiometry
				rstoi[eqn_step, exist_indx] += rstoi[eqn_step, reactant_step]
				jac_stoi[eqn_step, exist_indx] += -1*rstoi[eqn_step, reactant_step]
				# remove stoichiometry added above
				rstoi[eqn_step, reactant_step] = 1
				jac_stoi[eqn_step, reactant_step] = 0
				reactant_step -= 1 # ignore this duplicate
			else:
				rindx[eqn_step, reactant_step] = int(name_indx)
				y_arr[eqn_step, reactant_step] = int((eqn_step*max_no_reac)+reactant_step)
				y_rind = np.append(y_rind, int(name_indx))
				rr_arr = np.append(rr_arr, int(eqn_step))
		
			reactant_step += 1
					
		# number of reactants in this equation
		nreac[eqn_step] = int(reactant_step)
		
		# record 1D array of stoichiometries per equation
		rstoi_flat = np.append(rstoi_flat, rstoi[eqn_step, 0:int(reactant_step)])
		
		# right hand side of equations (gains)
		for product in products:

			if (re.findall(stoich_regex, product)[0] != ''):
				stoich_num = float(re.findall(stoich_regex, product)[0])
				name_only = re.sub(stoich_regex, '', product) # name with no stoich number

			elif (re.findall(stoich_regex, product)[0] == ''):
				stoich_num = 1.0
				name_only = product
			
			# store stoichiometry
			pstoi[eqn_step, product_step] = stoich_num
			jac_stoi[eqn_step, reactant_step+product_step] = 1*stoich_num
			if name_only not in spec_namelist: # if new component encountered
				spec_namelist.append(name_only)
				
				# convert MCM chemical names to SMILES
				# index where xml file name matches reaction component name
				if name_only in spec_name:
					name_indx = spec_name.index(name_only)
					name_SMILE = spec_smil[name_indx]
				else:
					sys.exit(str('Error: inside eqn_parser, chemical scheme name '+str(name_only)+' not found in xml file'))
				
				spec_list.append(name_SMILE) # list SMILE string of parsed species
				name_indx = comp_num # allocate index to this species
				# Generate pybel
				
				Pybel_object = pybel.readstring('smi', name_SMILE)
				# append to Pybel object list
				Pybel_objects.append(Pybel_object)
				
				comp_num += 1 # number of unique species
				
			

			else: # if it's a species already encountered
				# index of component already listed
				name_indx = spec_namelist.index(name_only)
				
			# store product index
			# check if index already present - i.e. component appears more than once
			if sum(pindx[eqn_step, 0:product_step]==int(name_indx))>0:
				# get existing index of this component
				exist_indx = (np.where(pindx[eqn_step, 0:product_step]==(int(name_indx))))[0]
				# add to existing stoichiometry
				pstoi[eqn_step, exist_indx] += pstoi[eqn_step, product_step]
				jac_stoi[eqn_step, reactant_step+exist_indx] += 1*pstoi[eqn_step, product_step]
				# remove stoichiometry added above
				pstoi[eqn_step, product_step] = 0
				jac_stoi[eqn_step, reactant_step+product_step] = 0
				product_step -= 1 # ignore this duplicate
			else:
				pindx[eqn_step, product_step] = int(name_indx)
				rr_arr_p = np.append(rr_arr_p, int(eqn_step))
				y_pind = np.append(y_pind, int(name_indx))

			product_step += 1
		
		# number of products in this equation
		nprod[eqn_step] = int(product_step)
		# record 1D array of stoichiometries per equation
		pstoi_flat = np.append(pstoi_flat, pstoi[eqn_step, 0:int(product_step)])
		
		# now that total number of components (reactants and products) 
		# in an equation is known, repeat the reactant indices over all 
		# components
		tot_comp = nreac[eqn_step]+nprod[eqn_step]
		for i in range(nreac[eqn_step]):
			jac_den_indx[eqn_step, i*tot_comp:(i+1)*tot_comp] = rindx[eqn_step, i]
			# also repeat the stoichiometries for every reactant
			if i>0:
				jac_stoi[eqn_step, i*tot_comp:(i+1)*tot_comp] = jac_stoi[eqn_step, 0:tot_comp] 
 		# number of Jacobian elements affected by this equation
		njac[eqn_step, 0] = tot_comp*nreac[eqn_step]
	
	# flatten index for arranging concentrations ready for reaction rate coefficient calculation
	y_arr = y_arr.flatten(order='C')
	y_arr = y_arr[y_arr != -9999] # remove fillers
	y_rind = y_rind.astype(int) # ensure integer type
	uni_y_rind = np.unique(y_rind) # unique index of reactants
	y_pind = y_pind.astype(int) # ensure integer type
	uni_y_pind = np.unique(y_pind) # unique index of products
	rr_arr = rr_arr.astype(int) # ensure integer type
	rr_arr_p = rr_arr_p.astype(int) # ensure integer type	
	# colptrs for sparse matrix of the change to reactants per equation
	reac_col = np.cumsum(nreac)-nreac
	# include final column
	reac_col = np.append(reac_col, reac_col[-1]+nreac[-1])
	# colptrs for sparse matrix of the change to products per equation
	prod_col = np.cumsum(nprod)-nprod
	# include final column
	prod_col = np.append(prod_col, prod_col[-1]+nprod[-1])

	print('Preparing Jacobian inputs')
	# ensure integer type
	jac_den_indx = jac_den_indx.astype(int)
	njac = njac.astype(int)
	# rows in Jacobian affected by this equation		 	
	rowvals = np.empty((0))
	# indices of rowvals representing each component being 
	# differentiated by in Jacobian, note, needs to start
	# as zero to represent first index of rowval
	colptrs = np.zeros((1))
	# track columns of Jacobian affected
	col_tr = 0
	# length of Jacobian
	len_jac = comp_num*(num_sb+1)

	# index for Jacobian - note can only be done after all equations checked as only then
	# is the total number of unique components known
	# equation loop
	for eqni in range(num_eqn):
		# total number of components in this equation
		tot_comp = nreac[eqni]+nprod[eqni]
		# combined index of reactants and products in this equation
		totindx = np.append(rindx[eqni, 0:nreac[eqni]], pindx[eqni, 0:nprod[eqni]])
		# reactant loop (equivalent to columns of 2D Jacobian)
		for i in range(nreac[eqni]):
			# index of flattened Jacobian being affected in 
			# this equation, note this doesn't need to 
			# account for water and seed material as below (when unique elements found) all
			# full Jacobian matrix elements that aren't indexed are
			# omitted to create the sparse Jacobian matrix  
			jac_indx[eqni, i*tot_comp:(i+1)*tot_comp] = rindx[eqni, i]*len_jac+(totindx)
			# check if rowvals already sufficiently long to contain this reactant  
			if rindx[eqni, i]<=(len(colptrs)-2):

				# start index of rowvals where it is represented
				r_st = int(colptrs[rindx[eqni, i]])
				# end index of rowvals where it is represented
				r_en = int(colptrs[rindx[eqni, i]+1])
				# check if the components affected have been 
				# represented in previous equations
				for i2 in totindx:
					# if this reactant has affected this component in 
					# a previous reaction then jac_indx (set above) will
					# sum its effect and no need to account for
					# the Jacobian index with rowvals and 
					# colptrs				
					if sum(rowvals[r_st:r_en] == i2)>0:
						continue
					else: # account for in Jacobian
						new_el = np.array((float(i2))).reshape(1)
						# index of where new index enters rowvals (to ensure
						# indices per Jacobian column are present in
						# ascending order)
						r_in = r_st+int(sum(rowvals[r_st:r_en]<new_el))
						rowvals = np.concatenate([rowvals[0:r_in], new_el, rowvals[r_in::]])
						colptrs[rindx[eqni, i]+1::] += 1
						r_en += 1
						
						
						
			else: # if not accounted for in previous reaction
				# rows in Jacobian affected by this reactant 
				# in this equation	 	
				rowvals = np.append(rowvals, np.unique(totindx))
				while (rindx[eqni, i]>col_tr):
					colptrs = np.append(colptrs, colptrs[-1])
					col_tr += 1
				# indices of rowvals
				colptrs = np.append(colptrs, colptrs[-1]+len(np.unique(totindx)))
				# track number of columns of Jacobian affected				
				col_tr += 1

	# because above jac_indx contains the indices of data assuming
	# a complete, rather than sparse matrix, now correct to affect
	# the sparse matrix by reducing indices where gaps occur in
	# the complete data matrix
	# get unique values in jac_indx
	uni = np.unique(jac_indx)
	# now reduce these indices based on the number of indices
	# omitted
	diff = 0 # prepare for number of empty points
	for i in range(1, len(uni)):
		diff += (uni[i]-uni[i-1])-1
		jac_indx[jac_indx==uni[i]] -= diff
	
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

	# particle influence on Jacobian part ------------------------------------------------------------
	# loop through jacobian index to check whether the centre 
	# diagonal for gas effect on gas components is due to be 
	# affected by gas-phase reactions or whether particle 
	# partitioning is the first effect.  Similarly accommodate 
	# for particle effect for: gas on particle effect, particle on gas effect
	# and particle on particle effect
	num_asb = num_sb-wall_on # number of actual particle size bins
	jac_part_indx = np.zeros((comp_num+2)*(num_asb+1)+((comp_num+2)*(num_asb*2))) # Jacobian index for particle effects
	
	part_cnt = 0 # count on jac_part_indx inputs
	
	if (num_asb>0):

		# loop through components in the gas-phase (add two 
		# to account for water and seed material)
		for compi in range(comp_num+2):
			# gas effect on gas part ----------------------------------------------------------
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
					jac_part_indx[part_cnt] = exist_indx
					
				else: # if diagonal not already affected, then include
					# get index
					new_indx = st_indx+(sum(rowvals[st_indx:en_indx]<compi))+1
					# modify indices for sparse Jacobian matrix
					jac_part_indx[part_cnt] = new_indx
					jac_indx[jac_indx>=new_indx] += 1
					new_el = np.array((compi)).reshape(1)
					rowvals = np.concatenate([rowvals[0:new_indx], 
						new_el, rowvals[new_indx::]])
					colptrs[compi+1::] += 1
					# increase en_indx to account for new row when
					# accommodating gas effect on wall below
					en_indx += 1
					
			else: # no rows yet attributed to this column, so need to include

				# modify indices for sparse Jacobian matrix
				jac_part_indx[part_cnt] = st_indx
				jac_indx[jac_indx>=st_indx] += 1
				new_el = np.array((compi)).reshape(1)
				rowvals = np.concatenate([rowvals[0:st_indx], new_el, rowvals[st_indx::]])
				colptrs[compi+1::] += 1
				# increase en_indx to account for new row when accommodating
				# gas effect on particle below
				en_indx += 1
				
			part_cnt += 1 # keep count on particle index
			
			# gas effect on particle part -----------------------------------------------------------
			# already know that the gas-phase photochemistry will not
			# have affected this part of the Jacobian, so just need to
			# modify sparse Jacobian inputs accordingly
			jac_part_indx[part_cnt:part_cnt+num_asb] = range(en_indx, en_indx+num_asb)
			jac_indx[jac_indx>=en_indx] += num_asb
			new_el = (np.array(range(comp_num+2+compi, 
				(comp_num+2)*(num_asb+1)+compi, (comp_num+2)))).reshape(-1)
			rowvals = np.concatenate([rowvals[0:en_indx], new_el, rowvals[en_indx::]])
			colptrs[compi+1::] += num_asb
			part_cnt += num_asb # keep count on wall index		
		
		
		# now, we must loop through components to include their particle effect on
		# gas phase and on particle components for the Jacobian, note if no wall this 
		# will include the final row in the final column of the Jacobian
		for sbi in range(num_asb): # size bin loop
			for compi in range(comp_num+2): # component loop
				# particle component effect on gas-phase
				jac_part_indx[part_cnt] = len(rowvals)
				new_el = np.array((compi))
				rowvals = np.append(rowvals, new_el)
				colptrs[(comp_num+2)*(sbi+1)+compi+1::] += 1
				part_cnt += 1
				# particle component effect on particle
				jac_part_indx[part_cnt] = len(rowvals)
				new_el = np.array(((comp_num+2)*(sbi+1)+compi))
				rowvals = np.append(rowvals, new_el)
				colptrs[(comp_num+2)*(sbi+1)+compi+1::] += 1
				part_cnt += 1
	# wall influence on Jacobian part ----------------------------------------------------------------

	# loop through jacobian index to check whether the centre 
	# diagonal for gas effect on gas components is due to be 
	# affected by gas-phase reactions or whether wall 
	# partitioning is the first effect.  Similarly accommodate 
	# for wall effect for: gas on wall effect, wall on gas effect
	# and wall on wall affect
	jac_wall_indx = np.zeros(((comp_num+2)*4))
	wall_cnt = 0 # count on jac_wall_indx inputs

	if wall_on>0:
		# loop through components in the gas-phase (add two 
		# to account for water and seed material)
		for compi in range(comp_num+2):
			# gas effect on gas part ----------------------------------------------------------
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
					jac_indx[jac_indx>=new_indx] += 1
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
				jac_indx[jac_indx>=st_indx] += 1
				new_el = np.array((compi)).reshape(1)
				rowvals = np.concatenate([rowvals[0:st_indx], new_el, rowvals[st_indx::]])
				colptrs[compi+1::] += 1
				# increase en_indx to account for new row when accommodating
				# gas effect on wall below
				en_indx += 1
				
			wall_cnt += 1 # keep count on wall index
			
			# gas effect on wall part -----------------------------------------------------------
			# already know that the gas-phase photochemistry will not
			# have affected this part of the Jacobian, so just need to
			# modify sparse Jacobian inputs accordingly
			jac_wall_indx[wall_cnt] = en_indx
			# bump gas effect and particle effect indices up by one to 
			# account for addition to Jacobian
			jac_indx[jac_indx>=en_indx] += 1
			if num_asb>0:
				jac_part_indx[jac_part_indx>=en_indx] += 1
			new_el = np.array(((comp_num+2)*num_sb+compi)).reshape(1)
			rowvals = np.concatenate([rowvals[0:en_indx], new_el, rowvals[en_indx::]])
			colptrs[compi+1::] += 1
			wall_cnt += 1 # keep count on wall index
		
		# wall effect on
		# gas phase and on wall components for the Jacobian, note this will include the
		# final row in the final column of the Jacobian
		for compi in range(comp_num+2):
			# wall component effect on itself in gas-phase
			jac_wall_indx[wall_cnt] = len(rowvals)
			new_el = np.array((compi)).reshape(1)
			rowvals = np.append(rowvals, new_el)
			colptrs[(comp_num+2)*num_sb+compi+1::] += 1
			wall_cnt += 1
			# wall component effect on itself on wall
			jac_wall_indx[wall_cnt] = len(rowvals)
			new_el = np.array(((comp_num+2)*num_sb+compi)).reshape(1)
			rowvals = np.append(rowvals, new_el)
			colptrs[(comp_num+2)*num_sb+compi+1::] += 1
			wall_cnt += 1
	# end of wall influence on Jacobian part ---------------------------------------------------
		 		
	if (num_sb == 0):

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
	jac_indx = jac_indx.astype(int)
	jac_part_indx = jac_part_indx.astype(int)
	jac_wall_indx = jac_wall_indx.astype(int)	
	rowvals = rowvals.astype(int)	
	colptrs = colptrs.astype(int)
	
	return(rindx, rstoi, pindx, pstoi, reac_coef, spec_namelist, spec_list, 
			Pybel_objects, nreac, nprod, comp_num, jac_stoi, 
			jac_den_indx, njac, jac_indx, 				
			rowvals, colptrs, y_arr, y_rind, uni_y_rind, y_pind, 
			uni_y_pind, reac_col, prod_col, rstoi_flat, pstoi_flat, 
			rr_arr, rr_arr_p, jac_wall_indx, jac_part_indx)
