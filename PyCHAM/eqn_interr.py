'''module to interrogate equations to withdraw essential information for solution'''
# code to extract the equation information for chemical reactions required for 
# their solution in PyCHAM 

import numpy as np
import re
import formatting
import pybel
import sys

def eqn_interr(num_eqn, eqn_list, aqeqn_list, chem_scheme_markers, spec_name, 
		spec_smil, num_sb, wall_on):
				
	# inputs: ----------------------------------------------------------------------------
	# num_eqn - number of equations
	# eqn_list - gas-phase equations in strings
	# aqeqn_list - aqueous-phase equations in strings
	# chem_scheme_markers - markers for separating sections of the chemical scheme
	# spec_name - name string of components in xml file (not SMILES)
	# spec_smil - SMILES from xml file
	# num_sb - number of size bins
	# wall_on - marker for whether to include wall partitioning
	# ------------------------------------------------------------------------------------
	
	# preparatory part ----------------------------------------------------
	# matrix to record indices of reactants (cols) in each equation (rows)
	rindx = np.zeros((num_eqn[0], 1)).astype(int)
	# matrix to arrange concentrations when reaction rate coefficient calculated
	y_arr = (np.ones((num_eqn[0], 1)).astype(int))*-9999
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
	pindx = np.zeros((num_eqn[0], 1)).astype(int)
	# matrix to record stoichiometries of reactants (cols) in each equation (rows)
	rstoi = np.ones((num_eqn[0], 1))
	jac_stoi = np.zeros((num_eqn[0], 1))
	# 1D array to record stoichiometries of reactants per equarion
	rstoi_flat = np.empty((0))
	# 1D array to record stoichiometries of products per equarion
	pstoi_flat = np.empty((0))
	# matrix to record stoichiometries of products (cols) in each equation (rows)
	pstoi = np.zeros((num_eqn[0], 1))
	# arrays to store number of reactants and products in gas-phase equations
	nreac = np.empty(num_eqn[0], dtype=np.int8)
	nprod = np.empty(num_eqn[0], dtype=np.int8)
	# list for equation reaction rate coefficients
	reac_coef = []
	# matrix containing index of components who are denominators in the
	# calculation of equation derivatives in the Jacobian
	jac_den_indx = np.zeros((num_eqn[0], 1))
	# total number of Jacobian elements per equation
	njac = np.zeros((num_eqn[0], 1))
	# indices of Jacobian to affect per equation (rows)
	jac_indx = np.zeros((num_eqn[0], 1))
	# a new list for the name strings of components presented in the scheme (not SMILES)
	comp_namelist = []
	# list for the SMILE strings of components present in the chemical scheme
	comp_list = []
	# list of Pybel objects of components in chemical scheme
	Pybel_objects = []
	comp_num = 0 # count the number of unique components in the chemical scheme
	# ---------------------------------------------------------------------

	max_no_reac = 0.0 # log maximum number of reactants in a reaction
	max_no_prod = 0.0 # log maximum number of products in a reaction

	# loop through gas-phase equations line by line and extract the required information
	for eqn_step in range(num_eqn[0]):
		
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
			rindx = np.append(rindx, (np.zeros((num_eqn[0], 1))).astype(int), axis=1)
			rstoi = np.append(rstoi, (np.ones((num_eqn[0], 1))), axis=1)
			y_arr = np.append(y_arr, (np.ones((num_eqn[0], 1))*-9999).astype(int), axis=1)
			y_arr_fixer = ((np.arange(0, num_eqn[0], dtype = 'int')).reshape(-1, 1)).repeat(max_no_reac, axis=1)
			y_arr[y_arr!=-9999] = y_arr[y_arr!=-9999]+y_arr_fixer[y_arr!=-9999] 
		while max_no_prod > np.minimum(pindx.shape[1], pstoi.shape[1]): 
			pindx = np.append(pindx, (np.zeros((num_eqn[0], 1))).astype(int), axis=1)
			pstoi = np.append(pstoi, (np.zeros((num_eqn[0], 1))), axis=1)
		while ((len(reactants)**2.0+len(reactants)*len(products))>jac_indx.shape[1]):
			jac_indx = np.append(jac_indx, (np.zeros((num_eqn[0], 1))), axis=1)
			jac_den_indx = np.append(jac_den_indx, (np.zeros((num_eqn[0], 1))), axis=1)		
			jac_stoi = np.append(jac_stoi, (np.zeros((num_eqn[0], 1))), axis=1)

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

			if name_only not in comp_namelist: # if new component encountered
				comp_namelist.append(name_only) # add to chemical scheme name list
			
				# convert MCM chemical names to SMILES
				if (name_only in spec_name):
					# index where xml file name matches reaction component name
					name_indx = spec_name.index(name_only)
					name_SMILE = spec_smil[name_indx] # SMILES of component
				else:
					print(str('Error: inside eqn_parser, chemical scheme name '+str(name_only)+' not found in xml file'))
					sys.exit()
			
				comp_list.append(name_SMILE) # list SMILE names
				name_indx = comp_num # allocate index to this species
				# generate pybel object
				Pybel_object = pybel.readstring('smi', name_SMILE)
				# append to Pybel object list
				Pybel_objects.append(Pybel_object)
				
				comp_num += 1 # number of unique species
				

			else: # if it's a species already encountered it will be in comp_list
				# existing index
				name_indx = comp_namelist.index(name_only)
			
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
			if name_only not in comp_namelist: # if new component encountered
				comp_namelist.append(name_only)
				
				# convert MCM chemical names to SMILES
				# index where xml file name matches reaction component name
				if name_only in spec_name:
					name_indx = spec_name.index(name_only)
					name_SMILE = spec_smil[name_indx]
				else:
					print('Error: inside eqn_interr, chemical scheme name '+str(name_only)+' not found in xml file')
					sys.exit()
				
				comp_list.append(name_SMILE) # list SMILE string of parsed species
				name_indx = comp_num # allocate index to this species
				# Generate pybel
				
				Pybel_object = pybel.readstring('smi', name_SMILE)
				# append to Pybel object list
				Pybel_objects.append(Pybel_object)
				
				comp_num += 1 # number of unique species
				
			

			else: # if it's a species already encountered
				# index of component already listed
				name_indx = comp_namelist.index(name_only)
				
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
			if (i > 0):
				jac_stoi[eqn_step, i*tot_comp:(i+1)*tot_comp] = jac_stoi[eqn_step, 0:tot_comp] 
 		# number of Jacobian elements affected by this equation
		njac[eqn_step, 0] = tot_comp*nreac[eqn_step]
	
	# flatten index for arranging concentrations ready for reaction rate coefficient calculation
	y_arr_g = y_arr.flatten(order='C')
	y_arr_g = y_arr[y_arr != -9999] # remove fillers
	y_rind_g = y_rind.astype(int) # ensure integer type
	uni_y_rind_g = (np.unique(y_rind)).astype(int) # unique index of reactants
	y_pind_g = y_pind.astype(int) # ensure integer type
	uni_y_pind_g = (np.unique(y_pind)).astype(int) # unique index of products
	rr_arr_g = rr_arr.astype(int) # ensure integer type
	rr_arr_p_g = rr_arr_p.astype(int) # ensure integer type	
	# colptrs for sparse matrix of the change to reactants per equation
	reac_col_g = np.cumsum(nreac)-nreac
	# include final column
	reac_col_g = np.append(reac_col_g, reac_col_g[-1]+nreac[-1])
	# colptrs for sparse matrix of the change to products per equation
	prod_col_g = np.cumsum(nprod)-nprod
	# include final column
	prod_col_g = np.append(prod_col_g, prod_col_g[-1]+nprod[-1])
	# tag other gas-phase arrays
	rindx_g = rindx
	pindx_g = pindx
	rstoi_g = rstoi
	pstoi_g = pstoi
	jac_stoi_g = jac_stoi
	rstoi_flat_g = rstoi_flat
	pstoi_flat_g = pstoi_flat
	nreac_g = nreac
	nprod_g = nprod
	reac_coef_g = reac_coef
	jac_den_indx_g = jac_den_indx.astype(int)
	njac_g = njac.astype(int)
	jac_indx_g = jac_indx
	jac_indx_g = jac_indx_g.astype(int)

	# repeat for aqueous-phase reactions ----------------------------------
	# preparatory part ----------------------------------------------------
	# matrix to record indices of reactants (cols) in each equation (rows)
	rindx = np.zeros((num_eqn[1], 1)).astype(int)
	# matrix to arrange concentrations when reaction rate coefficient calculated
	y_arr = (np.ones((num_eqn[1], 1)).astype(int))*-9999
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
	pindx = np.zeros((num_eqn[1], 1)).astype(int)
	# matrix to record stoichiometries of reactants (cols) in each equation (rows)
	rstoi = np.ones((num_eqn[1], 1))
	jac_stoi = np.zeros((num_eqn[1], 1))
	# 1D array to record stoichiometries of reactants per equarion
	rstoi_flat = np.empty((0))
	# 1D array to record stoichiometries of products per equarion
	pstoi_flat = np.empty((0))
	# matrix to record stoichiometries of products (cols) in each equation (rows)
	pstoi = np.zeros((num_eqn[1], 1))
	# arrays to store number of reactants and products of equations
	nreac = np.empty(num_eqn[1], dtype=np.int8)
	nprod = np.empty(num_eqn[1], dtype=np.int8)
	# list for equation reaction rate coefficients
	reac_coef = []
	# matrix containing index of components who are denominators in the
	# calculation of equation derivatives in the Jacobian
	jac_den_indx = np.zeros((num_eqn[1], 1))
	# total number of Jacobian elements per equation
	njac = np.zeros((num_eqn[1], 1))
	# indices of Jacobian to affect per equation (rows)
	jac_indx = np.zeros((num_eqn[1], 1))
	# ---------------------------------------------------------------------

	max_no_reac = 0.0 # log maximum number of reactants in a reaction
	max_no_prod = 0.0 # log maximum number of products in a reaction

	# loop through aqueous-phase equations line by line and extract the required information
	for eqn_step in range(num_eqn[1]):
		
		line = aqeqn_list[eqn_step] # extract this line
		
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
			rindx = np.append(rindx, (np.zeros((num_eqn[1], 1))).astype(int), axis=1)
			rstoi = np.append(rstoi, (np.ones((num_eqn[1], 1))), axis=1)
			y_arr = np.append(y_arr, (np.ones((num_eqn[1], 1))*-9999).astype(int), axis=1)
			y_arr_fixer = ((np.arange(0, num_eqn[1], dtype = 'int')).reshape(-1, 1)).repeat(max_no_reac, axis=1)
			y_arr[y_arr!=-9999] = y_arr[y_arr!=-9999]+y_arr_fixer[y_arr!=-9999] 
		while max_no_prod > np.minimum(pindx.shape[1], pstoi.shape[1]): 
			pindx = np.append(pindx, (np.zeros((num_eqn[1], 1))).astype(int), axis=1)
			pstoi = np.append(pstoi, (np.zeros((num_eqn[1], 1))), axis=1)
		while ((len(reactants)**2.0+len(reactants)*len(products))>jac_indx.shape[1]):
			jac_indx = np.append(jac_indx, (np.zeros((num_eqn[1], 1))), axis=1)
			jac_den_indx = np.append(jac_den_indx, (np.zeros((num_eqn[1], 1))), axis=1)		
			jac_stoi = np.append(jac_stoi, (np.zeros((num_eqn[1], 1))), axis=1)

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

			if name_only not in comp_namelist: # if new component encountered
				comp_namelist.append(name_only) # add to chemical scheme name list
			
				# convert MCM chemical names to SMILES
				if name_only in spec_name:
					# index where xml file name matches reaction component name
					name_indx = spec_name.index(name_only)
					name_SMILE = spec_smil[name_indx] # SMILES of component
				else:
					print(str('Error: inside eqn_parser, chemical scheme name '+str(name_only)+' not found in xml file'))
					sys.exit()
			
				comp_list.append(name_SMILE) # list SMILE names
				name_indx = comp_num # allocate index to this species
				# Generate pybel
				Pybel_object = pybel.readstring('smi', name_SMILE)
				# append to Pybel object list
				Pybel_objects.append(Pybel_object)
				
				comp_num += 1 # number of unique species
				

			else: # if it's a species already encountered it will be in comp_list
				# existing index
				name_indx = comp_namelist.index(name_only)
			
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
			if name_only not in comp_namelist: # if new component encountered
				comp_namelist.append(name_only)
				
				# convert MCM chemical names to SMILES
				# index where xml file name matches reaction component name
				if name_only in spec_name:
					name_indx = spec_name.index(name_only)
					name_SMILE = spec_smil[name_indx]
				else:
					print('Error: inside eqn_interr, chemical scheme name '+str(name_only)+' not found in xml file')
					sys.exit()
				
				comp_list.append(name_SMILE) # list SMILE string of parsed species
				name_indx = comp_num # allocate index to this species
				
				# generate pybel object
				Pybel_object = pybel.readstring('smi', name_SMILE)
				# append to Pybel object list
				Pybel_objects.append(Pybel_object)
				
				comp_num += 1 # number of unique species
				
			

			else: # if it's a species already encountered
				# index of component already listed
				name_indx = comp_namelist.index(name_only)
				
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
			if (i > 0):
				jac_stoi[eqn_step, i*tot_comp:(i+1)*tot_comp] = jac_stoi[eqn_step, 0:tot_comp] 
 		# number of Jacobian elements affected by this equation
		njac[eqn_step, 0] = tot_comp*nreac[eqn_step]
	
	# account for gas-phase in Jacobian denominator index
	jac_den_indx += (comp_num+2)

	# flatten index for arranging concentrations ready for reaction rate coefficient calculation
	y_arr_aq = y_arr.flatten(order='C')
	y_arr_aq = y_arr[y_arr != -9999] # remove fillers
	y_rind_aq = y_rind.astype(int) # ensure integer type
	uni_y_rind_aq = (np.unique(y_rind)).astype(int) # unique index of reactants
	y_pind_aq = y_pind.astype(int) # ensure integer type
	uni_y_pind_aq = (np.unique(y_pind)).astype(int) # unique index of products
	rr_arr_aq = rr_arr.astype(int) # ensure integer type
	rr_arr_p_aq = rr_arr_p.astype(int) # ensure integer type	
	# colptrs for sparse matrix of the change to reactants per equation
	reac_col_aq = np.cumsum(nreac)-nreac
	# colptrs for sparse matrix of the change to products per equation
	prod_col_aq = np.cumsum(nprod)-nprod
	if (len(reac_col_aq) > 0): # if aqueous-phase reaction present	
		# include final columns
		reac_col_aq = np.append(reac_col_aq, reac_col_aq[-1]+nreac[-1])
		prod_col_aq = np.append(prod_col_aq, prod_col_aq[-1]+nprod[-1])
	
	# tag other aqueous-phase arrays
	rindx_aq = rindx
	pindx_aq = pindx
	rstoi_aq = rstoi
	pstoi_aq = pstoi
	jac_stoi_aq = jac_stoi
	rstoi_flat_aq = rstoi_flat
	pstoi_flat_aq = pstoi_flat
	nreac_aq = nreac
	nprod_aq = nprod
	reac_coef_aq = reac_coef
	jac_den_indx_aq = jac_den_indx.astype(int)
	njac_aq = njac.astype(int)
	jac_indx_aq = jac_indx
	jac_indx_aq = jac_indx_aq.astype(int)


	return(rindx_g, rstoi_g, pindx_g, pstoi_g, reac_coef_g, 
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
			comp_num)
