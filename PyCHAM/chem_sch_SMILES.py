'''module to scan the chemical scheme to extract all unique component names'''
# based on the code originally contained in eqn_interr (as of PyCHAM version 3.0.1)
# the supplied chemical scheme is interrogated to extract the unique component names
# and their SMILES strings, which requires reference to the supplied xml file

import sch_interr
import xml_interr
import re

def chem_scheme_SMILES_extr(sch_name, xml_name, chem_scheme_markers):

	# inputs: ------------------------------------------
	# sch_name - path to chemical scheme
	# xml_name - path to xml name
	# chem_scheme_markers - markers for chemical scheme
	# ----------------------------------------------------
	
	# initiate with empty error message
	err_mess = ''

	f_open_eqn = open(sch_name, mode='r') # open the chemical scheme file
	# read the file and store everything into a list
	total_list_eqn = f_open_eqn.readlines()
	f_open_eqn.close() # close file

	# interrogate scheme to list equations
	[eqn_list, aqeqn_list, eqn_num, rrc, rrc_name, 
		RO2_names] = sch_interr.sch_interr(total_list_eqn, chem_scheme_markers)
	
	# interrogate xml to list all component names and SMILES
	[comp_smil, comp_name] = xml_interr.xml_interr(xml_name)
	
	comp_namelist = [] # list for chemical scheme names of components in the chemical scheme
	comp_list = [] # list for the SMILE strings of components present in the chemical scheme

	for eqn_step in range(eqn_num[0]): # loop through gas-phase reactions

		line = eqn_list[eqn_step] # extract this line
		
		# work out whether equation or reaction rate coefficient part comes first
		eqn_start = str('.*\\' +  chem_scheme_markers[10])
		rrc_start = str('.*\\' +  chem_scheme_markers[9])
		# get index of these markers, note span is the property of the match object that
		# gives the location of the marker
		eqn_start_indx = (re.match(eqn_start, line)).span()[1]
		rrc_start_indx = (re.match(rrc_start, line)).span()[1]

		if (eqn_start_indx>rrc_start_indx):
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
		# reactants with stoichiometry number and omit any photon
		reactants = [i for i in eqn_split[:eqmark_pos] if i != '+' and i != 'hv']
		# products with stoichiometry number
		products = [t for t in eqn_split[eqmark_pos+1:] if t != '+']
		
		stoich_regex = r"^\d*\.\d*|^\d*" # necessary for checking for stoichiometries
		
		for reactant in reactants: # left hand side of equations (losses)
		
			if (re.findall(stoich_regex, reactant)[0] != ''): # when stoichiometry included
				# name with no stoichiometry number
				name_only = re.sub(stoich_regex, '', reactant)
			elif (re.findall(stoich_regex, reactant)[0] == ''): # when stoichiometry excluded
				name_only = reactant
			
			if (name_only not in comp_namelist): # if new component encountered
			
				comp_namelist.append(name_only) # add to chemical scheme name list
				
				# convert MCM chemical names to SMILES
				if (name_only in comp_name):
					# index where xml file name matches reaction component name
					name_indx = comp_name.index(name_only)
					name_SMILE = comp_smil[name_indx] # SMILES of component
					comp_list.append(name_SMILE) # list SMILE names
				else:
					err_mess = str('Error - chemical scheme name '+str(name_only)+' not found in xml file')
					H2Oi = 0 # filler
					return(comp_namelist, comp_list, err_mess, H2Oi)
		
		for product in products: # right hand side of equations (gains)
			
			if (re.findall(stoich_regex, product)[0] != ''): # when stoichiometry included
				# name with no stoichiometry number
				name_only = re.sub(stoich_regex, '', product)
			elif (re.findall(stoich_regex, product)[0] == ''): # when stoichiometry excluded
				name_only = product
			
			if (name_only not in comp_namelist): # if new component encountered
			
				comp_namelist.append(name_only) # add to chemical scheme name list
			
				# convert MCM chemical names to SMILES
				# index where xml file name matches reaction component name
				if name_only in comp_name:
					name_indx = comp_name.index(name_only)
					name_SMILE = comp_smil[name_indx]
					comp_list.append(name_SMILE) # list SMILE names
				else:
					err_mess = str('Error - chemical scheme name '+str(name_only)+' not found in xml file')
					return(comp_namelist, comp_list, err_mess, H2Oi)
		
	
	# check for water presence in chemical scheme via its SMILE string
	# count on components
	indx = -1
	for single_chem in comp_list:
		indx += 1
		# ensure this is water rather than single oxygen (e.g. due to ozone photolysis 
		# (O is the MCM chemical scheme name for single oxygen))
		if (single_chem == 'O' and comp_namelist[indx] != 'O'):
				H2Oi = indx
				
	else: # if not in chemical scheme, water would be the next component appended to component list
		H2Oi = len(comp_list)
		comp_namelist.append('H2O')
	
	# if no error message but no equations identified then tell user the chemical scheme markers is probably wrong
	if (err_mess == '' and eqn_num[0] == 0):
		err_mess = 'Error: no gas-phase reactions seen, this could be due to the chemical scheme marker input (chem_scheme_markers in the model variables input) not corresponding to the chemical scheme file, please see README for more guidance.'
	
	
	return(comp_namelist, comp_list, err_mess, H2Oi)