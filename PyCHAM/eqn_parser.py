''' module for parsing the equation file, including converting MCM names of components into SMILES strings - see the README for how to format the equation and xml files'''

import os
import re
import collections
import sys
import datetime
import numpy as np
import pybel
import formatting
import xmltodict # for opening and converting xml files to python dictionaries
import ipdb
from water_calc import water_calc

# ----------Extraction of eqn info----------
# Extract the mechanism information
def extract_mechanism(filename, xmlname, TEMP, PInit, testf, RH, 
						start_sim_time, lat, lon, act_flux_path, DayOfYear, 
						chem_scheme_markers, photo_par_file):

	
	# inputs: ----------------------------------------------------------------------------
	# TEMP - experiment temperature (K)
	# testf - flag for operating in normal mode (0) or testing mode (1)
	# chem_scheme_markers - markers for different sections of the chemical scheme
	# photo_par_file - path (from PyCHAM home directory) to file containing photolysis
	#					information (absorption cross sections and quantum yields)
	# ------------------------------------------------------------------------------------
	
	if testf == 1: # for just testing mode
		return(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    
    # calculate gas-phase concentrations of M, N2 and O2 (molecules/cc (air))
	# 1.0e-6 converts from molecules/m3 to molecules/cc
	# R and Avogadro's constant set the same as in atmosphereFunctions.f90 of AtChem2
	M_val = (PInit/(8.3144621*TEMP)*6.02214129e+23)*1.0e-6
	# N2 and O2 given the same multiplication as in atmosphereFunctions.f90 of AtChem2
	N2_val = M_val*0.7809
	O2_val = M_val*0.2095
	# water concentration (C_H2O) (molecules/cc (air))
	[C_H2O, Psat_water, H2O_mw] = water_calc(TEMP, RH, 6.02214129e+23)
    
	# Open the file
	f_open_eqn = open(filename, mode='r')
	
	# read the file and store everything into a list
	# reaction rates have units /s
	total_list_eqn = f_open_eqn.readlines()
	f_open_eqn.close()
	if (f_open_eqn.closed == False):
		print('IOError')
		print('Eqn file not closed')
		sys.exit()
	
	naked_list_eqn = [] # empty list for equation reactions
	RO2_names = [] # empty list for peroxy radicals
	rrc = [] # empty list for reaction rate coefficients
	rrc_name = [] # empty list for reaction rate coefficient labels

	eqn_flag = 0 # don't collate reaction equations until seen
	pr_flag = 0 # don't collate peroxy radicals until seen
	rrc_flag = 0 # don't collate reaction rate coefficients until seen
	
	# obtain lists for reaction rate coefficients, peroxy radicals and equation reactions
	# using chemical scheme input file separators
	for line in total_list_eqn:
		
		line1 = line.strip() # remove bounding white space
		
		# --------------------------------------------------------------------------------
		# generic reaction rate coefficients part
		# look out for start of coefficients
		if line1 == str(chem_scheme_markers[3]):
			rrc_flag = 1 # ready to collate reaction rate coefficients
		matchObj = re.match(chem_scheme_markers[5], line1) # look out for end of coefficients
		if matchObj:
			rrc_flag = 0 # stop appending reaction equations
		# collate reaction rate coefficients
# 		if (rrc_flag==1 and re.match( r"(.*) = (.*)", line1)!=None): 
		if (rrc_flag==1 and re.match(str('(.*)' +' = ' + '(.*)'), line1)!=None): 
			# remove end characters
			line2 = line1.replace(str(chem_scheme_markers[4]), '')
			# remove all white space
			line2 = line2.replace(' ', '')
			# convert fortran-type scientific notation to python type
			line2 = formatting.SN_conversion(line2)
			# ensure rate coefficient is python readable
			line2 = formatting.convert_rate_mcm(line2)
			rrc.append(line2.strip())
			# get just name of reaction rate coefficient
			rrc_name.append((line2.split('=')[0]).strip())
		
		# --------------------------------------------------------------------------------
		# peroxy radical part
		# now start logging peroxy radicals
		if (re.match(chem_scheme_markers[6], line1) != None):
			pr_flag = 1
		if line1 == chem_scheme_markers[8]: # no longer log peroxy radicals
			pr_flag=0
		if (pr_flag==1):
			line2 = line1.split(chem_scheme_markers[7])
			
			for line3 in line2: # loop through elements in line
				if len(line3.split('='))>1: # in case of RO2 = ...
					line3 = (line3.split('='))[1]
				if len(line3.split(';'))>1: # in case of RO2 list finishing with ...;
					line3 = (line3.split(';'))[0]
				if (line3.strip() == '' or line3.strip() == chem_scheme_markers[6]):
					continue
				else:
					RO2_names.append(line3.strip())
		# --------------------------------------------------------------------------------
		# reaction equation part
		if line1 == chem_scheme_markers[0]: # look out for start of equations
			eqn_flag = 1 # ready to compile equations
		if (eqn_flag==1 and re.match(chem_scheme_markers[1], line1) != None):
			naked_list_eqn.append(line1) # begin storing reaction equations
			
		matchObj = re.match(chem_scheme_markers[2], line1) # look out for end of equations
		if matchObj:
			eqn_flag = 0 # stop appending reaction equations
			
		# --------------------------------------------------------------------------------
	
	# format the equation list

	# naked_list_eqn contains everything except the comments starting with //
	print('Now parsing the eqn info...\n')
	
	# first loop through equation list to concatenate any multiple lines that hold just
	# one equation, note this depends on a symbol representing the start and end of
	# and equation (inside MCM this is % for start and ; for end)
	for iline in range (0, len(naked_list_eqn)):
		
		# keep track of when we reach end of naked_list_eqn (required as naked_list_eqn
		# may shorten)
		if (iline+1)==len(naked_list_eqn):
			break
			
		# if one equation already on one line
		if (re.search(r'%', naked_list_eqn[iline])!=None and re.search(r';', naked_list_eqn[iline])!=None):
			continue
		# if spread across more than one line
		else:
			con_count = 1 # count number of lines that need concatenating
			iline2 = iline # index for further lines
			while con_count!=0:
				iline2 += 1 # move onto next line
				naked_list_eqn[iline] = str(naked_list_eqn[iline]+naked_list_eqn[iline2])
				if re.search(r';', naked_list_eqn[iline2])!=None:
					# remove concatenated line(s) from original list
					naked_list_eqn = naked_list_eqn[0:iline]+naked_list_eqn[iline2+1::]
					con_count = 0 # finish concatenating
	
	num_eqn = len(naked_list_eqn) # get number of equations
	
	
	# --open and initialise the xml file for converting chemical names to SMILES-----
	with open(xmlname) as fd:
		doc = xmltodict.parse(fd.read())

	a = doc['mechanism']['species_defs']['species']
	spec_numb = list(('0',) * len(a))
	spec_name = list(('0',) * len(a))
	spec_smil = list(('0',) * len(a))
	
	for i in range(len(a)):
		spec_numb[i] = a[i]['@species_number']
		spec_name[i] = a[i]['@species_name']
		if "smiles" in a[i]:
			spec_smil[i] = a[i]['smiles']
		elif spec_name[i][0]=='O' or spec_name[i][0]=='H':
			 spec_smil[i] = '['+spec_name[i]+']'
		else:
			 spec_smil[i] = spec_name[i] 
    
	species_step = 0 # log the number of unique species
	max_no_reac = 0.0 # log maximum number of reactants in a reaction
	max_no_prod = 0.0 # log maximum number of products in a reaction
	
	species_step = 0 # ready for equation loop
	
	# initialising lists
	
	# matrix to record indices of reactants (cols) in each equation (rows)
	rindx = np.zeros((num_eqn, 1)).astype(int)
	# matrix to record indices of products (cols) in each equation (rows)
	pindx = np.zeros((num_eqn, 1)).astype(int)
	# matrix to record stoichometries of reactants (cols) in each equation (rows)
	rstoi = np.zeros((num_eqn, 1))
	# matrix to record stoichometries of products (cols) in each equation (rows)
	pstoi = np.zeros((num_eqn, 1))
	# array to store number of reactants and products in an equation
	nreac = np.empty(num_eqn, dtype=np.int8)
	nprod = np.empty(num_eqn, dtype=np.int8)
	# list for equation reaction rate coefficients
	reac_coef = []
	# list for components' SMILE strings
	spec_list = []
	# list of Pybel objects
	Pybel_objects = []
	# a new list for the name strings of species presenting in the scheme (not SMILES)
	spec_namelist = []
	
	# Loop through the equations line by line and extract the information
	for eqn_step in range(num_eqn):
	
		line = naked_list_eqn[eqn_step] # extract this line
		
		# split the line into 2 parts: equation; rate coef 
		# (fac format doesnt have id for each equation)
		# extract the equation (in a string)
		eqn_regex = r"\:.*\;" # eqn starts with a : and end with a ;
		eqn = re.findall(eqn_regex, line)[0][1:-1].strip()
		
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
			rstoi = np.append(rstoi, (np.zeros((num_eqn, 1))), axis=1)
		while max_no_prod > np.minimum(pindx.shape[1], pstoi.shape[1]): 
			pindx = np.append(pindx, (np.zeros((num_eqn, 1))).astype(int), axis=1)
			pstoi = np.append(pstoi, (np.zeros((num_eqn, 1))), axis=1)

		# extract the reaction rate constant (in a string)
		rate_regex = r"\%.*\:" # rate coef starts with a % and end with a :
		# rate_ex: rate coefficient expression in a string
		rate_ex = re.findall(rate_regex, line)[0][1:-1].strip()
		# convert fortran-type scientific notation to python type
		rate_ex = formatting.SN_conversion(rate_ex)
		# convert the rate coefficient expressions into Python readable commands
		rate_ex = formatting.convert_rate_mcm(rate_ex)
		if (rate_ex.find('EXP') != -1):
			print(rate_ex)
			sys.exit()
		
		# store the reaction rate for this equation (/s once any inputs applied)
		reac_coef.append(rate_ex)
		
		# extract the stoichiometric number of the specii in current equation
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
			
			# store stoichometry
			rstoi[eqn_step, reactant_step] = stoich_num
			
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
				name_indx = species_step # allocate index to this species
				# Generate pybel
				Pybel_object = pybel.readstring('smi', name_SMILE)
				# append to Pybel object list
				Pybel_objects.append(Pybel_object)
				
				species_step += 1 # number of unique species
				

			else: # if it's a species already encountered it will be in spec_list
				# existing index
				name_indx = spec_namelist.index(name_only)
			
			# store reactant index
			# check if index already present - i.e. component appears more than once
			if sum(rindx[eqn_step, 0:reactant_step]==int(name_indx))>0:
				# get pre-existing index of this component
				exist_indx = np.where(rindx[eqn_step, 0:reactant_step]==(int(name_indx)))
				# add to pre-existing stoichometry
				rstoi[eqn_step, exist_indx] += rstoi[eqn_step, product_step]
				rstoi[eqn_step, reactant_step] = 0 # remove stoichometry added above
				reactant_step -= 1 # ignore this duplicate product
			else:
				rindx[eqn_step, reactant_step] = int(name_indx)

			reactant_step += 1
			
		# number of reactants in this equation
		nreac[eqn_step] = int(reactant_step)
		
		# right hand side of equations (gains)
		for product in products:

			if (re.findall(stoich_regex, product)[0] != ''):
				stoich_num = float(re.findall(stoich_regex, product)[0])
				name_only = re.sub(stoich_regex, '', product) # name with no stoich number

			elif (re.findall(stoich_regex, product)[0] == ''):
				stoich_num = 1.0
				name_only = product
			
			# store stoichometry
			pstoi[eqn_step, product_step] = stoich_num
			
			if name_only not in spec_namelist: # if new component encountered
				spec_namelist.append(name_only)
				
				# convert MCM chemical names to SMILES
				# index where xml file name matches reaction component name
				if name_only in spec_name:
					name_indx = spec_name.index(name_only)
					name_SMILE = spec_smil[name_indx]
				else:
					sys.exit(str('Error: inside eqn_parser, chemical scheme name '+str(spec_name)+' not found in xml file'))
				
				spec_list.append(name_SMILE) # list SMILE string of parsed species
				name_indx = species_step # allocate index to this species
				# Generate pybel
				
				Pybel_object = pybel.readstring('smi', name_SMILE)
				# append to Pybel object list
				Pybel_objects.append(Pybel_object)
				
				species_step += 1 # number of unique species
				
			

			else: # if it's a species already encountered
				# index of component already listed
				name_indx = spec_namelist.index(name_only)
				
			# store product index
			# check if index already present - i.e. component appears more than once
			if sum(pindx[eqn_step, 0:product_step]==int(name_indx))>0:
				exist_indx = np.where(pindx[eqn_step, 0:product_step]==(int(name_indx))) # get pre-existing index of this component
				# add to pre-existing stoichometry
				pstoi[eqn_step, exist_indx] += pstoi[eqn_step, product_step]
				pstoi[eqn_step, product_step] = 0 # remove stoichometry added above
				product_step -= 1 # ignore this duplicate product
			else:
				pindx[eqn_step, product_step] = int(name_indx)
			product_step += 1
		
		# number of products in this equation
		nprod[eqn_step] = int(product_step)
		
	if len(spec_list)!=len(spec_namelist):
		sys.exit('Error: inside eqn_parser, length of spec_list is different to length of spec_namelist and the SMILES in the former should align with the chemical scheme names in the latter')	
	
	# number of columns in rindx and pindx
	reacn = rindx.shape[1]
	prodn = pindx.shape[1]  
	
	# create a 2 column array, the first column with the RO2 list index of any RO2 species
	# that appear in the species list, the second column for its index in the species list
	RO2_indices = write_RO2_indices(spec_namelist, RO2_names)
	
	# automatically generate the Rate_coeffs module that will allow rate coefficients to
	# be calculated inside ode_gen module
	# now create reaction rate file (reaction rates are set up to have units /s)
	write_rate_file(reac_coef, rrc, rrc_name, M_val, N2_val, O2_val, TEMP, C_H2O, testf)

	# number of photolysis reactions, if this relevant
	cwd = os.getcwd() # address of current working directory
	if photo_par_file == str(cwd + '/PyCHAM/photofiles/MCMv3.2'):
		Jlen = 62 # for MCM (default name of photolysis parameters)
	else: # need to find out number of photolysis reactions
		# use Fortran indexing to be consistent with MCM photochemical reaction numbers
		Jlen = 1 
		# open file to read
		f = open(str(photo_par_file), 'r')
		for line in f: # loop through line
			if line.strip() == str('J_'+str(Jlen) + '_axs'):
				Jlen += 1

	# print the brief info for the simulation to the screen
	print('Briefing:')
	print('Total number of equations: %i' %(num_eqn))
	print('Total number of species found in chemical scheme file: %i\n' %(species_step))
	
# 	print(rindx[940,:], 'whoop',pindx[940,:])
# 	print(spec_namelist[312])
# 	print(spec_namelist[1])
# 	count = 0
# 	for i in spec_namelist:
# 		if i == "O3":
# 			print(count, i)
# 		count+=1
# 	ipdb.set_trace()
	# outputs: 
	
	# rindx  - matrix to record indices of reactants (cols) in each equation (rows)
	# pindx - indices of equation products (cols) in each equation (rows)
	# rstoi - matrix to record stoichometries of reactants (cols) in each equation (rows)
	# pstoi - matrix to record stoichometries of products (cols) in each equation (rows)
	# reac_coef - list for equation reaction rate coefficients
	# spec_list - list for components' SMILE strings
	# Pybel_objects - list of Pybel objects
	# species_step - number of species
	# num_eqn - number of equations
	# nreac - number of reactants in each equation
	# max_no_jaci - number of columns for Jacobian index matrix
	# nprod - number of products per equation
	# prodn - number of columns in pindx
	# reacn - rindx number of columns
	# M_val - gas-phase concentration of M (molecules/cc (air))
	# N2_val - gas-phase concentration of nitrogen (molecules/cc (air))
	# O2_val - gas-phase concentration of oxygen (molecules/cc (air))
	# spec_namelist - list of component names used in the chemical reaction file
	
	return (rindx, pindx, rstoi, pstoi, reac_coef, spec_list, Pybel_objects, num_eqn, 
			species_step, RO2_indices, nreac,
			nprod, prodn, reacn, M_val, N2_val, O2_val, C_H2O, 
			Psat_water, H2O_mw, spec_namelist, Jlen)



# This function generates a python script that calculate rate coef. numerically
# main part by Dave (/s)
def write_rate_file(reac_coef, rrc, rrc_name, M, N2, O2, TEMP, C_H2O, testf):
	if testf==0:
		f = open('PyCHAM/Rate_coeffs.py', mode='w')
	if testf==2:
		f = open('Rate_coeffs.py', mode='w')
	f.write('\'\'\'module for calculating gas-phase reaction rate coefficients, automatically generated by eqn_parser\'\'\'\n')
	f.write('\n')
	f.write('##################################################################################################### \n') # python will convert \n to os.linesep
	f.write('# Python function to hold expressions for calculating rate coefficients for a given equation number # \n') # python will convert \n to os.linesep
	f.write('#    Copyright (C) 2017  David Topping : david.topping@manchester.ac.uk                             # \n')
	f.write('#                                      : davetopp80@gmail.com                                       # \n')
	f.write('#    Personal website: davetoppingsci.com                                                           # \n')
	f.write('#                                                                                                   # \n')
	f.write('#                                                                                                   # \n')
	f.write('#                                                                                                   # \n')
	f.write('##################################################################################################### \n')    
	f.write('# Minor modified by XSX\n')
	f.write('# File Created at %s\n' %(datetime.datetime.now()))
	f.write('\n')
	f.write('import numpy\n')
	f.write('import PhotolysisRates\n')
	f.write('\n')

	# following part is the function (there should be an indent at the start of each line)
	# suggest using 1 Tab
	f.write('def evaluate_rates(RO2, H2O, TEMP, lightm, time, lat, lon, act_flux_path, DayOfYear, M, N2, O2, photo_par_file, Jlen):\n')
	f.write('\n')
	f.write('	# ------------------------------------------------------------------------')
	f.write('	# inputs:\n')
	f.write('	# M - third body concentration (molecules/cc (air))\n')
	f.write('	# N2 - nitrogen concentration (molecules/cc (air))\n')
	f.write('	# O2 - oxygen concentration (molecules/cc (air))\n')
	f.write('	# RO2: specified by the chemical scheme. eg: subset of MCM\n')
	f.write('	# H2O, TEMP: given by the user\n')
	f.write('	# lightm: given by the user and is 0 for lights off and 1 for on\n')
	f.write('	# reaction rate coefficients and their names parsed in eqn_parser.py \n')
	f.write('	# Jlen - number of photolysis reactions')
	f.write('\n')
	f.write('	# calculate reaction rates with given by chemical scheme\n')
	# code to calculate rate coefficients given by chemical scheme file
	for line in rrc:
		f.write('	%s \n' %line)
	f.write('\n')
	f.write('	# estimate and append photolysis rates\n')
	f.write('	J = PhotolysisRates.PhotolysisCalculation(time, lat, lon, TEMP, act_flux_path, DayOfYear, photo_par_file, Jlen)\n')
	f.write('\n')
	f.write('	if lightm == 0:\n')
	f.write('		J = [0]*len(J)\n')

	# calculate the rate coef. numerically for each equation
	f.write('	rate_values = numpy.zeros(%i)\n' %(len(reac_coef)))
	# BE NOTIFIED!!!: before writing the script, 'reac_coef' must be converted to 
	# python-compatible format
	f.write('	# reac_coef has been formatted so that python can recognize it\n')
	for eqn_key in range (len(reac_coef)):
		f.write('	rate_values[%s] = %s\n' %(eqn_key, reac_coef[eqn_key]))
	f.write('	\n')
	f.write('	return rate_values\n')
	f.close()

# function to automatically generate a module that is used to record the tendency
# of components (components specified by the user, their index given by rec_comp_index) 
# to change in response to box model 
# mechanisms - gets called inside ode_gen on each time step
def write_dydt_rec():

	f = open('PyCHAM/dydt_rec.py', mode='w')
	f.write('\'\'\'module for calculating and recording change tendency of components, automatically generated by eqn_parser\'\'\'\n')
	f.write('\n')
	f.write('# File Created at %s\n' %(datetime.datetime.now()))
	f.write('\n')
	f.write('import numpy as np \n')
	f.write('\n')
	# following part is the function (there should be an indent at the start of each line)
	# suggest 1 Tab
	f.write('def dydt_rec(y, rindx, rstoi, reac_coef, pindx, nprod, step, dydt_vst, nreac, num_sb, num_speci, pconc, core_diss, Psat, kelv_fac, kimt, kwgt, Cw, act_coeff):\n')
	f.write('	# loop through components to record the tendency of change \n')
	f.write('	for compi in dydt_vst.get(\'comp_index\'): \n')
	f.write('		# open relevant dictionary value \n')
	f.write('		dydt_rec = dydt_vst.get(compi) \n')
	f.write('		# keep count on relevant reactions \n')
	f.write('		reac_count = 0 \n')
	f.write('		# loop through relevant reactions \n')
	f.write('		for i in dydt_rec[0,0:-2]: # final two rows for particle- and wall-partitioning \n')
	f.write('			i = int(i) # ensure index is integer # this necessary because the dydt_rec array is float (the tendency to change records beneath its first row are float) \n')
	f.write('			# estimate gas-phase change tendency for every reaction involving this component \n')
	f.write('			gprate = ((y[rindx[i, 0:nreac[i]]]**rstoi[i, 0:nreac[i]]).prod())*reac_coef[i] \n')
	f.write('			# identify whether this component reacted or produced\n')
	f.write('			if sum(rindx[i, 0:nreac[i]]==compi)>0: \n')
	f.write('				dydt_rec[step+1, reac_count] -= ((gprate)) #*3600)/np.abs(y[compi]))*100.0 \n')
	f.write('			if sum(pindx[i, 0:nprod[i]]==compi)>0: \n')
	f.write('				dydt_rec[step+1, reac_count] += ((gprate)) #*3600)/np.abs(y[compi]))*100.0 \n')
	f.write('			reac_count += 1 \n')
	f.write('		# now estimate and record tendency to change due to particle- and wall-partitioning  \n')
	f.write('		# particle-partitioning \n')
	f.write('		for ibin in range(num_sb-1): # size bin loop\n')
	f.write('			Csit = y[num_speci*(ibin+1):num_speci*(ibin+2)]\n')
	f.write('			conc_sum = np.zeros((1)) \n')
	f.write('			if pconc>0.0: # if seed particles present \n')
	f.write('				conc_sum[0] = ((Csit[0:-1].sum())+Csit[-1]*core_diss)\n')
	f.write('			else: \n')
	f.write('				conc_sum[0] = Csit.sum() \n')
	f.write('			# prevent numerical error due to division by zero \n')
	f.write('			ish = conc_sum==0.0 \n')
	f.write('			conc_sum[ish] = 1.0e-40 \n')
	f.write('			# particle surface gas-phase concentration (molecules/cc (air)) \n')
	f.write('			Csit = (Csit[compi]/conc_sum)*Psat[compi, 0]*kelv_fac[ibin]*act_coeff[compi, 0] \n')
	f.write('			# partitioning rate (molecules/cc.s) \n')
	f.write('			dydt_all = kimt[compi, ibin]*(y[compi]-Csit) \n')
	f.write('			# gas-phase change (molecules/cc/s) \n')
	f.write('			dydt_rec[step+1, reac_count] -= dydt_all \n')
	f.write('		# wall-partitioning \n')
	f.write('		if (kwgt)>1.0e-10: \n')
	f.write('			# concentration at wall (molecules/cc (air)) \n')
	f.write('			Csit = y[num_speci*num_sb:num_speci*(num_sb+1)] \n')
	f.write('			Csit = (Psat[:,0]*(Csit/Cw)*act_coeff[compi, 0])\n')
	f.write('			dydt_all = (kwgt)*(y[compi]-Csit[compi]) \n')
	f.write('			# gas-phase change (molecules/cc/s) \n')
	f.write('			dydt_rec[step+1, reac_count+1] -= dydt_all \n')
	f.write('		\n')
	f.write('	return(dydt_vst) \n')
	f.close()				
						
	
# This function defines RO2 which is given by an MCM file
# this function is used when certain reaction rate coefficients are a function a RO2
# and is called on by the extract_mechanism function above,
# whilst the resulting RO2 index is used inside rate_valu_calc.py
def write_RO2_indices(smiles_array, RO2_names):
    
    # store the names of RO2 species which are present in the equation file
    # get a list of INDICES of RO2 that present in the equation file 
    # (or total species dict)
    # empty list for RO2_indices
    RO2_indices0 = []
    RO2_indices = []
    
    for name in RO2_names:
        
        if name in smiles_array:
            # get the RO2 index
            index0 = RO2_names.index(name)
            RO2_indices0.append(index0)
            # get the smiles_array index for this RO2 species
            index1 = smiles_array.index(name)
            RO2_indices.append(index1)
    
    # Ensure elements in RO2_indices are int (iterable)
    RO2_indices0 = (np.asarray(RO2_indices0, dtype=int)).reshape(-1, 1)
    RO2_indices = (np.asarray(RO2_indices, dtype=int)).reshape(-1, 1)
    RO2_indices = np.hstack((RO2_indices0, RO2_indices))
    
    return (RO2_indices)

# function to convert reaction rate coefficients to commands, and therefore quantities
# that will be used by Rate_coeffs.py
def write_rrc(rrc, rrc_name):

	f = open('PyCHAM/rrc_calc.py', mode='w')
	f.write('\'\'\'module for calculating reaction rate coefficients, automatically generated by eqn_parser\'\'\'\n')
	f.write('\n')
	f.write('# File Created at %s\n' %(datetime.datetime.now()))
	f.write('\n')
	f.write('import numpy as numpy \n')
	f.write('import PhotolysisRates \n')
	f.write('\n')
	# following part is the function (there should be an indent at the start of each line)
	# suggest 1 Tab
	f.write('def rrc_calc(rrc_name, TEMP, H2O, M, N2, O2, time, lat, lon, act_flux_path, DayOfYear):\n')
	f.write('	\n')
	f.write('	# reaction rate coefficients obtained from user\'s chemical scheme file \n')
	f.write('	rrc_constants = [] # empty list to hold rate coefficients\n')
	for line in rrc:
		f.write('	%s \n' %line)
	f.write('	\n')
	f.write('	for i in rrc_name: \n')
	f.write('		rrc_constants.append(locals()[i]) \n')
	f.write('	\n')
	f.write('	# estimate and append photolysis rates \n')
	f.write('	j = PhotolysisRates.PhotolysisCalculation(time, lat, lon, TEMP, act_flux_path, DayOfYear) \n')
	f.write('	rrc_constants.append(j) \n')
	f.write('	# append photolysis rate label \n')
	f.write('	rrc_name.append(\'J\') \n')
	f.write('	return(rrc_constants, rrc_name) \n')
	f.close()
	
	return