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

# ----------Extraction of eqn info----------
# Extract the mechanism information
def extract_mechanism(filename, xmlname, TEMP, PInit, Comp0, testf):

	# inputs:
	# testf - flag for operating in normal mode (0) or testing mode (1)
	
	if testf == 1:
		return(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    	
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

	# equation list without comments
	naked_list_eqn = formatting.remove_comments(total_list_eqn)
	
	# calculate gas-phase concentrations of M, N2 and O2 (molecules/cc (air))
	# 1.0e-6 converts from molecules/m3 to molecules/cc
	M_val = (PInit/(8.3144598*TEMP)*6.0221409e+23)*1.0e-6
	N2_val = M_val*0.79
	O2_val = M_val*0.2096
			
	# format the equation list
	
	# naked_list_eqn contains everything except the comments starting with //
	print('Now parsing the eqn info...\n')
	num_eqn = len(naked_list_eqn)
	
	
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
	
	# convert input component names for components present in gas phase at experiment
	# start from chemical scheme names to SMILES
	init_SMIL = []
	
	for species_step in range(len(Comp0)): 
		name_indx = spec_name.index(Comp0[species_step])
		init_SMIL.append(spec_smil[name_indx])
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
	nreac = np.zeros((num_eqn))
	nprod = np.zeros((num_eqn))
	# list for equation reaction rate coefficients
	reac_coef = []
	# list for species
	spec_list = []
	# list of Pybel objects
	Pybel_objects = []
	# a new list for the name strings of species presenting in the scheme (not SMILES)
	spec_namelist = []
	
	# Loop through the equations line by line and extract the information
	for eqn_step in range(num_eqn):
	
		line = naked_list_eqn[eqn_step]
		
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

		# extract the rate constant (in a string)
		rate_regex = r"\%.*\:" # rate coef starts with a : and end with a ；
		# rate_ex: rate coefficient expression in a string
		rate_ex = re.findall(rate_regex, line)[0][1:-1].strip()
		# convert fortran-type scientific notation to python type
		rate_ex = formatting.SN_conversion(rate_ex)
		# convert the rate coefficient expressions into Python readable commands
		rate_ex = formatting.convert_rate_mcm(rate_ex)
		
		# store the reaction rate for this equation (/s once any inputs applied)
		reac_coef.append(rate_ex)
		
		# extract the stoichiometric number of the specii in current equation
		reactant_step = 0
		product_step = 0
		stoich_regex = r"^\d*\.\d*|^\d*"
		numr = len(reactants) # number of reactants in this equation
		
		
		# left hand side of equations (losses)
		for reactant in reactants:
		
			if reactant not in spec_namelist:
				spec_namelist.append(reactant)
			if (re.findall(stoich_regex, reactant)[0] != ''):
				stoich_num = float(re.findall(stoich_regex, reactant)[0])
				name_only = re.sub(stoich_regex, '', reactant) # name with no stoich number
			elif (re.findall(stoich_regex, reactant)[0] == ''):
				stoich_num = 1.0
				name_only = reactant
			
			# store stoichometry
			rstoi[eqn_step, reactant_step] = stoich_num

			# convert MCM chemical names to SMILES
			# index where xml file MCM name matches MCM name
			if name_only in spec_name:
				
				name_indx = spec_name.index(name_only)
				name_only = spec_smil[name_indx]
			
			if (name_only not in spec_list):
				spec_list.append(name_only) # log parsed species
				name_indx = species_step # allocate index to this species
				# Generate pybel
				Pybel_object = pybel.readstring('smi', name_only)
				# append to Pybel object list
				Pybel_objects.append(Pybel_object)
				
				species_step += 1 # number of unique species
				

			else: # if it's a species already encountered
				# pre-defined number of species
				name_indx = spec_list.index(name_only)

			# store reactant index
			rindx[eqn_step, reactant_step] = int(name_indx)

			reactant_step += 1
			
		# number of reactants in this equation
		nreac[eqn_step] = reactant_step
		
		# right hand side of equations (gains)
		for product in products:
			if product not in spec_namelist:
				spec_namelist.append(product)
			if (re.findall(stoich_regex, product)[0] != ''):
				stoich_num = float(re.findall(stoich_regex, product)[0])
				name_only = re.sub(stoich_regex, '', product) # name with no stoich number
			elif (re.findall(stoich_regex, product)[0] == ''):
				stoich_num = 1.0
				name_only = product
			
			# store stoichometry
			pstoi[eqn_step, product_step] = stoich_num
						
			# convert MCM chemical names to SMILES
			# index where xml file MCM name matches MCM name
			if name_only in spec_name:
				
				name_indx = spec_name.index(name_only)
				name_only = spec_smil[name_indx]
			
			if (name_only not in spec_list):
				spec_list.append(name_only) # log parsed species
				name_indx = species_step # allocate index to this species
				# Generate pybel
				
				Pybel_object = pybel.readstring('smi', name_only)
				# append to Pybel object list
				Pybel_objects.append(Pybel_object)
				
				species_step += 1 # number of unique species
				

			else: # if it's a species already encountered
				# pre-defined number of species
				name_indx = spec_list.index(name_only)

			# store product index
			pindx[eqn_step, product_step] = int(name_indx)
			
			product_step += 1
			
		# number of products in this equation
		nprod[eqn_step] = product_step
		
	# number of columns in rindx and pindx
	reacn = rindx.shape[1]
	prodn = pindx.shape[1]  
	
	# create a 2 column array, the first column with the RO2 list index of any RO2 species
	# that appear in the species list, the second column for its index in the species list
	RO2_indices = write_RO2_indices(spec_namelist)
	
	# print the brief info for the simulation to the screen
	print('Briefing:')
	print('Total number of equations: %i' %(num_eqn))
	print('Total number of species: %i\n' %(species_step))

	# outputs: 
	
	# rindx  - matrix to record indices of reactants (cols) in each equation (rows)
	# pindx - indices of equation products (cols) in each equation (rows)
	# rstoi - matrix to record stoichometries of reactants (cols) in each equation (rows)
	# pstoi - matrix to record stoichometries of products (cols) in each equation (rows)
	# reac_coef - list for equation reaction rate coefficients
	# spec_list - list for species
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
	# init_SMIL - SMILE string for each component
	# spec_namelist - names of components as given in the equation file
	
	return (rindx, pindx, rstoi, pstoi, reac_coef, spec_list, Pybel_objects, num_eqn, 
			species_step, RO2_indices, nreac,
			nprod, prodn, reacn, M_val, N2_val, O2_val, init_SMIL, spec_namelist)



# This function generates a python script that calculate rate coef. numerically
# main part by Dave (/s)
def write_rate_file(filename, reac_coef, mcm_constants, MCMConstNameList, testf, M, N2, 
					O2):
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
    f.write('# File Created as %s\n' %(datetime.datetime.now()))
    f.write('\n')
    f.write('import numpy\n')
    f.write('\n')

    # following part is the function (there should be an indent at the start of each line)
    # suggest using 4 SPACES instead of 1 Tab
    f.write('def evaluate_rates(RO2, H2O, TEMP, lightm):\n')
    f.write('    # mcm_constants_dict: given by mcm_constants.py\n')
    f.write('    # RO2: specified by the chemical scheme. eg: subset of MCM\n')
    f.write('    # H2O, TEMP: given by the user\n')
    f.write('    # lightm: given by the user and is 0 for lights off and 1 for on\n')
    f.write('    # Creating reference to constant values used in rate expressions\n')
    f.write('    # mcm_constants_dict: given by MCM_constants.py\n')
    # retrive generic rate coef calculated by MCM_constants.py
    for rate_key in range (len(mcm_constants)):
        f.write('    %s = %s \n' %(MCMConstNameList[rate_key], mcm_constants[rate_key]))
    f.write('    \n')
    f.write('    if lightm == 0:\n')
    f.write('    	J = [0]*len(J)\n')
    f.write('    # Environmental Variables: M, O2, N2\n')
    f.write('    M = ' + str(M) + ' # 3rd body; number of molecules in per unit volume\n')
    f.write('    N2 = ' + str(N2) + ' # Nitrogen mass mixing ratio : 79%\n')
    f.write('    O2 = ' + str(O2) + ' # Oxygen mass mixing ratio : 20.96%\n')
    
    # calculate the rate coef. numerically for each equation
    f.write('    rate_values = numpy.zeros(%i)\n' %(len(reac_coef)))
    # BE NOTIFIED!!!: before writing the script, 'reac_coef' must be converted to 
    # python-compatible format
    f.write('    # reac_coef has been formatted so that python can recognize it\n')
    for eqn_key in range (len(reac_coef)):
        f.write('    rate_values[%s] = %s\n' %(eqn_key, reac_coef[eqn_key]))
    
    f.write('\n')
    f.write('    return rate_values\n')
    f.close()

# This function defines RO2 which is given by an MCM file
# this function is used when certain rate coef. is a function a RO2
# species_dict2array = {pybel_object: species_num}

def write_RO2_indices(smiles_array):
    # the list below is an universal list or a full list of RO2. 
    # Only part of it might be used in real calculation.
#     RO2_names = ['NBUTOLAO2','HO3C4O2','BU1ENO3O2','C43NO34O2','BZBIPERO2','CH3O2','C2H5O2','HOCH2CH2O2',
#     'ETHENO3O2','C6H5C2H4O2','EBZBIPERO2','ISOPAO2','ISOPBO2','ISOPCO2','ISOPDO2',
#     'NISOPO2','CH3CO3','C2H5CO3','NC3H7O2','IC3H7O2','HYPROPO2','IPROPOLO2','PRONO3BO2',
#     'PRONO3AO2','C6H5CH2O2','TLBIPERO2','CH3COCH2O2','BUT2OLO2','C42NO33O2','IC4H9O2','TC4H9O2','IBUTOLBO2','TBUTOLO2',
#     'MPRANO3O2','MPRBNO3O2','IPEAO2','IPEBO2','IPECO2','MXYBIPERO2','MXYLO2','MEKAO2',
#     'MEKCO2','MEKBO2','NC4H9O2','SC4H9O2','HEXAO2','HEXBO2','HEXCO2','PEAO2','PEBO2','PECO2','OXYBIPERO2','OXYLO2',
#     'PXYBIPERO2','PXYLO2','MPRKAO2','CO2C54O2','HO2C5O2','DIEKAO2','DIEKBO2','BZEMUCO2',
#     'BZEMUCCO3','C5DIALO2','PHENO2','NPHENO2','EBZMUCO2','EBZMUCCO3','C715CO2O2','EBENZOLO2','NEBNZOLO2','HCOCO3','HMVKAO2',
#     'HMVKBO2','MVKO2','MACO3','MACRO2','TLEMUCO2','TLEMUCCO3','C615CO2O2','CRESO2',
#     'NCRESO2','MXYMUCO2','MXYMUCCO3','C726CO5O2','MXYOLO2','NMXYOLO2','OXYMUCO2',
#     'OXYMUCCO3','MC6CO2O2','OXYOLO2','NOXYOLO2','PXYMUCO2','PXYMUCCO3','C6M5CO2O2','PXYOLO2','NPXYOLO2','HO3C3CO3',
#     'CO3C4NO3O2','MALDIALCO3','EPXDLCO3','C3DIALO2','MALDIALO2','HOCH2CO3',
#     'NO3CH2CO3','C6H5CH2CO3','C6DCARBBO2','C58O2','HC4ACO3','HC4CCO3','C57O2',
#     'C59O2','NC4CO3','C510O2','HO1C3O2','CH3CHOHCO3','PRNO3CO3','C6H5CO3','C5CO14O2',
#     'IBUTOLCO2','IBUTALBO2','IBUTALCO2','IPRCO3','IPRHOCO3','MPRBNO3CO3','M2BUOL2O2',
#     'HM2C43O2','BUT2CO3','C52O2','ME2BUOLO2','H2M3C4O2','MIPKAO2','MIPKBO2','ME2BU2OLO2',
#     'PROL11MO2','HO2M2C4O2','C3MCODBCO3','MXYLCO3','EPXMDLCO3','C3MDIALO2','HO1CO3C4O2','CO2C3CO3','BIACETO2',
#     'NBUTOLBO2','BUTALO2','C3H7CO3','HO1C4O2','C5H11CO3','HO1C6O2','HEX2ONAO2','HEX2ONBO2','HEX2ONCO2','HO2C6O2','HO3C6O2',
#     'HEX3ONDO2','HEX3ONCO2','HEX3ONBO2','HEX3ONAO2','C4CHOBO2','C4H9CO3','HO1C5O2','PE2ENEBO2','HO3C5O2','OXYLCO3','EPXM2DLCO3',
#     'C4MCO2O2','PXYLCO3','CO23C54O2','HO2CO4C5O2','CO24C53O2','HO2C4CO3','HO2C4O2','HOCO3C54O2','CO3C4CO3','BZFUO2','NBZFUO2','HCOCOHCO3','EBFUO2','BUTALAO2',
#     'NEBFUO2','C6DICARBO2','C7CO3OHO2','MVKOHBO2','MVKOHAO2','CO2H3CO3','ACO3','TLFUO2','NTLFUO2','C5DICARBO2','MC3CODBCO3','C4M2ALOHO2','C4MCODBCO3','C5MCO2OHO2',
#     'MXYFUO2','C23O3MO2','NMXYFUO2','PXYFUO2','MCOCOMOXO2','NPXYFUO2','MC5CO2OHO2','MC4CODBCO3','OXYFUO2','C6OTKETO2',
#     'NOXYFUO2','DMKOHO2','C4CO2O2','C51O2','HCOCH2O2','PBZQO2','NBZQO2','NCATECO2','NNCATECO2','CO3H4CO3','PEBQO2','NPEBQO2','ENCATECO2','ENNCATECO2','HOC2H4CO3',
#     'MECOACETO2','PTLQO2','NPTLQO2','MNCATECO2','MNNCATECO2','HOIPRCO3','IBUDIALCO3','PROPALO2','HOIBUTCO3','HO2C43CO3','C56O2',
#     'C53O2','C41CO3','PROL1MCO3','H2M2C3CO3','C54O2','CHOMOHCO3','MXYQO2','NMXYQO2','MXNCATECO2','MXNNCATCO2','HO2C3CO3','HOC3H6CO3','C63O2','CO2HOC61O2','CO24C6O2',
#     'CO25C6O2','C61O2','CO23C65O2','HO3C5CO3','C6HO1CO3O2','C3COCCO3','PEN2ONE1O2','C6CO34O2','C6CO3OH5O2','HO3C4CO3','HO13C5O2','TMB1FUO2','NTMB1FUO2','OXYQO2',
#     'NOXYQO2','OXNCATECO2','OXNNCATCO2','PXYQO2','NPXYQO2','PXNCATECO2','PXNNCATCO2',
#     'CO2C4CO3','HO13C4O2','MALANHYO2','DNPHENO2','NDNPHENO2','DNEBNZLO2','NDNEBNZLO2','H13CO2CO3','DNCRESO2','NDNCRESO2',\
#     'HOBUT2CO3','ACCOCOMEO2','MMALANHYO2','DNMXYOLO2','NDNMXYOLO2','CO3C5CO3','C6O4KETO2','DNOXYOLO2','NDNOXYOLO2',
#     'TL4OHNO2O2','DNPXYOLO2','NDNPXYOLO2','CO23C4CO3','HCOCH2CO3','C5CO2OHCO3','ECO3CO3','C7OHCO2CO3','ACCOMECO3',
#     'CH3COCO3','C6CO2OHCO3','C42CO3','H13C43CO3','C4COMOHCO3','C23O3MCO3','C23O3CCO3',
#     'C7CO2OHCO3','C62O2','C5M2OHOCO3','C6MOHCOCO3','HO13C3CO3','C4CO2DBCO3','C7CO2DBCO3','C5CO2DBCO3','C5CO234O2',
#     'C5CO34CO3','C4DBM2CO3','C5DBCO2CO3','C5CO23O2','EMALANHYO2','APINOOA','APINOOB','ELVOCp']
    
    # the following RO2 names are for 'apinene_HOM_full.eqn'
    # the following RO2 names are for 'apinene_HOM_full.eqn'
    RO2_names = [
    'NAPINAO2',
    'NAPINBO2',
    'APINAO2',
    'APINBO2',
    'APINCO2',
    'C107O2',
    'C109O2',
    'C96O2',
    'NC101O2',
    'C96CO3',
    'C720O2',
    'PINALO2',
    'C108O2',
    'C89CO3',
    'C920CO3',
    'C920O2',
    'C97O2',
    'C85CO3',
    'C85O2',
    'CH3COCH2O2',
    'CH3CO3',
    'CH3O2',
    'C719O2',
    'NC102O2',
    'C106O2',
    'C717O2',
    'C811CO3',
    'C89O2',
    'C921O2',
    'C98O2',
    'C86O2',
    'CO235C6CO3',
    'CHOC3COCO3',
    'NC71O2',
    'C811O2',
    'CHOC3COO2',
    'H3C25C6CO3',
    'H3C25C6O2',
    'CO235C6O2',
    'C716O2',
    'C810O2',
    'C922O2',
    'C614O2',
    'C511O2',
    'C812O2',
    'C721CO3',
    'C721O2',
    'HCOCH2CO3',
    'BIACETO2',
    'HOCH2CO3',
    'H3C2C4CO3',
    'HMVKAO2',
    'CO23C4CO3',
    'C312COCO3',
    'CHOCOCH2O2',
    'NC72O2',
    'C514O2',
    'HCOCH2O2',
    'C621O2',
    'C813O2',
    'C722O2',
    'HCOCO3',
    'CO2H3CO3',
    'NC61CO3',
    'C516O2',
    'C44O2',
    'H1C23C4CO3',
    'H1C23C4O2',
    'APINOOA',
    'APINOOB',
    'R_HOM_O2'
    ]
    
    # store the names of RO2 species which are present in the equation file
    # get a list of INDICES of RO2 that present in the equation file (or total species dict)
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
