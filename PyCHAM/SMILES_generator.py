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
'''function to generate xml contents from chemical schemes'''
# because some chemical schemes don't come with xml files, but
# do have identifiers for components that aids SMILES generation

def xml_cont_gen(): # define function

	# import dependencies
	import re # for parsing
	import numpy as np # for math and matrix operations
	import xmltodict # for opening and converting xml files to python dictionaries
	import openbabel.pybel as pybel # for converting SMILE strings to pybel objects

	# -----------------------------------------------------------------------------------
	# tell code what method to take for estimating SMILES
	# if SMILES to be formed from number following atom letter in chemical scheme name 
	# of component set to = 1
	# if SMILES to be formed from a combination of precursors to the ROOR accretion file 
	# (RO2+RO2 = ROOR (ROOR has two oxygens less than the sum of the precursor RO2 oxygens))
	# then set to = 2
	number_of_atom_stated = 2
	# -----------------------------------------------------------------------------------

	# open chemical scheme file
	chem_sch = open('C:\\Users\\Psymo\\Desktop\\PyCHAM\\PyCHAM\\PyCHAM\\input\\CH4_AP_BZ_wautox_wbz\\AP_BZ_MCM_PRAMAP_autoAPRAMBZ_scheme.dat', mode='r')
	# read the file and store everything into a list
	total_list_eqn = chem_sch.readlines()
	chem_sch.close() # close file

	# markers to isolate sections of chemical scheme
	chem_sch_mrk = ['{', 'RO2', '+', 'C(ind_', ')', '', '&', '', '', ':', '}', ';', ''] #autoPRAM format  
	# PRAM format
	#['%', 'RO2', '+', '', '', ';', '+', ';', '', '%', ':', ';']
	
	eqn_list = [] # empty list for equations

	comp_list = [] # empty list for components

	# loop through equations
	for line in total_list_eqn:
		line1 = line.strip() # remove bounding white space

		# gas-phase reaction equation part
		# ^ means occurs at start of line and, first \ means 
		# second \ can be interpreted 
		# and second \ ensures recognition of marker
		marker = str('^\\' +  chem_sch_mrk[0])
		
		# first check is whether equation start marker is present
		if (re.match(marker, line1) != None):
			
			# second check is whether markers for starting reaction rate coefficients
			# part, and markers for end of equation lines, are present
			eqn_markers = [str('.*\\' +  chem_sch_mrk[9]), str('.*\\' +  chem_sch_mrk[11])]
			if (re.match(eqn_markers[0], line1) != None and 
				re.match(eqn_markers[1], line1) != None):
		
				eqn_list.append(line1) # store reaction equations
	
	eqn_list_firsts = [] # empty list for equations where components first occur

	for line in eqn_list: # loop through reactions
	
		# start of equation preparation part --------------------------

		# work out whether equation or reaction rate coefficient part comes first
		eqn_start = str('.*\\' +  chem_sch_mrk[10])
		rrc_start = str('.*\\' +  chem_sch_mrk[9])
		# get index of these markers, note span is the property of the match object that
		# gives the location of the marker
		eqn_start_indx = (re.match(eqn_start, line)).span()[1]
		rrc_start_indx = (re.match(rrc_start, line)).span()[1]
		
		if (eqn_start_indx > rrc_start_indx):
			eqn_sec = 1 # equation is second part
		else:
			eqn_sec = 0 # equation is first part
		
		# split the line into 2 parts: equation and rate coefficient
		# . means match with anything except a new line character., when followed by a * 
		# means match zero or more times (so now we match with all characters in the line
		# except for new line characters, so final part is stating the character(s) we 
		# are specifically looking for, \\ ensures the marker is recognised
		if eqn_sec == 1:
			eqn_markers = str('\\' +  chem_sch_mrk[10]+ '.*\\' +  chem_sch_mrk[11])
		else: # end of equation part is start of reaction rate coefficient part
			eqn_markers = str('\\' +  chem_sch_mrk[10]+ '.*\\' +  chem_sch_mrk[9])

		# extract the equation as a string ([0] extracts the equation section and 
		# [1:-1] removes the bounding markers)
		eqn = re.findall(eqn_markers, line)[0][1:-1].strip()
		
		# ensure there are spaces either side of the = sign
		if eqn[eqn.index('=')-1] != ' ':
			eqn = str(eqn[0:eqn.index('=')] + ' ' + eqn[eqn.index('=')::])
		if (len(eqn)-1>eqn.index('=')): # note that some equations do not contain reactants
			if eqn[eqn.index('=')+1] != ' ':
				eqn = str(eqn[0:eqn.index('=')+1] + ' ' + eqn[eqn.index('=')+1::])
		
		# ensure there are spaces either side of + signs
		# number of +s
		pcnt = eqn.count('+')
		plusindx = 0 # initiating index
		
		for i in range(pcnt): # loop through + signs
		
			plusindx = eqn[plusindx::].index('+') + plusindx
		
			if (eqn[plusindx-1] != ' '):
				eqn = str(eqn[0:plusindx] + ' ' + eqn[plusindx::])
				plusindx += 1 # move index up
			if (eqn[plusindx+1] != ' '):
				eqn = str(eqn[0:plusindx+1] + ' ' + eqn[plusindx+1::])
		
		
			# set new plusindx to search from
			plusindx += 1

		# end of equation preparation part ----------------------------
	
		# loop through spaces (and any tabs) in this equation
		for i in eqn.split():
			
			if i == '+' or i == '=': # continue if just symbols
				continue
			else: # check for unique components in this reaction
				if i in comp_list:
					continue
				else:
					comp_list.append(i)
					eqn_list_firsts.append(eqn) # list for equations where components first occur			
								
	# now that we have names of all components in equation file, 
	# check for their presence in xml file
	xml_name = 'C:\\Users\\Psymo\\Desktop\\PyCHAM\\PyCHAM\\PyCHAM\\input\\CH4_AP_BZ_wautox_wbz\\MCM_PRAM_xml.xml'
	
	with open(xml_name) as fd:
		doc = xmltodict.parse(fd.read())
	
	a = doc['mechanism']['species_defs']['species']
	
	fd.close()
	
	# remove last two lines from xml, note these are 
	# replaced at end of this function
	with open(xml_name, 'r') as fr:
		lines = fr.readlines()
		
		with open(xml_name, 'w') as fw:
			for line in lines:
				if '</species_defs>' in line: 
					continue
				if '</mechanism>' in line: 
					continue
				fw.write(line)
		
		fw.close()
	fr.close()

	# prepare arrays to fill	
	comp_numb = list(('0',) * len(a))
	comp_name = list(('0',) * len(a))
	comp_smil = list(('0',) * len(a))
	
	for i in range(len(a)):
		comp_numb[i] = a[i]['@species_number']
		comp_name[i] = a[i]['@species_name']
		if (comp_name[i] == ''): # nothing to register here
			continue
	
		if ("smiles" in a[i]):
			comp_smil[i] = a[i]['smiles']
		else: # if no SMILE string explicitly given
	
			succ = 0 # assume no success
			if (comp_name[i] == 'O3'):
				comp_smil[i] = '[O-][O+]=O'
				succ = 1
			if (comp_name[i] == 'NO2'):
				comp_smil[i] = '[N+](=O)[O-]'
				succ = 1
			if (comp_name[i] == 'NO3'):
				comp_smil[i] = '[N+](=O)([O-])[O]'
				succ = 1
			if (succ == 0):
				try: # first try assuming that SMILE string is represented by component name
					comp_smil[i] = comp_name[i]
					Pybel_object = pybel.readstring('smi', comp_smil[i])
				except:
					err_mess_new = str('Error: a smiles string was not found for component ' + str(comp_name[i]) + ' in the xml file, nor could its name be interpreted as a SMILE string')
					break
		
	comp_numi = -1 # keep count on components
	# loop through components in chemical scheme file to 
	# check on presence in xml
	for comp in comp_list:
	
		comp_numi += 1 # keep count on components	

		if comp in comp_name: # if present then continue
			continue
		
		if comp == 'CARENE':
				SMILES = 'CC1=CCC2CC1C2(C)C' # just bung in the same SMILES as APINENE
				# write lines for xml file
				lines = '<species species_name=\"' + comp + '\" species_number=\"n\">\n' + '<smiles>' + SMILES + '</smiles>\n' + '</species>\n'
	
				# write this line to bottom of xml file
				with open(xml_name, "a") as f:
					f.write(lines)
				f.close()
				continue
		
		else: # if absent from xml, then action

			if (number_of_atom_stated == 1):

				# number of carbon, hydrogen and oxygen
				Cn = 0; Hn = 0; On = 0; Nn = 0
	
				# starting flags for counting atoms
				C_flag = 0; H_flag = 0; O_flag = 0; N_flag = 0

				# loop through string elements
				for stri in comp:
				
					if stri == 'C':
						if (C_flag == 1): # if finished counting on carbons
							Cn += float(count_now)
							C_flag = 0
						if (H_flag == 1): # if finished counting on hydrogens
							Hn += float(count_now)
							H_flag = 0
						if (O_flag == 1): # if finished counting on oxygens
							On += float(count_now)
							O_flag = 0
						if (N_flag == 1): # if finished counting on nitrogens
							if count_now == '':
								Nn += 1
							else:
								Nn += float(count_now)
							N_flag = 0

					
						C_flag = 1 # know what we're counting
						count_now = '' # start count

						continue

					if (stri == 'H'):
						if (C_flag == 1): # if finished counting on carbons
							Cn += float(count_now)
							C_flag = 0
						if (H_flag == 1): # if finished counting on hydrogens
							Hn += float(count_now)
							H_flag = 0
						if (O_flag == 1): # if finished counting on oxygens
							On += float(count_now)
							O_flag = 0
						if (N_flag == 1): # if finished counting on nitrogens
							if count_now == '':
								Nn += 1
							else:
								Nn += float(count_now)
							N_flag = 0
					
						H_flag = 1 # know what we're counting
						count_now = '' # start count

						continue
	

					if (stri == 'O'):
						if (C_flag == 1): # if finished counting on carbons
							Cn += float(count_now)
							C_flag = 0
						if (H_flag == 1): # if finished counting on hydrogens
							Hn += float(count_now)
							H_flag = 0
						if (O_flag == 1): # if finished counting on oxygens
							On += float(count_now)
							O_flag = 0
						if (N_flag == 1): # if finished counting on nitrogens
							if count_now == '':
								Nn += 1
							else:
								Nn += float(count_now)
							N_flag = 0

						O_flag = 1 # know what we're counting
						count_now = '' # start count

						continue

					if (stri == 'N'):
						if (C_flag == 1): # if finished counting on carbons
							Cn += float(count_now)
							C_flag = 0
						if (H_flag == 1): # if finished counting on hydrogens
							Hn += float(count_now)
							H_flag = 0
						if (O_flag == 1): # if finished counting on oxygens
							On += float(count_now)
							O_flag = 0
						if (N_flag == 1): # if finished counting on nitrogens
							if count_now == '':
								Nn += 1
							else:
								Nn += float(count_now)
							N_flag = 0

						N_flag = 1 # know what we're counting
						count_now = '' # start count

						continue

					try:
						test = float(stri)
						count_now = str(count_now+stri)
					
					except:
						continue


				try:
					# get any final numbers accounted for
					if C_flag == 1 or H_flag == 1 or O_flag == 1 or N_flag == 1:
						if (C_flag == 1): # if finished counting on carbons
							Cn += float(count_now)
							C_flag = 0
						if (H_flag == 1): # if finished counting on hydrogens
							Hn += float(count_now)
							H_flag = 0
						if (N_flag == 1): # if finished counting on nitrogens
							if count_now == '':
								Nn += 1
							else:
								Nn += float(count_now)
							N_flag = 0
						if (O_flag == 1): # if finished counting on oxygens
							if count_now == '':
								On += 1
							else:
								On += float(count_now)
							O_flag = 0


					# assemble SMILES string
					SMILES = ''
					for i in range(int(Cn)):
						SMILES = str(SMILES + 'C')
					# comment out H if H not accepted in SMILES
					#for i in range(int(Hn)):
					#	SMILES = str(SMILES + 'H')
					for i in range(int(On)):
						SMILES = str(SMILES + 'O')
					for i in range(int(Nn)):
						SMILES = str(SMILES + 'N')
			
					# write lines for xml file
					lines = '<species species_name=\"' + comp + '\" species_number=\"n\">\n' + '<smiles>' + SMILES + '</smiles>\n' + '</species>\n'
				
					# write this line to bottom of xml file
					with open(xml_name, "a") as f:
						f.write(lines)
					f.close()
				
				except: # move onto next component
					continue


			# if SMILES to be formed from a combination of precursors to the ROOR accretion file 
			# (RO2+RO2 = ROOR (ROOR has two oxygens less than the sum of the precursor RO2 oxygens))
			if (number_of_atom_stated == 2):
				
				# combination of reactant SMILE strings
				reac_SMILE = ''
				print(comp)
				print(comp_numi)
				print(len(eqn_list_firsts))
				if comp not in eqn_list_firsts[comp_numi].split():
					print(str('Error, component ' + str(comp) + ' not found in equation: ' + str(eqn_list_firsts[comp_numi])))

				# get the names of the precursors
				for i in eqn_list_firsts[comp_numi].split():

					if i == '+': # skip past plus symbol
						continue

					if i == '=': # finished looking at reactants
						break # stop looking through equation parts


					# get SMILE string of reactant
					reac_SMILE = str(reac_SMILE + comp_smil[comp_name.index(i)])

				# assuming the component in question is a product of RO2 + RO2 = ROOR, then remove two oxygens
				for oxygen_count in range (2):
					# get index of first oxygen
					oxy_indx = reac_SMILE.index('O')
			
					oxy_indx_st  = oxy_indx # starting index
					oxy_indx_fi  = oxy_indx # finishing index
					
					# check indices either side
					if (reac_SMILE[oxy_indx-1:oxy_indx] == '[') or (reac_SMILE[oxy_indx-1:oxy_indx] == '(') or (reac_SMILE[oxy_indx-1:oxy_indx] == '='):
						if reac_SMILE[oxy_indx+1] == '(' or reac_SMILE[oxy_indx+1] == '[' or reac_SMILE[oxy_indx+1] == 'O' or reac_SMILE[oxy_indx+1:oxy_indx+3] == '=(' or reac_SMILE[oxy_indx+1:oxy_indx+3] == '=[':							# do not want to remove punctuation
							oxy_indx_st = oxy_indx_st # starting index
						else:
							oxy_indx_st = oxy_indx-1 # starting index
					if (reac_SMILE[oxy_indx-2:oxy_indx] == '[=') or (reac_SMILE[oxy_indx-2:oxy_indx] == '(='):
						oxy_indx_st = oxy_indx-2 # starting index

					if (reac_SMILE[oxy_indx+1] == ']') or (reac_SMILE[oxy_indx+1] == ')') or (reac_SMILE[oxy_indx+1] == '='):
						oxy_indx_fi = oxy_indx+1 # starting index
					if (reac_SMILE[oxy_indx+1:oxy_indx+3] == '=]') or (reac_SMILE[oxy_indx+1:oxy_indx+3] == '=)'):
						oxy_indx_fi = oxy_indx+2 # starting index
					
					reac_SMILE = str(reac_SMILE[0:oxy_indx_st] + reac_SMILE[oxy_indx_fi+1::])

				# write lines for xml file
				lines = '<species species_name=\"' + comp + '\" species_number=\"n\">\n' + '<smiles>' + reac_SMILE + '</smiles>\n' + '</species>\n'
				
				# write this line to bottom of xml file
				with open(xml_name, "a") as f:
					f.write(lines)
				f.close()
	
	# finish the xml file correctly
	with open(xml_name, "a") as f:
		f.write("</species_defs>\n")
		f.write("</mechanism>\n")
	f.close()			

	return()

xml_cont_gen() # call function