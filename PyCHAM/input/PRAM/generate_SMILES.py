'''function to generate xml contents from chemical schemes'''
# because some chemical schemes don't come with xml files, but
# do have identifiers for species that aids SMILES generation

def xml_cont_gen(): # define function

	# import dependencies
	import re # for parsing
	import numpy as np # for math and matrix operations
	import xmltodict # for opening and converting xml files to python dictionaries
	import pybel # for converting SMILE strings to pybel objects

	# open scheme
	chem_sch = open('/Users/Simon_OMeara/Desktop/PRAM_KPP.txt', mode='r') # open the chemical scheme file
	# read the file and store everything into a list
	total_list_eqn = chem_sch.readlines()
	chem_sch.close() # close file

	# markers to isolate sections of chemical scheme based on PRAM format
	chem_sch_mrk = ['%', 'RO2', '+', '', '', ';', '+', ';', '', '%', ':', ';']
	
	eqn_list = [] # # empty list for equations

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
								
	# now that we have names of all components in equation file, 
	# check for their presence in xml file
	xml_name = '/Users/Simon_OMeara/Desktop/PRAM_xml.xml'
	
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
		
	
	# loop through components in chemical scheme file to 
	# check on presence in xml
	for comp in comp_list:
		if comp in comp_name: # if present then continue
			continue
		print(comp)
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
	
	# ensure to finish the xml file correctly
	with open(xml_name, "a") as f:
		f.write("</species_defs>\n")
		f.write("</mechanism>\n")
	f.close()			

	return()

xml_cont_gen() # call function