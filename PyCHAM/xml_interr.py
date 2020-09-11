'''interrogate the xml file'''
# opens and extracts the component names and associated SMILE strings
# from the xml file

import xmltodict # for opening and converting xml files to python dictionaries

# define function
def xml_interr(xml_name):

	# inputs: --------------------------------------------------------
	# xml_name - name of xml file
	# ----------------------------------------------------------------

	with open(xml_name) as fd:
		doc = xmltodict.parse(fd.read())

	a = doc['mechanism']['species_defs']['species']
	
	# prepare arrays to fill	
	comp_numb = list(('0',) * len(a))
	comp_name = list(('0',) * len(a))
	comp_smil = list(('0',) * len(a))
	
	for i in range(len(a)):
		comp_numb[i] = a[i]['@species_number']
		comp_name[i] = a[i]['@species_name']
		if "smiles" in a[i]:
			comp_smil[i] = a[i]['smiles']
		elif comp_name[i][0]=='O' or comp_name[i][0]=='H':
			 comp_smil[i] = '['+comp_name[i]+']'
		else:
			 comp_smil[i] = comp_name[i]
 
	return(comp_smil, comp_name)
