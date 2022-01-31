##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2022 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
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
