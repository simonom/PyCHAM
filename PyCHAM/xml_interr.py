##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2024 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
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
import openbabel.pybel as pybel # for converting SMILE strings to pybel objects

# define function
def xml_interr(xml_name):

	# inputs: --------------------------------------------------------
	# xml_name - name of xml file
	# ----------------------------------------------------------------

	# start with no error message
	err_mess_new = ''

	with open(xml_name) as fd:
		try:
			doc = xmltodict.parse(fd.read())
		except:
			err_mess_new = 'Error: xml file could not be interpreted, please check file'
			return(err_mess_new, [], [])

	a = doc['mechanism']['species_defs']['species']
	
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
	return(err_mess_new, comp_smil, comp_name)