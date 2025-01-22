##########################################################################################
#                                                                                        #
#    Copyright (C) 2018-2024 Simon O'Meara : simon.omeara@manchester.ac.uk               #
#                                                                                        #
#    All Rights Reserved.                                                                #
#    This file is part of PyCHAM                                                         #
#                                                                                        #
#    PyCHAM is free software: you can redistribute it and/or modify it under             #
#    the terms of the GNU General Public License as published by the Free Software       #
#    Foundation, either version 3 of the License, or (at your option) any later          #
#    version.                                                                            #
#                                                                                        #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT               #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              #
#    details.                                                                            #
#                                                                                        #
#    You should have received a copy of the GNU General Public License along with        #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                #
#                                                                                        #
##########################################################################################
'''module to write the equations for mass transfer coefficient to surfaces'''
# writes the mass_trans_coeff_eq module based on the user inputs

import datetime

# define function
def write(mtc_str, self):

	# inputs: ----------------------------------------------
	# mtc_str - string from user inputs containing the equation for
	#	mass transfer coefficient of component(s) to surface(s)
	# self - reference to the PyCHAM object
	# --------------------------------------------------------
	
	# convert Dig string to DStar_org
	mtc_str = mtc_str.replace('D_ig', 'DStar_org')

	# create new  file - will contain module for both deliquescence and efflorescence
	f = open(str(self.PyCHAM_path + '/PyCHAM/mass_trans_coeff_eq.py'), mode='w')

	# loop through string elements to find any reference diffusion coefficients
	for stri in range(len(mtc_str)):
		# in case we want the diffusion coefficient of a specific component
		if (mtc_str[stri:stri+2] == 'D_' and  mtc_str[stri:stri+4] != 'D_ig'):
			# index where reference component finishes
			fini = [] # prepare a list to hold inidices
			if (')' in mtc_str[stri::]):
				fini.append(mtc_str[stri::].index(')'))
			if (']' in mtc_str[stri::]):
				fini.append(mtc_str[stri::].index(']'))
			if ('(' in mtc_str[stri::]):
				fini.append(mtc_str[stri::].index('('))
			if ('[' in mtc_str[stri::]):
				fini.append(mtc_str[stri::].index('['))
			if ('*' in mtc_str[stri::]):
				fini.append(mtc_str[stri::].index('*'))
			if ('/' in mtc_str[stri::]):
				fini.append(mtc_str[stri::].index('/'))
			if ('-' in mtc_str[stri::]):
				fini.append(mtc_str[stri::].index('-'))
			if ('+' in mtc_str[stri::]):
				fini.append(mtc_str[stri::].index('+'))
			ref_comp = mtc_str[stri+2:stri+min(fini)]
			# the replacement string
			rep_str = ('DStar_org[self.comp_namelist.index(\'%s\')]' %ref_comp)
			
			# replace diffusion coefficient term of reference component in
			# user-provided string with string referencing the diffusion
			# coefficient array (cm2/s)
			mtc_str = mtc_str.replace(mtc_str[stri:stri+min(fini)], 
				rep_str)
			
	# loop through string elements to find any reference diffusion coefficients
	for stri in range(len(mtc_str)):

		# in case we want the diffusion coefficient of all components
		if (mtc_str[stri:stri+2] == 'D_' and  mtc_str[stri:stri+4] == 'D_ig'):

			# replace diffusion coefficient term of all components in
			# user-provided string with string referencing the diffusion
			# coefficient array (cm2/s)
			mtc_str = mtc_str.replace(mtc_str[stri:stri+4], 
				'DStar_org[:]')

	
	f.write('##########################################################################\n')
	f.write('#                                                                                        #\n')
	f.write('#    Copyright (C) 2018-2024 Simon O\'Meara : simon.omeara@manchester.ac.uk               #\n')
	f.write('#                                                                                        #\n')
	f.write('#    All Rights Reserved.                                                                #\n')
	f.write('#    This file is part of PyCHAM                                                         #\n')
	f.write('#                                                                                        #\n')
	f.write('#    PyCHAM is free software: you can redistribute it and/or modify it under             #\n')
	f.write('#    the terms of the GNU General Public License as published by the Free Software       #\n')
	f.write('#    Foundation, either version 3 of the License, or (at your option) any later          #\n')
	f.write('#    version.                                                                            #\n')
	f.write('#                                                                                        #\n')
	f.write('#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT               #\n')
	f.write('#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       #\n')
	f.write('#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              #\n')
	f.write('#    details.                                                                            #\n')
	f.write('#                                                                                        #\n')
	f.write('#    You should have received a copy of the GNU General Public License along with        #\n')
	f.write('#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                #\n')
	f.write('#                                                                                        #\n')
	f.write('##################################################################\n')
	f.write('\'\'\'solution of the mass transfer coefficient for components to surfaces\'\'\'\n')
	f.write('# module to estimate mass transfer coefficient for components to surface, generated in mod_var_read and called from both partit_var_prep and cham_up\n')
	f.write('# File Created at %s\n' %(datetime.datetime.now()))	
	f.write('\n')
	f.write('# function for mass transfer coefficient\n')
	f.write('def mtc(DStar_org, TEMP, num_comp, kwn, self):\n')
	f.write('	\n')
	f.write('	# inputs: -----------------\n')
	f.write('	# DStar_org - gas-phase diffusion coefficients of components (cm2/s)\n')
	f.write('	# TEMP - temperature now (K)\n')
	f.write('	# num_comp - number of components\n')
	f.write('	# kwn - mass transfer coefficients of components to wall\n')
	f.write('	# self - the PyCHAM object\n')
	f.write('	# ---------------------------\n')
	f.write('	\n')
	f.write('	# mass transfer coefficient to surface(s) (/s)\n')
	f.write('	kwn = %s\n' %mtc_str)
	f.write('	kwn = kwn.reshape(self.wall_on, num_comp) # correct shape\n')
	f.write('	\n')
	f.write('	return(kwn)\n')

	f.close() # close file
	
	return()
