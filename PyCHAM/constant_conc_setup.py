##########################################################################################
#                                                                                        #
#    Copyright (C) 2018-2026 Simon O'Meara : simon.omeara@manchester.ac.uk               #
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
'''module to setup the arrays for constant concentration components'''
# based on user input

import numpy as np

def constant_conc_setup(erf, err_mess, self):

	# empty array for storing components with constant concentration
	self.con_C_indx = np.zeros((self.const_comp.shape[0], 
		self.const_comp.shape[1])).astype('int')

	# components with constant concentration
	for i in range(self.const_comp.shape[0]):
		for it in range(self.const_comp.shape[1]):
			try:
				# index of where constant concentration components occur in list 
				# of components
				self.con_C_indx[i, it] = self.comp_namelist.index(
					self.const_comp[i, it])
				
			except:
				# if a component doesn't appear at a given time then provide
				# a marker for no component
				if (self.const_comp[i, it] == ''):
					self.con_C_indx[i, it] = -1e6
					continue	
				
				# if water then we know it will be the next 
				# component to
				# be appended to the component list
				if (self.const_comp[i, it] == 'H2O'):
					self.con_C_indx[i] = len(
						self.comp_namelist)
				else: # if not water
					erf = 1 # raise error
					err_mess = str('Error: constant ' +
					'concentration ' +
					'component with name ' + 
					str(self.const_comp[i]) +  
					'has not been identified in the ' +
					'chemical scheme, ' +
					'please check it is present and the ' +
					'chemical scheme markers are correct')
	
	return(erf, err_mess, self) # end of constant_conc_setup function
