########################################################################
#								       #
# Copyright (C) 2018-2024					       #
# Simon O'Meara : simon.omeara@manchester.ac.uk			       #
#								       #
# All Rights Reserved.                                                 #
# This file is part of PyCHAM                                          #
#                                                                      #
# PyCHAM is free software: you can redistribute it and/or modify it    #
# under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or    #
# (at  your option) any later version.                                 #
#                                                                      #
# PyCHAM is distributed in the hope that it will be useful, but        #
# WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# General Public License for more details.                             #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with PyCHAM.  If not, see <http://www.gnu.org/licenses/>.      #
#                                                                      #
########################################################################
'''determines number of photochemical reactions'''
# based on the photochemistry file path, the number of reactions is
# inferred

import os

# define function
def photo_num(photo_file):

	# inputs: -------------------------------------------------------------
	# photo_file - path to file containing abosrption cross-sections and
	# quantum yields
	# ---------------------------------------------------------------------
	# number of photolysis reactions, if this relevant
	# if using the MCM photofiles
	if (photo_file[-26::] == '/PyCHAM/photofiles/MCMv3.2'):
		Jlen = 62 # for MCM (default name of photolysis parameters)
	else: # need to find out number of photolysis reactions
		# use Fortran indexing to be consistent with MCM photochemical reaction numbers
		Jlen = 1 
		# open file to read
		f = open(str(photo_file), 'r')
		for line in f: # loop through line
			if line.strip() == str('J_'+str(Jlen) + '_axs'):
				Jlen += 1
	return(Jlen)
