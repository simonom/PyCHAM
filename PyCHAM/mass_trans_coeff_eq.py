##########################################################################
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
##################################################################
'''solution of the mass transfer coefficient for components to surfaces'''
# module to estimate mass transfer coefficient for components to surface, generated in mod_var_read and called from both partit_var_prep and cham_up
# File Created at 2025-02-14 16:21:21.602842

# function for mass transfer coefficient
def mtc(DStar_org, TEMP, num_comp, kwn, self):
	
	# inputs: -----------------
	# DStar_org - gas-phase diffusion coefficients of components (cm2/s)
	# TEMP - temperature now (K)
	# num_comp - number of components
	# kwn - mass transfer coefficients of components to wall
	# self - the PyCHAM object
	# ---------------------------
	
	# mass transfer coefficient to surface(s) (/s)
	kwn = (1./75.)*(DStar_org/DStar_org[self.comp_namelist.index('C10H16O8iso1')])
	kwn = kwn.reshape(self.wall_on, num_comp) # correct shape
	
	return(kwn)
