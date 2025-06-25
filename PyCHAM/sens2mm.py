####################################################################
#                                                                                        #
#    Copyright (C) 2018-2025 Simon O'Meara : simon.omeara@manchester.ac.uk               #
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
'''solving the sensitivity (Hz/ppt) of instrument to molar mass (g/mol)'''
# module to estimate the sensitivity of an instrument to the molar mass of components, for example a Chemical Ionisiation Mass Spectrometer
# File Created at 2025-06-09 15:17:44.678180

import numpy as np

# function for sensitivity
def sens2mm(caller, y_MM, Cn):
	
	# inputs: -----------------
	# caller - flag for the calling function
	# y_MM - molar mass (g/mol) of components in question
	# Cn - carbon number
	# ---------------------------
	
	fac_per_comp = np.ones((len(y_MM))) # sensitivity (Hz/ppt) per molar mass (g/mol) 
	fac_per_comp = np.array((fac_per_comp)).reshape(-1) # reshape 
	
	if (caller == 3): # called on to plot sensitivity to molar mass
		import matplotlib.pyplot as plt 
		plt.ion()
		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
		ax0.plot(y_MM, fac_per_comp)
		ax0.set_title('Sensitivity of instrument to molar mass')
		ax0.set_ylabel('Sensitivity (fraction (0-1))', size = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		ax0.set_xlabel('Molar Mass ($\mathrm{g\,mol^{-1}}$)', fontsize=14)
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	
	return(fac_per_comp)