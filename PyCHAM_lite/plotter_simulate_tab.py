#########################################################################								       #
# Copyright (C) 2018-2024					       #
# Simon O'Meara : simon.omeara@manchester.ac.uk			       ##								       #
# All Rights Reserved.                                                 #
# This file is part of PyCHAM                                          #
#                                                                      #
# PyCHAM is free software: you can redistribute it and/or modify it    ## under the terms of the GNU General Public License as published by    #
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
'''plot checks on PyCHAM inputs from the simulate tab of the GUI'''
# with the exception of the number size distribution plot 
# (plotter_nsd module), this module plots the requested checks on
# PyCHAM inputs as selected from the PyCHAM GUI simulate tab

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
# for customised colormap
from matplotlib.colors import LinearSegmentedColormap
# set colormap tick labels to standard notation
import matplotlib.ticker as ticker 
import numpy as np
import math

def plotter_taf(self): # define function for photolysis stuff

	import os
	import pickle # for storing inputs

	# inputs: --------------------------
	# self - self-reference to PyCHAM
	# ----------------------------------
	
	self.testf = 4 # modify test flag value
		
	# now run program up to the plot
	from middle import middle # prepare to communicate with main program
		
	note_messf = 0 # cancel note message flag
		
	for prog in middle(self): # call on modules to solve problem
			
		# note that for checking on photolysis rates, middle is called here,
		# then inside middle var_checker is called
		if (isinstance(prog, str)): # check if it's a message
			mess = prog
			if (mess[0:4] == 'Stop'): # if it's an error message
				return()
	
	return()

def plotter_gpdc(self): # define function for gas-phase diffusion coefficients

	import os
	import pickle # for storing inputs

	# inputs: --------------------------
	# self - self-reference to PyCHAM
	# ----------------------------------

	self.testf = 2 # modify test flag value
		
	# now run program up to the gas-phase diffusion coefficient plot
	from middle import middle # prepare to communicate with main program
		
	note_messf = 0 # cancel note message flag
		
	for prog in middle(self): # call on modules to solve problem
			
		
		if (isinstance(prog, str)): # check if it's a message
			mess = prog
			if (mess[0:4] == 'Stop'): # if it's an error message
				return()
	
	return()

def plotter_gpmts(self): # define function for gas-phase mean thermal speeds

	import os
	import pickle # for storing inputs

	# inputs: --------------------------
	# self - self-reference to PyCHAM
	# ----------------------------------

	self.testf = 3 # modify test flag value
		
	# now run program up to the gas-phase diffusion coefficient plot
	from middle import middle # prepare to communicate with main program
		
	note_messf = 0 # cancel note message flag
		
	for prog in middle(self): # call on modules to solve problem
			
		
		if (isinstance(prog, str)): # check if it's a message
			mess = prog
			if (mess[0:4] == 'Stop'): # if it's an error message
				return()

def plotter_mm(self): # define function for molar masses

	import os
	import pickle # for storing inputs

	# inputs: --------------------------
	# self - reference to PyCHAM program
	# ----------------------------------

	self.testf = 3.1 # modify test flag value
		
	# now run program up to the gas-phase diffusion coefficient plot
	from middle import middle # prepare to communicate with main program
		
	note_messf = 0 # cancel note message flag
		
	for prog in middle(self): # call on modules to solve problem
			
		
		if (isinstance(prog, str)): # check if it's a message
			mess = prog
			if (mess[0:4] == 'Stop'): # if it's an error message
				return()
				
def plotter_vp(self): # define function for vapour pressures

	# inputs: --------------------------
	# self - reference to PyCHAM object
	# ----------------------------------
	
	self.testf = 3.2 # modify test flag value
		
	# now run program up to the vapour pressure plot
	# prepare to communicate with main program
	from middle import middle 
		
	note_messf = 0 # cancel note message flag
		
	for prog in middle(self): # call on modules to solve problem
		
		if (isinstance(prog, str)): # check if it's a message
			mess = prog
			# if it's an error message
			if (mess[0:4] == 'Stop'): 
				return()

def plotter_nucfunc(self):

	# inputs: --------------------------
	# self - reference to PyCHAM object
	# ----------------------------------

	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	plt.ion()
	if (self.nucv3 > 0):
		nsum1 = self.nucv1*np.exp(self.nucv2*(np.exp(\
			-self.update_stp/self.nucv3)))
		if (math.isinf(nsum1)):
			ax0.text(x=0.5, y=0.5, s=str('nucleation '
			'parameters give infinite value, please modify'
			', \nand note that an example of nucleation '
			'parameters is available \nin '
			'PyCHAM/input/gas-phase_ex'), size=12, 
			color='black')

		else:
			# run over all times
			trange = np.arange(0, self.tot_time, 
				self.update_stp)
			nsum1 = self.nucv1*np.exp(self.nucv2*(np.exp(\
			-trange/self.nucv3)))
			ax0.plot(trange/3.6e3, nsum1)
			ax0.set_ylabel('Total number concentration '
				'of newly \nnucleated particles '
				'(particles/cm3)', fontsize = 14)	
			ax0.set_xlabel('Time through simulation '
				'(hours)', fontsize = 14)
		
			ax0.yaxis.set_tick_params(labelsize = 14, 
				direction = 'in')
			ax0.xaxis.set_tick_params(labelsize = 14, 
				direction = 'in')

	else:
	
		ax0.text(x=0.5, y=0.5, s=str('nucv3 is 0, therefore '
			'duration of nucleation is zero \nand no '
			'nucleation occurs'), size=12, color='black')

	
	plt.show()
