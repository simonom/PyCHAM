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
'''module to produce plot to check on variables provided as input'''
# using the variable checker buttons from the simulate tab in the GUI,
# this module runs the necessary operations to produce some of the 
# requested plots

import numpy as np # for math and array handling

def var_checker(tot_time, Jlen, update_stp, err_mess, erf, self): # define function

	# inputs: -----------------------------------------------------------------------
	# self.testf - flag for whether in checking mode and what to check
	# self.light_stat - status of light (on=1 or off=0)
	# self.light_time - times through experiment light status corresponds to
	# self.TEMP - temperatures inside chamber (K)
	# self.tempt - times through experiment (s) temperatures correspond to
	# tot_time - total time to run experiment for (s)
	# self.af_path - path to file stating actinic flux
	# self.DayOfYear - day number of the year (1-365)
	# self.photo_file - name of file with with estimates for photolysis absorption
	# 	cross-sections and quantum yields
	# Jlen - number of photolysis reactions
	# self.tf - sunlight transmission factor
	# tstep - suggested interval between integrations (s)
	# err_mess - any existing error mesages
	# erf -any existing error flags
	# self.tf_UVC - transmission factor for 254 nm wavelength light (0-1)
	# ---------------------------------------------------------------------------------

	if (self.testf == 4): # checking on estimated photolysis rates throughout simulation
		
		# following the method used in rate_coeffs
		import photolysisRates
		
		sumt = 0 # time through experiment (s)
		self.sumt = 0.

		# ensure numpy arrays
		self.tempt = np.array(self.tempt)
		self.TEMP = np.array(self.TEMP)
		
		while (sumt < tot_time-tot_time*1.e-10):

			# identify relevant light status
			lindx = np.sum(self.light_time >= sumt)
		
			if (self.light_stat[lindx-1] == 0): # if lights off
				if (sumt == 0): # initiate results array (time in rows, photolysis channels in columns)
					Jres = (np.zeros(Jlen)).reshape(1, -1)
				else:
					Jres = np.concatenate((Jres, (np.zeros(Jlen)).reshape(1, -1)), axis = 0)
					
			else: # if lights on
				TEMPn = self.TEMP[np.sum(self.tempt >= sumt)-1] # temperature inside chamber now (K)
			
				# estimate and append photolysis rates
				J = photolysisRates.PhotolysisCalculation(TEMPn, 
					Jlen, sumt, self)
				
				if (sumt == 0): # initiate results array (time in rows, photolysis channels in columns)
					Jres = J.reshape(1, -1)
				else:
					Jres = np.concatenate((Jres, J.reshape(1, -1)), axis = 0)
					
			sumt += update_stp # time through experiment (s)
			self.sumt += update_stp # time through experiment (s)

		# ignore first column of Jres as this left empty by photolysisRates
		Jres = Jres[:, 1::]
		Jlen = Jres.shape[1] # total number of equations
		import matplotlib.pyplot as plt
		from matplotlib.colors import BoundaryNorm
		from matplotlib.ticker import MaxNLocator
		from matplotlib.colors import LinearSegmentedColormap # for customised colormap
		import matplotlib.ticker as ticker # set colormap tick labels to standard notation
		import math as math
		plt.ion() # allow plotting
		
		# times to plot against (hours through experiment)
		thr = (np.arange(0.0, tot_time, update_stp))/3600.0
		
		# number of plots
		nplot = math.ceil((Jlen)/10.)
		
		for spi in range(nplot):
			# prepare plot
			fig, (ax0) = plt.subplots(1, 1, figsize = (7, 3))
			
			# photolysis reaction numbers being plotted now
			if ((spi+1)*10 <= Jlen): # if within bounds
				prn = (np.arange(spi*10, (spi+1)*10)).astype(int)
			else: # will be reaching final reaction
				prn = (np.arange(spi*10, Jlen)).astype(int)
			
			# loop through each of these photolysis reactions
			for pri in prn:

				# plot photolysis rates (/s), note reaction number label starts from 1 
				# (not pythonic counting)
				ax0.plot(thr, Jres[:, pri], label = str('Photolysis Reaction # ' + str(pri+1)))
			
			ax0.set_xlabel(r'Time through experiment (hr)', fontsize = 14)
			ax0.set_ylabel('Photolysis rate ($\mathrm{s^{-1}}$)', fontsize = 14)
			# set location of x ticks
			ax0.set_title('Reaction numbers in legend start from 1 and \ntherefore align directly with MCM photolysis numbers', fontsize = 14)
			ax0.legend()
		
		err_mess = 'Stop'
		erf = 1

	return(err_mess, erf)