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
'''plots results for the wall-phase (from gas-wall partitioning) 
temporal profiles of specified components'''
# simulation results are represented graphically

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import os
import retr_out
import numpy as np
import scipy.constants as si

def plotter(caller, comp_names_to_plot, self):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# self.dir_path - path to folder containing results files to plot
	# comp_names_to_plot - chemical scheme names of components to plot
	# self - reference to GUI
	# --------------------------------------------------------------------------

	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = self.ro_obj.yrec
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = self.ro_obj.Nrec_wet
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec = self.ro_obj.rad
	space_mode = self.ro_obj.spacing
	group_indx = self.ro_obj.gi
	
	# number of actual particle size bins
	num_asb = (num_sb-wall_on)

	if (caller == 0):
		plt.ion() # show results to screen and turn on interactive mode

	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))

	if (comp_names_to_plot): # if component names specified
	
		# concentration plot ---------------------------------------------	
		for i in range(len(comp_names_to_plot)):
			
			# start by assuming it's not a group to be plotted
			group_flag = 0

			if comp_names_to_plot[i].strip() == 'H2O':
				indx_plot = np.array((H2Oi)).reshape(1)
				group_flag = 1
			if (comp_names_to_plot[i].strip() == 'HOMRO2'):
				indx_plot = (np.array((group_indx['HOMRO2'])))
				group_flag = 1
			if (group_flag == 0):
				# will work if provided components were in simulation 
				# chemical scheme
				try:
					# get index of this specified component, removing 
					# any white space
					indx_plot = np.array((comp_names.index(
						comp_names_to_plot[i].strip()))).reshape(1)
				except:
					self.l203a.setText(str('Component ' + 	
					comp_names_to_plot[i] + ' not found in chemical ' +
					'scheme used for this simulation'))
					# set border around error message
					if (self.bd_pl == 1):
						self.l203a.setStyleSheet(0., 
							'2px dashed red', 0., 0.)
						self.bd_pl = 2
					else:
						self.l203a.setStyleSheet(0., 
							'2px solid red', 0., 0.)
						self.bd_pl = 1
					
					plt.ioff() # turn off interactive mode
					plt.close() # close figure window
					return()
			
			if (wall_on > 0):
				# wall-phase concentration (from gas-wall partitioning 
				# and production on wall) (# molecules/cm3)
				# (summed over walls)
				# prepare to hold abundance
				conc = np.zeros((len(timehr)))
				# loop over components in this group (or just one loop if a 
				# single component)

				for ip in range(len(indx_plot)):
					conc_ip = np.sum(yrec[:, ((num_asb+1)*num_comp)
					+indx_plot[ip]::num_comp], axis = 1)

					# convert abundance from molecules/cm3 to ug/m3
					conc += ((conc_ip/si.N_A)*y_MW[0, indx_plot[ip]])*1.e12
			
			else:
				self.l203a.setText(str('Wall not considered in ' +
					 'this simulation'))
				
				# set border around error message
				if (self.bd_pl == 1):
					self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
					self.bd_pl = 2
				if (self.bd_pl >= 2):
					self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
					self.bd_pl = 1

				plt.ioff() # turn off interactive mode
				plt.close() # close figure window
				return()
				
			# plot this component
			ax0.plot(timehr, conc, '+', linewidth = 4., 
				label = str(str(comp_names_to_plot[i]+
				' (wall (from gas-wall partitioning))')))

		ax0.set_ylabel(r'Concentration ($\rm{\mu}$g$\,$m$\rm{^{-3}}$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)

		# end of wall-phase concentration ----------------------------
	

	# display
	if (caller == 2):
		plt.show()	

	return()