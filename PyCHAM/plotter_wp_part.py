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
'''plots results for the wall-phase (from particle partitioning to wall) temporal profiles of specified components'''
# simulation results are represented graphically

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import os
import numpy as np
import scipy.constants as si

def plotter(caller, dir_path, comp_names_to_plot, self):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# comp_names_to_plot - chemical scheme names of components to plot
	# self - reference to GUI
	# --------------------------------------------------------------------------
	
	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = np.zeros((self.ro_obj.yrec.shape[0], self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], self.ro_obj.Nrec_wet.shape[1]))
	Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
	Ndry = np.zeros((self.ro_obj.Nrec_dry.shape[0], self.ro_obj.Nrec_dry.shape[1]))
	Ndry[:, :] = self.ro_obj.Nrec_dry[:, :]
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec= np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	space_mode = self.ro_obj.spacing
	Cfac = self.ro_obj.cfac
	group_indx = self.ro_obj.gi
	y_MV = (np.array((self.ro_obj.comp_MV))).reshape(1, -1)
	yrec_p2w = self.ro_obj.part_to_wall
	PsatPa = self.ro_obj.vpPa
	O_to_C = self.ro_obj.O_to_C
	H2Oi = self.ro_obj.H2O_ind
	
	# number of actual particle size bins
	num_asb = (num_sb-wall_on)
	

	if (caller == 0):
		plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))

	if (comp_names_to_plot): # if component names specified
	
		# plotting section ---------------------------------------------	
		for i in range(len(comp_names_to_plot)):
			
			if comp_names_to_plot[i].strip() == 'H2O':
				indx_plot = H2Oi
			else:
				try: # will work if provided components were in simulation chemical scheme
					# get index of this specified component, removing any white space
					indx_plot = comp_names.index(comp_names_to_plot[i].strip())
				except:
					self.l203a.setText(str('Component ' + comp_names_to_plot[i] + ' not found in chemical scheme used for this simulation'))
					# set border around error message
					if (self.bd_pl == 1):
						self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
						self.bd_pl = 2
					else:
						self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
						self.bd_pl = 1
					
					plt.ioff() # turn off interactive mode
					plt.close() # close figure window
					return()
			
			if (wall_on > 0):
				# total concentration on wall (from particle deposition to wall) (molecules/cm3)
				conc = (yrec_p2w[:, indx_plot::num_comp]).sum(axis = 1)
				
			else:
				self.l203a.setText(str('Wall not considered in this simulation'))
				
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

				
			# concentration in ug/m3
			conc = ((conc/si.N_A)*y_MW[0, indx_plot])*1.e12
			# plot this component
			ax0.plot(timehr, conc, '+', linewidth = 4., label = str(str(comp_names[indx_plot ]+' (wall (from particle deposition to wall))')))

		ax0.set_ylabel(r'Concentration ($\rm{\mu}$g$\,$m$\rm{^{-3}}$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)

		# end of gas-phase concentration sub-plot ---------------------------------------
	

	# display
	if (caller == 2):
		plt.show()	

	return()