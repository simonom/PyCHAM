##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2023 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
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
'''plot the temporal profile of chamber environment (temperature, pressure, relative humidity)'''
# plots temperature, pressure and relative humidity with time through simulation

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import os
import retr_out
import numpy as np
import scipy.constants as si

# define function
def plotter(caller, dir_path, self):

	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# self - reference to GUI
	# --------------------------------------------------------------------------	

	wall_on = self.ro_obj.wf
	yrec = np.zeros((self.ro_obj.yrec.shape[0], self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], self.ro_obj.Nrec_wet.shape[1]))
	Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
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
	cham_env = self.ro_obj.env_cond

	# check whether chamber environment variables were saved and therefore
	# retrieved
	if (cham_env == []):
		mess = str('Please note, no chamber environmental variables were found, perhaps the simulation was completed in a PyCHAM version predating this functionality')
		self.l203a.setText(mess)
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

	# prepare figure
	plt.ion() # display figure in interactive mode
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	# ensure that all axes can be seen
	plt.subplots_adjust(left=0.1, right=0.8)
	
	# parasite axis part ---------------------------------------------------
	par1 = ax0.twinx() # first parasite axis
	par2 = ax0.twinx() # second parasite axis
	# Offset the right spine of par2.  The ticks and label have already been
	# placed on the right by twinx above.
	par2.spines["right"].set_position(("axes", 1.11))
	# Having been created by twinx, par2 has its frame off, so the line of its
	# detached spine is invisible.  First, activate the frame but make the patch
	# and spines invisible.
	make_patch_spines_invisible(par2)
	# second, show the right spine
	par2.spines['right'].set_visible(True)
	# -------------------------------------------------------------------------
	
	# plot temporal profiles
	# temperature on original vertical axis
	p0, = ax0.plot(timehr, cham_env[:, 0], 'k', label = 'temperature (K)')
	ax0.set_ylabel('Temperature (K)', size=16)
	ax0.yaxis.label.set_color('black')
	ax0.tick_params(axis='y', colors='black')
	ax0.spines['left'].set_color('black')
	ax0.yaxis.set_tick_params(direction = 'in', which = 'both')
	
	# pressure on first parasite axis
	p1, = par1.plot(timehr, cham_env[:, 1], '--m', label = 'pressure (Pa)')
	par1.set_ylabel('Pressure (Pa)', rotation=270, size=16, labelpad = 15)
	par1.yaxis.label.set_color('magenta')
	par1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e')) # set tick format for vertical axis
	par1.tick_params(axis='y', colors='magenta')
	par1.spines['right'].set_color('magenta')
	par1.yaxis.set_tick_params(direction = 'in', which = 'both')
	
	# relative humidity on second parasite axis
	p2, = par2.plot(timehr, cham_env[:, 2], '-.b', label = 'relative humidity (fraction)')
	par2.set_ylabel('Relative Humidity (0-1)', rotation=270, size=16, labelpad = 15)
	par2.yaxis.label.set_color('blue')
	par2.tick_params(axis='y', colors='blue')
	par2.spines['right'].set_color('blue')
	par2.yaxis.set_tick_params(direction = 'in', which = 'both')
	
	
	
	ax0.xaxis.set_tick_params(direction = 'in', which = 'both')
		
	ax0.set_xlabel('Time through experiment (hours)', size=16)
	plt.legend(fontsize=14, handles=[p0, p1, p2], loc=4)
	

	return()

# function to makes spines invisible
def make_patch_spines_invisible(ax):
	ax.set_frame_on(True)
	ax.patch.set_visible(False)
	for sp in ax.spines.values():
		sp.set_visible(False)
