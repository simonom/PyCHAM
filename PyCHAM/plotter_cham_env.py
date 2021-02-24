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

	# retrieve results
	(_, _, _, _, _, _, _, timehr, _, 
		_, _, _, _, _, _, _, 
		_, _, _, _, _, _, _, _, cham_env) = retr_out.retr_out(dir_path)

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