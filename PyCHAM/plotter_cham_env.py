'''plot the temporal profile of chamber environment (temperature, pressure, relative humidity)'''
# plots temperature, pressure and relative humidity with time through simulation

import matplotlib.pyplot as plt
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

	# prepare figure
	plt.ion() # display figure in interactive mode
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
	# plot temporal profiles
	ax0.plot(timehr, cham_env[:, 0], label = 'temperature (K)')
	ax0.plot(timehr, cham_env[:, 1], label = 'pressure (Pa)')
	ax0.plot(timehr, cham_env[:, 2], label = 'relative humidity (fraction)')
	ax0.yaxis.set_tick_params(direction = 'in')
		
	ax0.set_xlabel('Time through experiment (hours)')
	ax0.set_ylabel('Temperature (K))')

	return()