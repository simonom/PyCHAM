'''module to produce plot to check on variables provided as input'''
# using the variable checker buttons from the simulate tab in the GUI,
# this module runs the necessary operations to produce some of the 
# requested plots

import numpy as np # for math and array handling

def var_checker(testf, light_stat, light_time, daytime, lat, lon, temp, tempt, tot_time, 
	act_flux_path, DayOfYear, photo_par_file, Jlen, tf, tf_UVC, update_stp, 
	err_mess, erf): # define function

	# inputs: -----------------------------------------------------------------------
	# testf - flag for whether in checking mode and what to check
	# light_stat - status of light (on=1 or off=0)
	# light_time - times through experiment light status corresponds to
	# daytime - time of day experiment starts (for natural light photolysis)
	# lat - latitude
	# lon - longitude
	# temp - temperatures inside chamber (K)
	# tempt - times through experiment (s) temperatures correspond to
	# tot_time - total time to run experiment for (s)
	# act_flux_path - path to file stating actinic flux
	# DayOfYear - day number of the year (1-365)
	# photo_par_file - name of file with with estimates for photolysis absorption
	# 	cross-sections and quantum yields
	# Jlen - number of photolysis reactions
	# tf - sunlight transmission factor
	# tf_UVC - transmission factor for 254 nm wavelength light (0-1)
	# tstep - suggested interval between integrations (s)
	# err_mess - any existing error mesages
	# erf -any existing error flags
	# ---------------------------------------------------------------------------------

	if (testf == 4): # checking on estimated photolysis rates throughout simulation
		
		# following the method used in rate_coeffs
		import photolysisRates
		
		sumt = 0 # time through experiment (s)
		
		# ensure numpy arrays
		tempt = np.array(tempt)
		temp = np.array(temp)
		
		while (sumt < tot_time):

			# identify relevant light status
			lindx = np.sum(light_time >= sumt)
		
			if (light_stat[lindx-1] == 0): # if lights off
				if (sumt == 0): # initiate results array (time in rows, photolysis channels in columns)
					Jres = (np.zeros(Jlen)).reshape(1, -1)
				else:
					Jres = np.concatenate((Jres, (np.zeros(Jlen)).reshape(1, -1)), axis = 0)
					
			else: # if lights on
				TEMP = temp[np.sum(tempt >= sumt)-1] # temperature inside chamber now (K)
			
				# estimate and append photolysis rates
				J = photolysisRates.PhotolysisCalculation(daytime+sumt, lat, lon, TEMP, 
					act_flux_path, DayOfYear, photo_par_file, Jlen, tf, sumt, tf_UVC)
			
				if (sumt == 0): # initiate results array (time in rows, photolysis channels in columns)
					Jres = J.reshape(1, -1)
				else:
					Jres = np.concatenate((Jres, J.reshape(1, -1)), axis = 0)
					
			sumt += update_stp # time through experiment
		
		import matplotlib.pyplot as plt
		from matplotlib.colors import BoundaryNorm
		from matplotlib.ticker import MaxNLocator
		from matplotlib.colors import LinearSegmentedColormap # for customised colormap
		import matplotlib.ticker as ticker # set colormap tick labels to standard notation
		import math as math
		plt.ion() # allow plotting
		
		# times to plot against (hours through experiment)
		thr = (np.arange(0.0, sumt, update_stp))/3600.0
		
		# number of plots
		nplot = math.ceil(Jlen/10.)
		
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
			ax0.set_title('Note that reaction numbers start from 1', fontsize = 14)
			ax0.legend()
		
		err_mess = 'Stop'
		erf = 1

	return(err_mess, erf)