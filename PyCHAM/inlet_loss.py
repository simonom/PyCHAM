'''module for inlet loss correction'''
# when particles pass through an inlet system this module accounts for size-dependent loss

import numpy as np
import datetime
import importlib

def inlet_loss(call, Nwet, xn, yp, loss_func_str, losst, num_comp):

	# inputs: -----------------------------------------------
	# call - the calling module flag
	# Nwet - number concentration of particles entering the inlet (# particles/cm3)
	# xn - radii of particles, corrected for relative humidity in inlet (um)
	# yp - concentrations of components in the particle phase, corrected for 
	#	relative humidity in inlet (molecules/cm3)
	# loss_func_str - string stating particle loss rate (fraction/s) 
	#	as a function of particle size (um)
	# losst - time of passage through inlet (s)
	# num_comp - number of components
	# ---------------------------------------------------------

	# create new  file - will contain the inlet loss rate as a function of particle size
	f = open('PyCHAM/inlet_loss_func.py', mode='w')
	
	f.write('\'\'\'solving the size-dependent loss rate of particles in instrument inlet\'\'\'\n')
	f.write('# module to estimate loss rate (fraction /s) of particles during passage through an instrument inlet \n')
	f.write('# File Created at %s\n' %(datetime.datetime.now()))	
	f.write('\n')
	f.write('import numpy as np\n')
	f.write('\n')
	f.write('# function for loss rate\n')
	f.write('def inlet_loss_func(Dp_all):\n')
	f.write('	\n')
	f.write('	# inputs: -------------------\n')
	f.write('	# Dp_all - diameter (um) of size bins (columns) per time step (rows)\n')
	f.write('	# -----------------------------\n')
	f.write('	\n')
	f.write('	sd_lrate = np.zeros((Dp_all.shape[0], Dp_all.shape[1]))\n')
	f.write('	try: # in case loss function string is acceptable\n')
	f.write('		# estimate loss rate (fraction/s)\n')
	f.write('		for it in range(Dp_all.shape[0]):\n')
	f.write('			Dp = Dp_all[it, :]\n')
	f.write('			sd_lrate[it, :] = %s\n' %(loss_func_str))
	f.write('	except: # in case of issue\n')
	f.write('		sd_lrate = \'Error, function of loss rate of particles during passage through inlet failed, please revise\'\n')
	f.write('	\n')
	f.write('	return(sd_lrate)')
	f.close()
	
	import inlet_loss_func
	importlib.reload(inlet_loss_func)
		

	# size-dependent rates of particles lost over entire inlet passage (fraction /s)
	# for all size bins (columns) and times (rows)
	sd_lrate = inlet_loss_func.inlet_loss_func(xn*2.)
	
	if (call == 3): # if button pressed to plot loss rate as a function of diameter

		if isinstance(sd_lrate, str):
			if (sd_lrate[0:5] == 'Error'):
				return(sd_lrate)		

		import matplotlib.pyplot as plt
		plt.ion()
		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
		# plot weighting dependency against response time
		if ((sd_lrate == 0).any()):
			ax0.semilogx((xn*2.)[0, :], sd_lrate[0, :], '-k')
		else:
			ax0.loglog((xn*2.)[0, :], sd_lrate[0, :], '-k')

		# plot details
		ax0.set_title('Loss rate of particles during passage through instrument inlet')
		ax0.set_ylabel('Loss rate (fraction of particles (0-1) s$^{-1}$)', size = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		ax0.set_xlabel('Particle diameter ($\mu\,$m)', fontsize=14)
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		plt.show()
		return([])

	# integrate over entire time in inlet (fraction)
	sd_lrate = sd_lrate*losst
	
	Nwet -= Nwet*sd_lrate # correct for inlet loss (# particles/cm3)
	# repeat over components and correct for loss of components (molecules/cm3)
	yp -= yp*np.repeat(sd_lrate, num_comp, axis=1)

	return(Nwet, yp)