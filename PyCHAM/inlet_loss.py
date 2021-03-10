'''module for inlet loss correction'''
# when particles pass through an inlet system this module accounts for size-dependent loss

import numpy as np
import datetime
import importlib

def inlet_loss(Nwet, xn, yp, loss_func_str, losst, num_comp):

	# inputs: -----------------------------------------------
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
	f.write('	# estimate loss rate (fraction/s)\n')
	f.write('	for it in range(Dp_all.shape[0]):\n')
	f.write('		Dp = Dp_all[it, :]\n')
	f.write('		sd_lrate[it, :] = %s\n' %(loss_func_str))
	f.write('	\n')
	f.write('	return(sd_lrate)')
	f.close()
	
	import inlet_loss_func
	importlib.reload(inlet_loss_func)

	# size-dependent rates of particles lost over entire inlet passage (fraction /s)
	# for all size bins (columns) and times (rows)
	sd_lrate = inlet_loss_func.inlet_loss_func(xn*2.)
	
	# integrate over entire time in inlet (fraction)
	sd_lrate = sd_lrate*losst
	
	Nwet -= Nwet*sd_lrate # correct for inlet loss (# particles/cm3)
	# repeat over components and correct for loss of components (molecules/cm3)
	yp -= yp*np.repeat(sd_lrate, num_comp, axis=1)

	return(Nwet, yp)