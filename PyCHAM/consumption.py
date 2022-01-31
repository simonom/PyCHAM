'''estimating consumption of a component over simulation'''
# a module to calculate the consumed mass concentration of a 
# component introduced artificially to chamber, so not produced
# by means other than through injection

import scipy.constants as si # scientific constants
import retr_out # retrieving information
import numpy as np # for arithmetic

def cons(comp_chem_schem_name, dir_path, self, caller):

	# inputs: -------------------------------
	# comp_chem_schem_name - chemical scheme name of component
	# dir_path - path to results
	# self - reference to GUI
	# caller - flag for calling function
	# ---------------------------------------

	# retrieve results
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_mw, Nwet, comp_names, y_MV, _, wall_on, space_mode, indx_plot, 
		comp0, _, PsatPa, OC, H2Oi, seedi, _, _, _, tot_in_res, 
		_) = retr_out.retr_out(dir_path)

	try:
		# get index of component of interest
		compi = comp_names.index(comp_chem_schem_name)
	except:
		self.l203a.setText(str('Error - could not find component ' + comp_chem_schem_name + ' in the simulated system.  Please check whether the name matches that in the chemical scheme.'))		
		# set border around error message
		if (self.bd_pl == 1):
			self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
			self.bd_pl = 2
		else:
			self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
			self.bd_pl = 1
		return() # return now

	# get consumption
	try:
		indx_int = (np.where(tot_in_res[0, :].astype('int') == compi))[0][0]
	
	except:
		self.l203a.setText(str('Error - could not find a consumption value for component ' + comp_chem_schem_name + '.  Please check whether the name matches that in the chemical scheme and that it had gas-phase influx through the model variables file.'))		
		# set border around message
		if (self.bd_pl == 1):
			self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
			self.bd_pl = 2
		else:
			self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
			self.bd_pl = 1
		return() # return now

	# indices for supplied time
	indxt = (timehr>=self.tmin)*(timehr<=self.tmax)
	
	# cumulative influxed over all times (ug/m3)
	tot_in = tot_in_res[1::, indx_int]
	
	# influxes over interested time period
	tot_in = tot_in[indxt]

	# gas-phase concentration (ppb) over all times
	yrecn = yrec[:, compi]

	# gas-phase concentration (ppb) over interested time period (ppb)
	yrecn = yrecn[indxt]

	# change in gas-phase concentration (ug/m3)
	yrecn = (((yrecn[-1]-yrecn[0])*Cfac[-1])/si.N_A)*y_mw[compi]*1.e12

	# total influx over this time (ug/m3)
	tot_in = tot_in[-2]-tot_in[0]

	# total consumed (ug/m3)
	cons = tot_in-yrecn
	
	if (caller == 0): # call from the consumption button
	
		self.l203a.setText(str('Consumption of ' + comp_chem_schem_name + ': ' + str(cons) + ' ' + u'\u03BC' + 'g/m' + u'\u00B3'))
		# set border around message
		if (self.bd_pl == 1):
			self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
			self.bd_pl = 2
		else:
			self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
			self.bd_pl = 1
		return() # return now

	if (caller == 1): # call from the yield button

		# concentrations of components in particle phase at end of simulation (# molecules/cm3)
		SOA = yrec[-2, num_comp:num_comp*(num_sb-wall_on+1)]
		# remove seed and water in all size bins
		SOA[seedi[0]::num_comp] = 0.
		SOA[H2Oi::num_comp] = 0.
		
		# convert from # molecules/cm3 to ug/m3
		SOA = ((SOA/si.N_A)*np.tile(y_mw, (num_sb-wall_on)))*1.e12

		# sum for total (ug/m3)
		SOA = np.sum(SOA)

		yld = SOA/cons

	
		self.l203a.setText(str('Yield of ' + comp_chem_schem_name + ': ' + str(yld)))
		# set border around message
		if (self.bd_pl == 1):
			self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
			self.bd_pl = 2
		else:
			self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
			self.bd_pl = 1
		return() # return now

	
	
	return()