'''unit test for calculating accommodation coefficients'''
# Checks that the automatically generated accommodation coefficient
# calculation module functions correctly.  Assumes calling from the
# PyCHAM home folder.

import numpy as np
import os
import pickle # used for storing values
import sys
import shutil

def test_accom_coeff_calc(): # define function

	print('Running test for accom_coeff_calc.py, the starting gas- and particle-phase concentrations of O3 and APINENE have been set equal, as have their vapour pressures, but their accommodation coefficients are different, therefore, the resulting time profile of their gas-phase concentrations should be different')

	# inputs: ------------
	# --------------------

	# ensure modules can be seen 
	# (assumes calling function is in the home folder)
	sys.path.append(str(os.getcwd() + '/PyCHAM'))


	# state model variables
	
	# names ---------------------------------------------------------------------------------
	# name of folder to save results to
	sav_nam = 'test_accom_coeff_calc_res'
	sch_name = os.getcwd()+'/PyCHAM/unit_tests/input/test_scheme_nonreac' # chemical scheme file path
	# markers to isolate sections of chemical scheme based on MCM KPP format
	chem_sch_mark = ['%', 'RO2', '+', '', '', ';', '+', ';', '', '%', ':', ';']
	xml_name = os.getcwd()+'/PyCHAM/input/Example_Run_xml.xml' # xml file path

	# times ----------------------------------------------------------------------------------
	# time interval between updates to integration inputs (s)
	update_stp = 1.
	tot_time = 1. # total time to integrate over (s)
	save_step = 1. # time interval between saving results (s)
	
	# initial atmosphere -----------------------------------------------------------------
	temp = np.array((298.15)).reshape(1) # temperature of experiment (K)
	tempt = np.array((0.0)).reshape(1) # time that temperatures reached (s)	
	RH = 0.65 # humidity of experiment (fraction of 1)
	Press = 9.8e4 # air pressure during experiment (Pa)
	
	# particle section ----------------------------------------------------------------
	num_sb = 1 # number of particle size bins
	# initial concentration of particles (# particles/cc (air))
	pconc = np.ones((1, 1))*1.e4
	pconct = np.zeros((1, 1)) # times of particle injection (s)
	seed_mw = 132.14 # molecular weight of seed material (g/mol)
	core_diss = 1.0 # dissociation constant of seed material
	seed_dens = 1.0 # density of seed material (g/cc)
	seed_name = 'core' # name of component forming seed material
	lowsize = 0. # smallest size bin boundary (radius) (um)
	uppsize = 5.e-1 # largest size bin boundary (radius) (um)
	space_mode = 'lin' # treatment for spacing between size bins
	std = np.zeros((1, 1)) # standard deviation for particle number size distribution
	mean_rad = np.zeros((1, 1)) # mean radius for particle number size distribution (um)
	new_partr = 2.e-7 # radius of newly nucleated particles (cm)
	# nucleation parameters
	nucv1 = 0.
	nucv2 = 0.
	nucv3 = 0.
	nuc_comp = ['core'] # chemical scheme name of nucleating component
	# marker to say whether or not to adapt integration time interval 
	# and initial condition update to nucleation
	nuc_ad = 1

	# component inputs ------------------------------------------------------------
	# chemical scheme name of components present initially,
	comp0 = np.array(('O3', 'APINENE'))
	# initial concentrations (ppb),
	y0 = np.array((30., 30.))	
	con_infl_nam = [] # chemical scheme names of components with continuous influx
	con_infl_t = np.array(())
	# influx rate of components with continuous influx (ppb/s)
	con_infl_C = np.array(())
	# times of component influx (s)
	const_infl_t = []	
	# chemical scheme name of components with constant concentration	
	const_comp = []
	# name(s) of component(s) injected instantaneously after start of experiment
	Compt = []
	# times at which instantaneous injection of component(s) occur after 
	# experiment start (s)
	injectt = []
	# concentration(s) (ppb) of component(s) injected instantaneously after 
	# experiment start
	Ct = []

	# lights -------------------------------------------------------------------------------
	light_stat = [0] # light status
	light_time = np.zeros((1, 1)) # time that light status attained (s)
	daytime = 0. # time of day experiment starts (s)
	lat = 0. # latitude of experiment (degrees)
	lon = 0. # longitude of experiment (degrees)
	af_path = 'no' # path to customised (non-MCM) actinic flux file
	# path to file containing absorption cross-sections and quantum yields
	photo_path = str(os.getcwd() + '/PyCHAM/photofiles/' + 'MCMv3.2')
	dayOfYear = 1 # number of days since 31st December that experiment conducted
	tf = 1. # transmission factor for natural sunlight (0-1 fraction)
	# marker to say whether or not to adapt integration time interval 
	# and initial condition update to changing natural light intensity
	light_ad = 1

	# deposition of particles and vapours to wall ------------------------------------------
	wall_on = 1 # marker for whether to consider wall (0 for no, 1 for yes)
	Cw = 0. # effective absorbing mass of wall (g/m3 (air))
	kw = 0. # gas-wall mass transfer coefficient (/s)

	inflectDp = 0. # diameter of deposition function inflection
	pwl_xpre = 0. # gradient before inflection
	pwl_xpro = 0. # gradient after inflection
	inflectk = 0. # rate at inflection
	chamSA = 42. # chamber surface area (m2) 
	Rader = -1 # flag for deposition to wall treatment (0 for customised, 1 for Rader and McMurry (1985))
	p_char = 0. # average number of charges per particle (/particle)
	e_field = 0. # average electric field inside chamber (g.m/A.s3)

	# miscellaneous ------------------------------------------------------------------------
	# chemical scheme name of components to track the change tendencies of	
	dydt_trak = []
	# chemical scheme names of components with vapour pressures manually set
	vol_comp = ['O3', 'APINENE']
	# manually assigned vapour pressures (Pa)
	volP = [1.e-4, 1.e-4]
	# names of components (corresponding to chemical scheme name) with 
	# activity coefficient stated in act_user
	act_comp = []
	# user-specified activity coefficients of components with names given 
	# in act_comp
	act_user = []
	# names of components with user-defined accommodation coefficients	
	accom_comp = ['O3', 'APINENE']
	# user-defined accommodation coefficients
	accom_val = ['2.0', '6.09e-08/radius']
	uman_up = 0 # marker for whether to update the UManSysProp folder

	int_tol = [1.e-3, 1.e-4]
	coag_on = 1 # whether to model coagulation
	dil_fac = 0. # dilution factor

	# prepare for pickling
	list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, num_sb, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, core_diss, seed_dens, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac]

	
	# path to store for variables
	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
		
	with open(input_by_sim, 'wb') as f: # the file to be used for pickling
		pickle.dump(list_vars,f) # pickle
		f.close() # close

	import middle # the main call to program
	middle.middle()

	# call to plot
	print('Plotting and saving standard results graph')
	import plotter
	plotter.plotter(2) # plot results

	# delete results folder		
	shutil.rmtree(str(os.getcwd() + '/PyCHAM/output/test_scheme_nonreac'))	

	return()
	

test_accom_coeff_calc() # call function
