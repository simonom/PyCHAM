'''default model variables'''
# the default model variables inputs

import numpy as np
import os
import pickle # used for storing values

def def_mod_var(caller): # define function

	# inputs: -----------------------------------------------
	# caller - mark for calling function
	# -------------------------------------------------------

	# names ---------------------------------------------------------------------------------
	# name of folder to save results to
	sav_nam = 'default_res_name'
	sch_name = os.getcwd()+'/PyCHAM/input/example_scheme.txt'
	# markers to isolate sections of chemical scheme based on MCM KPP format
	chem_sch_mark = ['%', 'RO2', '+', '', '', ';', '+', ';', '', '%', ':', ';']
	xml_name = os.getcwd()+'/PyCHAM/input/example_xml.xml' # xml file path

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
	siz_stru = 0 # size structure (0 for moving-centre, 1 for full-moving)
	num_sb = 0 # number of particle size bins
	# whether particle number concentrations expressed by modes (0) or explicitly	
	pmode = 0
	# initial concentration of particles (# particles/cc (air))
	pconc = np.zeros((1, 1))
	pconct = np.zeros((1, 1)) # times of particle injection (s)
	seed_mw = 132.14 # molecular weight of seed material (g/mol)
	seed_diss = [1.] # dissociation constant of seed material
	seed_dens = 1. # density of seed material (g/cc)
	seed_name = ['core'] # name of component forming seed material
	seedVr = [1.] # volume ratio of seed components
	lowsize = 0. # smallest size bin boundary (radius) (um)
	uppsize = 5.e-1 # largest size bin boundary (radius) (um)
	space_mode = 'lin' # treatment for spacing between size bins
	std = np.ones((1, 1))*1.2 # standard deviation for particle number size distribution
	mean_rad = np.ones((1, 1))*-1.e6 # mean radius for particle number size distribution (um)
	new_partr = 2.e-7 # radius of newly nucleated particles (cm)
	# nucleation parameters
	nucv1 = 0.
	nucv2 = 0.
	nucv3 = 0.
	nuc_comp = ['core'] # chemical scheme name of nucleating component
	# marker to say whether or not to adapt integration time interval 
	# and initial condition update to nucleation
	nuc_ad = 1
	ser_H2O = 1 # whether to serialise gas-particle partitioning of water

	# component inputs ------------------------------------------------------------
	# chemical scheme name of components present initially,
	comp0 = np.array(())
	# initial concentrations (ppb),
	y0 = np.array(())	
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
	Ct = np.zeros((0, 0))

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
	# chemical scheme names of components with densities manually set
	dens_comp = []
	# manually assigned densities (g/cc)
	dens = []
	# chemical scheme names of components with vapour pressures manually set
	vol_comp = []
	# manually assigned vapour pressures (Pa)
	volP = []
	# names of components (corresponding to chemical scheme name) with 
	# activity coefficient stated in act_user
	act_comp = []
	# user-specified activity coefficients of components with names given 
	# in act_comp
	act_user = []
	# names of components with user-defined accommodation coefficients	
	accom_comp = []
	# user-defined accommodation coefficients
	accom_val = []
	# the gas-particle partitioning cutoff (Pa)
	partit_cutoff = []

	if (caller == 0): # called from PyCHAM
		uman_up = 0 # marker for whether to update the UManSysProp folder
	if (caller == 1): # called from Travis
		uman_up = 1
	int_tol = [1.e-3, 1.e-4]
	coag_on = 1 # whether to model coagulation
	dil_fac = 0. # dilution factor

	# prepare for pickling
	list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O]

	
	# path to store for variables
	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
		
	with open(input_by_sim, 'wb') as f: # the file to be used for pickling
		pickle.dump(list_vars,f) # pickle
		f.close() # close


	return(sav_nam, sch_name, chem_sch_mark, xml_name, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O)
