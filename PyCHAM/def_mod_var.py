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
'''default model variables'''
# the default model variables, loaded by PyCHAM on starting

import numpy as np
import os
import pickle # used for storing values

def def_mod_var(caller, self): # define function

	# inputs: -----------------------------------------------
	# caller - mark for calling function
	# self - reference to program
	# -------------------------------------------------------

	# default input files ------------------------------
	# default chemical scheme
	self.sch_name = os.getcwd()+'/PyCHAM/input/gas-phase_ex/ex_chem_scheme.txt'
	self.xml_name = os.getcwd()+'/PyCHAM/input/gas-phase_ex/ex_xml.xml' # xml file path
	self.inname = 'Default' # model variables file name
	
	# general ---------------------------------------------------------------------------------
	# name of folder to save results to
	sav_nam = 'default_res_name'
	# markers to isolate sections of chemical scheme based on MCM KPP format
	self.chem_sch_mrk = ['{', 'RO2', '+', 'C(ind_', ')','' , '&', '' , '', ':', '}', ';', '']
	# time interval between updates to integration inputs (s)
	self.update_stp = 1.
	self.tot_time = 1. # total time to integrate over (s)
	self.save_step = 1. # time interval between saving results (s)
	if (caller == 0): # called from PyCHAM
		uman_up = 0 # marker for whether to update the UManSysProp folder
	if (caller == 1): # called from Travis
		uman_up = 1
	int_tol = [1.e-3, 1.e-4] # integration tolerances (absolute first, relative second)
	self.testf = 0 # whether in testing mode or not
	self.pars_skip = 0 # whether to skip chemical scheme parsing
	# whether to remove influxes of components that aren't seen in chemical scheme
	self.remove_influx_not_in_scheme = 0
	
	# chamber environment -----------------------------------------------------------------
	self.TEMP = np.array((298.15)).reshape(1) # temperature of experiment (K)
	self.tempt = np.array((0.0)).reshape(1) # time that temperatures reached (s)	
	self.RH = np.array(([0.65])) # humidity of experiment (fraction of 1)
	self.RHt = np.array(([0])) # time through simulation (s) RH reached
	Press = 9.8e4 # air pressure during experiment (Pa)
	self.dil_fac = np.zeros(1) # dilution factor (volume fraction per second)
	self.dil_fact = np.zeros(1) # dilution factor times through experiment (s)
	
	# particle section ----------------------------------------------------------------
	siz_stru = 0 # size structure (0 for moving-centre, 1 for full-moving)
	num_sb = 0 # number of particle size bins
	# whether particle number concentrations expressed by modes (0) or explicitly	
	pmode = 0
	# concentration of particles (# particles/cm3 (air))
	pconc = np.zeros((1, 1))
	pconct = np.zeros((1, 1)) # times of particle injection (s)
	pcont = np.zeros((1, 1), dtype=int) # whether particle injection instantaneous or continuous
	seed_mw = (np.array((132.14))).reshape(1) # molecular weight of seed material (g/mol)
	seed_diss = [1.] # dissociation constant of seed material
	seed_dens = np.ones((1)) # density of seed material (g/cc)
	seed_name = ['core'] # name of component forming seed material
	seedx = np.ones((1, 1)) # mole fraction of dry seed components
	self.Vwat_inc = 1 # seed particle number size distribution does include water volume
	self.seed_eq_wat = 1 # allow water equilibration with seed particle before experiment starts
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
	coag_on = 1 # whether to model coagulation
	# flag for particle-phase history with respect to water (0 for dry and therefore 
	# on the deliquescence curve, 1 for wet and therefore on the efflorescence curve)
	wat_hist = 1
	# string describing the deliquescence relative humidity (fraction 0-1) dependence on 
	# temperature
	drh_str = str('0.*TEMP')
	# string describing the efflorescence relative humidity (fraction 0-1) dependence on 
	# temperature
	erh_str = str('0.*TEMP')
	# fraction of total gas-particle partitioning coefficient below which the 
	# partitioning coefficient is set to zero, e.g. because surface area of a size bin
	# is relatively very small
	z_prt_coeff = 1.e-9

	# gas inputs ------------------------------------------------------------
	# chemical scheme name of components present initially
	comp0 = np.array(())
	# initial concentrations (ppb)
	y0 = np.array(())	
	self.con_infl_nam = [] # chemical scheme names of components with continuous influx
	# influx rate of components with continuous influx (ppb/s)
	self.con_infl_C = np.empty(0)
	# times of component influx (s)
	self.con_infl_t = np.empty(0)
	# chemical scheme name of components with constant concentration	
	self.const_comp = []
	# Chemical scheme names of components injected instantaneously after start of experiment
	Compt = []
	# times at which instantaneous injection of component(s) occur after 
	# experiment start (s)
	injectt = np.empty(0)
	# concentration(s) (ppb) of component(s) injected instantaneously after 
	# experiment start
	Ct = np.zeros((0, 0))
	# the gas-particle partitioning cutoff (Pa)
	partit_cutoff = []

	# lights -------------------------------------------------------------------------------
	self.light_stat = np.zeros((1), dtype='int') # light status
	self.light_stat_now = 0 # ready for light checking when chem_sch_SMILES called
	self.light_time = np.zeros((1, 1)) # time that light status attained (s)
	self.daytime = 0. # time of day experiment starts (s)
	self.lat = 0. # latitude of experiment (degrees)
	self.lon = 0. # longitude of experiment (degrees)
	self.af_path = 'no' # path to customised (non-MCM) actinic flux file
	# path to file containing absorption cross-sections and quantum yields
	self.photo_path = str(os.getcwd() + '/PyCHAM/photofiles/' + 'MCMv3.2')
	self.dayOfYear = 1 # number of days since 31st December that experiment conducted
	self.tf = 1. # transmission factor for natural sunlight (0-1 fraction)
	self.tf_range = 0 # flag for no wavelength-dependency of transmission factor
	# marker to say whether or not to adapt integration time interval 
	# and initial condition update to changing natural light intensity
	self.light_ad = 1
	self.tf_UVC = np.ones((1, 1)) # transmission factor for 254 nm wavelength artificial light
	self.tf_UVCt = np.zeros((1, 1)) # time (s) for transmission factor for 254 nm wavelength artificial light
	# if natural sunlight paraemeters have been set in a previous run then delete
	if hasattr(self, 'secx'):
		delattr(self, 'secx')
	if hasattr(self, 'cosx'):
		delattr(self, 'cosx')	


	# deposition of particles and vapours to wall ------------------------------------------
	self.wall_on = 1 # marker for whether to consider wall (0 for no, 1 for yes)
	self.Cw = np.zeros((1)) # effective absorbing mass of wall (g/m3 (air))
	self.kw = np.zeros((1, 1)) # gas-wall mass transfer coefficient (/s)

	inflectDp = 1.e-6 # diameter of deposition function inflection
	pwl_xpre = 0. # gradient before inflection
	pwl_xpro = 0. # gradient after inflection
	inflectk = 0. # rate at inflection
	chamSA = 42. # chamber surface area (m2)
	chamV = 18. # chamber volume (m3, set to MAC value)
	Rader = -1 # flag for deposition to wall treatment (0 for customised, 1 for Rader and McMurry (1985))
	p_char = 0. # average number of charges per particle (/particle)
	e_field = 0. # average electric field inside chamber (g.m/A.s3)

	# specific component properties ---------------------------------------------------------
	# chemical scheme name of components to track the change tendencies of	
	self.dydt_trak = []
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

	# fixing to observations ---------------------------------------------
	self.obs_file = [] # path to observations file
	# components to fix to observed
	self.obs_comp_i = []
	# --------------------------------------------------------------------

	# prepare for pickling
	list_vars = [sav_nam, comp0, y0, Press, 
			siz_stru, num_sb, pmode, pconc, pconct, lowsize, 
			uppsize, space_mode, std, mean_rad, 
			Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, 
			seedx, dens_comp, dens, vol_comp, 
			volP, act_comp, act_user, accom_comp, accom_val, uman_up, 
			int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, 
			coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, 
			Rader, p_char, e_field, partit_cutoff, ser_H2O, 
			wat_hist, drh_str, erh_str, pcont, 
			z_prt_coeff, chamV]

	
	# path to store for variables
	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
		
	with open(input_by_sim, 'wb') as f: # the file to be used for pickling
		pickle.dump(list_vars,f) # pickle
		f.close() # close


	return(sav_nam, comp0, y0, Press, siz_stru, num_sb, 
		pmode, pconc, pconct, lowsize, uppsize, 
		space_mode, std, mean_rad, Compt, 
		injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedx,  
		dens_comp, dens, vol_comp, volP, act_comp, act_user, 
		accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
		nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, 
		pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, 
		partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str, pcont, 
		z_prt_coeff, chamV, self)