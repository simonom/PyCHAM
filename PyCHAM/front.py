'''module called on by __main__ to initiate the box model'''

import numpy as np
from ode_gen import ode_gen
import matplotlib.pyplot as plt
import os
import eqn_parser # to parse the .scm file
from init_conc_func import init_conc_func
from pp_intro import pp_intro
from kimt_prep import kimt_prep
from saving import saving
import time # timing how long operations take
import user_input as ui
import pickle # for storing inputs
from volat_calc import volat_calc

def run(testf):
	
	# inputs:
	# testf - test flag, 1 for test mode (called by test_front.py), 2 for test mode 
	# (called by test_PyCHAM.py), 0 for normal mode
	
	if testf==2:
		print('"Run Model" button works fine')
		return()
	if testf==1:
		print('calling user_input.py')
	# module to ask, receive and return required inputs
	[fname, num_sb, lowersize, uppersize, end_sim_time, resfname, tstep_len, 
	tmax, TEMP, PInit, RH, lat, lon, start_sim_time, save_step, Cw, 
	ChamR, nucv1, nucv2, nucv3, nuc_comp, inflectDp, pwl_xpre,  
	pwl_xpro, inflectk, xmlname, init_conc, Comp0, Rader, voli, volP, 
	pconc, std, loc, scale, core_diss, light_stat, light_time, kgwt, testm] = ui.run(0, testf)
	
	if testm == 1:
		print('PyCHAM calls front fine, now returning to PyCHAM.py')
		print('Please select the "Plot Results" button to ensure this works fine')
		return()
	
	if testf==1:
		print('user_input.py called and returned fine')
		print('calling eqn_parser.extract_mechanism')
	
	# obtain gas-phase reaction info
	[rindx, pindx, rstoi, pstoi, reac_coef, spec_list, Pybel_objects, num_eqn, num_speci, 
		RO2_indices, nreac, nprod, prodn, 
		reacn, M_val, N2_val, O2_val, 
		init_SMIL] = eqn_parser.extract_mechanism(fname, xmlname, 
							TEMP, PInit, Comp0, testf)
	
	if testf==1:
		print('eqn_parser.extract_mechanism called and returned fine')
		print('calling init_conc_func')
	# set up initial gas-phase concentration array
	[y, H2Oi, Psat_water, y_mw, num_speci, 
	Cfactor, y_indx_plot, corei] = init_conc_func(num_speci, init_SMIL, spec_list, 
							init_conc, TEMP, RH, M_val, N2_val, O2_val, reac_coef, fname, 
							PInit, start_sim_time, lat, lon, Pybel_objects, testf, pconc)
	
	
	if testf==1:
		print('init_conc_func called and returned fine')
		print('calling kimt_prep')
	# set up partitioning variables
	[DStar_org, mfp, accom_coeff, therm_sp, surfT, Cw, kgwt] = kimt_prep(y_mw, TEMP, 
														num_speci, testf, Cw, kgwt)
	
	# volatility (molecules/cc (air)) and density (rho, kg/m3) of components
	if testf==1:
		print('kimt_prep called and returned fine')
		print('calling volat_calc.py')
	[Psat, y_dens, Psat_Pa] = volat_calc(spec_list, Pybel_objects, TEMP, H2Oi, num_speci,  
								Psat_water, voli, volP, testf, corei)
	
	if testf==1:
		print('volat_calc called and returned fine')
		print('calling wall_prep')
	
	if testf==1:
		print('wall_prep called and returned fine')
		print('calling pp_intro')
	# set up particle phase part
	[y, N_perbin, x, Varr, Vbou, rad0, Vol0, rbou, 
							new_partr, MV, num_sb, nuc_comp] = pp_intro(y, 
							num_speci, spec_list, Pybel_objects, TEMP, H2Oi, 
							mfp, accom_coeff, y_mw, surfT, DStar_org, 
							RH, num_sb, lowersize, uppersize, pconc, tmax, nuc_comp, 
							voli, volP, testf, std, loc, scale, 
							therm_sp, Cw, y_dens, Psat, core_diss, kgwt)
							
	
	t1 = time.clock() # get wall clock time before call to solver
	if testf==1:
		print('pp_intro called and returned fine')
		print('calling ode_gen')
	
	# call on ode function
	[t_out, y_mat, Nresult, x2] = ode_gen(tstep_len, y, num_speci, num_eqn, rindx, pindx, 
				rstoi, pstoi, H2Oi, TEMP, RO2_indices, 
				num_sb, Psat, mfp, accom_coeff, surfT, y_dens, N_perbin,
				DStar_org, y_mw, x, core_diss, Varr, Vbou, RH, rad0, Vol0,
				end_sim_time, pconc, save_step, 
				rbou, therm_sp, Cw, light_time, light_stat,
				nreac, nprod, prodn,
				reacn, new_partr, MV, nucv1, nucv2, nucv3, inflectDp, pwl_xpre, 
				pwl_xpro, inflectk, nuc_comp, ChamR, Rader, PInit, testf, kgwt)
				
	
	t2 = time.clock() # get wall clock time after call to solver
	if testf==0: # in normal mode
		print('time taken=')
		print(t2-t1)
	
		# make new pickle dump to store the indices and names of interesting gas-phase 
		# components along with initial pickle variables
		list_vars = [fname, resfname, y_indx_plot, Comp0]
	if testf==1:
		print('ode_gen called and returned fine')
		print('dumping variables in pickle file')
		# dummy list of variables to dump
		list_vars = ['fnametest','resfnametest', 0, 0]
		with open('test_var_store.pkl','wb') as f:
			pickle.dump(list_vars,f)
	
	if testf==00:
		with open('PyCHAM/var_store.pkl','wb') as f:
			pickle.dump(list_vars,f) 
	
	if testf==1:
		print('dumped successfully')
	# save data
	output_by_sim = saving(fname, y_mat, t_out, Nresult, x2, num_sb, y_mw, num_speci, 
							resfname, rbou, Cfactor, MV, testf)
	if testf==1:
		print('saving called and returned successfully')
	return()