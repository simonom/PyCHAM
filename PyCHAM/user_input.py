'''module collects and passes user input information to front.py'''

import importlib
import numpy as np
import sys
import os


def run(source, testf):
	
	# inputs:
	# source - flag for whether front (0) or plotting (1) is calling
	# testf - flag for whether operating in normal mode (0) or test (1 or 2)
	if testf==1: # testing mode
		# return dummies to continue test
		if source==0:
			return(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
		
	
	import pickle
	
	if testf==2: # testing mode
		var_store_name = 'test_var_store.pkl'
	else:
		var_store_name = 'PyCHAM/var_store.pkl'
	
	with open(var_store_name,'rb') as pk:
		# read in variables from gui
		if source == 0: # when called from front.py
			[fname, num_sb, lowersize, uppersize, end_sim_time, resfname, tstep_len, 
			TEMP, PInit, RH, lat, lon, dt_start, Cw, save_step, ChamSA, 
			nucv1, nucv2, nucv3, nuc_comp, inflectDp, pwl_xpre, pwl_xpro, 
			inflectk, Rader, xmlname, C0, Comp0, voli, volP, pconc, 
			std, loc, scale, core_diss, light_stat, light_time, kgwt] = pickle.load(pk)	

			
			# convert chamber surface area (m2) to spherical equivalent radius (m)
			# (below eq. 2 in Charan (2018))
			ChamR = (ChamSA/(4.0*np.pi))**0.5
		if source == 1:	# when called from res_plot_super.py
			[fname, resfname, y_indx_plot, Comp0] = pickle.load(pk)
				
					
	output_root = 'output'
	filename = os.path.basename(fname)
	filename = os.path.splitext(filename)[0]
	
	pathname = os.path.dirname(sys.argv[0])        
	dir_path = os.path.abspath(pathname)    
	
	# one folder for one simulation - only relevant if called from front 
	# (rather than res_plot_super)
	output_by_sim = os.path.join(dir_path, 'output', filename, resfname)
	
	
	if os.path.isdir(output_by_sim)==True and source==0:
		sys.exit('Results file name (' +output_by_sim+ ') already exists, please use an alternative')
	
	if source == 0:
		return(fname, num_sb, lowersize, uppersize, end_sim_time, resfname, 
		tstep_len, tstep_len, TEMP, PInit, RH, lat, lon, dt_start, save_step, 
		Cw, ChamR, nucv1, nucv2, nucv3, nuc_comp, inflectDp, 
		pwl_xpre, pwl_xpro, inflectk, xmlname, C0, Comp0, Rader, voli, volP, pconc, std, 
		loc, scale, core_diss, light_stat, light_time, kgwt, 0)
		
	if source == 1:
		return(fname, resfname, y_indx_plot, Comp0)
		
if __name__ == "__main__":
    run()