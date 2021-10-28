'''plot checks on PyCHAM inputs from the simulate tab of the GUI'''
# with the exception of the number size distribution plot 
# (plotter_nsd module), this module plots the requested checks on
# PyCHAM inputs as selected from the PyCHAM GUI simulate tab

def plotter_taf(self): # define function for total actinic flux

	import os
	import pickle # for storing inputs

	# inputs: --------------------------
	# self - self-reference to PyCHAM
	# ----------------------------------

	# prepare by opening existing inputs, ready for modification
	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
	with open(input_by_sim, 'rb') as pk:
		[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
		tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
		Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
		save_step, const_comp, Compt, injectt, Ct, seed_name,
		seed_mw, seed_diss, seed_dens, seedx,
		light_stat, light_time, daytime, lat, lon, af_path, 
		dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, 
		dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
		accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
		nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
		inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, 
		ser_H2O, wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, 
		z_prt_coeff, tf_UVC, testf, chamV] = pickle.load(pk)
		pk.close() # close pickle file
	
	testf = 4 # modify test flag value

	# pickle with new testf
	list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
		tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on, 
		Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, 
		uppsize, space_mode, std, mean_rad, save_step, const_comp, 
		Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, 
		seedx, light_stat, light_time, daytime, lat, lon, af_path, 
		dayOfYear, photo_path, tf, light_ad, con_infl_nam, 
		con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, 
		volP, act_comp, act_user, accom_comp, accom_val, uman_up, 
		int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, 
		coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, 
		Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, 
		wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, 
		z_prt_coeff, tf_UVC, testf, chamV]

	with open(input_by_sim, 'wb') as pk:
		pickle.dump(list_vars, pk) # pickle
		pk.close() # close
		
	# now run program up to the plot
	from middle import middle # prepare to communicate with main program
		
	note_messf = 0 # cancel note message flag
		
	for prog in middle(): # call on modules to solve problem
			
		
		if (isinstance(prog, str)): # check if it's a message
			mess = prog
			if (mess[0:4] == 'Stop'): # if it's an error message
				return()
	
	return()

def plotter_gpdc(self): # define function for gas-phase diffusion coefficients

	import os
	import pickle # for storing inputs

	# inputs: --------------------------
	# self - self-reference to PyCHAM
	# ----------------------------------

	# prepare by opening existing inputs, ready for modification
	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
	with open(input_by_sim, 'rb') as pk:
		[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
		tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
		Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
		save_step, const_comp, Compt, injectt, Ct, seed_name,
		seed_mw, seed_diss, seed_dens, seedx,
		light_stat, light_time, daytime, lat, lon, af_path, 
		dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, 
		dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
		accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
		nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
		inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, 
		ser_H2O, wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, 
		z_prt_coeff, tf_UVC, testf, chamV] = pickle.load(pk)
		pk.close() # close pickle file
	
	testf = 2 # modify test flag value

	# pickle with new testf
	list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
		tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on, 
		Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, 
		uppsize, space_mode, std, mean_rad, save_step, const_comp, 
		Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, 
		seedx, light_stat, light_time, daytime, lat, lon, af_path, 
		dayOfYear, photo_path, tf, light_ad, con_infl_nam, 
		con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, 
		volP, act_comp, act_user, accom_comp, accom_val, uman_up, 
		int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, 
		coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, 
		Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, 
		wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, 
		z_prt_coeff, tf_UVC, testf, chamV]

	with open(input_by_sim, 'wb') as pk:
		pickle.dump(list_vars, pk) # pickle
		pk.close() # close
		
	# now run program up to the gas-phase diffusion coefficient plot
	from middle import middle # prepare to communicate with main program
		
	note_messf = 0 # cancel note message flag
		
	for prog in middle(): # call on modules to solve problem
			
		
		if (isinstance(prog, str)): # check if it's a message
			mess = prog
			if (mess[0:4] == 'Stop'): # if it's an error message
				return()
	
	return()

def plotter_gpmts(self): # define function for gas-phase mean thermal speeds

	import os
	import pickle # for storing inputs

	# inputs: --------------------------
	# self - self-reference to PyCHAM
	# ----------------------------------

	# change the test flag to value 3:

	# prepare by opening existing inputs, ready for modification
	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
	with open(input_by_sim, 'rb') as pk:
		[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
		tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
		Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
		save_step, const_comp, Compt, injectt, Ct, seed_name,
		seed_mw, seed_diss, seed_dens, seedx,
		light_stat, light_time, daytime, lat, lon, af_path, 
		dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, 
		dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
		accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
		nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
		inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, 
		ser_H2O, wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, 
		z_prt_coeff, tf_UVC, testf, chamV] = pickle.load(pk)
		pk.close() # close pickle file
	
	testf = 3 # modify test flag value

	# pickle with new testf
	list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
		tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on, 
		Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, 
		uppsize, space_mode, std, mean_rad, save_step, const_comp, 
		Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, 
		seedx, light_stat, light_time, daytime, lat, lon, af_path, 
		dayOfYear, photo_path, tf, light_ad, con_infl_nam, 
		con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, 
		volP, act_comp, act_user, accom_comp, accom_val, uman_up, 
		int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, 
		coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, 
		Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, 
		wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, 
		z_prt_coeff, tf_UVC, testf, chamV]

	with open(input_by_sim, 'wb') as pk:
		pickle.dump(list_vars, pk) # pickle
		pk.close() # close
		
	# now run program up to the gas-phase diffusion coefficient plot
	from middle import middle # prepare to communicate with main program
		
	note_messf = 0 # cancel note message flag
		
	for prog in middle(): # call on modules to solve problem
			
		
		if (isinstance(prog, str)): # check if it's a message
			mess = prog
			if (mess[0:4] == 'Stop'): # if it's an error message
				return()