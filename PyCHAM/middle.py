'''calls modules to setup simulation'''
# the modules necessary to setup a simulation are called here

import eqn_pars
import user_input as ui
import ode_updater
import init_conc
import prop_calc
import partit_var_prep
import pp_intro
import time
import save
import pickle
import os

# define function
def middle():
	
	# get required inputs
	[sav_nam, sch_name, chem_sch_mrk, xml_name, update_stp, tot_time, 
		comp0, y0, temp, tempt, RH, Pnow, wall_on,
		Cw, kw, siz_str, num_sb, pconc, pconct,
		lowsize, uppsize, space_mode, std, mean_rad,
		save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, 
		core_diss, core_dens, light_stat, light_time, daytime, lat, 
		lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, 
		const_infl_t, con_infl_C, dydt_trak, vol_comp, volP, 
		act_comp, act_user, accom_comp, accom_coeff_user, uman_up, 
		int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, 
		inflectDp, pwl_xpre, pwl_xpro, inflectk, ChamR, Rader, p_char, 
		e_field, dil_fac] = ui.share(0)

	# parse the chemical scheme equation file to convert equations
	# into usable code
	[rindx, pindx, rstoi, pstoi, nreac, nprod, comp_num, 
	jac_stoi, njac, jac_den_indx, 
	jac_indx, RO2_indx, comp_list, Pybel_objects, eqn_num, 
	comp_namelist, Jlen, y_arr, y_rind, uni_y_rind, y_pind, uni_y_pind, 
	reac_col, prod_col, rstoi_flat, 
	pstoi_flat, rr_arr, rr_arr_p, rowvals, colptrs, jac_wall_indx, 
	jac_part_indx] = eqn_pars.extr_mech(sch_name, 
		chem_sch_mrk, xml_name, photo_path, con_infl_nam, int_tol, wall_on, 
		(num_sb+wall_on), const_comp)
	
	# set initial concentrations (molecules/cc)
	[y, H2Oi, y_mw, num_comp, Cfactor, indx_plot, corei, dydt_vst, comp_namelist, 
	inj_indx, core_diss, Psat_water, 
	nuci, con_infl_C, nrec_steps] = init_conc.init_conc(comp_num, comp0, y0, temp[0], RH, 
	Pnow, Pybel_objects, 0, pconc, dydt_trak, tot_time, save_step, rindx, 
	pindx, eqn_num[0], nreac, nprod, 
	comp_namelist, Compt, seed_name, 
	seed_mw, core_diss, nuc_comp, con_infl_C)
	
	# dump new pickle file ready for plotting script to use
	list_vars = [sav_nam, sch_name, indx_plot, comp0]
	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
	with open(input_by_sim, 'wb') as f: # the file to be used for pickling
		pickle.dump(list_vars, f) # pickle
		f.close() # close

	# get component properties
	[Psat, y_dens, Psat_Pa] = prop_calc.prop_calc(comp_list, Pybel_objects, temp[0], H2Oi, 
		num_comp, Psat_water, vol_comp, volP, 0, corei, pconc,
		uman_up, core_dens, comp_namelist, 0, nuci, nuc_comp, num_sb)
	
	# prepare for the calcuation of partitioning variables
	[DStar_org, mfp, accom_coeff, therm_sp, surfT, Cw, act_coeff, 
		R_gas, NA] = partit_var_prep.prep(y_mw, 
		temp[0], num_comp, 0, Cw, act_comp, act_user, accom_comp, 
		accom_coeff_user, comp_namelist, num_sb, num_sb)
	count = 0	
	
	# prepare particle phase and wall
	[y, N_perbin, x, Varr, Vbou, rad0, Vol0, rbou, MV, num_sb, nuc_comp, 
	rbou00, upper_bin_rad_amp, np_sum] = pp_intro.pp_intro(y, num_comp, Pybel_objects, temp[0], H2Oi, 
		mfp, accom_coeff, y_mw, surfT, DStar_org, RH, siz_str, num_sb, lowsize, 
		uppsize, pconc, pconct, nuc_comp, 0, std, mean_rad, 
		therm_sp, Cw, y_dens, Psat, core_diss, kw, space_mode, corei,
		comp_namelist, act_coeff, wall_on)
	
	print('Calling integration routine, starting timer')
	st_time = time.time()
	
	# solve problem
	[trec, yrec, dydt_vst, Cfactor_vst, Nres_dry, Nres_wet, x2] = ode_updater.ode_updater(update_stp, 
		tot_time, save_step, y, rindx, 
		pindx, rstoi, pstoi, nreac, nprod, jac_stoi, njac, 
		jac_den_indx, jac_indx, RO2_indx, H2Oi, temp, tempt, 
		Pnow, light_stat, light_time, daytime, lat, lon, af_path, 
		dayOfYear, photo_path, Jlen, con_infl_C, nrec_steps, 
		dydt_vst, num_sb, num_comp, corei, core_diss, Psat, 
		mfp, therm_sp,  
		accom_coeff, y_mw, surfT, R_gas, NA, y_dens, DStar_org, 
		x, Varr, act_coeff, Cw, kw, Cfactor, tf, light_ad, y_arr, 
		y_rind, 
		uni_y_rind, y_pind, uni_y_pind, reac_col, prod_col, 
		rstoi_flat, pstoi_flat, rr_arr, rr_arr_p, rowvals, 
		colptrs, wall_on, jac_wall_indx, jac_part_indx, Vbou, 
		N_perbin, Vol0, rad0, np_sum, new_partr, nucv1, nucv2, 
		nucv3, nuc_comp, nuc_ad, RH, coag_on, inflectDp, pwl_xpre, 
		pwl_xpro, inflectk, ChamR, Rader, p_char, e_field, 
		injectt, inj_indx, Ct, pconc, pconct, mean_rad, lowsize, 
		uppsize, std, rbou, const_infl_t)
	
	time_taken = time.time()-st_time
	print('Simulation complete, wall clock time elapsed since first call to solver: ', time_taken, ' s')		

	print('Saving results')
	# save results
	save.saving(sch_name, yrec, Nres_dry, Nres_wet, trec, sav_nam, 
		dydt_vst, num_comp, Cfactor_vst, 0, 
		num_sb, comp_namelist, dydt_trak, y_mw, MV, time_taken, 
		seed_name, x2, rbou, wall_on, space_mode, rbou00, upper_bin_rad_amp)
	print('Saved results')
	
	
	return()
