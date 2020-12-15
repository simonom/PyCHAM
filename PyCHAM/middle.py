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
		Cw, kw, siz_str, num_sb, pmode, pconc, pconct,
		lowsize, uppsize, space_mode, std, mean_rad,
		save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, 
		core_diss, core_dens, seedVr, light_stat, light_time, daytime, lat, 
		lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, 
		const_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, 
		volP, act_comp, act_user, accom_comp, accom_coeff_user, uman_up, 
		int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, 
		inflectDp, pwl_xpre, pwl_xpro, inflectk, ChamR, Rader, p_char, 
		e_field, dil_fac, partit_cutoff, ser_H2O] = ui.share(0)
	
	# parse the chemical scheme equation file to convert equations
	# into usable code
	[rindx_g, pindx_g, rstoi_g, pstoi_g, nreac_g, nprod_g, 
	jac_stoi_g, njac_g, jac_den_indx_g, 
	jac_indx_g, y_arr_g, y_rind_g, uni_y_rind_g, y_pind_g, uni_y_pind_g, 
	reac_col_g, prod_col_g, rstoi_flat_g, 
	pstoi_flat_g, rr_arr_g, rr_arr_p_g, rowvals, colptrs, jac_wall_indx, 
	jac_part_indx, comp_num, RO2_indx, comp_list, Pybel_objects, eqn_num, 
	comp_namelist, Jlen, 
	rindx_aq, rstoi_aq, pindx_aq, pstoi_aq, reac_coef_aq, 
	nreac_aq, nprod_aq, jac_stoi_aq, 
	jac_den_indx_aq, njac_aq, jac_indx_aq, 				
	y_arr_aq, y_rind_aq, uni_y_rind_aq, y_pind_aq, 
	uni_y_pind_aq, reac_col_aq, prod_col_aq, rstoi_flat_aq, pstoi_flat_aq, 
	rr_arr_aq, rr_arr_p_aq, comp_xmlname, comp_smil] = eqn_pars.extr_mech(sch_name, 
	chem_sch_mrk, xml_name, photo_path, con_infl_nam, int_tol, wall_on, 
	(num_sb+wall_on), const_comp)
	
	# set initial concentrations (molecules/cc)
	[y, H2Oi, y_mw, num_comp, Cfactor, indx_plot, corei, dydt_vst, comp_namelist, 
	inj_indx, core_diss, Psat_water, 
	nuci, nrec_steps, seedi] = init_conc.init_conc(comp_num, comp0, y0, temp[0], RH, 
	Pnow, Pybel_objects, 0, pconc, dydt_trak, tot_time, save_step, rindx_g, 
	pindx_g, eqn_num[0], nreac_g, nprod_g, 
	comp_namelist, Compt, seed_name,
	seed_mw, core_diss, nuc_comp, comp_xmlname, comp_smil)
	
	
	# dump new pickle file ready for plotting script to use
	list_vars = [sav_nam, sch_name, indx_plot, comp0]
	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
	with open(input_by_sim, 'wb') as f: # the file to be used for pickling
		pickle.dump(list_vars, f) # pickle
		f.close() # close

	# get component properties
	[Psat, y_dens, Psat_Pa] = prop_calc.prop_calc(comp_list, Pybel_objects, temp[0], H2Oi, 
		num_comp, Psat_water, vol_comp, volP, 0, corei, pconc,
		uman_up, core_dens, comp_namelist, 0, nuci, nuc_comp, num_sb, dens_comp, dens,
		seed_name)

	# prepare for the calcuation of partitioning variables
	[mfp, accom_coeff, therm_sp, surfT, Cw, act_coeff, 
		R_gas, NA, diff_vol, Dstar_org] = partit_var_prep.prep(y_mw, 
		temp[0], num_comp, 0, Cw, act_comp, act_user, accom_comp, 
		accom_coeff_user, comp_namelist, num_sb, num_sb, Pnow, 
		Pybel_objects, comp_smil)
	count = 0
	
	# prepare particle phase and wall
	[y, N_perbin, x, Varr, Vbou, rad0, Vol0, rbou, MV, num_sb, nuc_comp, 
	rbou00, ub_rad_amp, np_sum] = pp_intro.pp_intro(y, num_comp, Pybel_objects, temp[0],
	 H2Oi, mfp, accom_coeff, y_mw, surfT, RH, siz_str, num_sb, lowsize, 
		uppsize, pmode, pconc, pconct, nuc_comp, 0, std, mean_rad, 
		therm_sp, y_dens, Psat, core_diss, kw, space_mode, seedVr,
		comp_namelist, act_coeff, wall_on, partit_cutoff, Pnow, seedi)
	
	print('Calling integration routine, starting timer')
	st_time = time.time()
	
	# solve problem
	[trec, yrec, dydt_vst, Cfactor_vst, Nres_dry, Nres_wet, x2, rbou_rec] = ode_updater.ode_updater(update_stp, 
		tot_time, save_step, y, rindx_g, 
		pindx_g, rstoi_g, pstoi_g, nreac_g, nprod_g, jac_stoi_g, njac_g, 
		jac_den_indx_g, jac_indx_g, RO2_indx, H2Oi, temp, tempt, 
		Pnow, light_stat, light_time, daytime, lat, lon, af_path, 
		dayOfYear, photo_path, Jlen, con_infl_C, nrec_steps, 
		dydt_vst, siz_str, num_sb, num_comp, seedi, seed_name, seedVr, 
		core_diss, Psat, mfp, therm_sp,  
		accom_coeff, y_mw, surfT, R_gas, NA, y_dens, 
		x, Varr, act_coeff, Cw, kw, Cfactor, tf, light_ad, y_arr_g,
		y_rind_g, 
		uni_y_rind_g, y_pind_g, uni_y_pind_g, reac_col_g, prod_col_g, 
		rstoi_flat_g, pstoi_flat_g, rr_arr_g, rr_arr_p_g, rowvals, 
		colptrs, wall_on, jac_wall_indx, jac_part_indx, Vbou, 
		N_perbin, Vol0, rad0, np_sum, new_partr, nucv1, nucv2, 
		nucv3, nuc_comp, nuc_ad, RH, coag_on, inflectDp, pwl_xpre, 
		pwl_xpro, inflectk, ChamR, Rader, p_char, e_field, 
		injectt, inj_indx, Ct, pmode, pconc, pconct, mean_rad, lowsize, 
		uppsize, std, rbou, const_infl_t, MV,
		rindx_aq, 
		pindx_aq, rstoi_aq, pstoi_aq, nreac_aq, nprod_aq, jac_stoi_aq, njac_aq, 
		jac_den_indx_aq, jac_indx_aq, y_arr_aq,
		y_rind_aq, 
		uni_y_rind_aq, y_pind_aq, uni_y_pind_aq, reac_col_aq, prod_col_aq, 
		rstoi_flat_aq, pstoi_flat_aq, rr_arr_aq, rr_arr_p_aq, eqn_num,
		partit_cutoff, diff_vol, Dstar_org, corei, ser_H2O)
	
	time_taken = time.time()-st_time
	print('Simulation complete, wall clock time elapsed since first call to solver: ', time_taken, ' s')		

	print('Saving results')
	# save results
	save.saving(sch_name, yrec, Nres_dry, Nres_wet, trec, sav_nam, 
		dydt_vst, num_comp, Cfactor_vst, 0, 
		num_sb, comp_namelist, dydt_trak, y_mw, MV, time_taken, 
		seed_name, x2, rbou_rec, wall_on, space_mode, rbou00, ub_rad_amp)
	print('Saved results')
	
	
	return()
