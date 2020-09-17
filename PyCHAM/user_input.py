'''opens and returns user-defined variables from pickle file'''
# acts as the conduit between the user input process and the 
# software collecting the inputs

#import ui_check
import pickle
import os
import ui_check

# define function
def share(source):

	# inputs: ---------------------------------------------------------------------
	# source - marker for calling script
	# -----------------------------------------------------------------------------

	if (source == 0): # when called from PyCHAM script

		# path to store for variables
		input_by_sim = str(os.getcwd() + '/var_store.pkl')
		with open(input_by_sim, 'rb') as pk:
			[sav_nam, sch_name, chem_sch_mark, xml_name, update_stp, 
			tot_time, comp0, y0, temp, tempt, RH, Press, wall_on,
			Cw, kw, num_sb, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
			save_step, const_comp, Compt, injectt, Ct, seed_name,
			seed_mw, core_diss, seed_dens,
			light_stat, light_time, daytime, lat, lon, af_path, 
			dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, 
			dydt_trak, vol_comp, volP, act_comp, act_user, 
			accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
			nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
			inflectk, chamSA, Rader, p_char, e_field, dil_fac] = pickle.load(pk)
			pk.close()
			
			os.remove(input_by_sim) # delete pickle file

			# check on inputs
			[wall_on, pconc, lowsize, std, mean_rad, new_partr, chamR, chem_sch_mark, 
			af_path, int_tol, update_stp, tot_time] = ui_check.ui_check(sav_nam, sch_name, wall_on, 0, 				num_sb, pconc, pconct, lowsize, std, mean_rad, new_partr, chamSA, chem_sch_mark, af_path, 
			int_tol, update_stp, tot_time, RH, uman_up)
	
		return(sav_nam, sch_name, chem_sch_mark, xml_name, update_stp, 
			tot_time, comp0, y0, temp, tempt, RH, Press, wall_on,
			Cw, kw, num_sb, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
			save_step, const_comp, Compt, injectt, Ct, seed_name,
			seed_mw, core_diss, seed_dens,
			light_stat, light_time, daytime, lat, lon, af_path, 
			dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, 
			dydt_trak, vol_comp, volP, act_comp, act_user, 
			accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
			nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
			inflectk, chamR, Rader, p_char, e_field, dil_fac)

	if (source == 1): # when called from plotting script
		with open(input_by_sim, 'rb') as pk:
			[sav_name, sch_name, indx_plot, Comp0] = pickle.load(pk)
			pk.close()

			os.remove(input_by_sim) # delete pickle file

		return(sav_name, sch_name, indx_plot, Comp0)
