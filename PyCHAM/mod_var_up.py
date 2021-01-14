'''module to update displayed model variables in PyCHAM GUI'''

import os
import pickle

def mod_var_up(self):

	# inputs: ------------------------
	# self - reference to GUI
	# ---------------------------------
	
	# open model variables
	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
		
	with open(input_by_sim, 'rb') as pk:
		[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
		tot_time, comp0, y0, temp, tempt, RH, Press, wall_on,
		Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
		save_step, const_comp, Compt, injectt, Ct, seed_name,
		seed_mw, seed_diss, seed_dens, seedVr,
		light_stat, light_time, daytime, lat, lon, af_path, 
		dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, 
		dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
		accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
		nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
		inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O] = pickle.load(pk)
	pk.close()

	
	# contents of model variables update
	self.l9a.setText(sav_nam);
	self.l10a.setText((str(chem_sch_mark)).replace('\'', '').replace(' ', '')[1:-1])
	self.l11a.setText((str(tot_time)).replace('\'', '').replace(' ', ''))
	self.l12a.setText((str(update_stp)).replace('\'', '').replace(' ', ''))
	self.l13a.setText((str(save_step)).replace('\'', '').replace(' ', ''))#
	self.l14a.setText((str(temp)).replace('\'', '').replace(' ', '')[1:-1])
	self.l15a.setText((str(tempt)).replace('\'', '').replace(' ', '')[1:-1])
	self.l16a.setText((str(RH)).replace('\'', '').replace(' ', ''))
	self.l17a.setText((str(Press)).replace('\'', '').replace(' ', ''))
	self.l18a.setText((str(siz_stru)).replace('\'', '').replace(' ', ''))
	self.l19a.setText((str(num_sb)).replace('\'', '').replace(' ', ''))
	self.l20a.setText((str(pconc)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l21a.setText((str(pconct)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l22a.setText((str(seed_mw)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l23a.setText((str(seed_diss)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l24a.setText((str(seed_dens)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l25a.setText((str(seed_name)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l26a.setText((str(seedVr)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l27a.setText((str(lowsize)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l28a.setText((str(uppsize)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l29a.setText((str(space_mode)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l30a.setText((str(std)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l31a.setText((str(mean_rad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l32a.setText((str(new_partr)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l33a.setText((str(nucv1)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l34a.setText((str(nucv2)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l35a.setText((str(nucv3)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l36a.setText((str(nuc_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l37a.setText((str(nuc_ad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l38a.setText((str(ser_H2O)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l40a.setText((str(comp0)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l41a.setText((str(y0)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l42a.setText((str(con_infl_nam)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l43a.setText((str(con_infl_C)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l44a.setText((str(con_infl_t)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l45a.setText((str(const_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l46a.setText((str(Compt)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l47a.setText((str(injectt)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l48a.setText((str(Ct)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		
	self.show()
	return()