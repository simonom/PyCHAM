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
		tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
		Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
		save_step, const_comp, Compt, injectt, Ct, seed_name,
		seed_mw, seed_diss, seed_dens, seedVr,
		light_stat, light_time, daytime, lat, lon, af_path, 
		dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, 
		dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
		accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
		nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
		inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, 
		wat_hist, drh_str, erh_str, pcont] = pickle.load(pk)
	pk.close()

	
	# contents of model variables update
	self.l9a.setText(sav_nam);
	self.l10a.setText((str(chem_sch_mark)).replace('\'', '').replace(' ', '')[1:-1])
	self.l11a.setText((str(tot_time)).replace('\'', '').replace(' ', ''))
	self.l12a.setText((str(update_stp)).replace('\'', '').replace(' ', ''))
	self.l13a.setText((str(save_step)).replace('\'', '').replace(' ', ''))#
	self.l13_1a.setText((str(uman_up)).replace('\'', '').replace(' ', ''))
	self.l13_2a.setText((str(int_tol)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	
	
	self.l14a.setText((str(temp)).replace('\'', '').replace(' ', '')[1:-1])
	self.l15a.setText((str(tempt)).replace('\'', '').replace(' ', '')[1:-1])
	self.l16a.setText((str(RH.tolist())).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l16c.setText((str(RHt.tolist())).replace('\'', '').replace(' ', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l17a.setText((str(Press)).replace('\'', '').replace(' ', ''))
	self.l17_1a.setText((str(dil_fac)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
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
	self.l38_1a.setText((str(coag_on)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l38_2a.setText((str(wat_hist)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l38_3a.setText((str(drh_str)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l38_4a.setText((str(erh_str)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l38_5a.setText((str(pcont)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	
	self.l40a.setText((str(comp0)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l41a.setText((str(y0)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l42a.setText((str(con_infl_nam)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l43a.setText((str(con_infl_C)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l44a.setText((str(con_infl_t.tolist())).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l45a.setText((str(const_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l46a.setText((str(Compt)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l47a.setText((str(injectt.tolist())).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l48a.setText((str(Ct)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l49a.setText((str(light_stat.tolist())).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l50a.setText((str(light_time.tolist())).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l51a.setText((str(daytime)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l52a.setText((str(lat)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l53a.setText((str(lat)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l54a.setText((str(af_path)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l55a.setText((str(photo_path)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l56a.setText((str(dayOfYear)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l57a.setText((str(tf)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l58a.setText((str(light_ad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l59a.setText((str(wall_on)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l60a.setText((str(Cw)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l61a.setText((str(kw)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l62a.setText((str(inflectDp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l63a.setText((str(pwl_xpre)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l64a.setText((str(pwl_xpro)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l65a.setText((str(inflectk)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l66a.setText((str(chamSA)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l67a.setText((str(Rader)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l68a.setText((str(p_char)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l69a.setText((str(e_field)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l70a.setText((str(dydt_trak)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l71a.setText((str(dens_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l72a.setText((str(dens)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l73a.setText((str(vol_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l74a.setText((str(volP)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l75a.setText((str(act_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l76a.setText((str(act_user)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l77a.setText((str(accom_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	self.l78a.setText((str(accom_val)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
	
	self.show()
	return()