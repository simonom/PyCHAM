##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2022 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
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
'''opens and returns user-defined variables from pickle file'''
# acts as the conduit between the user input process and the 
# software collecting the inputs

import pickle
import os
import numpy as np

# define function
def share(self):

	# inputs: ---------------------------------------------------------------------
	# -------------------------------------------------------------------------------

	# path to store for variables
	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
	with open(input_by_sim, 'rb') as pk:
		[sav_nam, update_stp, 
		tot_time, comp0, y0, RH, RHt, Press,
		Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
		save_step, Compt, injectt, Ct, seed_name,
		seed_mw, seed_diss, seed_dens, seedx,
		con_infl_t, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
		accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
		nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
		inflectk, chamSA, Rader, p_char, e_field, partit_cutoff, ser_H2O, 
		wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, z_prt_coeff, 
		chamV] = pickle.load(pk)
		pk.close()
		
		
		# convert chamber surface area (m2) to spherical equivalent radius (m)
		# (below eq. 2 in Charan et al. 2018, doi.org/10.1080/02786826.2018.1474167)
		chamR = (chamSA/(4.*np.pi))**0.5

		# in addition to other variables already saved to self, 
		# ensure that every variable passed to middle is also available in self.
		# this allows saving of all initial states of variables
		self.sav_nam_orig = sav_nam
		self.update_stp_orig = update_stp
		self.tot_time_orig = tot_time
		self.comp0_orig = comp0
		self.y0_orig = y0
		self.RH_orig = RH
		self.RHt_orig = RHt
		self.Press_orig = Press
		self.Cw_orig = Cw
		self.kw_orig = kw
		self.siz_stru_orig = siz_stru
		self.num_sb_orig = num_sb
		self.pmode_orig = pmode
		self.pconc_orig = pconc
		self.pconct_orig = pconct
		self.lowsize_orig = lowsize
		self.uppsize_orig = uppsize
		self.space_mode_orig = space_mode
		self.std_orig = std
		self.mean_rad_orig = mean_rad
		self.save_step_orig = save_step
		self.Compt_orig = Compt
		self.injectt_orig = injectt
		self.Ct_orig = Ct
		self.seed_name_orig = seed_name
		self.seed_mw_orig = seed_mw
		self.seed_diss_orig = seed_diss
		self.seed_dens_orig = seed_dens
		self.seedx_orig = seedx
		self.con_infl_t_orig = con_infl_t
		self.dens_comp_orig = dens_comp
		self.dens_orig = dens
		self.vol_comp_orig = vol_comp
		self.volP_orig = volP
		self.act_comp_orig = act_comp
		self.act_user_orig = act_user
		self.accom_comp_orig = accom_comp
		self.accom_val_orig = accom_val
		self.uman_up_orig = uman_up
		self.int_tol_orig = int_tol
		self.new_partr_orig = new_partr
		self.nucv1_orig = nucv1
		self.nucv2_orig = nucv2
		self.nucv3_orig = nucv3
		self.nuc_comp_orig = nuc_comp
		self.nuc_ad_orig = nuc_ad
		self.coag_on_orig = coag_on
		self.inflectDp_orig = inflectDp
		self.pwl_xpre_orig = pwl_xpre
		self.pwl_xpro_orig = pwl_xpro
		self.inflectk_orig = inflectk
		self.chamR_orig = chamR
		self.Rader_orig = Rader
		self.p_char_orig = p_char
		self.e_field_orig = e_field
		self.partit_cutoff_orig = partit_cutoff
		self.ser_H2O_orig = ser_H2O
		self.wat_hist_orig = wat_hist
		self.drh_str_orig = drh_str
		self.erh_str_orig = erh_str
		self.pcont_orig = pcont
		self.Vwat_inc_orig = Vwat_inc
		self.seed_eq_wat_orig = seed_eq_wat
		self.z_prt_coeff_orig = z_prt_coeff
		self.chamSA_orig = chamSA
		self.chamV_orig = chamV
		self.light_stat_orig = self.light_stat
		self.light_time_orig = self.light_time
		self.daytime_orig = self.daytime
		self.lat_orig = self.lat
		self.lon_orig = self.lon
		self.af_path_orig = self.af_path
		self.dayOfYear_orig = self.dayOfYear
		self.photo_path_orig = self.photo_path
		self.con_infl_C_orig = self.con_infl_C
		self.tf_orig = self.tf
		self.light_ad_orig = self.light_ad
		self.sch_name_orig = self.sch_name
		self.dydt_trak_orig = self.dydt_trak
		self.inname_orig = self.inname
		self.dil_fac_orig = self.dil_fac
		self.tf_UVC_orig = self.tf_UVC

	return(sav_nam, update_stp, 
		tot_time, comp0, y0, RH, RHt, Press,
		Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, 
		std, mean_rad, save_step, Compt, injectt, Ct, seed_name,
		seed_mw, seed_diss, seed_dens, seedx,
		con_infl_t, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
		accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
		nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
		inflectk, chamR, Rader, p_char, e_field, partit_cutoff, ser_H2O, 
		wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, z_prt_coeff, 
		chamSA, chamV)