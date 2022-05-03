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
def share():

	# inputs: ---------------------------------------------------------------------
	# -------------------------------------------------------------------------------

	# path to store for variables
	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
	with open(input_by_sim, 'rb') as pk:
		[sav_nam, chem_sch_mark, update_stp, 
		tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
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

	return(sav_nam, chem_sch_mark, update_stp, 
		tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
		Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, 
		std, mean_rad, save_step, Compt, injectt, Ct, seed_name,
		seed_mw, seed_diss, seed_dens, seedx,
		con_infl_t, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
		accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
		nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
		inflectk, chamR, Rader, p_char, e_field, partit_cutoff, ser_H2O, 
		wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, z_prt_coeff, 
		chamSA, chamV)