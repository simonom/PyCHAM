##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2023 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
#                                                                                       											 #
#    All Rights Reserved.                                                                									 #
#    This file is part of PyCHAM                                                         									 #
#                                                                                        											 #
#    PyCHAM is free software: you can redistribute it and/or modify it under              						 #
#    the terms of the GNU General Public License as published by the Free Software       					 #
#    Foundation, either version 3 of the License, or (at your option) any later          						 #
#    version.                                                                            										 #
#                                                                                        											 #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT         						 #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       			 #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              				 #
#    details.                                                                            										 #
#                                                                                        											 #
#    You should have received a copy of the GNU General Public License along with        					 #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                 							 #
#                                                                                        											 #
##########################################################################################
'''calls modules to setup simulation'''
# the modules necessary to setup a simulation are called here

import eqn_pars
import eqn_pars_skipper
import user_input as ui
import init_conc
import prop_calc
import partit_var_prep
import pp_intro
import time
import save
import os
import ode_updater
import tot_in # preparing record of gas-phase influxes

def middle(self): # define function
	
	# inputs: -------------------------------------------------------
	# self - reference to program
	# ---------------------------------------------------------------
	st_time = time.time()
	# get required inputs
	[sav_nam, comp0, y0, Pnow,
		siz_str, num_sb, pmode, pconc, pconct,
		lowsize, uppsize, space_mode, std, mean_rad,
		Compt, injectt, Ct, seed_name, seed_mw, 
		core_diss, seed_dens, seedx, 
		dens_comp, dens, vol_comp, 
		volP, act_comp, act_user, accom_comp, accom_coeff_user, uman_up, 
		int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, 
		inflectDp, pwl_xpre, pwl_xpro, inflectk, ChamR, Rader, p_char, 
		e_field, partit_cutoff, ser_H2O, wat_hist, drh_str, 
		erh_str, pcont, z_prt_coeff, 
		chamSA, chamV] = ui.share(self)
	
	if (self.pars_skip == 0): # if not skipping parsing of chemical scheme
		# parse the chemical scheme equation file to convert equations
		# into usable code
		[rowvals, colptrs, jac_wall_indx, 
		jac_part_indx, jac_extr_indx, comp_num, rel_SMILES, 
		Pybel_objects, Jlen, comp_xmlname, comp_smil, erf, err_mess, 
		self] = eqn_pars.extr_mech(int_tol, (num_sb+self.wall_on), drh_str, erh_str, sav_nam,
		pcont, self)

	if (self.pars_skip == 1): # if skipping parsing of chemical scheme
		[rowvals, colptrs, jac_wall_indx, 
		jac_part_indx, jac_extr_indx, comp_num, rel_SMILES, 
		Pybel_objects, Jlen, comp_xmlname, comp_smil, erf, 
		err_mess] = eqn_pars_skipper.eqn_pars_skipper(self)

	# if needed, then run operations to produce variable checker plot 
	# from the simulate tab
	if (self.testf == 4):
		import var_checker
		[err_mess, erf] = var_checker.var_checker(Jlen, err_mess, erf, self)
	
	# if error raised, then tell GUI to display and to stop program
	if (erf == 1):
		yield err_mess
	
	# set initial concentrations (# molecules/cm3)
	[y, H2Oi, y_mw, num_comp, Cfactor, indx_plot, corei, 
	inj_indx, core_diss, Psat_water, 
	nuci, nrec_steps, erf, err_mess, NOi, HO2i, NO3i, self, 
	rel_SMILES] = init_conc.init_conc(comp_num, 
	comp0, y0, Pnow, Pybel_objects, 0, pconc, self.eqn_num[0], Compt, seed_name,
	seed_mw, core_diss, nuc_comp, comp_xmlname, comp_smil, rel_SMILES, self)
	# if error raised, then tell GUI to display it and to stop programme
	if (erf == 1):
		yield err_mess
	
	tempt_cnt = 0 # count on chamber temperatures
	
	if (self.pars_skip == 0): # if not skipping component properties
		# get component properties
		[self, err_mess, erf] = prop_calc.prop_calc(rel_SMILES, Pybel_objects, 
			H2Oi, num_comp, Psat_water, vol_comp, volP, 0, corei, pconc,
			uman_up, seed_dens, 0, nuci, nuc_comp, dens_comp, dens,
			seed_name, y_mw, tempt_cnt, self)
	
	# if error raised, then tell GUI to display and stop program
	if (erf == 1):
		yield err_mess
	
	# prepare for the calculation of partitioning variables
	[mfp, accom_coeff, therm_sp, surfT, act_coeff, 
		R_gas, NA, diff_vol, Dstar_org, err_mess, self] = partit_var_prep.prep(y_mw, 
		self.TEMP[0], num_comp, act_comp, act_user, accom_comp, 
		accom_coeff_user, num_sb+self.wall_on, num_sb, Pnow, 
		Pybel_objects, comp_smil, self)
	
	if (err_mess != ''): # if error raised or in testing mode then stop
		yield err_mess
	
	# prepare particle phase
	[y, N_perbin, x, Varr, Vbou, rad0, Vol0, rbou, MV, num_sb, nuc_comp, 
	rbou00, ub_rad_amp, np_sum] = pp_intro.pp_intro(y, num_comp, Pybel_objects, self.TEMP[0],
	 H2Oi, mfp, accom_coeff, y_mw, surfT, siz_str, num_sb, lowsize, 
		uppsize, pmode, pconc, pconct, nuc_comp, 0, std, mean_rad, 
		therm_sp, core_diss, space_mode, seedx,
		act_coeff, partit_cutoff, Pnow, 
		pcont, seed_mw, R_gas, self)
	
	# estimate total inputs of emitted components (ug/m3)
	[tot_in_res, Compti, tot_in_res_indx] = tot_in.tot_in(y0, Cfactor, comp0, y_mw, Compt, self)
	
	# solve problem
	for prog in ode_updater.ode_updater(y, H2Oi, 
		Pnow, Jlen, nrec_steps, 
		siz_str, num_sb, num_comp, seed_name, seedx, 
		core_diss, mfp, therm_sp,  
		accom_coeff, y_mw, surfT, R_gas, NA, 
		x, Varr, act_coeff, Cfactor, rowvals, 
		colptrs, jac_wall_indx, jac_part_indx, jac_extr_indx, Vbou, 
		N_perbin, Vol0, rad0, np_sum, new_partr, nucv1, nucv2, 
		nucv3, nuci, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, 
		pwl_xpro, inflectk, ChamR, Rader, p_char, e_field, 
		injectt, inj_indx, Ct, pmode, pconc, pconct, mean_rad, lowsize, 
		uppsize, std, rbou, MV,
		partit_cutoff, diff_vol, Dstar_org, corei, ser_H2O, 
		sav_nam, space_mode, 
		rbou00, ub_rad_amp, indx_plot, comp0, rel_SMILES,
		wat_hist, Pybel_objects, pcont, NOi, 
		HO2i, NO3i, z_prt_coeff, tot_in_res,
		Compti, tot_in_res_indx, chamSA, chamV, tempt_cnt, self, vol_comp, volP):

		yield prog # update progress bar	

	
	return()
