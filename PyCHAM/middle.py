########################################################################
#                                                                      #
# Copyright (C) 2018-2025                                              #
# Simon O'Meara : simon.omeara@manchester.ac.uk                        #
#                                                                      #
# All Rights Reserved.                                                 #
# This file is part of PyCHAM                                          #
#                                                                      #
# PyCHAM is free software: you can redistribute it and/or modify it    #
# under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or    #
# (at  your option) any later version.                                 #
#                                                                      #
# PyCHAM is distributed in the hope that it will be useful, but        #
# WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# General Public License for more details.                             #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with PyCHAM.  If not, see <http://www.gnu.org/licenses/>.      #
#                                                                      #
########################################################################
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
import ode_updater_su
import tot_in # preparing record of gas-phase influxes

def middle(self): # define function
	
	# inputs: -----------------------------------------------------
	# self - reference to program
	# -------------------------------------------------------------
	st_time = time.time()
	# get required inputs
	[y0, siz_str, num_sb,
		lowsize, uppsize, std,
		Compt, injectt, Ct, seed_mw, seed_dens, 
		dens_comp, dens, vol_comp, 
		volP, act_comp, act_user, accom_comp, 
		accom_coeff_user, uman_up, 
		int_tol, new_partr, coag_on, 
		inflectDp, pwl_xpre, pwl_xpro, inflectk, ChamR, 
		Rader, p_char, 
		e_field, ser_H2O, wat_hist, drh_str, 
		erh_str, z_prt_coeff] = ui.share(self)
	
	# if not skipping parsing of chemical scheme	
	if (self.pars_skip == 0 or self.pars_skip == 2 or self.pars_skip == 4): 
		# parse the chemical scheme equation file to convert 
		# equations into usable code
		[rowvals, colptrs, comp_num, 
		Jlen, erf, err_mess, 
		self] = eqn_pars.extr_mech(int_tol, 
			(num_sb+self.wall_on), 
			drh_str, erh_str, self)
	
	
	# if skipping parsing of chemical scheme
	if (self.pars_skip == 1 or self.pars_skip == 3):
		[rowvals, colptrs, comp_num, Jlen, 
		erf, err_mess] = eqn_pars_skipper.eqn_pars_skipper(self)

	# if needed, then run operations to produce variable checker 
	# plot from the simulate tab
	if (self.testf == 4):
		import var_checker
		[err_mess, erf] = var_checker.var_checker(Jlen, 
			err_mess, erf, self)
		
	# if error raised, then tell GUI to display and to stop program
	if (erf == 1):
		yield err_mess
	
	# set initial concentrations (# molecules/cm3)
	[y, H2Oi, y_mw, num_comp, Cfactor, indx_plot, corei, 
	inj_indx, nuci, nrec_steps, erf, err_mess, NOi, HO2i, NO3i, y0,
	self] = init_conc.init_conc(comp_num, 
	y0, 0, self.eqn_num[0], Compt,
	seed_mw, self)
	
	
	# if error raised, then tell GUI to display it and to stop 
	# programme
	if (erf == 1):
		yield err_mess
	
	tempt_cnt = 0 # count on chamber temperatures
	
	# if not skipping component properties
	if (self.pars_skip == 0 or self.pars_skip == 2 or self.pars_skip == 3 or self.pars_skip == 4): 
		# get component properties
		[self, err_mess, erf] = prop_calc.prop_calc(H2Oi, 
			num_comp, vol_comp, volP, 0, corei,
			uman_up, seed_dens, 0, nuci, dens_comp, dens,
			y_mw, tempt_cnt, self)
	
	# if error raised, then tell GUI to display and stop program
	if (erf == 1):
		yield err_mess
	
	# prepare for the calculation of partitioning variables
	[mfp, accom_coeff, therm_sp, surfT, act_coeff, 
		R_gas, NA, diff_vol, Dstar_org, err_mess, 
		self] = partit_var_prep.prep(y_mw, 
		self.TEMP[0], num_comp, act_comp, act_user, accom_comp, 
		accom_coeff_user, num_sb+self.wall_on, num_sb, self)
	
	# if error raised or in testing mode then stop
	if (err_mess != ''):
		yield err_mess
	
	# prepare particle phase
	[y, N_perbin, x, Varr, Vbou, rad0, Vol0, rbou, MV, num_sb, 
	rbou00, ub_rad_amp, np_sum] = pp_intro.pp_intro(y, num_comp, 
	self.TEMP[0], H2Oi, mfp, accom_coeff, y_mw, surfT, siz_str, 
	num_sb, lowsize, uppsize, 0, std, 
	therm_sp, act_coeff, seed_mw, R_gas, self)
	
	# estimate total inputs of emitted components (ug/m^3)
	[tot_in_res, Compti, tot_in_res_indx] = tot_in.tot_in(y0, 
		Cfactor, y_mw, Compt, self)
	
	# in case user has specified to spin-up simulation
	if (self.spin_up > 0):

		# spin-up problem
		[y, N_perbin, x, Varr, rbou, 
		Vbou] = ode_updater_su.ode_updater_su(y, H2Oi, 
		Jlen, nrec_steps, 
		siz_str, num_sb, num_comp, 
		mfp, therm_sp,  
		accom_coeff, y_mw, surfT, R_gas, NA, 
		x, Varr, act_coeff, Cfactor, rowvals, colptrs, Vbou, 
		N_perbin, Vol0, rad0, np_sum, new_partr, 
		nuci, coag_on, inflectDp, pwl_xpre, 
		pwl_xpro, inflectk, ChamR, Rader, p_char, e_field, 
		injectt, inj_indx, Ct, 
		lowsize, uppsize, std, rbou, MV,
		diff_vol, Dstar_org, corei, ser_H2O, 
		rbou00, ub_rad_amp, indx_plot,
		wat_hist, NOi, 
		HO2i, NO3i, z_prt_coeff, tot_in_res,
		Compti, tot_in_res_indx, tempt_cnt, self,
		vol_comp, volP)

	# solve problem
	for prog in ode_updater.ode_updater(y, H2Oi, 
		Jlen, nrec_steps, 
		siz_str, num_sb, num_comp, mfp, therm_sp,  
		accom_coeff, y_mw, surfT, R_gas, NA, 
		x, Varr, act_coeff, Cfactor, rowvals, colptrs, Vbou, 
		N_perbin, Vol0, rad0, np_sum, new_partr, 
		nuci, coag_on, inflectDp, pwl_xpre, 
		pwl_xpro, inflectk, ChamR, Rader, p_char, e_field, 
		injectt, inj_indx, Ct, 
		lowsize, uppsize, std, rbou, MV,
		diff_vol, Dstar_org, corei, ser_H2O, 
		rbou00, ub_rad_amp, indx_plot,
		wat_hist, NOi, 
		HO2i, NO3i, z_prt_coeff, tot_in_res,
		Compti, tot_in_res_indx, tempt_cnt, self, 
		vol_comp, volP):

		yield prog # update progress bar	

	return()
