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
'''module to record required results'''
# keeps track of model variables

import numpy as np
import scipy.constants as si
import importlib

def rec(save_cnt, trec, yrec, Cfactor_vst, y, sumt,
	rindx, rstoi, rrc, pindx, pstoi, nprod, 
	nreac, num_sb, num_comp, pconc, core_diss, Psat, kelv_fac, 
	kimt, kw, Cw, act_coeff, Cfactor, Nres_dry, Nres_wet, x2, x,
	MV, H2Oi, Vbou, rbou, rbou_rec, 
	yrec_p2w, C_p2w, cham_env, temp_now, Pnow, tot_in_res, 
	tot_in_res_ft, self):
	
	# inputs: ------------------------------------------------------------
	# save_cnt - count on saving steps
	# trec - the time through simulation record (s)
	# yrec - concentration record (# molecules/cm3 (air))
	# Cfactor_vst - record of the conversion factor for # molecules/cm3 
	# to ppb
	# y - concentrations (molecules/cc (air))
	# sumt - cumulative time through simulation (s)
	# rindx - indices of reactants
	# rstoi - stoichiometries of reactants
	# rrc - reaction rate coefficients
	# pindx - indices of products
	# pstoi - stoichiometries of products
	# nprod - number of products
	# nreac - number of reactants
	# num_sb - number of size bins
	# num_comp - number of components
	# pconc - particle concentrations now (# particles/cm3 (air))
	# core_diss - dissociation constant of seed
	# Psat - pure component saturation vapour pressure 
	#	(molecules/cm3 (air))
	# kelv_fac - Kelvin factor
	# kimt - partitioning coefficient (/s)
	# kw - gas-wall mass transfer coefficient (/s)
	# Cw - effective absorbing mass of wall (molecules/cc (air))
	# act_coeff - activity coefficients
	# Cfactor - conversion factor from molecules/cc to ppb
	# Nres_dry - number size distributions excluding water (#/cc (air))
	# Nres_wet - number size distributions including water (#/cc (air))
	# x2 - radii of single particle (um)
	# x - radii of single particles per size bin now (um)
	# MV - molar volume of components (cm3/mol)
	# H2Oi - index for water
	# Vbou - volume boundaries per size bin (um3)
	# rbou - size bin radius boundaries (um)
	# self.wall_on - marker for whether wall turned on
	# rbou_rec - size bin radius boundary record (um)
	# yrec_p2w - record of concentration of components on 
	#	the wall due to particle-wall deposition (molecules/cc)
	# C_p2w - concentration of components on the wall due to 
	#	particle-wall deposition, stacked by component first then by
	#	size bin (molecules/cc)
	# cham_env - chamber environmental conditions (temperature (K), 
	# pressure (Pa) and relative humdity
	# temp_now - chamber temperature (K)
	# Pnow - chamber pressure (Pa)
	# tot_in_res - cumulative influx of injected components (ug/m3)
	# tot_in_res_ft - record of continuous influx of injected 
	#	components (ug/m3)
	# self - reference to program
	# --------------------------------------------------------------------

	trec[save_cnt] = sumt # track recording times (s)
	yrec[save_cnt, :] = y # track concentrations (# molecules/cm3/s)
	# track conversion factor for gas phase (ppb/# molecules/cm3)
	Cfactor_vst[save_cnt, 0] = Cfactor
	
	# single particle radius (um) at size bin centre 
	# including contribution of water
	if ((num_sb-self.wall_on) > 0): # if particle size bins present
		x2[save_cnt, :] = x
		# single particle radius boundaries (um) including contribution of water
		rbou_rec[save_cnt, :] = rbou
		
		# record component concentrations on the wall due to particle loss to wall
		# (molecules/cc (air))
		yrec_p2w[save_cnt, :] = C_p2w
	
	# estimate particle number size distributions ----------------------------------
	Vnew = np.zeros((num_sb-self.wall_on))
	ish = (pconc[:, 0] > 1.e-10)
	
	if (ish.sum() > 0):
		# rearrange particle concentrations into size bins in rows, components in columns
		Cn = y[num_comp:num_comp*(num_sb-self.wall_on+1)].reshape(num_sb-self.wall_on, num_comp)
		# new volume of single particle per size bin (um3) excluding volume of water	
		Vnew[ish] = (np.sum((Cn[ish, :]/(si.N_A*pconc[ish]))*MV[:, 0]*1.e12, 1)-
				((Cn[ish, H2Oi]/(si.N_A*pconc[ish, 0]))*MV[H2Oi, 0]*1.e12))
		# loop through size bins to find number of particles in each 
		# (# particle/cm3 (air))
		for Ni in range(0, (num_sb-self.wall_on)):
			ish = (Vnew>=Vbou[Ni])*(Vnew<Vbou[Ni+1])
			Nres_dry[save_cnt, Ni] = pconc[ish, 0].sum()
		Nres_wet[save_cnt, :] = pconc[:, 0] # record with water
	
	# end of number size distribution part ----------------------------------------
	
	# chamber environmental conditions ----------------------------------
	
	cham_env[save_cnt, 0] = temp_now # temperature (K)
	cham_env[save_cnt, 1] = Pnow # pressure (Pa)
	cham_env[save_cnt, 2] = y[H2Oi]/Psat[0, H2Oi] # relative humidity (fraction (0-1))
	# --------------------------------------------------------------------------------

	# cumulative influx of injected components (ug/m3)
	tot_in_res_ft[save_cnt+1, :] += tot_in_res

	save_cntf = 1 # flag for increasing number of recordings 

	return(trec, yrec, Cfactor_vst, save_cntf, Nres_dry, Nres_wet, x2, rbou_rec, 
		yrec_p2w, cham_env, tot_in_res_ft)
