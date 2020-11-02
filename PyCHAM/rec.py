'''module to record required results'''
# keeps track of model variables

import numpy as np
import scipy.constants as si

def rec(save_cnt, trec, yrec, dydt_vst, Cfactor_vst, y, sumt,
	rindx, rstoi, rrc, pindx, pstoi, nprod, 
	nreac, num_sb, num_comp, pconc, core_diss, Psat, kelv_fac, 
	kimt, kw, Cw, act_coeff, Cfactor, Nres_dry, Nres_wet, x2, x,
	MV, H2Oi, Vbou, rbou, wall_on, rbou_rec, seedi):
	
	# inputs: ------------------------------------------------------------
	# save_cnt - count on saving steps
	# trec - the time through simulation record (s)
	# yrec - concentration record (molecules/cc (air))
	# dydt_vst - record of change tendencies (molecules/cc/s)
	# Cfactor_vst - record of the conversion factor for moleculec/ss 
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
	# pconc - particle concentrations now (#/cc (air))
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
	# MV - molar volume of components (um3/mol)
	# H2Oi - index for water
	# Vbou - volume boundaries per size bin (um3)
	# rbou - size bin radius boundaries (um)
	# wall_on - marker for whether wall turned on
	# rbou_rec - size bin radius boundary record (um)
	# seedi - index of seed component
	# --------------------------------------------------------------------

	trec[save_cnt] = sumt # track recording times (s)
	yrec[save_cnt, :] = y # track concentrations (molecules/cc/s)
	# track conversion factor (ppb/molecules/cc)
	Cfactor_vst[save_cnt, 0] = Cfactor
	
	# single particle radius (um) at size bin centre 
	# including contribution of water
	if (num_sb-wall_on > 0):
		x2[save_cnt, :] = x
		# single particle radius boundaries (um) including contribution of water
		rbou_rec[save_cnt, :] = rbou 
	
	# estimate particle number size distributions ----------------------------------
	Vnew = np.zeros((num_sb-wall_on))
	ish = pconc[:, 0]>1.0e-10
	
	if ish.sum()>0:
		# rearrange particle concentrations into size bins in rows, components in columns
		Cn = y[num_comp:num_comp*(num_sb-wall_on+1)].reshape(num_sb-wall_on, num_comp)
		# new volume of single particle per size bin (um3) excluding volume of water	
		Vnew[ish] = (np.sum((Cn[ish, :]/(si.N_A*pconc[ish]))*MV[:, 0]*1.e12, 1)-
				((Cn[ish, H2Oi]/(si.N_A*pconc[ish, 0]))*MV[H2Oi]*1.e12))
		# loop through size bins to find number of particles in each 
		# (# particle/cc (air))
		for Ni in range(0, (num_sb-wall_on)):
			ish = (Vnew>=Vbou[Ni])*(Vnew<Vbou[Ni+1])
			Nres_dry[save_cnt, Ni] = pconc[ish, 0].sum()
		Nres_wet[save_cnt, :] = pconc[:, 0] # record with water
	
	# end of number size distribution part ----------------------------------------


	if len(dydt_vst)>0:
		# record any change tendencies of specified components
		import dydt_rec
		dydt_vst = dydt_rec.dydt_rec(y, rindx, rstoi, rrc, pindx, pstoi, nprod, save_cnt, 
					dydt_vst, nreac, num_sb, num_comp, pconc, core_diss, Psat, kelv_fac, 
					kimt, kw, Cw, act_coeff, seedi)
	save_cnt += 1 # track number of recordings 

	return(trec, yrec, dydt_vst, Cfactor_vst, save_cnt, Nres_dry, Nres_wet, x2, rbou_rec)
