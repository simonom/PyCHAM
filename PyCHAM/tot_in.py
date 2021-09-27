'''function to estimate total inputs of components'''
# for estimating the total input of components to the simulation
# based on the user-supplied model variables

import scipy.constants as si # for scientific constants
import numpy as np # for arithmetic

def tot_in(init_conc, Cfac, comp0, comp_namelist, y_mw, con_infl_nam,
		const_infl_t, tot_time, con_infl_C, Compt): # define function

	# inputs: ----------------------------------------------
	# init_conc - initial concentrations (ppb)
	# Cfac - factor to convert ppb to # molecules/cm3
	# comp0 - chemical scheme names of components present 
	#	at simulation start
	# comp_namelist - list of all chemical scheme names
	# y_mw - molar mass of components (g/mol)
	# con_infl_nam - names of components injected after 
	# 	experiment start
	# const_infl_t - times at which injections occur (s)
	# tot_time - total experiment time (s)
	# con_infl_C - concentrations injected after component start (ppb)
	# Compt - name of components injected instantaneously 
	#	after start of experiment
	# ------------------------------------------------------

	# indices of components injected
	tot_in_res_indx = []
	# total concentration of components injected
	tot_in_res_con = []

	ccnt = 0 # count on components

	for cnam in comp0: # loop through components present initially

		
		ci = comp_namelist.index(cnam) # index within all components
		tot_in_res_indx.append(ci) # remember component index

		# initial input ug/m3, note *1e12 converts g/cm3 to ug/m3
		Czero = ((init_conc[ccnt]*Cfac)/si.N_A)*y_mw[ci]*1.e12
		tot_in_res_con.append(Czero[0])
		
		ccnt += 1 # count on components

	# instantaneous injection -------------------------------------------
	Compti = [] # indices in record for components injected instantaneously
	# index for this record of components injected instantaneously after experiment start
	# loop through chemical scheme names of components injected instantaneously
	for cnam in Compt:
		if comp_namelist.index(cnam) in tot_in_res_indx:
			Compti.append(tot_in_res_indx.index(comp_namelist.index(cnam)))
			continue
		else:
			Compti.append(len(tot_in_res_nam))
			tot_in_res_indx.append(comp_namelist.index(cnam))
			tot_in_res_con.append(0.)

	# continuous injection ---------------------------------------------
	cont_inf_reci = [] # index inside record
	cont_inf_i = [] # index inside all concentrations
	# loop through components injected after experiment start
	for cnam in con_infl_nam: # loop through components after experiment start
		# index where this component ordered among all components
		cont_inf_i.append(comp_namelist.index(cnam))

		if comp_namelist.index(cnam) in tot_in_res_indx:
			cont_inf_reci.append(tot_in_res_indx.index(comp_namelist.index(cnam)))
			continue
		else:
			cont_inf_reci.append(len(tot_in_res_indx))
			tot_in_res_indx.append(comp_namelist.index(cnam))
			tot_in_res_con.append(0.)
	
	# array to hold all information (component indices in first 
	# column, influxes in second column)
	tot_in_res = np.array((tot_in_res_con))
	
	return(tot_in_res, Compti, cont_inf_reci, cont_inf_i, tot_in_res_indx)