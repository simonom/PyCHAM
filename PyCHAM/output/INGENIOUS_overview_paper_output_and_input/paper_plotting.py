'''code to plot EAC2024 abstract and poster results'''

# set the path to the folder where results are saved
dir_path = '/Users/user/Documents/GitHub/PyCHAM/PyCHAM/output/full_chem_scheme'

def abs_and_pos_plot(dir_path): # define function

	# import dependencies start ---------------------
	import matplotlib.pyplot as plt
	import matplotlib as mpl
	from matplotlib.lines import Line2D
	import matplotlib.ticker as ticker
	from matplotlib.gridspec import GridSpec
	import numpy as np
	import scipy.constants as si
	import os
	import sys
	import platform
	import scipy
	# import dependencies end -----------------------

	import retr_out
	
	# names of primary organic components (both indoor and outdoor (pri_org))
	pri_org_comps = ['pri_org', 'PNNCATCOOH', 'C1213OH', 'NLIMOOH', 'CO25C6CHO']

	# prepare results array, with simulations in rows and in columns:
	# arithemtic mean (over time) ozone concentration, arithmetic mean over
	# time terpene concentrations, secondary organic particulate matter of
	# indoor and outdoor origin, 
	# secondary particulate matter of indoor origin only, 
	# organic particulate matter of outdoor origin, primary
	# organic particulate matter, total organic particulate matter, total
	# inorganic particulate matter and total particulate matter (excluding)
	# water
	res_arr = np.zeros((0, 9))

	# prepare results array to hold fraction of organic particulate mass
	# comprised of: secondary of indoor origin (column 0), organic of
	# outdoor origin (column 1) and primary organic of indoor origin
	# (column 3)
	res_arr_frac_org_pm = np.zeros((0, 3))

	# prepare results array to hold total PM2.5 mass
	# comprised of: secondary of all source origin (column 0),
	# secondary of just indoor origin (column 1), organic of
	# outdoor origin (column 2) and primary organic of indoor origin
	# (column 3) and organic aerosol (sum of previous three) (column 4)
	# and inorganic aerosol of both outdoor and indoor origin (column 5)
	# and total PM2.5 (all components) (column 6)
	res_arr_pm2p5 = np.zeros((0, 7))

	# prepare results array to hold fraction of total PM2.5 mass
	# comprised of: secondary of indoor and outdoor origin (column 0),
	# organic of indoor origin (column 1),
	# organic of outdoor origin (column 2) and primary organic of indoor origin
	# (column 3) and organic aerosol (sum of previous three) (column 4)
	# and inorganic aerosol of both outdoor and indoor origin (column 5)
	res_arr_frac_pm2p5 = np.zeros((0, 6))

	# prepare to hold ventilation rates of each simulation (acr (/hr))
	vent_rate = []

	# prepare to hold simulation names
	sim_names = []

	res_num = 0 # count on simulations

	# prepare to hold (for overview paper): indoor-generated inorganic, indoor-
	# generated primary organics, indoor-generated secondary organics, 
	# outdoor-generated inorganic, outdoor-generated primary organics and 
	# outdoor-generated secondary organics (components in rows and ventilation
	# rates (high and low in columns)), comparing households with same activity
	# but two different ventilation rates
	overv_res = np.zeros((6, 2))
	# prepare to hold time series for the overview paper too (sources in rows, 
	# scenarios in columns and times in third dimension)
	overv_res_ts = np.zeros((6, 2, 145))

	# labels for overview plot components
	overv_labels = ['outdoor-generated\nsecondary organics', 
		'outdoor-generated\nprimary organics', 'outdoor-generated\ninorganic', 
		'indoor-generated\nsecondary organics', 'indoor-generated\nprimary organics', 
		'indoor-generated\ninorganic']
	# bar colours for overview plot
	overv_colours = [(1., 0., 0.), (0., 0., 1.), (0., 1., 0.), (0.5, 0., 0.), 
		(0., 0., 0.5), (0., 0.5, 0.)]

	# loop through results, note that using np.sort
	# means that the first simulation we see of a given
	# air change rate and light transmission factor
	# is the one with no indoor activity
	for resi in np.sort(os.listdir(dir_path)):

		# importing results start -------------------------------
		# if path join needed
		if (dir_path[-1]!= '/' and dir_path[-1] != '\\'):
			res_path = str(dir_path + '/')
		# complete path to this item in this folder
		res_path = str(res_path + resi)

		# create the self object
		self = self_def(res_path)


		# import results
		try:
			# 
			for prog in retr_out.retr_out(self):
				prog = prog

			# include rows for the result of this simulation
			res_arr = np.concatenate((res_arr, 
				np.zeros((1, res_arr.shape[1]))), axis=0)
			res_arr_frac_org_pm = np.concatenate((res_arr_frac_org_pm, 
				np.zeros((1, res_arr_frac_org_pm.shape[1]))), axis=0)
			res_arr_pm2p5 = np.concatenate((res_arr_pm2p5, 
				np.zeros((1, res_arr_pm2p5.shape[1]))), axis=0)
			res_arr_frac_pm2p5 = np.concatenate((res_arr_frac_pm2p5, 
				np.zeros((1, res_arr_frac_pm2p5.shape[1]))), axis=0)
			# store simulation name
			sim_names.append(resi)
		except:
			continue

		# importing results end ---------------------------------
	
		# get required variables from self start ----------------
		wall_on = self.ro_obj.wf
		yrec = np.zeros((self.ro_obj.yrec.shape[0], 
			self.ro_obj.yrec.shape[1]))
		yrec[:, :] = self.ro_obj.yrec[:, :]
		num_comp = self.ro_obj.nc
		# number of size bins (including number of surfaces)
		num_sb = self.ro_obj.nsb
		Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], 
		self.ro_obj.Nrec_wet.shape[1]))
		Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
		Ndry = np.zeros((self.ro_obj.Nrec_dry.shape[0], 
			self.ro_obj.Nrec_dry.shape[1]))
		Ndry[:, :] = self.ro_obj.Nrec_dry[:, :]
		timehr = self.ro_obj.thr
		comp_names = self.ro_obj.names_of_comp
		rel_SMILES = self.ro_obj.rSMILES
		y_MM = self.ro_obj.comp_MW
		H2Oi = self.ro_obj.H2O_ind
		seedi = self.ro_obj.seed_ind
		rbou_rec = np.zeros((self.ro_obj.rad.shape[0], 
			self.ro_obj.rad.shape[1]))
		rbou_rec[:, :] = self.ro_obj.rad[:, :]
		# hydrogen to carbon ratios of components
		HtoC_ratio = self.ro_obj.HyC
		# oxygen to carbon ratio of components
		OtoC_ratio = self.ro_obj.O_to_C
		# total concentration of components influxed 
		# throughout simulation (ug/m3), with component index 
		# (relative to all components) in the first row and its
		# cumulative injected concentrations in following rows
		tot_influx = self.ro_obj.total_influx

		# get required variables from self end ------------------

		# store air changes /hour
		vent_rate.append(float(resi[resi.index('_')+1:resi.index('acr')]))

		# get number of particle size bins
		num_asb = num_sb-wall_on
		
		# need to get particulate matter mass of each component 
		# per simulation time step
		# first, isolate just particle-phase concentrations 
		# (# molecules/cm3) (time in rows, components and 
		# particle size bins in columns)
		yp_molec = yrec[:, num_comp:-(num_comp*wall_on)]

		# prepare to hold particle concentrations of components
		# in ug/m3
		yp = np.zeros((yp_molec.shape[0], yp_molec.shape[1]))

		# need to convert number concentration (molecules/cm3) 
		# of components to mass concentration (ug/m3)
		# first spread molar mass over particle size bins
		y_MM = np.tile(np.array((y_MM)).reshape(1, -1), (1, num_asb))

		# convert to mol/cm3
		yp[:, :] = yp_molec[:, :]/si.N_A

		# convert to g/cm3
		yp = yp*y_MM

		# convert to ug/m3
		yp = yp*1.e12

		# get indices of organics
		org_indx = []
		# get indices of inorganics
		inorg_indx = []
		# get indices of secondary organics
		sec_org_indx = []
		# get indices of VOCs with 10 or more carbons
		tenormoC_VOC_indx = []
		# get indices of black\elemental carbon
		bc_indx = []
		# get indices of primary organics
		pri_org_indx = []

		compi = 0 # count components
		for SMILESi in rel_SMILES:
			# if an organic component
			if (SMILESi.count('c') + SMILESi.count('C')>0
			and comp_names[compi] != 'bc' and comp_names[compi] != 'bcin'):
				# remember organic component index
				org_indx.append(compi)
				# check it's not one of the primary
				# organic components
				if comp_names[compi] not in pri_org_comps:
					# remember secondary organic
					# component index
					sec_org_indx.append(compi)
				# if it is a primary organic component
				if comp_names[compi] in pri_org_comps:
					# remember index
					pri_org_indx.append(compi)
				
				
			else: # if not an organic component
				# if not water
				if (comp_names[compi] != 'H2O'):
					# remember inorganic component index
					inorg_indx.append(compi)
				# black/elemental carbon
				if ('bc' in comp_names[compi]):
					bc_indx.append(compi)

			if (SMILESi.count('c') + SMILESi.count('C')>=10
			and SMILESi.count('o') + SMILESi.count('O')==0):
				tenormoC_VOC_indx.append(compi)	
			
			compi += 1 # count components
		
		# total influx mass concentration of ozone (ug/m3)
		# by end of simulation
		res_arr[-1, 0] = tot_influx[-1, tot_influx[0,:]==comp_names.index('O3')]

		# get indices of tot_influx where VOCs with 10 or more carbons are
		tot_influx_bigVOC_indx = []
		for bigVOCi in tenormoC_VOC_indx:
			tot_influx_bigVOC_indx.append(tot_influx[0,:].tolist().index(bigVOCi))
			
		# total influx mass concentration of VOCs with 10 or more carbons (ug/m3)
		# by end of simulation
		res_arr[-1, 1] = sum(tot_influx[-1, tot_influx_bigVOC_indx])

		# sum particle-phase components over particle size bins (ug/m3)
		yp_sb_tot = np.zeros((yp.shape[0], num_comp))
		for sbi in range(num_asb):
			yp_sb_tot[:, :] +=  yp[:, sbi*num_comp:(sbi+1)*num_comp]	


		##################################################################
		# in this section we need to subtract the secondary organics in
		# the particle phase of outdoor origin, using the simulation
		# results for the given air change rate and light intensity
		# but without activity

		# check if this is the simulation of a given air change rate
		# and light transmission factor without indoor activity
		if '0mops_0fry_0sho_0wst' in resi:

			# note that using 0:num_comp index below considers
			# components in just the smaller (PM2.5) particle size bin

			# if it is, then get the secondary organics of
			# outdoor origin (ug/m3)
			soopm = np.sum(yp_sb_tot[:, sec_org_indx], axis=1)
			soopm_pm2p5 = np.sum(yp[:, 0:num_comp][:, sec_org_indx], axis=1)

			# and get the primary organics of outdoor origin (ug/m3)
			poopm_pm2p5 = np.sum(yp[:, 0:num_comp][:, pri_org_indx], axis=1)

			# and get the inorganics of outdoor origin (ug/m3)
			ioopm_pm2p5 = (np.sum(yp[:, 0:num_comp][:, inorg_indx], axis=1)+
				np.sum(yp[:, 0:num_comp][:, bc_indx], axis=1))

		# sum the secondary organic particulate matter over
		# components and size bins (ug/m3)
		sopm_allsource = np.sum(yp_sb_tot[:, sec_org_indx], axis=1)
		sopm_indsource = np.sum(yp_sb_tot[:, sec_org_indx], axis=1)-soopm

		# and over just components in pm2.5
		sopm_pm2p5_allsource = np.sum(yp[:, 0:num_comp][:, sec_org_indx], axis=1)
		sopm_pm2p5_indsource = (np.sum(yp[:, 0:num_comp][:, sec_org_indx], axis=1)
			-soopm_pm2p5)

		# arithmetic mean of the secondary organic particulate matter
		# from indoors and outdoors
		# over the simulated time (ug/m3)
		res_arr[-1, 2] = np.sum(sopm_allsource, axis=0)/yp.shape[0]

		# arithmetic mean of the secondary organic particulate matter
		# from indoors only
		# over the simulated time (ug/m3)
		res_arr[-1, 3] = np.sum(sopm_indsource, axis=0)/yp.shape[0]

		# arithmetic mean of the secondary organic particulate matter
		# from indoors and outdoors
		# over the simulated time (ug/m3) in PM2.5
		res_arr_pm2p5[-1, 0] = np.sum(sopm_pm2p5_allsource, axis=0)/yp.shape[0]
		
		# arithmetic mean of the secondary organic particulate matter
		# from indoors only
		# over the simulated time (ug/m3) in PM2.5
		res_arr_pm2p5[-1, 1] = np.sum(sopm_pm2p5_indsource, axis=0)/yp.shape[0]

		##############################################################

		# indices of primary organic particulate matter of outside origin
		oopm_indx = [comp_names.index('pri_org')]

		# sum the organic particulate matter of outside origin 
		# (primary and secondary) over
		# components and size bins (ug/m3)
		oopm = np.sum(yp_sb_tot[:, oopm_indx], axis=1)+soopm
		# and over just components in pm2.5
		oopm_pm2p5 = np.sum(yp[:, 0:num_comp][:, oopm_indx], axis=1)+soopm_pm2p5

		# arithmetic mean of the organic particulate matter of
		# outside origin over the simulated time (ug/m3)
		res_arr[-1, 4] = np.sum(oopm, axis=0)/yp.shape[0]

		# arithmetic mean of the organic particulate matter of
		# outside origin over the simulated time (ug/m3) in PM2.5
		res_arr_pm2p5[-1, 2] = np.sum(oopm_pm2p5, axis=0)/yp.shape[0]


		# indices of indoor generated primary organic particulate matter	
		popm_indx = [comp_names.index('PNNCATCOOH'), 
			comp_names.index('C1213OH'), comp_names.index('NLIMOOH'), 
			comp_names.index('CO25C6CHO')]

		# sum the primary organic particulate matter of inside origin over
		# components and size bins (ug/m3)
		popm = np.sum(yp_sb_tot[:, popm_indx], axis=1)

		# and over just components in pm2.5
		popm_pm2p5 = np.sum(yp[:, 0:num_comp][:, popm_indx], axis=1)

		# arithmetic mean of the primary organic particulate matter of
		# inside origin over the simulated time (ug/m3)
		res_arr[-1, 5] = np.sum(popm, axis=0)/yp.shape[0]

		# arithmetic mean of the primary organic particulate matter of
		# over the simulated time (ug/m3) in PM2.5
		res_arr_pm2p5[-1, 3] = np.sum(popm_pm2p5, axis=0)/yp.shape[0]

		# get total mass concentration of organic particulate matter (ug/m3)
		# arithmetic mean over time (indoor secondary plus 
		# (outdoor primary and secondary) plus
		# indoor primary)
		res_arr[-1, 6] = (res_arr[-1, 3]+res_arr[-1, 4] +res_arr[-1, 5])

		# get total mass concentration of organic particulate matter (ug/m3)
		# (secondary from all sources and primary from all sources)
		res_arr_pm2p5[-1, 4] = (res_arr_pm2p5[-1, 0]+res_arr_pm2p5[-1, 3])	

		# sum the inorganic particulate matter of inside and 
		# outside origin over components and size bins (ug/m3)
		iopm = np.sum(yp_sb_tot[:, inorg_indx], axis=1)

		# and over just components in pm2.5
		iopm_pm2p5 = np.sum(yp[:, 0:num_comp][:, inorg_indx], axis=1)

		# arithmetic mean of the inorganic particulate matter of inside and 
		# outside origin over the simulated time (ug/m3)
		res_arr[-1, 7] = np.sum(iopm, axis=0)/yp.shape[0]

		# arithmetic mean of the inorganic particulate matter of inside and
		# outside origin over the simulated time (ug/m3) in PM2.5
		res_arr_pm2p5[-1, 5] = np.sum(iopm_pm2p5, axis=0)/yp.shape[0]

		# total particulate matter summed over organic and inorganic
		# components, but not water, and over size bins and arithmetic
		# mean average over time
		res_arr[-1, 8] = res_arr[-1, 6]+res_arr[-1, 7]
		
		res_num += 1 # count on simulations

		# fraction of organic mass comprised of 
		# secondary organic particulate matter 
		# from indoors and outdoors (0-1)
		res_arr_frac_org_pm[-1, 0] = res_arr[-1, 2]/res_arr[-1, 6]

		# fraction of organic mass comprised of 
		# organic particulate matter of outdoor origin
		res_arr_frac_org_pm[-1, 1] = res_arr[-1, 4]/res_arr[-1, 6]

		# fraction of organic mass comprised of 
		# primary organic particulate matter of indoor origin
		res_arr_frac_org_pm[-1, 2] = res_arr[-1, 5]/res_arr[-1, 6]

		# fraction of particulate matter components comprising
		# total mass concentration of PM2.5
		# first, get total (all components excluding water) 
		# PM2.5 mass concentrations
		res_arr_pm2p5[-1, 6] = (res_arr_pm2p5[-1, 4]+res_arr_pm2p5[-1, 5])
		# SOPM from indoor and outdoor sources
		res_arr_frac_pm2p5[-1, 0] = res_arr_pm2p5[-1, 0]/res_arr_pm2p5[-1, 6]
		# SOPM from indoor sources
		res_arr_frac_pm2p5[-1, 1] = res_arr_pm2p5[-1, 1]/res_arr_pm2p5[-1, 6]
		# OOOPM
		res_arr_frac_pm2p5[-1, 2] = res_arr_pm2p5[-1, 2]/res_arr_pm2p5[-1, 6]
		# POPM
		res_arr_frac_pm2p5[-1, 3] = res_arr_pm2p5[-1, 3]/res_arr_pm2p5[-1, 6]
		# OPM
		res_arr_frac_pm2p5[-1, 4] = res_arr_pm2p5[-1, 4]/res_arr_pm2p5[-1, 6]
		# IOPM
		res_arr_frac_pm2p5[-1, 5] = res_arr_pm2p5[-1, 5]/res_arr_pm2p5[-1, 6]

		print('resi: ', resi)
	
		# for the simulation with total PM2.5 over the WHO guideline and
		# a maximum relative concentration of indoor SOPM, first (here) get the
		# concentration of outdoor secondaries under high ventilation but
		# no indoor activity, use this to get the corresponding 
		# plot of contributors to indoor organics
		if 'sgr_3.2acr_0.0tf_0mops_0fry_0sho_0wst' in resi:

			# observed CIMS fractions for household M16 on Rikki's UKRI-MET 2024
			# presentation - low candle and cooking signal -> indicating
			# low indoor activity
			obs_cims_fracs = np.array((0.17, 0.48, 0.08, 0.02))
			# scale to equal 1
			obs_cims_fracs = obs_cims_fracs*(1./np.sum(obs_cims_fracs))

			# prepare to hold the oxidative potential of components
			# integrated over their PM2.5 mass concentration
			op_res = np.zeros((5, len(timehr)))

			# if it is, then get the oxidative potential of
			# secondary organics of
			# outdoor origin in PM2.5, e.g. Figure 6 of 
			# https://doi.org/10.5194/acp-22-1793-2022
			# but also: doi.org/10.5194/acp-18-9617-2018
			# and doi.org/10.1021/es505577w
			op_res[3, :] = np.sum(90.*yp[:, 0:num_comp][:, sec_org_indx], axis=1)

			# get the molecular concentration of primary organics of
			# outdoor origin (molecules/cm3)
			# primary outdoor organics (molecules/cm3)
			for sbi in range(num_asb): # loop through particle size bins
				cims_sources[1] += np.sum(np.sum(yp_molec[:, 
					sbi*num_comp:((sbi+1)*num_comp)][:, pri_org_indx]))

			# prepare figure for contribution to organics indoors based on
			# molecular concentration
			fig, (ax0) = plt.subplots(1, 1, figsize=(10, 5))

			# secondary organics (outdoors and indoors) (molecules/cm3)
			for sbi in range(num_asb): # loop through particle size bins
				cims_sources[0] += np.sum(np.sum(yp_molec[:, 
					sbi*num_comp:((sbi+1)*num_comp)][:, sec_org_indx]))

			# primary outdoor organics (molecules/cm3)
			for sbi in range(num_asb): # loop through particle size bins
				cims_sources[2] += np.sum(np.sum(yp_molec[:, 
					sbi*num_comp:((sbi+1)*num_comp)][:, pri_org_indx]))

			# now subtract the primary outdoor organics from the total
			# of primary organics to leave the primary indoor organics
			cims_sources[2] -= cims_sources[1]

			# convert molecular concentrations into fractions (0-1)
			cims_fracs = cims_sources[:]/np.sum(cims_sources)

			# prepare bar colours for simulated results
			bar_colours = ['red', (1., 0.61, 0), 'blue']
			# prepare bar colours for observed results
			bar_colours = ['red', (1., 0.61, 0), 'blue', (0.5, 0.5, 0.5)]

			for compi in range(len(cims_sources)):
				# stacked barplot of particulate matter 
				# molecular concentration fraction sources 
				# summed over 
				# components (per source) and summed over time
			
				# get vertical axis value to start this part of the
				# bar from
				bott_val = np.sum(cims_fracs[0:compi])
				# get height of this part of the bar
				height_now = cims_fracs[compi]

				ax0.bar(0., height_now, 
				label = cims_sources_labels[compi], 
				bottom = bott_val, color = bar_colours[compi], width = 0.25)

				# bar starting height from observations
				bott_val_obs = np.sum(obs_cims_fracs[0:compi])
				# height now from observations
				height_now_obs = obs_cims_fracs[compi]

				# observations for household with low cooking signal (M16)
				# on Rikki's source apportionment plot in UKRI-MET 2024
				# presentation
				ax0.bar(0.26, height_now_obs, 
				bottom = bott_val_obs, color = bar_colours[compi], 
				width = 0.25)

		
		# note that for overview paper (September 2024, Manchester email
		# correspondence with Nic Carslaw), we compare two household
		# simulations with identical indoor activity but different	
		# ventilation rates: 'sgr_0.1acr_1.0tf_1mops_1fry_1sho_0wst' (low
		# ventilation) and 'sgr_3.2acr_1.0tf_1mops_1fry_1sho_0wst' (high
		# ventilation)
		if 'sgr_0.1acr_1.0tf_1mops_1fry_1sho_0wst' in resi:
			
			# prepare to hold (for overview paper): indoor-generated inorganic, 
			# indoor-generated primary organics, indoor-generated secondary
			# organics, outdoor-generated inorganic, outdoor-generated primary
			# organics and outdoor-generated secondary organics (components in 
			# rows and ventilation rates (high and low in columns)), comparing 
			# households with same activity but two different ventilation rates

			# indoor-generated inorganic (ug/m3)
			overv_res[5, 0] = (np.mean((np.sum(yp[:, 0:num_comp][:, inorg_indx],
				 axis=1)+np.sum(yp[:, 0:num_comp][:, bc_indx], axis=1))-
				ioopm_pm2p5))
			overv_res_ts[5, 0, :] = ((np.sum(yp[:, 0:num_comp][:, inorg_indx],
				 axis=1)+np.sum(yp[:, 0:num_comp][:, bc_indx], axis=1))-
				ioopm_pm2p5)

			# indoor-generated primary organics (ug/m3)
			overv_res[4, 0] = (np.mean((np.sum(yp[:, 0:num_comp][:, pri_org_indx],
				 axis=1))-poopm_pm2p5))
			overv_res_ts[4, 0, :] = ((np.sum(yp[:, 0:num_comp][:, pri_org_indx],
				 axis=1))-poopm_pm2p5)

			# indoor-generated secondary organics (ug/m3)
			overv_res[3, 0] = (np.mean((np.sum(yp[:, 0:num_comp][:, sec_org_indx], 
				axis=1))-soopm_pm2p5))
			overv_res_ts[3, 0, :] = ((np.sum(yp[:, 0:num_comp][:, sec_org_indx], 
				axis=1))-soopm_pm2p5)

			# outdoor-generated inorganic (ug/m3)
			overv_res[2, 0] = np.mean(ioopm_pm2p5)
			overv_res_ts[2, 0, :] = (ioopm_pm2p5)

			# outdoor-generated primary organics (ug/m3)
			overv_res[1, 0] += np.mean(poopm_pm2p5)
			overv_res_ts[1, 0, :] = (poopm_pm2p5)

			# outdoor-generated secondary organics (ug/m3)
			overv_res[0, 0] = np.mean(soopm_pm2p5)
			overv_res_ts[0, 0, :] = (soopm_pm2p5)

		if 'sgr_3.2acr_1.0tf_1mops_1fry_1sho_0wst' in resi:
			
			# prepare to hold (for overview paper): indoor-generated inorganic, 
			# indoor-generated primary organics, indoor-generated secondary
			# organics, outdoor-generated inorganic, outdoor-generated primary
			# organics and outdoor-generated secondary organics (components in 
			# rows and ventilation rates (high and low in columns)), comparing 
			# households with same activity but two different ventilation rates
			
			# indoor-generated inorganic (ug/m3)
			overv_res[5, 1] = (np.mean((np.sum(yp[:, 0:num_comp][:, inorg_indx],
				 axis=1) + np.sum(yp[:, 0:num_comp][:, bc_indx], axis=1))-
				ioopm_pm2p5))
			overv_res_ts[5, 1, :] = ((np.sum(yp[:, 0:num_comp][:, inorg_indx],
				 axis=1)+np.sum(yp[:, 0:num_comp][:, bc_indx], axis=1))-
				ioopm_pm2p5)

			# indoor-generated primary organics (ug/m3)
			overv_res[4, 1] = (np.mean((np.sum(yp[:, 0:num_comp][:, pri_org_indx],
				 axis=1))-poopm_pm2p5))
			overv_res_ts[4, 1, :] = ((np.sum(yp[:, 0:num_comp][:, pri_org_indx],
				 axis=1))-poopm_pm2p5)

			# indoor-generated secondary organics (ug/m3)
			overv_res[3, 1] = (np.mean((np.sum(yp[:, 0:num_comp][:, sec_org_indx],
				 axis=1))-soopm_pm2p5))
			overv_res_ts[3, 1, :] = ((np.sum(yp[:, 0:num_comp][:, sec_org_indx],
				 axis=1))-soopm_pm2p5)

			# outdoor-generated inorganic (ug/m3)
			overv_res[2, 1] += np.mean(ioopm_pm2p5)
			overv_res_ts[2, 1, :] = (ioopm_pm2p5)

			# outdoor-generated primary organics (ug/m3)
			overv_res[1, 1] += np.mean(poopm_pm2p5)
			overv_res_ts[1, 1, :] = (poopm_pm2p5)

			# outdoor-generated secondary organics (ug/m3)
			overv_res[0, 1] = np.mean(soopm_pm2p5)
			overv_res_ts[0, 1, :] = (soopm_pm2p5)

			# get the mass fractions from each component
			overv_res_ts_frac = overv_res_ts[:, :, :]/np.sum(overv_res_ts, axis=0)

			# make the plot for overview paper (September 2024 and Manchester
			# email correspondence with Nic)
			# prepare figure for contribution to organics indoors based on
			# molecular concentration
			fig = plt.figure(figsize=(16, 6))
			gs = GridSpec(3, 3, figure=fig, width_ratios = [1, 3, 3], 
				height_ratios = [1, 1, 0.2])
			ax0 = fig.add_subplot(gs[0:2, 0])
			ax1 = fig.add_subplot(gs[0, 1])
			ax2 = fig.add_subplot(gs[1, 1])
			ax3 = fig.add_subplot(gs[0, 2])
			ax4 = fig.add_subplot(gs[1, 2])

			# parasite axis on time series plots
			#par1 = ax1.twinx() # parasite axis
			#par2 = ax2.twinx() # parasite axis
			
			for ci in range(overv_res.shape[0]): # loop through components

				height_now_lov = overv_res[ci, 0]
				bott_val_lov = np.sum(overv_res[0:ci, 0])

				height_now_hiv = overv_res[ci, 1]
				bott_val_hiv = np.sum(overv_res[0:ci, 1])
	
				# low ventilation result
				ax0.bar(0, height_now_lov, 
					bottom = bott_val_lov, label = overv_labels[ci],
					color = overv_colours[ci])
				
				# high ventilation result
				ax0.bar(1, height_now_hiv, 
					bottom = bott_val_hiv, color = overv_colours[ci])

			# time series for first scenario
			ax1.stackplot(timehr, overv_res_ts[:, 0, :], colors = overv_colours,
				labels=overv_labels)

			# time series of mass fractions for first scenario
			# loop through components
			for compi in range(overv_res_ts_frac.shape[0]):
				ax3.plot(timehr, overv_res_ts_frac[compi, 0, :], 
					color = overv_colours[compi],
					label = overv_labels[compi])

			# time series for second scenario
			ax2.stackplot(timehr, overv_res_ts[:, 1, :], colors = overv_colours,
				labels=overv_labels)

			# time series of mass fractions for first scenario
			# loop through components
			for compi in range(overv_res_ts_frac.shape[0]):
				ax4.plot(timehr, overv_res_ts_frac[compi, 1, :], 
					color = overv_colours[compi],
					label = overv_labels[compi])
			
			ax0.set_xticks([0., 1.], [str('ACR=0.1\n' + r'hour$\mathrm{^{-1}}$'), 
				str('ACR=3.2\n' + r'hour$\mathrm{^{-1}}$')])

			# set titles
			ax1.set_title(r'ACR=0.1 hour$\mathrm{^{-1}}$')
			ax2.set_title(r'ACR=3.2 hour$\mathrm{^{-1}}$')
			ax3.set_title(r'ACR=0.1 hour$\mathrm{^{-1}}$')
			ax4.set_title(r'ACR=3.2 hour$\mathrm{^{-1}}$')

			ax1.set_yscale('log') # set y-axis to log10
			ax2.set_yscale('log') # set y-axis to log10

			ax0.yaxis.set_tick_params(labelsize=14, direction = 'in', which='both')
			ax0.xaxis.set_tick_params(labelsize=14, direction = 'in', which='both')
		
			ax1.yaxis.set_tick_params(labelsize=14, direction = 'in', which='both')
			ax1.xaxis.set_tick_params(labelsize=14, direction = 'in', which='both')
			ax1.set_ylim(bottom=1.e-2, top = 7.e3)
			ax3.yaxis.set_tick_params(labelsize=14, direction = 'in', which='both')
			ax3.xaxis.set_tick_params(labelsize=14, direction = 'in', which='both')
			ax4.yaxis.set_tick_params(labelsize=14, direction = 'in', which='both')
			ax3.set_ylim(bottom=0.0, top = 1.3)

			# share the ax2 x axis with ax1
			ax1.sharex(ax2)
			ax3.sharex(ax4)

			ax2.yaxis.set_tick_params(labelsize=14, direction = 'in', which='both')
			ax2.xaxis.set_tick_params(labelsize=14, direction = 'in', which='both')
			ax2.set_ylim(bottom=2.e-1)
			ax4.yaxis.set_tick_params(labelsize=14, direction = 'in', which='both')
			ax4.xaxis.set_tick_params(labelsize=14, direction = 'in', which='both')
			#ax4.set_ylim(bottom=-0.4, top = 1.)

			locmin = ticker.FixedLocator((2.e-2, 3.e-2, 4.e-2, 5.e-2, 6.e-2, 7.e-2, 8.e-2, 9.e-2, 1.e-1, 2.e-1, 3.e-1, 4.e-1, 5.e-1, 6.e-1, 7.e-1, 8.e-1, 9.e-1, 2.e0, 3.e0, 4.e0, 5.e0, 6.e0, 7.e0, 8.e0, 9.e0, 1.e1, 2.e1, 3.e1, 4.e1, 5.e1, 6.e1, 7.e1, 8.e1, 9.e1, 2.e2))
			ax1.yaxis.set_minor_locator(locmin)
			#ax1.yaxis.set_minor_formatter(ticker.NullFormatter())

			ax0.set_ylabel(str(
				r'[PM$\mathrm{_{2.5}}$] $\mathrm{(\mu g \, m^{-3})}$'), 
				fontsize = 14)

			ax1.set_ylabel(str(
				r'[PM$\mathrm{_{2.5}}$] $\mathrm{(\mu g \, m^{-3})}$'), 
				fontsize = 14)
			
			ax2.set_ylabel(str(r'PM$\mathrm{_{2.5}}$ mass' +'\nfraction (0-1)'), 
				fontsize = 14)
			ax1.set_xlim(left=-1, right = 25)

			ax2.set_ylabel(str(
				r'[PM$\mathrm{_{2.5}}$] $\mathrm{(\mu g \, m^{-3})}$'), 
				fontsize = 14)
			ax3.set_ylabel(str(r'PM$\mathrm{_{2.5}}$ mass' + '\nfraction (0-1)'), 
				fontsize = 14)

			
			ax4.set_ylabel(str(r'PM$\mathrm{_{2.5}}$ mass' + '\nfraction (0-1)'), 
				fontsize = 14)
			ax2.set_xlim(left=-1, right = 25)

			ax2.set_xlabel(str('Time through day (hour)'), fontsize = 14)

			ax4.set_xlabel(str('Time through day (hour)'), fontsize = 14)

			ax4.set_xlim(left=-1, right = 25)

			ax0.legend(reverse=True, bbox_to_anchor=(7.5, -0.15), ncol=6)
			ax0.text(-0.5, 4., 'a)', fontsize = 14)
			ax1.text(-1.3, 1.4e4, 'b)', fontsize = 14)
			# activities start -------------------
			ax1.text(3., 1.5e3, 'mop', fontsize = 14)
			r0 = plt.Rectangle((6., 1.5e3), width=0.25, height=3.e3, 
				fill=None, alpha=1)
			ax1.add_patch(r0)
			ax1.text(4., 2.4e2, 'fry', fontsize = 14)
			r1 = plt.Rectangle((6., 3.e2), width=0.33, height=6.e2, 
				fill=None, alpha=1)
			ax1.add_patch(r1)
			ax1.text(7.5, 2.4e2, 'shower', fontsize = 14)
			r2 = plt.Rectangle((7., 3.e2), width=0.33, height=6.e2, 
				fill=None, alpha=1)
			ax1.add_patch(r2)

			ax3.text(3., 1.15, 'mop', fontsize = 14)
			r3 = plt.Rectangle((6., 1.15), width=0.25, height=0.1, 
				fill=None, alpha=1)
			ax3.add_patch(r3)
			ax3.text(4., 1.02, 'fry', fontsize = 14)
			r4 = plt.Rectangle((6., 1.02), width=0.33, height=0.1, 
				fill=None, alpha=1)
			ax3.add_patch(r4)
			ax3.text(7.5, 1.02, 'shower', fontsize = 14)
			r5 = plt.Rectangle((7., 1.02), width=0.33, height=0.1, 
				fill=None, alpha=1)
			ax3.add_patch(r5)
			# activities end -------------------

			ax2.text(-1.3, 2.8e2, 'c)', fontsize = 14)
			ax3.text(-1.4, 1.35, 'd)', fontsize = 14)
			ax4.text(-1.4, 1.05, 'e)', fontsize = 14)
			plt.setp(ax1.get_xticklabels(), visible=False)
			plt.setp(ax3.get_xticklabels(), visible=False)
			plt.tight_layout()
			plt.show()


	return()

# function to setup self
def self_def(dir_path_value):

	class testobj(object):
		pass

	self = testobj()
	self.dir_path = dir_path_value

	return(self)

abs_and_pos_plot(dir_path) # call function