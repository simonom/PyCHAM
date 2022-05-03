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
'''plots results for the particle-phase stuff'''
# simulation results are represented graphically:
# plots results for the total particle-phase concentration temporal profiles of specified components
# also the temporal profile or non-seed and non-water particle-phase

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import os
import retr_out
import numpy as np
import scipy.constants as si # for scientific constants

def plotter(caller, dir_path, comp_names_to_plot, self):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0 for individual components, 3 for 
	#	total excluding seed and water, 4 for top contributors to 
	#	particle-phase, 5 for particle surface area concentration, 
	#	6 for seed surface area concentration, 7 for top contributors
	#	to particle-phase excluding seed and water) 
	#	or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# comp_names_to_plot - chemical scheme names of components to plot
	# self - reference to GUI
	# --------------------------------------------------------------------------

	# chamber condition ---------------------------------------------------------
	# retrieve results
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, rel_SMILES, 
		y_MW, Nwet, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, H2Oi, seedi, _, _, group_indx, _, 
		ro_obj) = retr_out.retr_out(dir_path)
	
	# number of actual particle size bins
	num_asb = (num_sb-wall_on)

	if (caller == 0 or caller >= 3):
		plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize = (14, 7))

	if (comp_names_to_plot): # if component names specified
	
		# concentration plot ---------------------------------------------	
		for i in range(len(comp_names_to_plot)):
			
			if (comp_names_to_plot[i].strip() == 'H2O'):
				indx_plot = [H2Oi]
				indx_plot = np.array((indx_plot))
			if (comp_names_to_plot[i].strip() == 'RO2'):
				indx_plot = (np.array((group_indx['RO2i'])))
			if (comp_names_to_plot[i].strip() == 'RO'):
				indx_plot = (np.array((group_indx['ROi'])))
				
			if (comp_names_to_plot[i].strip() != 'H2O' and comp_names_to_plot[i].strip() != 'RO2' and comp_names_to_plot[i].strip() != 'RO'):
				try: # will work if provided components were in simulation chemical scheme
					# get index of this specified component, removing any white space
					indx_plot = [comp_names.index(comp_names_to_plot[i].strip())]
					indx_plot = np.array((indx_plot))
				except:
					self.l203a.setText(str('Component ' + comp_names_to_plot[i] + ' not found in chemical scheme used for this simulation'))
					# set border around error message
					if (self.bd_pl == 1):
						self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
						self.bd_pl = 2
					else:
						self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
						self.bd_pl = 1
					
					plt.ioff() # turn off interactive mode
					plt.close() # close figure window
					return()
			
			import scipy.constants as si # for scientific constants
			
			# particle-phase concentrations of all components (# molecules/cm3)
			if (wall_on == 1): # wall on
				ppc = yrec[:, num_comp:-num_comp]
			if (wall_on == 0): # wall off
				ppc = yrec[:, num_comp::]
			
			# particle-phase concentration of this component over all size bins (# molecules/cm3)
			conc = np.zeros((ppc.shape[0], num_sb-wall_on))
			for indxn in indx_plot: # loop through the indices
				concf = ppc[:, indxn::num_comp]
				# particle-phase concentration (ug/m3)
				conc[:, :] += ((concf/si.N_A)*y_MW[indxn])*1.e12
			
			if (comp_names_to_plot[i].strip() != 'RO2' and comp_names_to_plot[i].strip() != 'RO'): # if not a sum
				# plot this component
				ax0.plot(timehr, conc.sum(axis = 1), '+', linewidth = 4., label = str(str(comp_names[indx_plot[0]])+' (particle-phase)'))
			
			if (comp_names_to_plot[i].strip() == 'RO2'): # if is the sum of organic peroxy radicals
				ax0.plot(timehr, conc.sum(axis = 1), '-+', linewidth = 4., label = str(r'$\Sigma$RO2 (particle-phase)'))
			
			if (comp_names_to_plot[i].strip() == 'RO'): # if is the sum of organic alkoxy radicals
				ax0.plot(timehr, conc.sum(axis = 1), '-+', linewidth = 4., label = str(r'$\Sigma$RO (particle-phase)'))
				
		ax0.set_ylabel(r'Concentration ($\rm{\mu}$g$\,$m$\rm{^{-3}}$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)

		# end of gas-phase concentration sub-plot ---------------------------------------
	
	# if called by button to plot temporal profile of total particle-phase concentration 
	# excluding water and seed
	if (caller == 3):
		import scipy.constants as si
		# particle-phase concentrations of all components (# molecules/cm3)
		if (wall_on == 1): # wall on
			ppc = yrec[:, num_comp:-num_comp]
		if (wall_on == 0): # wall off
			ppc = yrec[:, num_comp::]

		# zero water and seed
		ppc[:, H2Oi::num_comp] = 0.
		for i in seedi: # loop through seed components
			ppc[:, i::num_comp] = 0.
		# tile molar weights over size bins and times
		y_mwt = np.tile(np.array((y_MW)).reshape(1, -1), (1, num_sb-wall_on))
		y_mwt = np.tile(y_mwt, (ppc.shape[0], 1))
		# convert from # molecules/cm3 to ug/m3
		ppc = (ppc/si.N_A)*y_mwt*1.e12
		# sum over components and size bins
		ppc = np.sum(ppc, axis=1)
		
		# plot
		ax0.plot(timehr, ppc, '+', linewidth = 4., label = 'total particle-phase excluding seed and water')
		ax0.set_ylabel(r'Concentration ($\rm{\mu}$g$\,$m$\rm{^{-3}}$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)

	# if called by button to plot top contributors to particle-phase
	if (caller == 4):
		import scipy.constants as si
		# particle-phase concentrations of all components (# molecules/cm3)
		if (wall_on == 1): # wall on
			ppc = yrec[:, num_comp:-num_comp]
		if (wall_on == 0): # wall off
			ppc = yrec[:, num_comp::]

		# tile molar weights over size bins and times
		y_mwt = np.tile(np.array((y_MW)).reshape(1, -1), (1, (num_sb-wall_on)))
		y_mwt = np.tile(y_mwt, (ppc.shape[0], 1))

		# convert to ug/m3 from # molecules/cm3
		ppc = ((ppc/si.N_A)*y_mwt)*1.e12

		# sum particle-phase concentration per component over time (ug/m3)
		ppc_t = (np.sum(ppc, axis=0)).reshape(1, -1)
		
		# sum particle-phase concentrations over size bins (ug/m3), 
		# but keeping components separate
		ppc_t = np.sum(ppc_t.reshape(num_sb-wall_on, num_comp), axis=0)
		
		# convert to mass contributions
		ppc_t = ((ppc_t/np.sum(ppc_t))*100.).reshape(-1, 1)
		
		# rank in ascending order
		ppc_ts = np.flip(np.sort(np.squeeze(ppc_t)))
		
		# take just top number as supplied by user
		ppc_ts = ppc_ts[0:self.e300r]
		
		# sum concentrations over size bins and components
		ppc_sbc = np.sum(ppc, axis=1)

		# loop through to plot
		for i in range(self.e300r):
			
			# get index
			indx = np.where(ppc_t == ppc_ts[i])[0][0]
			# get name
			namei = comp_names[indx]
			# sum for this component over size bins (ug/m3)
			ppci = np.sum(ppc[:, indx::num_comp], axis=1)
			# get contribution (%)
			ppci = (ppci[ppc_sbc>0.]/ppc_sbc[ppc_sbc>0.])*100.

			ax0.plot(timehr[ppc_sbc>0.], ppci, '-+', linewidth = 4., label = namei)
	
		ax0.set_ylabel(r'Contribution to particle-phase mass concentration (%)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)
	
	# if called by button to plot top contributors to particle-phase excluding seed and water
	if (caller == 7):

		import scipy.constants as si		

		# particle-phase concentrations of all components (# molecules/cm3)
		if (wall_on == 1): # wall on
			ppc = yrec[:, num_comp:-num_comp]
		if (wall_on == 0): # wall off
			ppc = yrec[:, num_comp::]

		for i in seedi: # loop through seed components
			ppc[:, i::num_comp] = 0. # zero seed

		ppc[:, H2Oi::num_comp] = 0. # zero water

		# tile molar weights over size bins and times
		y_mwt = np.tile(np.array((y_MW)).reshape(1, -1), (1, (num_sb-wall_on)))
		y_mwt = np.tile(y_mwt, (ppc.shape[0], 1))

		# convert to ug/m3 from # molecules/cm3
		ppc = ((ppc/si.N_A)*y_mwt)*1.e12

		# sum particle-phase concentration per component over time (ug/m3)
		ppc_t = (np.sum(ppc, axis=0)).reshape(1, -1)
		
		# sum particle-phase concentrations over size bins (ug/m3), 
		# but keeping components separate
		ppc_t = np.sum(ppc_t.reshape(num_sb-wall_on, num_comp), axis=0)
		
		# convert to mass contributions
		ppc_t = ((ppc_t/np.sum(ppc_t))*100.).reshape(-1, 1)
		
		# rank in ascending order
		ppc_ts = np.flip(np.sort(np.squeeze(ppc_t)))
		
		# take just top number as supplied by user
		ppc_ts = ppc_ts[0:self.e300r]
		
		# sum concentrations over size bins and components
		ppc_sbc = np.sum(ppc, axis=1)

		# loop through to plot
		for i in range(self.e300r):
			
			# get index
			indx = np.where(ppc_t == ppc_ts[i])[0][0]
			# get name
			namei = comp_names[indx]
			# sum for this component over size bins (ug/m3)
			ppci = np.sum(ppc[:, indx::num_comp], axis=1)
			# get contribution (%)
			ppci = (ppci[ppc_sbc>0.]/ppc_sbc[ppc_sbc>0.])*100.

			ax0.plot(timehr[ppc_sbc>0.], ppci, '-+', linewidth = 4., label = namei)
	
		ax0.set_ylabel(r'Contribution to particle-phase mass concentration (%)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)
	
	# if called by button to plot total particle surface area concentration (m2/m3)
	if (caller == 5):
		# get surface area of single particles per size bin 
		# at all times (m2), note radius converted from um 
		# to m
		asp = 4.*np.pi*(x*1.e-6)**2.
		# integrate over all particles in a size bin, 
		# note conversion from /cm3 to /m3 (m2/m3)
		asp = asp*(Nwet*1.e6)	
		# sum over size bins (m2/m3)
		asp = np.sum(asp, axis=1)
		# plot
		ax0.plot(timehr, asp, '-+', linewidth = 4.)
		ax0.set_ylabel(r'Total particle-phase surface area concentration ($\rm{m^{2}\,m^{-3}}$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')

	# if called by button to plot seed particle surface area concentration (m2/m3)
	if (caller == 6):

		import scipy.constants as si

		# particle-phase concentrations of all components (# molecules/cm3)
		if (wall_on == 1): # wall on
			ppc = yrec[:, num_comp:-num_comp]
		if (wall_on == 0): # wall off
			ppc = yrec[:, num_comp::]
		
		# empty results array for seed particle volume concentrations per size bin
		ppcs = np.zeros((ppc.shape[0], (num_sb-wall_on)))

		for i in seedi: # loop through seed indices
			# get just seed component particle-phase volume concentration 
			# (cm3/cm3)
			ppcs[:, :] += (ppc[:, i::num_comp]/si.N_A)*(np.array((y_MV))[i])
		
		# convert total volume to volume per particle (cm3)
		ppcs[Nwet>0] = ppcs[Nwet>0]/Nwet[Nwet>0]
		
		# convert volume to radius (cm)
		ppcs = ((3.*ppcs)/(4.*np.pi))**(1./3.)
		
		# convert cm to m
		ppcs = ppcs*1.e-2
		# get surface area over all particles in a size bin (m2/m3), 
		# note conversion from m2/cm3 to m2/m3
		ppcs = (4.*np.pi*ppcs**2.)*Nwet*1.e6
		# sum over all size bins (m2/m3)
		ppcs = np.sum(ppcs, axis=1)
		
		# plot
		ax0.plot(timehr, ppcs, '-+', linewidth = 4.)
		ax0.set_ylabel(r'Seed particle-phase surface area concentration ($\rm{m^{2}\,m^{-3}}$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')

	if (caller == 8): # if called by button to plot generational contribution to particle-phase mass

		import sch_interr # for interpeting chemical scheme
		import re # for parsing chemical scheme
		import scipy.constants as si

		gen_num = [] # empty results list
		ci = 0 # component count
		sch_name = ro_obj.sp
		inname = ro_obj.vp

		sch_name = str(dir_path + '/inputs/api_iso_ch4_mechHOM_scheme.txt')
		inname = str(dir_path + '/inputs/api_iso_ch4_mechHOM_model_var.txt')

		f_open_eqn = open(sch_name, mode='r') # open the chemical scheme file
		# read the file and store everything into a list
		total_list_eqn = f_open_eqn.readlines()
		f_open_eqn.close() # close file

		inputs = open(inname, mode= 'r' ) # open model variables file
		in_list = inputs.readlines() # read file and store everything into a list
		inputs.close() # close file
		for i in range(len(in_list)): # loop through supplied model variables to interpret

			# ----------------------------------------------------
			# if commented out continue to next line
			if (in_list[i][0] == '#'):
				continue
			key, value = in_list[i].split('=') # split values from keys
			# model variable name - a string with bounding white space removed
			key = key.strip()
			# ----------------------------------------------------

			if key == 'chem_scheme_markers' and (value.strip()): # formatting for chemical scheme
				chem_sch_mrk = [str(i).strip() for i in (value.split(','))]

		# interrogate scheme to list equations
		[eqn_list, aqeqn_list, eqn_num, rrc, rrc_name, 
			RO2_names] = sch_interr.sch_interr(total_list_eqn, chem_sch_mrk)	

		# loop through components to identify their generation
		for compi in comp_names[0:-2]: # don't include core and water
			
			comp_fin = 0 # flag for when component generation number found
			gen_num.append(0) # assume zero generation by default
			#print(ci, compi, rel_SMILES[ci], len(comp_names), len(rel_SMILES))
			# if an organic molecule
			if ('C' in rel_SMILES[ci] or 'c' in rel_SMILES[ci]):
				
				# loop through reactions to identify where this first occurs,
				# assuming that if first occurrence is as a reactant 
				# it is zero generation and if as a product is >zero generation
				for eqn_step in range(len(eqn_list)):

					line = eqn_list[eqn_step] # extract this line
					# work out whether equation or reaction rate coefficient part comes first
					eqn_start = str('.*\\' +  chem_sch_mrk[10])
					rrc_start = str('.*\\' +  chem_sch_mrk[9])
					# get index of these markers, note span is the property of the match object that
					# gives the location of the marker
					eqn_start_indx = (re.match(eqn_start, line)).span()[1]
					rrc_start_indx = (re.match(rrc_start, line)).span()[1]
		
					if (eqn_start_indx>rrc_start_indx):
						eqn_sec = 1 # equation is second part
					else:
						eqn_sec = 0 # equation is first part
		
					# split the line into 2 parts: equation and rate coefficient
					# . means match with anything except a new line character., when followed by a * 
					# means match zero or more times (so now we match with all characters in the line
					# except for new line characters, so final part is stating the character(s) we 
					# are specifically looking for, \\ ensures the marker is recognised
					if eqn_sec == 1:
						eqn_markers = str('\\' +  chem_sch_mrk[10]+ '.*\\' +  chem_sch_mrk[11])
					else: # end of equation part is start of reaction rate coefficient part
						eqn_markers = str('\\' +  chem_sch_mrk[10]+ '.*\\' +  chem_sch_mrk[9])

					# extract the equation as a string ([0] extracts the equation section and 
					# [1:-1] removes the bounding markers)
					eqn = re.findall(eqn_markers, line)[0][1:-1].strip()
		
					eqn_split = eqn.split()
					eqmark_pos = eqn_split.index('=')
					# reactants with stoichiometry number and omit any photon
					reactants = [i for i in eqn_split[:eqmark_pos] if i != '+' and i != 'hv']
					# products with stoichiometry number
					products = [t for t in eqn_split[eqmark_pos+1:] if t != '+']
					
					for ri in reactants:
						# note that no spaces or other punctuation included around 
						# the component name as extracted from the equation
						if compi in ri and len(compi) == len(ri):
							# assuming that components appearing first as a reactant 
							# must be zero generation
							gen_num[ci] = 0
							# finished with this component, break out of reactant loop
							comp_fin = 1
							break

					for pi in products:
						# note that no spaces or other punctuation included around 
						# the component name as extracted from the equation
						if compi in pi and len(compi) == len(pi):
							# check which reactant has earliest generation
							rcheck = []
							for ri in reactants:
								# get its index
								rindx = comp_names.index(ri)
								# if an organic, note if inorgnic then ignore
								if ('C' in rel_SMILES[rindx] or 'c' in rel_SMILES[rindx]):
									if (rindx < ci): # if generation number of recatant already known
										rcheck.append(gen_num[rindx])
							# identify earliest generation reactant
							egen = np.min(rcheck)
							# check whether this species has a charge, i.e. is open shell, note that according 
							# to DAYLIGHT, the square brackets in a SMILES string deontes when an atom has a 
							# valence other than normal, i.e. is in an excited (radical) state: 
							# https://daylight.com/dayhtml/doc/theory/theory.smiles.html
							if ('[' in rel_SMILES[ci] or ']' in rel_SMILES[ci]): # if open-shell
								gen_num[ci] = egen
							else: # if closed-shell
								gen_num[ci] = egen+1
							# finished with this component, break out of reactant loop
							comp_fin = 1
							break

					# finished with this component, break out of equation loop, 
					# move onto next component
					if (comp_fin == 1):
						break
				

			ci += 1 # component count
		# convert list to numpy array
		gen_num = np.array(gen_num)
		# append two zeros for core and water
		gen_num = np.concatenate((gen_num, np.zeros(2)), axis=0)
		# in case we want to see generation number per component -------------
		#ax0.plot(np.arange(len(comp_names)), gen_num, '+')
		#ax0.set_xticks(np.arange(len(comp_names)))
		#ax0.set_xticklabels(comp_names, rotation = 90)
		# --------------------------------------------------------------------
		
		# empty results matrix for contribution (%) to particle-phase mass per 
		# generation per time step
		res = np.zeros((len(timehr), int(np.max(gen_num))))
		# particle-phase concentrations of all components (# molecules/cm3)
		if (wall_on == 1): # wall on
			ppc = yrec[:, num_comp:-num_comp]
		if (wall_on == 0): # wall off
			ppc = yrec[:, num_comp::]

		y_MW = np.tile(y_MW, (len(timehr), num_sb-wall_on))

		# convert particle-phase concentrations into mass 
		# concentration (ug/m3) from # molecules/cm3		
		ppc[:, :] += ((ppc/si.N_A)*y_MW)*1.e12
		# zero water and seed contribution
		ppc[:, int(H2Oi)::int(num_comp)] = 0.
		for seei in seedi:
			ppc[:, int(seei)::int(num_comp)] = 0.
		# sum concentrations over all size bins 
		# and components per time step (ug/m3)
		ppc_sum = np.sum(ppc, axis= 1)

		# tile generation number over size bins
		gen_num = np.tile(gen_num, (num_sb-wall_on))
		# loop through generations
		for gi in range(int(np.max(gen_num))):
			# get percentage contribution from this generation
			res[:, gi] = (np.sum(ppc[:, gen_num == gi], axis = 1)/ppc_sum)*100.
			ax0.plot(timehr, res[:, gi], label = str('Generation ' + str(gi)))
		ax0.set_ylabel(r'% Contribution to organic particle-phase mass concentration ($\%$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)

		return()

	# display
	if (caller == 2):
		plt.show()	

	return()