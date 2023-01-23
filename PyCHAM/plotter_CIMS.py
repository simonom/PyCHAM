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
'''plots a replication of mass spectrum as reported by a chemical ionisation mass spectrometer (CIMS)'''
# simulation results are represented graphically

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import retr_out
import numpy as np
import scipy.constants as si
import importlib

def plotter_CIMS(self, res_in, tn, iont, sens_func):
	
	# inputs: -----------------
	# self - reference to PyCHAM class
	# res_in - inputs for resolution of molar mass to charge ratio (g/mol/charge)
	# tn - time through experiment to plot at (s)
	# iont - type of ionisation
	# sens_func - sensitivity to molar mass function
	# ---------------------------
	
	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = np.zeros((self.ro_obj.yrec.shape[0], self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], self.ro_obj.Nrec_wet.shape[1]))
	Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
	Ndry = np.zeros((self.ro_obj.Nrec_dry.shape[0], self.ro_obj.Nrec_dry.shape[1]))
	Ndry[:, :] = self.ro_obj.Nrec_dry[:, :]
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec= np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	space_mode = self.ro_obj.spacing
	Cfac = self.ro_obj.cfac
	group_indx = self.ro_obj.gi
	y_MV = (np.array((self.ro_obj.comp_MV))).reshape(1, -1)
	yrec_p2w = self.ro_obj.part_to_wall
	PsatPa = self.ro_obj.vpPa
	O_to_C = self.ro_obj.O_to_C
	H2Oi = self.ro_obj.H2O_ind
	
	# convert to 2D numpy array
	y_MW = np.array((y_MW)).reshape(-1, 1)
	
	# get index of time wanted
	ti = (np.where(np.abs(timehr-tn/3600.) == np.min(np.abs(timehr-tn/3600.))))[0][0]
	
	# convert yrec from 1D to 2D with times in rows, then select time wanted
	yrec = (yrec.reshape(len(timehr), num_comp*(num_sb+1)))[ti, :]
	
	# get gas-phase concentrations (ppt, note starting concentration is ppb)
	gp = yrec[0:num_comp]*1.e3  # conversion to ug/m3 (if wanted): /Cfac[ti]/si.N_A*y_MW[:, 0]*1.e12
	
	# get particle-phase concentrations (molecules/cm3)
	pp = yrec[num_comp:num_comp*(num_sb+1-wall_on)]
	
	# sum each component over size bins (molecules/cm3)
	pp = np.sum(pp.reshape(num_sb-wall_on, num_comp), axis=0)
	
	# convert to ppt
	pp = (pp/si.N_A)*Cfac[ti]*1.e6 # or convert to ug/m3: *y_MW[:, 0]*1.e12
	
	# correct for sensitivity to molar mass
	fac_per_comp = write_sens2mm(0, sens_func, y_MW)
	
	gp = gp*fac_per_comp[:]
	pp = pp*fac_per_comp[:]

	# if ionisation source molar mass to be added 
	# (e.g. because not corrected for in measurment software), then add
	if (int(iont[1]) == 1):
		if (iont[0] == 'I'): # (https://pubchem.ncbi.nlm.nih.gov/compound/Iodide-ion)
			y_MW += 126.9045
		if (iont[0] == 'N'): # (https://pubchem.ncbi.nlm.nih.gov/compound/nitrate)
			y_MW += 62.005
	
	# remove water
	gp = np.append(gp[0:H2Oi], gp[H2Oi+1::])
	pp = np.append(pp[0:H2Oi], pp[H2Oi+1::])
	y_MW = np.append(y_MW[0:H2Oi, 0], y_MW[H2Oi+1::, 0])

	# account for mass to charge resolution
	[pdf, comp_indx, comp_prob, mm_all] = write_mzres(1, res_in, y_MW)
	gpres = np.zeros((len(comp_indx)))
	ppres = np.zeros((len(comp_indx)))
	
	for pdfi in range(len(comp_indx)): # loop through resolution intervals
		gpres[pdfi] = np.sum(gp[comp_indx[pdfi]]*comp_prob[pdfi])
		ppres[pdfi] = np.sum(pp[comp_indx[pdfi]]*comp_prob[pdfi])
		
	plt.ion() # disply plot in interactive mode

	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
	ax0.semilogy(mm_all, gpres, '+m', markersize = 14, markeredgewidth = 5,  label = str('gas-phase'))
	ax0.semilogy(mm_all, ppres, 'xb', markersize = 14, markeredgewidth = 5, label = str('particle-phase'))
	
	ax0.set_title(str('Mass spectrum at ' + str(timehr[ti]) + ' hours'), fontsize = 14)
	ax0.set_xlabel(r'Mass/charge (Th)', fontsize = 14)
	ax0.set_ylabel(r'Concentration (ppt)', fontsize = 14)
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	ax0.legend(fontsize = 14)
	
	return()

# function for transforming and saving output in CIMS format
def write_CIMS_output(self):

	# import dependencies
	import scipy.stats as st
	import numpy as np
	import os

	# ------------------------------------------------------------
	# inputs: self - reference to PyCHAM class
	# ------------------------------------------------------------
	
	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = np.zeros((self.ro_obj.yrec.shape[0], self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], self.ro_obj.Nrec_wet.shape[1]))
	Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
	Ndry = np.zeros((self.ro_obj.Nrec_dry.shape[0], self.ro_obj.Nrec_dry.shape[1]))
	Ndry[:, :] = self.ro_obj.Nrec_dry[:, :]
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec= np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	space_mode = self.ro_obj.spacing
	Cfac = self.ro_obj.cfac
	group_indx = self.ro_obj.gi
	y_MV = (np.array((self.ro_obj.comp_MV))).reshape(1, -1)
	yrec_p2w = self.ro_obj.part_to_wall
	PsatPa = self.ro_obj.vpPa
	O_to_C = self.ro_obj.O_to_C
	H2Oi = self.ro_obj.H2O_ind

	# ensure numpy array
	y_MW = np.array((y_MW)).reshape(1, -1)

	with open(os.path.join(self.dir_path, 'MCMnames.txt'), "w") as output:
		for row in comp_names[0:-2]:
			
			output.write(str(row) + '\n')

	# preparing code as requested by Lukas in meeting on 19/01/2023
	# txt file of C H O N (4 integers) (columns) per species (rows)
	CHON_res = np.zeros((len(comp_names), 4))

	# txt file of concentration (columns) and rate of change of 
	# concentration (columns) for each component (rows)
	conc_res = np.zeros((len(comp_names), 2))

	# text file of NO, HO2, RO2 pool concentrations at steady 
	# state (# molecules/cm3)
	concs_res = np.zeros((3, 1))

	import pybel # required dependency

	# indices of RO2 components
	indx_plot = (np.array((group_indx['RO2i'])))

	# loop through components to calculate their C H O N and 
	# record,  but exclude water and core
	for icomp in range(len(comp_names)-2):
		# get the pybel object for this component
		Pybel_object = pybel.readstring('smi', rel_SMILES[icomp])	
		
		# get number of each atom
		if 'C' in Pybel_object.formula:
			try:
				numiC = Pybel_object.formula.index('C')
				try:
					numC = float(Pybel_object.formula[numiC+1])
					try:
						numC = float(Pybel_object.formula[numiC+1:numiC+3])
					except:
						numC = float(Pybel_object.formula[numiC+1])
				except:
					numC = 1
			except:
				numiC = 'str'

			if numiC != 'str':
				CHON_res[icomp, 0] += numC

		# get number of each atom
		if 'c' in Pybel_object.formula:
			try:
				numiC = Pybel_object.formula.index('c')
				try:
					numC = float(Pybel_object.formula[numiC+1])
					try:
						numC = float(Pybel_object.formula[numiC+1:numiC+3])
					except:
						numC = float(Pybel_object.formula[numiC+1])
				except:
					numC = 1
			except:
				numiC = 'str'

			if numiC != 'str':
				CHON_res[icomp, 0] += numC

		# get number of each atom
		if 'H' in Pybel_object.formula:
			try:
				numiC = Pybel_object.formula.index('H')
				try:
					numC = float(Pybel_object.formula[numiC+1])
					try:
						numC = float(Pybel_object.formula[numiC+1:numiC+3])
					except:
						numC = float(Pybel_object.formula[numiC+1])
				except:
					numC = 1
			except:
				numiC = 'str'

			if numiC != 'str':
				CHON_res[icomp, 1] += numC

		# get number of each atom
		if 'h' in Pybel_object.formula:
			try:
				numiC = Pybel_object.formula.index('h')
				try:
					numC = float(Pybel_object.formula[numiC+1])
					try:
						numC = float(Pybel_object.formula[numiC+1:numiC+3])
					except:
						numC = float(Pybel_object.formula[numiC+1])
				except:
					numC = 1
			except:
				numiC = 'str'

			if numiC != 'str':
				CHON_res[icomp, 1] += numC

		# get number of each atom
		if 'O' in Pybel_object.formula:
			try:
				numiC = Pybel_object.formula.index('O')
				try:
					numC = float(Pybel_object.formula[numiC+1])
					try:
						numC = float(Pybel_object.formula[numiC+1:numiC+3])
					except:
						numC = float(Pybel_object.formula[numiC+1])
				except:
					numC = 1
			except:
				numiC = 'str'

			if numiC != 'str':
				CHON_res[icomp, 2] += numC

		# get number of each atom
		if 'o' in Pybel_object.formula:
			try:
				numiC = Pybel_object.formula.index('o')
				try:
					numC = float(Pybel_object.formula[numiC+1])
					try:
						numC = float(Pybel_object.formula[numiC+1:numiC+3])
					except:
						numC = float(Pybel_object.formula[numiC+1])
				except:
					numC = 1
			except:
				numiC = 'str'

			if numiC != 'str':
				CHON_res[icomp, 2] += numC

		# get number of each atom
		if 'N' in Pybel_object.formula:
			try:
				numiC = Pybel_object.formula.index('N')
				try:
					numC = float(Pybel_object.formula[numiC+1])
					try:
						numC = float(Pybel_object.formula[numiC+1:numiC+3])
					except:
						numC = float(Pybel_object.formula[numiC+1])
				except:
					numC = 1
			except:
				numiC = 'str'

			if numiC != 'str':
				CHON_res[icomp, 3] += numC

		# get number of each atom
		if 'N' in Pybel_object.formula:
			try:
				numiC = Pybel_object.formula.index('n')
				try:
					numC = float(Pybel_object.formula[numiC+1])
					try:
						numC = float(Pybel_object.formula[numiC+1:numiC+3])
					except:
						numC = float(Pybel_object.formula[numiC+1])
				except:
					numC = 1
			except:
				numiC = 'str'

			if numiC != 'str':
				CHON_res[icomp, 3] += numC
				
	
		# concentration of this component at penultimate time step (# molecules/cm3)
		conc_res[icomp, 0] = yrec[-1, icomp]*Cfac[-1]
		# rate of change of this component at penultimate time step (# molecules/cm3/s)
		conc_res[icomp, 1] = (yrec[-1, icomp]-yrec[-2, icomp])*Cfac[-1]
		
		if comp_names[icomp] == 'NO':
			concs_res[0] = yrec[-1, icomp]*Cfac[-1]
		if comp_names[icomp] == 'HO2':
			concs_res[1] = yrec[-1, icomp]*Cfac[-1]
		if icomp in indx_plot: # RO2 pool
			concs_res[2] += yrec[-1, icomp]*Cfac[-1]
		
	
	# save
	np.savetxt(os.path.join(self.dir_path, 'CHON.txt'), CHON_res, delimiter=' ') 		
	np.savetxt(os.path.join(self.dir_path, 'Cin1stcol_dCDtin2ndcol.txt'), conc_res, delimiter=' ')
	np.savetxt(os.path.join(self.dir_path, 'CNO1strow_CHO22ndrow_CRO23rdrow.txt'), concs_res, delimiter=' ') 

	if self.iont[1] == 1: # if we need to add ioniser molar mass onto molar mass
		if self.iont[0] == 'I': # if iodide
			y_MW += 127.
		if self.iont[0] == 'N': # if nitrate
			y_MW += 62.
		if self.iont[0] == 'B': # if bromide
			y_MW += 119.

	# prepare results matrix (for CIMS output)
	CIMS_res = np.zeros((len(timehr), int(max(y_MW[0, :])/self.resol_in[0])))
	
	# m:z bin centres for instrument
	mzbins = np.arange(0, int(max(y_MW[0, :])), self.resol_in[0])

	# prepare results matrix (for CIMS output)
	CIMS_res = np.zeros((len(timehr)+1, len(mzbins)))

	# record m:z values in first row
	CIMS_res[0, :] = mzbins[:]

	# empty array to hold the fractional contribution to each m:z of each 
	# component, omitting water and core
	fc = np.zeros((num_comp-2, len(mzbins)))

	# loop through components to get m:z indices and fractional
	# contributions, omitting water and core
	for imm in range(len(y_MW[0, 0:-2])):
		
		# get the probability distribution function range for 
		# this molar mass centre (0-1 over entire molar mass range)
		pdf = st.norm.pdf(mzbins, y_MW[0, imm], self.resol_in[1])

		# ensure fractional contributions equal 1
		fc[imm, :] = pdf/sum(pdf)

	# loop through times
	for it in range(1, len(timehr)+1):
		
		# multiply concentrations of all gas-phase components (ppb) now by
		# their fractional contribution to each m:z bin (components vary by
		# row)
		fcn = (yrec[it-1, 0:(num_comp-2)].reshape(-1, 1))*fc
		
		# sum over components (ppb)
		fcn = (np.sum(fcn, axis=0)).reshape(1, -1)

		# remember in results (ppb)
		CIMS_res[it, :] = fcn[0, :] 

	# save CIMS-style output to csv file
	# saving both gas, particle and wall concentrations of components
	np.savetxt(os.path.join(self.dir_path, 'simulat_res_conv_to_CIMS_output'), CIMS_res, delimiter=',', header=str('time changes with rows which correspond to the time output file, m:z ratios is columns, first row is m:z values')) 		
	 
	# uncomment below to run check ----------------------------------------------------
	# check correct results saved (by comparing with output from the 'CIMS observations' button 
	#check_plot = np.loadtxt(os.path.join(self.dir_path, 'simulat_res_conv_to_CIMS_output'), delimiter=',', skiprows=1, dtype='str')
	#import matplotlib.pyplot as plt
	#plt.ion() # disply plot in interactive mode
	#fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7)) # prepare plot
	#y_MW = check_plot[0, :].astype('float')
	# chose the y-axis, note that choice of row index determines time through simulation
	#peaks = check_plot[-1, :].astype('float')
	#ax0.plot(y_MW, peaks)
	# uncomment above to run check ----------------------------------------------------


	return()

# function for plotting sensitivity to components
def write_sens2mm(caller, sens_func, y_MW):

	import datetime

	# inputs: --------------------------
	# caller - flag for the calling function
	# sens_func - the sensitivity (Hz/ppt) function to molar mass (g/mol)
	# y_MW - molar mass of components (g/mol)
	# ------------------------------------

	# create new  file
	f = open('PyCHAM/sens2mm.py', mode='w')
	f.write('##########################################################################################\n')
	f.write('#                                                                                        											 #\n')
	f.write('#    Copyright (C) 2018-2022 Simon O\'Meara : simon.omeara@manchester.ac.uk                  				 #\n')
	f.write('#                                                                                       											 #\n')
	f.write('#    All Rights Reserved.                                                                									 #\n')
	f.write('#    This file is part of PyCHAM                                                         									 #\n')
	f.write('#                                                                                        											 #\n')
	f.write('#    PyCHAM is free software: you can redistribute it and/or modify it under              						 #\n')
	f.write('#    the terms of the GNU General Public License as published by the Free Software       					 #\n')
	f.write('#    Foundation, either version 3 of the License, or (at your option) any later          						 #\n')
	f.write('#    version.                                                                            										 #\n')
	f.write('#                                                                                        											 #\n')
	f.write('#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT                						 #\n')
	f.write('#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       			 #\n')
	f.write('#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              				 #\n')
	f.write('#    details.                                                                            										 #\n')
	f.write('#                                                                                        											 #\n')
	f.write('#    You should have received a copy of the GNU General Public License along with        					 #\n')
	f.write('#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                 							 #\n')
	f.write('#                                                                                        											 #\n')
	f.write('##########################################################################################\n')
	f.write('\'\'\'solving the sensitivity (Hz/ppt) of instrument to molar mass (g/mol)\'\'\'\n')
	f.write('# module to estimate the sensitivity of an instrument to the molar mass of components, for example a Chemical Ionisiation Mass Spectrometer\n')
	f.write('# File Created at %s\n' %(datetime.datetime.now()))	
	f.write('\n')
	f.write('import numpy as np\n')
	f.write('\n')
	f.write('# function for sensitivity\n')
	f.write('def sens2mm(caller, y_MW):\n')
	f.write('	\n')
	f.write('	# inputs: -----------------\n')
	f.write('	# caller - flag for the calling function\n')
	f.write('	# y_MW - molar mass (g/mol) of components in question\n')
	f.write('	# ---------------------------\n')
	f.write('	\n')
	f.write('	fac_per_comp = %s # sensitivity (Hz/ppt) per molar mass (g/mol) \n' %(sens_func))
	f.write('	fac_per_comp = np.array((fac_per_comp)).reshape(-1) # reshape \n')
	f.write('	if (len(fac_per_comp) == 1): # if just a single value then tile across components \n')
	f.write('		fac_per_comp = np.tile(fac_per_comp, len(y_MW)) # if just a single value then tile across components \n')
	f.write('	\n')
	f.write('	if (caller == 3): # called on to plot sensitivity to molar mass\n')
	f.write('		import matplotlib.pyplot as plt \n')
	f.write('		plt.ion()\n')
	f.write('		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))\n')
	f.write('		ax0.plot(y_MW, fac_per_comp)\n')
	f.write('		ax0.set_title(\'Sensitivity of instrument to molar mass\')\n')
	f.write('		ax0.set_ylabel(\'Sensitivity (fraction (0-1))\', size = 14)\n')
	f.write('		ax0.yaxis.set_tick_params(labelsize = 14, direction = \'in\', which = \'both\')\n')
	f.write('		ax0.set_xlabel(\'Molar Mass ($\mathrm{g\,mol^{-1}}$)\', fontsize=14)\n')
	f.write('		ax0.xaxis.set_tick_params(labelsize = 14, direction = \'in\', which = \'both\')\n')
	f.write('	\n')
	f.write('	return(fac_per_comp)')
	f.close() # close file
	
	# get sensitivity for each component
	import sens2mm
	importlib.reload(sens2mm)
	fac_per_mass = sens2mm.sens2mm(caller, y_MW)

	return(fac_per_mass)

# function for plotting mass:charge resolution
def write_mzres(caller, res_in, y_mw):

	import datetime

	# inputs: --------------------------
	# caller - flag for the calling function
	# res_in - inputs for the mass:charge resolution (start point and width descriptors)
	# y_mw - the molar mass range to describe (g/mol)
	# ------------------------------------

	# create new  file
	f = open('PyCHAM/mzres.py', mode='w')
	f.write('##########################################################################################\n')
	f.write('#                                                                                        											 #\n')
	f.write('#    Copyright (C) 2018-2022 Simon O\'Meara : simon.omeara@manchester.ac.uk                  				 #\n')
	f.write('#                                                                                       											 #\n')
	f.write('#    All Rights Reserved.                                                                									 #\n')
	f.write('#    This file is part of PyCHAM                                                         									 #\n')
	f.write('#                                                                                        											 #\n')
	f.write('#    PyCHAM is free software: you can redistribute it and/or modify it under              						 #\n')
	f.write('#    the terms of the GNU General Public License as published by the Free Software       					 #\n')
	f.write('#    Foundation, either version 3 of the License, or (at your option) any later          						 #\n')
	f.write('#    version.                                                                            										 #\n')
	f.write('#                                                                                        											 #\n')
	f.write('#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT                						 #\n')
	f.write('#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       			 #\n')
	f.write('#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              				 #\n')
	f.write('#    details.                                                                            										 #\n')
	f.write('#                                                                                        											 #\n')
	f.write('#    You should have received a copy of the GNU General Public License along with        					 #\n')
	f.write('#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                 							 #\n')
	f.write('#                                                                                        											 #\n')
	f.write('##########################################################################################\n')
	f.write('\'\'\'solving probability density function of mass:charge resolution\'\'\'\n')
	f.write('# module to estimate the probability density function that is demonstrative of an instrument\'s mass:charge resolution, for example a Chemical Ionisiation Mass Spectrometer\n')
	f.write('# File Created at %s\n' %(datetime.datetime.now()))	
	f.write('\n')
	f.write('import numpy as np\n')
	f.write('import scipy.stats as st\n')
	f.write('\n')
	f.write('# function for sensitivity\n')
	f.write('def mzres(caller, res_in, y_mw):\n')
	f.write('	\n')
	f.write('	# inputs: -----------------\n')
	f.write('	# caller - flag for the calling function\n')
	f.write('	# res_in - inputs for the mass:charge resolution (start point and width descriptors)\n')
	f.write('	# y_mw - molar mass (g/mol) of components in question\n')
	f.write('	# ---------------------------\n')
	f.write('	\n')
	f.write('	if (caller == 3): # called on to plot\n')
	f.write('		import matplotlib.pyplot as plt \n')
	f.write('		plt.ion()\n')
	f.write('		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))\n')
	f.write('	\n')
	f.write('	y_mw = np.array((y_mw)) # ensure numpy array rather than list\n')
	f.write('	comp_indx = [] # empty list to hold results\n')
	f.write('	comp_prob = [] # empty list to hold results\n')
	f.write('	maxmm = np.max(y_mw) + res_in[0]\n')
	f.write('	mm_acc = 0. + %s # count on accumulated molar mass (g/mol) \n' %(res_in[0]))
	f.write('	# get maximum probability possible\n')
	f.write('	pdfm = st.norm.pdf(mm_acc, mm_acc, %s)\n' %(res_in[1]))
	f.write('	# loop through until upper end of molar mass reached \n')
	f.write('	while (mm_acc < maxmm):\n')
	f.write('		pdf = st.norm.pdf(y_mw, mm_acc, %s)\n' %(res_in[1]))
	f.write('		try: # in case a maximum can be identified\n')
	f.write('			pdf = pdf/pdfm # ensure that probability at distribution peak is 1\n')
	f.write('			# minimum and maximum molar masses covered significantly by this resolution interval\n')
	f.write('			mm = [np.min(y_mw[pdf>1.e-2]), np.max(y_mw[pdf>1.e-2])]\n')
	
	f.write('			if (caller == 1): # if called from non-test module\n')
	f.write('				if (len(y_mw[pdf > 1.e-2]) > 0): # if components do contribute to this interval\n')
	f.write('					# store indices of components contributing to this mass:charge interval\n')
	f.write('					ci = (np.where((y_mw >= mm[0])*(y_mw <= mm[1]) == 1))[0][:]\n')
	f.write('					comp_indx.append(ci)\n')
	f.write('					# store probability of contribution to this resolution interval\n')
	f.write('					comp_prob.append(pdf[ci])\n')
	f.write('				else: # if components do not contribute to this interval\n')
	f.write('					comp_indx.append([])\n')
	f.write('					comp_prob.append([])\n')
	f.write('		except: # no maximum, so no components contributing to this interval\n')
	f.write('			if (caller == 1):\n')
	f.write('				comp_indx.append([])\n')
	f.write('				comp_prob.append([])\n')
	f.write('		if (caller == 3): # called on to plot\n')
	f.write('			ax0.plot(y_mw[pdf>1.e-2], pdf[pdf>1.e-2])\n')
	f.write('		mm_acc += res_in[0] # keep count on accumulated molar mass \n')
	f.write('	\n')
	f.write('	if (caller == 3): # called on to plot\n')
	f.write('		ax0.set_title(\'Sensitivity of instrument due to mass:charge resolution\')\n')
	f.write('		ax0.set_ylabel(\'Probability of inclusion in resolution interval (fraction (0-1))\', size = 14)\n')
	f.write('		ax0.yaxis.set_tick_params(labelsize = 14, direction = \'in\', which = \'both\')\n')
	f.write('		ax0.set_xlabel(\'Molar Mass ($\mathrm{g\,mol^{-1}}$)\', fontsize=14)\n')
	f.write('		ax0.xaxis.set_tick_params(labelsize = 14, direction = \'in\', which = \'both\')\n')
	f.write('	\n')
	f.write('	# remember the range of molar masses representing mass:charge resolution\n')
	f.write('	mm_all = np.arange((0. + res_in[0]), (mm_acc), res_in[0]) \n')
	f.write('	return(pdf, comp_indx, comp_prob, mm_all)')
	f.close() # close file
	
	# get resolution result for mass spectrum
	import mzres
	importlib.reload(mzres)
	[pdf, comp_indx, comp_prob, mm_all] = mzres.mzres(caller, res_in, y_mw)

	return(pdf, comp_indx, comp_prob, mm_all)