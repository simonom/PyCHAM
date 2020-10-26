'''code to plot results from AtChem2 and compare against PyCHAM results'''
# introduction
# aim is to calculate deviation of PyCHAM photochemistry output from AtChem2
# use the AtChem2_apinene_scheme.txt in the Results folder of the GMD_paper for chemical 
# scheme input
# use the xml file from the PyCHAM inputs folder
# use the Photo_chem_inputs_hiNOx.txt and Photo_chem_inputs_loNOx.txt of the 
# GMD_paper/Results folder, for the high and low NOx simulations, respectively, where for 
# the low NOx case, initial concentration of NO2 is set to 0 ppb whilst for the high NOx 
# case it's set to 9.8 ppb
# important that exactly the same chemical scheme as used in the AtChem2 simulation is
# used, this saved as AtChem2_apinene_scheme.txt in the same folder as this file
# for the simulations to compare with AtChem2 use:
# update_step = 90.0
# recording_time_step = 90.0
# whereas, for the temporal resolution sensitivity simulations use the appropriate
# resolution for update_step and set recording_time_step to be the same

# PyCHAM inputs are saved in GMD_paper_plotting_scripts/Photo_chem_inputs_hiNOx.txt

# to run AtChem2:

# prepare conda environment with python 3.6:
# conda create -n AtChem2 python=3.6 numpy scipy
# activate:
# conda activate AtChem2
# clone the AtChem2 repository to wanted folder:
# git clone https://github.com/AtChem/AtChem2.git
# cd into this new folder
# create new folder inside new repository to hold dependencies:
# mkdir atchem-lib
# check all dependencies listed in the wiki (https://github.com/AtChem/AtChem2/wiki) are available;
# fortran compiler:
# which gfortran
# python:
# python -V
# cmake:
# cmake --version
# check on BLAS and LAPACK from inside python:
# python
# import numpy as np
# np.__config__.show()
# quit()
# state the path for LAPACK libararies inside the file tools/install/install_cvode.sh (line beginning LAPACK_LIBS=)
# note that despite the last direction the installation of Sundials on my linux stated that BLAS and LAPACK tests failed
# so that they wouldn't be able to support, however this didn't prevent AtChem2 functioning
# in the folder tools/install, install dependencies, where the third argument states the location of fortran compiler:
# ./install_cvode.sh ~/Documents/AtChem2/AtChem2/atchem-lib /usr/bin/gfortran
# ./install_openlibm.sh ~/Documents/AtChem2/AtChem2/atchem-lib /usr/bin/gfortran
# copy AtChem2/tools/install/Makefile.skel to the home directory and rename Makefile
# inside the renames Makefile set the paths to cvode and openlibm (on linux this was):
# CVODELIB     = /home/simonom/Documents/AtChem2/AtChem2/atchem-lib/cvode/lib
# OPENLIBMDIR  = /home/simonom/Documents/AtChem2/AtChem2/atchem-lib/openlibm-0.4.1
# to install and compile AtChem2:
# build/build_atchem2.sh mcm/mechanism_test.fac
# if successful, this will produce a file call atchem2, which is an executable
# to run:
# ./atchem2

# to compare with PyCHAM, ensure both AtChem2 and PyCHAM use the chemical scheme: AtChem2_apinene_scheme - for AtChem2 this should have a .fac file extension and be stored in the AtChem2 mcm folder
# for AtChem2 model variables are stored in the model/configuration folder.  For both the high and low NOx case set APINENE and O3 concentrations inside the initialConcentrations.config to 5.2e11, whilst for high NOx set NO2 to 2.4123e11 and for low NOx set NO2 to 0.
# inside outputRates.config state HCHO, inside outputSpecies.config set to:
# APINENE
# O3
# NO2
# NO
# OH
# HO2
# CH3O2
# O
# HCHO
# check that model.parameters and environmentVariables.config agree with PyCHAM

# to build:
# build/build_atchem2.sh mcm/AtChem2_apinene_scheme.fac
# then to run:
# ./atchem2
# then copy the model/output/speciesConcentrations.output and lossRates.output and productionRates.output
# to the required folder 

# assumes calling from the PyCHAM home folder

import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os

# get current working directory
cwd = os.getcwd()
# ensure modules can be seen 
# (assumes calling from the home folder)
sys.path.append(str(os.getcwd() + '/PyCHAM'))
import retr_out

# ----------------------------------------------------------------------------------------
# AtChem2 part
# open saved files
Atfname = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/photo_chem_data/AtChem2_APINENE/hiNOx/speciesConcentrations.output')
inputs = open(Atfname, mode='r') # open results
# read the file and store everything into a list
in_list = inputs.readlines()
inputs.close() # close file
for i in range(len(in_list)):
	if i==0: # obtain headers
		comp_names = [] # empty list
		comp_nam = (in_list[i].split(' ')) # list of component names
		for ii in comp_nam:
			if len(ii.strip())==0: # omit white space
				continue
			else:
				comp_names.append(ii.strip())

		# create empty array to hold results
		gconc = np.empty((len(in_list)-1, len(comp_names)))

	else:
		compC0 = (in_list[i].split(' ')) # list of component concentrations
		comp_num = 0 # number of components at this time
		for ii in compC0:
			if len(ii.strip())==0: # omit white space
				continue
			else:
				gconc[i-1, comp_num] = ii.strip()
				comp_num += 1

# check that first column is time, otherwise throw error and exit
if str(comp_names[0]) != 't' and str(comp_names[0]) != 'time':
	sys.exit('Error, first column is not titled t or time, instead it is called: ' + str(comp_names[0]))
else: # convert to hours from s
	gconc[1::, 0] = gconc[1::, 0]/3600.
# ----------------------------------------------------------------------------------------
# PyCHAM part
# file name
Pyfname = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/photo_chem_data/PyCHAM_APINENE/hiNOx/PyCHAM_comp_hiNOx')

# required outputs
(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, xfm, t_array, PyCHAM_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(Pyfname)

# convert from ppb to molecules/cc (air)
yrec = yrec*((np.array((Cfac))).reshape(-1,1))

# ----------------------------------------------------------------------------------------
# comparative statistics

# if AtChem values don't synchronise with PyCHAM, then interpolate, note we assume the 
# time arrays have the same units
if (len(gconc[:, 0]) != len(t_array)):
	interp_flag = 1
elif np.abs(sum(gconc[:, 0]-t_array))>1.0e-3:
	interp_flag = 1
else:
	interp_flag = 0
	
if interp_flag == 1:
	gconc_int = np.empty((len(t_array), len(comp_names)))
	# adopt PyCHAM times
	gconc_int[:, 0] = t_array

	# loop through the component names from the AtChem array
	for i in range( 1, len(gconc[1::])): # loop through AtChem components:
		# interpolate to PyCHAM times
		gconc_int[:, i] = np.interp(gconc_int[:, 0], gconc[:, 0], gconc[:, i])

else:
	gconc_int = gconc

# empty array for fractional deviation
frac_dev = np.zeros((len(t_array), len(gconc[1::])))
# loop through the component names from the AtChem array
for i in range(1, comp_num): # loop through AtChem components
	
	# find index of corresponding component in PyCHAM results
	ind = PyCHAM_names.index(comp_names[i])
	nz_ind = np.where(gconc_int[:, i]!=0)
	frac_dev[nz_ind, i-1] = ((yrec[nz_ind, ind]-gconc_int[nz_ind, i])/np.max(gconc_int[nz_ind, i]))*100.0



# ----------------------------------------------------------------------------------------
# make plot with gas-phase concentration deviations shown
compnum = 0 # count on components
fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(8,6), sharex=True)
for i in comp_names[1::]: # loop through components (excluding time in column 0)
	if str(i)=="CH3O2" or str(i)=="HO2" or str(i)=="O" or str(i)=="CO":
		compnum += 1
		continue
	if (str(i) == "OH"): # force color for OH
		ax0.plot(t_array, frac_dev[:, compnum], 'm', label=str(i))
		compnum += 1		
		continue	
	if (str(i) == "HCHO"): # force color for HCHO
		ax0.plot(t_array, frac_dev[:, compnum], 'k', label=str(i))
		compnum += 1		
		continue
	
	ax0.plot(t_array, frac_dev[:, compnum], label=str(i))	
	compnum += 1

# plt.title('PyCHAM-AtChem2 Fractional Deviation of Gas-phase Concentration')
ax0.set_ylabel(r'Deviation (%)', fontsize=14)
ax0.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
ax0.xaxis.set_tick_params(direction = 'in', which = 'both')
ax0.legend(fontsize = 12, loc = 'lower right')
ax0.text(x=-1.8, y=0.7, s='(a)', size=14)
ax0.yaxis.set_tick_params(direction = 'in', which = 'both')


# ----------------------------------------------------------------------------------------
# temporal profile of gas-phase concentrations (ppb) (used for EAC abstract)

# fig, (ax4) = plt.subplots(1, 1, figsize=(8,6))
# comp_num = 1
# for i in comp_names[1::]: # loop through Atchem2 components (excluding time in column 0)
# 	if i=='HO2' or i=='O':
# 		comp_num += 1 # Atchem2 index
# 		continue
# 	ind = PyCHAM_names.index(i) # PyCHAM index for this component
# 	# temporal (24 hour clock) profiles of gas-phase concentration (ppb), note gas-phase
# 	# concentration conversion from molecules/cc (air) to ppb (air) 
# 	ax4.semilogy(gconc[:,0]/3600.0, gconc[:,comp_num]/Cfactor, 
# 						label=str('Atchem2 '+str(i)))
# 	ax4.semilogy(t_array/3600.0, y[:, ind], '--', 
# 						label=str('PyCHAM '+str(i)))
# 	comp_num += 1 # Atchem2 index
# ax4.set_ylabel(r'Gas-phase concentration (ppb)', fontsize=12)
# ax4.set_xlabel(r'Time of day', fontsize=12)
# ax4.yaxis.set_tick_params(size=12)
# ax4.xaxis.set_tick_params(size=12)
# ax4.set_ylim([1.0e-10, 2.0e2])
# ax4.legend(fontsize=10)
# plt.show()

#-----------------------------------------------------------------------------------------
# AtChem2 part
# open saved files
Atfname = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/photo_chem_data/AtChem2_APINENE/loNOx/speciesConcentrations.output')

inputs = open(Atfname, mode='r') # open results
# read the file and store everything into a list
in_list = inputs.readlines()
inputs.close() # close file
for i in range(len(in_list)):
	if i==0: # obtain headers
		comp_names = [] # empty list
		comp_nam = (in_list[i].split(' ')) # list of component names
		for ii in comp_nam:
			if len(ii.strip())==0: # omit white space
				continue
			else:
				comp_names.append(ii.strip())

		# create empty array to hold results
		gconc = np.empty((len(in_list)-1, len(comp_names)))

	else:
		compC0 = (in_list[i].split(' ')) # list of component concentrations
		comp_num = 0 # number of components at this time
		for ii in compC0:
			if len(ii.strip())==0: # omit white space
				continue
			else:
				gconc[i-1, comp_num] = ii.strip()
				comp_num += 1

# check that first column is time, otherwise throw error and exit
if str(comp_names[0]) != 't' and str(comp_names[0]) != 'time':
	sys.exit('Error, first column is not titled t or time, instead it is called: ' + str(comp_names[0]))
else: # convert to hours from s
	gconc[1::, 0] = gconc[1::, 0]/3600.
# ----------------------------------------------------------------------------------------
# PyCHAM part
# file name
Pyfname = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/photo_chem_data/PyCHAM_APINENE/loNOx/PyCHAM_comp_loNOx')

# required outputs
(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array2, PyCHAM_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(Pyfname)

# convert from ppb to molecules/cc (air)
y = y*((np.array((Cfac))).reshape(-1, 1))

# ----------------------------------------------------------------------------------------
# comparative statistics
# empty fractional deviation matrix
frac_dev2 = np.empty((len(t_array2), len(comp_names)-1))

# if AtChem values don't synchronise with PyCHAM, then interpolate, note we assume the 
# time arrays have the same units
if len(gconc[:,0])!=len(t_array2):
	interp_flag = 1
elif np.abs(sum(gconc[:,0]-t_array2))>1.0e-3:
	interp_flag = 1
else:
	interp_flag = 0

if interp_flag == 1:
	gconc_int = np.empty((len(t_array2), len(comp_names)))
	# adopt PyCHAM times
	gconc_int[:,0] = t_array2

	# loop through the component names from the AtChem array
	for i in range( 1, len(comp_names[1::])): # loop through AtChem components:
		# interpolate to PyCHAM times
		gconc_int[:, i] = np.interp(gconc_int[:, 0], gconc[:, 0], gconc[:, i])

else:
	gconc_int = gconc

# empty array for fractional deviation
frac_dev2 = np.zeros((len(t_array2), len(comp_names[1::])))
# loop through the component names from the AtChem array
for i in range(1, len(comp_names[1::])): # loop through AtChem components:
	# find index of corresponding component in PyCHAM results
	ind = PyCHAM_names.index(comp_names[i])
	nz_ind = np.where(gconc_int[:, i]!=0)
	frac_dev2[nz_ind, i-1] = ((y[nz_ind, ind]-gconc_int[nz_ind, i])/np.max(gconc_int[nz_ind, i]))*100.0



# ----------------------------------------------------------------------------------------
# make plot with all gas-phase concentrations shown
compnum = 0 # count on components

for i in comp_names[1::]: # loop through components (excluding time in column 0)
	if str(i)=="CH3O2" or str(i)=="HO2" or str(i)=="O" or str(i)=="CO" or str(i)=="NO2" or str(i)=="NO":
		compnum += 1
		continue
	
	if (str(i) == "OH"): # force color for OH
		ax1.plot(t_array2, frac_dev2[:, compnum], 'm', label=str(i))
		compnum += 1		
		continue	
	if (str(i) == "HCHO"): # force color for HCHO
		ax1.plot(t_array2, frac_dev2[:, compnum], 'k', label=str(i))
		compnum += 1		
		continue
	
	ax1.plot(t_array2, frac_dev2[:,compnum], label=str(i))	
	compnum += 1

# plt.title('PyCHAM-AtChem2 Fractional Deviation of Gas-phase Concentration')
ax1.set_ylabel(r'Deviation (%)', fontsize=14)
ax1.legend(fontsize=12)
ax1.set_xlabel(r'Time of day (hours)', fontsize=14)
ax1.text(x=-1.8, y=0.07, s='(b)', size=14)
ax1.xaxis.set_tick_params(direction = 'in', which = 'both')
ax1.yaxis.set_tick_params(direction = 'in', which = 'both')
fig.savefig('fig03.png')
plt.show()

# ----------------------------------------------------------------------------------------
# comparison of reaction rates recorded by PyCHAM and AtChem2

# first, the high NOx case

# open AtChem2 loss reaction results
Atfname = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/photo_chem_data/AtChem2_APINENE/hiNOx/lossRates.output')

inputs = open(Atfname, mode= 'r') # open results
# read the file and store everything into a list
in_list = inputs.readlines()
inputs.close() # close file

# empty list for AtChem2 reaction rates (molecules/cc (air).s)
# first row for time, second for reaction number and third for rate
Atres = np.zeros((0, 3))

for i in range(1, len(in_list)): # loop through lines in list, skipping heading row
	Atres = np.append(Atres, np.zeros((1, 3)), axis=0) # create empty row to hold results
	
	# loop through the line elements to check if reaction number is the desired
	elem_num = 0
	
	for ii in (in_list[i].split(' ')): # loop through elements in line
		if ii == ' ' or ii == '' : # ignore white space
			continue
		else:
			elem_num += 1 # record number of inputs in row
			
			if elem_num == 1: # remember time (s)
				Atres[i-1, 0] = float(ii)
			if elem_num == 4: # record reaction number
				Atres[i-1, 1] = float(ii)
			if (elem_num==5): # record reaction rate (molecules/cc (air).s)
				Atres[i-1, 2] = float(ii)
					
num_lossr = int(len(np.unique(Atres[:, 1]))) # number of loss reactions
				
# open AtChem2 production reaction results
Atfname = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/photo_chem_data/AtChem2_APINENE/hiNOx/productionRates.output')

inputs = open(Atfname, mode='r') # open results
# read the file and store everything into a list
in_list = inputs.readlines()
inputs.close() # close file

# empty list for AtChem2 reaction rates (molecules/cc (air).s)
# first row for time, second for reaction number and third for rate
Atres_prod = np.zeros((0, 3))

for i in range(1, len(in_list)): # loop through lines in list, skipping heading row
	Atres_prod = np.append(Atres_prod, np.zeros((1, 3)), axis=0) # empty column to hold results
	
	# loop through the line elements to check if reaction number is the desired
	elem_num = 0
	
	for ii in (in_list[i].split(' ')): # loop through elements in line
		if ii == ' ' or ii == '' : # ignore white space
			continue
		else:
			elem_num += 1 # record number of inputs in row
			
			if elem_num == 1: # remember time (s)
				Atres_prod[i-1, 0] = float(ii)
			if elem_num == 4: # record reaction number
				Atres_prod[i-1, 1] = float(ii)
			if (elem_num == 5): # record reaction rate (molecules/cc (air).s)
				Atres_prod[i-1, 2] = float(ii)

num_prodr = int(len(np.unique(Atres_prod[:, 1]))) # number of production reactions

# open the PyCHAM tracking results
Pyfname = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/photo_chem_data/PyCHAM_APINENE/hiNOx/PyCHAM_comp_hiNOx')

fname = str(Pyfname+'/HCHO_rate_of_change')
dydt = np.loadtxt(fname, delimiter = ',', skiprows = 1) # skiprows = 1 omits header

fname = str(Pyfname+'/time')
Pyt = np.loadtxt(fname, delimiter = ',', skiprows = 1) # skiprows = 1 omits header

PyCHAM_reac_num = dydt[0, :] # reaction numbers stored in PyCHAM results

Pyres = np.zeros((Atres.shape[0], Atres.shape[1]))

# for loss reactions
Pyres[:, 0] = Atres[:, 0] # times (s)
# reaction numbers, subtract one to account for Fortran indexing
Pyres[:, 1] = Atres[:, 1]-1

# PyCHAM indices for loss reactions
Ati = 0 # count on AtChem  times

for i in range(len(Pyt)):
	
	if Pyt[i] == Atres[Ati*num_lossr, 0]:
		
		for ii in range(4):
			# -1 from Atres index to account for the Fortran indexing
			PyCHAM_ind = np.where(PyCHAM_reac_num==Atres[ii, 1]-1) 
			# record reaction rates
			Pyres[Ati*num_lossr+ii, 2] = dydt[i, PyCHAM_ind]
		Ati += 1
		
# repeat for production reactions
Pyres_prod = np.zeros((Atres_prod.shape[0], Atres_prod.shape[1]))		
Pyres_prod[:, 0] = Atres_prod[:, 0] # times (s)
# reaction numbers, subtract one to account for Fortran indexing
Pyres_prod[:, 1] = Atres_prod[:, 1]-1

# PyCHAM indices for production reactions
Ati = 0 # count on AtChem  times
for i in range(len(Pyt)): # loop through PyCHAM times
	
	if Pyt[i] == Atres_prod[Ati*num_prodr, 0]: # if times match
		
		for ii in range(num_prodr):
			
			# -1 from Atres index to account for the Fortran indexing
			PyCHAM_ind = np.where(PyCHAM_reac_num == (Atres_prod[ii, 1]-1))
			
			# record production reaction rates
			Pyres_prod[Ati*num_prodr+ii, 2] = dydt[i, PyCHAM_ind]
		Ati += 1 # keep count on AtChem times

# ----------------------------------------------------------------------------------------
# comparative statistics

# empty array for fractional deviation, first row for times, second for equation number
# and third for fractional deviation
frac_dev_loss = np.zeros((Atres[:, 0].shape[0], 1))
# loop through loss reactions
for i in range(num_lossr):
	# % difference, multiply AtChem result by -1 to represent loss
	frac_dev_loss[i::num_lossr, 0] = ((Pyres[i::num_lossr, 2]-(-1.0*Atres[i::num_lossr, 2]))/
									np.max(Atres[i::num_lossr, 2]))*100.0

# empty array for fractional deviation, first row for times, second for equation number
# and third for fractional deviation
frac_dev_prod = np.zeros((Atres_prod[:, 0].shape[0], 1))
# loop through loss reactions
for i in range(num_prodr):
	# % difference
	frac_dev_prod[i::num_prodr, 0] = ((Pyres_prod[i::num_prodr, 2]-(Atres_prod[i::num_prodr, 2]))/
				np.max(Atres_prod[i::num_prodr, 2]))*100.0


# identify the two loss and two production reactions with greatest deviation
max_dev1 = 0.
max_dev0 = 0.
loss_reac_i = np.zeros((2)).astype(int)
for i in range(num_lossr): # loop through loss reactions
	# maximum deviation for this reaction	
	max_dev2 = np.max(np.abs(frac_dev_loss[i::num_lossr, 0]))
	
	if (max_dev2>max_dev0): # highest maximum seen so far
		max_dev0 = max_dev2 # track maximum seen
		loss_reac_i[0] = int(i) # record index of highest maximum
	if (max_dev2<max_dev0 and max_dev2>max_dev1): # second highest maximum seen so far
		max_dev1 = max_dev2
		loss_reac_i[1] = int(i) # record second highest maximum

# identify the two production reactions with greatest deviation
max_dev1 = 0.
max_dev0 = 0.
prod_reac_i = np.zeros((2)).astype(int)
for i in range(num_prodr): # loop through production reactions
	# maximum deviation for this reaction
	max_dev2 = np.max(np.abs(frac_dev_prod[i::num_prodr, 0]))
	
	if (max_dev2>max_dev0): # highest maximum seen so far
		max_dev0 = max_dev2 # track maximum seen
		prod_reac_i[0] = int(i) # record index of highest maximum
	if (max_dev2<max_dev0 and max_dev2>max_dev1): # second highest maximum seen so far
		max_dev1 = max_dev2
		prod_reac_i[1] = int(i) # record second highest maximum

# print the reaction numbers where disagreement is greatest (Python indexing, starting at 
# 0, so add 1 when comparing with line numbers in equation text file)
print('for the high NOx case (Python indexing, starting at 0, so add 1 when comparing with line numbers in equation text file): ')
print('equation numbers producing two highest deviations for loss: ', Pyres[loss_reac_i[:], 1])
print('equation numbers producing two highest deviations for production: ', Pyres_prod[prod_reac_i[:], 1])
print('now manually go find the equations these numbers relate to and add to plot labels')

# plot results for hihg NOx case
fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10,6), sharex=True) # prepare plot
ax0.plot(Pyres[loss_reac_i[0]::num_lossr, 0]/3600., frac_dev_loss[loss_reac_i[0]::num_lossr, 0], linewidth = 3, label=r'$\mathrm{HCHO\rightarrow CO + 2HO2}$')
ax0.plot(Pyres[loss_reac_i[1]::num_lossr, 0]/3600., frac_dev_loss[loss_reac_i[1]::num_lossr, 0], linewidth = 3, label=r'$\mathrm{HCHO\rightarrow H2 + CO}$')
ax0.plot(Pyres_prod[prod_reac_i[0]::num_prodr, 0]/3600., frac_dev_prod[prod_reac_i[0]::num_prodr, 0], linewidth = 3, label=r'$\mathrm{C516O\rightarrow CO13C3CO2H + HCHO + HO2}$')
ax0.plot(Pyres_prod[prod_reac_i[1]::num_prodr, 0]/3600., frac_dev_prod[int(prod_reac_i[1])::num_prodr, 0], linewidth = 3, label=r'$\mathrm{C6140\rightarrow CO23C4CHO + HCHO + HO2}$')
ax0.text(x=-0.5, y=1.5, s='(a)', size=16)
ax0.set_ylabel(r'Deviation (%)', fontsize=16)
ax0.legend(fontsize=12, loc = 'lower right')
ax0.xaxis.set_tick_params(direction = 'in', which = 'both')
ax0.yaxis.set_tick_params(direction = 'in', which = 'both')
ax0.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))

# ----------------------------------------------------------------------------------------
# Repeat for low NOx case

# open AtChem2 loss reaction results
Atfname = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/photo_chem_data/AtChem2_APINENE/loNOx/lossRates.output')

inputs = open(Atfname, mode='r') # open results
# read the file and store everything into a list
in_list = inputs.readlines()
inputs.close() # close file

# empty list for AtChem2 reaction rates (molecules/cc (air).s)
# first row for time, second for reaction number and third for rate
Atres = np.zeros((0, 3))

for i in range(1, len(in_list)): # loop through lines in list, skipping heading row
	Atres = np.append(Atres, np.zeros((1, 3)), axis=0) # create empty row to hold results
	
	# loop through the line elements to check if reaction number is the desired
	elem_num = 0
	
	for ii in (in_list[i].split(' ')): # loop through elements in line
		if ii == ' ' or ii == '' : # ignore white space
			continue
		else:
			elem_num += 1 # record number of inputs in row
			
			if elem_num == 1: # remember time (s)
				Atres[i-1, 0] = float(ii)
			if elem_num == 4: # record reaction number
				Atres[i-1, 1] = float(ii)
			if (elem_num==5): # record reaction rate (molecules/cc (air).s)
				Atres[i-1, 2] = float(ii)

num_lossr = int(len(np.unique(Atres[:, 1]))) # number of loss reactions					
					
# open AtChem2 production reaction results
Atfname = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/photo_chem_data/AtChem2_APINENE/loNOx/productionRates.output')

inputs = open(Atfname, mode='r') # open results
# read the file and store everything into a list
in_list = inputs.readlines()
inputs.close() # close file

# empty list for AtChem2 reaction rates (molecules/cc (air).s)
# first row for time, second for reaction number and third for rate
Atres_prod = np.zeros((0, 3))

for i in range(1, len(in_list)): # loop through lines in list, skipping heading row
	Atres_prod = np.append(Atres_prod, np.zeros((1, 3)), axis=0) # create empty row to hold results
	
	# loop through the line elements to check if reaction number is the desired
	elem_num = 0
	
	for ii in (in_list[i].split(' ')): # loop through elements in line
		if ii == ' ' or ii == '' : # ignore white space
			continue
		else:
			elem_num += 1 # record number of inputs in row
			
			if elem_num == 1: # remember time (s)
				Atres_prod[i-1, 0] = float(ii)
			if elem_num == 4: # record reaction number
				Atres_prod[i-1, 1] = float(ii)
			if (elem_num==5): # record reaction rate (molecules/cc (air).s)
				Atres_prod[i-1, 2] = float(ii)

num_prodr = int(len(np.unique(Atres_prod[:, 1]))) # number of production reactions

# open the PyCHAM tracking results
Pyfname = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/photo_chem_data/PyCHAM_APINENE/loNOx/PyCHAM_comp_loNOx')

fname = str(Pyfname+'/HCHO_rate_of_change')
dydt = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header

fname = str(Pyfname+'/time')
Pyt = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header

PyCHAM_reac_num = dydt[0, :] # reaction numbers stored in PyCHAM results
Pyres = np.zeros((Atres.shape[0], Atres.shape[1]))

# for loss reactions
Pyres[:, 0] = Atres[:, 0] # times (s)
# reaction numbers, -1 to account for Fortran indexing
Pyres[:, 1] = Atres[:, 1]-1 

# PyCHAM indices for loss reactions
Ati = 0 # count on AtChem  times

for i in range(len(Pyt)): # loop through PyCHAM times
	
	if Pyt[i] == Atres[Ati*num_lossr, 0]:
		
		for ii in range(num_lossr): # loop through HCHO loss reactions
			# -1 from Atres index to account for the Fortran indexing
			PyCHAM_ind = np.where(PyCHAM_reac_num == Atres[ii, 1]-1) 
			# record reaction rates
			Pyres[Ati*num_lossr+ii, 2] = dydt[i, PyCHAM_ind]
		Ati += 1 # keep count on AtChem times
		
# repeat for production reactions
Pyres_prod = np.zeros((Atres_prod.shape[0], Atres_prod.shape[1]))		
Pyres_prod[:, 0] = Atres_prod[:, 0] # times (s)
# reaction numbers, -1 to account for Fortran indexing
Pyres_prod[:, 1] = Atres_prod[:, 1]-1

# PyCHAM indices for production reactions
Ati = 0 # count on AtChem  times
for i in range(len(Pyt)):
	if Pyt[i] == Atres_prod[Ati*num_prodr, 0]:
		
		for ii in range(num_prodr):
			# -1 from Atres index to account for the Fortran indexing
			PyCHAM_ind = np.where(PyCHAM_reac_num == Atres_prod[ii, 1]-1) 
			# record reaction rates
			Pyres_prod[Ati*num_prodr+ii, 2] = dydt[i, PyCHAM_ind]
		Ati += 1

# ----------------------------------------------------------------------------------------
# comparative statistics

# empty array for fractional deviation, first row for times, second for equation number
# and third for fractional deviation
frac_dev_loss = np.zeros((Atres[:, 0].shape[0], 1))

for i in range(num_lossr): # loop through loss reactions
	
	# if maximum rate doesn't exceed 1 molecule/cc.s (air), then ignore
	if (np.max(Atres[i::num_lossr, 2])<1.):
		frac_dev_loss[i::num_lossr, 0] = 0.
		continue
	
	# % difference, multiply AtChem result by -1 to represent loss
	frac_dev_loss[i::num_lossr, 0] = ((Pyres[i::num_lossr, 2]-(-1.0*Atres[i::num_lossr, 2]))/
									np.max(Atres[i::num_lossr, 2]))*100.0

							
# empty array for fractional deviation, first row for times, second for equation number
# and third for fractional deviation
frac_dev_prod = np.zeros((Atres_prod[:, 0].shape[0], 1))
# loop through loss reactions
for i in range(num_prodr):
	# if maximum rate doesn't exceed 1 molecule/cc.s (air), then ignore
	if np.max(Atres_prod[i::num_prodr, 2])<1.0:
		continue
	# % difference
	frac_dev_prod[i::num_prodr, 0] = ((Pyres_prod[i::num_prodr, 2]-(Atres_prod[i::num_prodr, 2]))/
									np.max(Atres_prod[i::num_prodr, 2]))*100.0

# identify the two loss and two production reactions with greatest deviation
max_dev1 = 0.
max_dev0 = 0.
loss_reac_i = np.zeros((2)).astype(int)
for i in range(num_lossr): # loop through loss reactions
	# maximum deviation for this reaction	
	max_dev2 = np.max(np.abs(frac_dev_loss[i::num_lossr, 0]))
	
	if (max_dev2 > max_dev0): # highest maximum seen so far
		max_dev0 = max_dev2 # track maximum seen
		loss_reac_i[0] = int(i) # record index of highest maximum
	if (max_dev2 < max_dev0 and max_dev2 > max_dev1): # second highest maximum seen so far
		max_dev1 = max_dev2
		loss_reac_i[1] = int(i) # record second highest maximum

# identify the two loss and two production reactions with greatest deviation
max_dev1 = 0.0
max_dev0 = 0.0
prod_reac_i = np.zeros((2)).astype(int)
for i in range(num_prodr): # loop through loss reactions
	# maximum deviation for this reaction	
	max_dev2 = np.max(np.abs(frac_dev_prod[i::num_prodr, 0]))
	
	if (max_dev2 > max_dev0): # highest maximum seen so far
		max_dev0 = max_dev2 # track maximum seen
		prod_reac_i[0] = int(i) # record index of highest maximum
	if (max_dev2 < max_dev0 and max_dev2 > max_dev1): # second highest maximum seen so far
		max_dev1 = max_dev2
		prod_reac_i[1] = int(i) # record second highest maximum

# print the reaction numbers where disagreement is greatest (Python indexing, starting at 
# 0, so add 1 when comparing with line numbers in equation text file)
print('for the low NOx case (Python indexing, starting at 0, so add 1 when comparing with line numbers in equation text file): ')
print('equation numbers producing two highest deviations for loss: ', Pyres[loss_reac_i[:], 1])
print('equation numbers producing two highest deviations for production: ', Pyres_prod[prod_reac_i[:], 1])
print('now manually go find the equations these numbers relate to and add to plot labels')

# plot results for low NOx case
ax1.plot(Pyres[loss_reac_i[0]::num_lossr, 0]/3600.0, frac_dev_loss[loss_reac_i[0]::num_lossr, 0], linewidth=3, label=r'$\mathrm{HCHO\rightarrow CO + 2HO2}$')
ax1.plot(Pyres[loss_reac_i[1]::num_lossr, 0]/3600.0, frac_dev_loss[loss_reac_i[1]::num_lossr, 0], linewidth=3, label=r'$\mathrm{HCHO\rightarrow H2 + CO}$')
ax1.plot(Pyres_prod[prod_reac_i[0]::num_prodr, 0]/3600., frac_dev_prod[prod_reac_i[0]::num_prodr, 0], linewidth=3, label=r'$\mathrm{HOCH2CO3\rightarrow HCHO + HO2}$')
ax1.plot(Pyres_prod[prod_reac_i[1]::num_prodr, 0]/3600.0, frac_dev_prod[prod_reac_i[1]::num_prodr, 0], linewidth=3, label=r'$\mathrm{C621O\rightarrow HCHO + H1C23C4CHO + HO2}$')


ax1.text(x=-0.5, y=3.7, s='(b)', size=16)
ax1.set_xlabel(r'Time of day (hours)', fontsize=16)
ax1.set_ylabel(r'Deviation (%)', fontsize=16)
ax1.legend(fontsize=12, loc = [0.6, 0.6,])
ax1.yaxis.set_tick_params(direction = 'in', which = 'both')
ax1.xaxis.set_tick_params(direction = 'in', which = 'both')
ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))


fig.savefig('fig04.png')
plt.show()