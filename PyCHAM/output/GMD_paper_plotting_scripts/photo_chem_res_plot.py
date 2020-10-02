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

import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os

# get current working directory
cwd = os.getcwd()

# ----------------------------------------------------------------------------------------
# AtChem2 part
# open saved files
Atfname = str(cwd + '/photo_chem_data/AtChem2_APINENE/hiNOx/speciesConcentrations.output')
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
# ----------------------------------------------------------------------------------------
# PyCHAM2 part
# file name
Pyfname = str(cwd + '/photo_chem_data/PyCHAM_APINENE/hiNOx/AtChem2_comp_hiNOx')

# name of file where experiment constants saved (number of size bins and whether wall 
# included)
fname = str(Pyfname+'/model_and_component_constants')

const_in = open(fname)
const = {} # prepare to create dictionary
for line in const_in.readlines():

	# convert to python list
	dlist = []
	for i in line.split(',')[1::]:
		if str(line.split(',')[0]) == 'number_of_size_bins':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'number_of_components':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'molecular_weights_g/mol_corresponding_to_component_names' or  str(line.split(',')[0]) == 'molecular_volumes_cm3/mol':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			dlist.append(float(i))
		if str(line.split(',')[0]) == 'component_names':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			i = i.strip('\'')
			dlist.append(str(i))
		if str(line.split(',')[0]) == 'factor_for_multiplying_ppb_to_get_molec/cm3':
			dlist.append(float(i))
			
	const[str(line.split(',')[0])] = dlist

num_sb = int((const['number_of_size_bins'])[0]) # number of size bins
num_speci = int((const['number_of_components'])[0]) # number of species
y_mw = const['molecular_weights_g/mol_corresponding_to_component_names']
y_MV = const['molecular_volumes_cm3/mol']
PyCHAM_names = const['component_names']
# conversion factor to change gas-phase concentrations from molecules/cc 
# (air) into ppb
Cfactor = float((const['factor_for_multiplying_ppb_to_get_molec/cm3'])[0])

# name of file where concentration (molecules/cc (air)) results saved
fname = str(Pyfname+'/concentrations_all_components_all_times_gas_particle_wall')
y = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header)
# convert gas-phase results from ppb to molecules/cm3 (air), for consistency with AtChem2
# results
y[:, 0:num_speci] = y[:, 0:num_speci]*Cfactor

# withdraw times
fname = str(Pyfname+'/time')
t_array = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header

# ----------------------------------------------------------------------------------------
# comparative statistics

frac_dev = np.empty((len(t_array), len(comp_names)-1)) # empty fractional deviation matrix

# if AtChem values don't synchronise with PyCHAM, then interpolate, note we assume the 
# time arrays have the same units
if len(gconc[:,0])!=len(t_array):
	interp_flag = 1
elif np.abs(sum(gconc[:,0]-t_array))>1.0e-3:
	interp_flag = 1
else:
	interp_flag = 0
	
if interp_flag == 1:
	gconc_int = np.empty((len(t_array), len(comp_names)))
	# adopt PyCHAM times
	gconc_int[:, 0] = t_array

	# loop through the component names from the AtChem array
	for i in range( 1, len(comp_names[1::])): # loop through AtChem components:
		# interpolate to PyCHAM times
		gconc_int[:, i] = np.interp(gconc_int[:, 0], gconc[:, 0], gconc[:, i])

else:
	gconc_int = gconc

# empty array for fractional deviation
frac_dev = np.zeros((len(t_array), len(comp_names[1::])))
# loop through the component names from the AtChem array
for i in range(1, len(comp_names[1::])): # loop through AtChem components:
	# find index of corresponding component in PyCHAM results
	ind = PyCHAM_names.index(comp_names[i])
	nz_ind = np.where(gconc_int[:, i]!=0)
	frac_dev[nz_ind, i-1] = ((y[nz_ind, ind]-
							gconc_int[nz_ind, i])/np.max(gconc_int[nz_ind, i]))*100.0



# ----------------------------------------------------------------------------------------
# make plot with gas-phase concentration deviations shown
compnum = 0 # count on components
fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(8,6), sharex=True)
for i in comp_names[1::]: # loop through components (excluding time in column 0)
	if str(i)=="CH3O2" or str(i)=="HO2" or str(i)=="O" or str(i)=="CO":
		compnum += 1
		continue
	ax0.plot(t_array/3600.0, frac_dev[:, compnum], label=str(i))
	compnum += 1
# plt.title('PyCHAM-AtChem2 Fractional Deviation of Gas-phase Concentration')
ax0.set_ylabel(r'Deviation (%)', fontsize=14)
ax0.yaxis.set_tick_params(size=14)
ax0.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
ax0.xaxis.set_tick_params(size=14)
ax0.legend(fontsize=12)
ax0.text(x=-1.8, y=8.6, s='(a)', size=14)
ax0.yaxis.set_tick_params(size=13)

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
Atfname = str(cwd + '/photo_chem_data/AtChem2_APINENE/loNOx/speciesConcentrations.output')

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

# ----------------------------------------------------------------------------------------
# PyCHAM2 part
# file name
Pyfname = str(cwd + '/photo_chem_data/PyCHAM_APINENE/loNOx/AtChem2_comp_loNOx')


# name of file where experiment constants saved (number of size bins and whether wall 
# included)
# name of file where experiment constants saved (number of size bins and whether wall 
# included)
fname = str(Pyfname+'/model_and_component_constants')

const_in = open(fname)
const = {} # prepare to create dictionary
for line in const_in.readlines():

	# convert to python list
	dlist = []
	for i in line.split(',')[1::]:
		if str(line.split(',')[0]) == 'number_of_size_bins':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'number_of_components':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'molecular_weights_g/mol_corresponding_to_component_names' or  str(line.split(',')[0]) == 'molecular_volumes_cm3/mol':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			dlist.append(float(i))
		if str(line.split(',')[0]) == 'component_names':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			i = i.strip('\'')
			dlist.append(str(i))
		if str(line.split(',')[0]) == 'factor_for_multiplying_ppb_to_get_molec/cm3':
			dlist.append(float(i))
			
	const[str(line.split(',')[0])] = dlist

num_sb = int((const['number_of_size_bins'])[0]) # number of size bins
num_speci = int((const['number_of_components'])[0]) # number of species
y_mw = const['molecular_weights_g/mol_corresponding_to_component_names']
y_MV = const['molecular_volumes_cm3/mol']
PyCHAM_names = const['component_names']
# conversion factor to change gas-phase concentrations from molecules/cc 
# (air) into ppb
Cfactor = float((const['factor_for_multiplying_ppb_to_get_molec/cm3'])[0])

# name of file where concentration (molecules/cc (air)) results saved
fname = str(Pyfname+'/concentrations_all_components_all_times_gas_particle_wall')
y = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header)
# convert gas-phase results from ppb to molecules/cm3 (air), for consistency with AtChem2
# results
y[:, 0:num_speci] = y[:, 0:num_speci]*Cfactor

# withdraw times
fname = str(Pyfname+'/time')
t_array2 = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header)
# add starting time of day (s)
t_array2 = t_array2+gconc[0, 0]

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
	frac_dev2[nz_ind, i-1] = ((y[nz_ind, ind]-
							gconc_int[nz_ind, i])/np.max(gconc_int[nz_ind, i]))*100.0



# ----------------------------------------------------------------------------------------
# make plot with all gas-phase concentrations shown
compnum = 0 # count on components

for i in comp_names[1::]: # loop through components (excluding time in column 0)
	if str(i)=="CH3O2" or str(i)=="HO2" or str(i)=="O" or str(i)=="CO" or str(i)=="NO2" or str(i)=="NO":
		compnum += 1
		continue
	ax1.plot(t_array2/3600.0, frac_dev2[:,compnum], label=str(i))
	compnum += 1
# plt.title('PyCHAM-AtChem2 Fractional Deviation of Gas-phase Concentration')
ax1.set_ylabel(r'Deviation (%)', fontsize=14)
ax1.yaxis.set_tick_params(size=14)
ax1.xaxis.set_tick_params(size=14)
ax1.legend(fontsize=12)
ax1.set_xlabel(r'Time of day (hours)', fontsize=14)
ax1.text(x=-1.8, y=0.18, s='(b)', size=14)
ax1.xaxis.set_tick_params(size=13)
ax1.yaxis.set_tick_params(size=13)
fig.savefig('fig03.png')
plt.show()

# ----------------------------------------------------------------------------------------
# comparison of reaction rates recorded by PyCHAM and AtChem2

# first, the high NOx case

# open AtChem2 loss reaction results
Atfname = str(cwd + '/photo_chem_data/AtChem2_APINENE/hiNOx/lossRates.output')

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
					
					
# open AtChem2 production reaction results
Atfname = str(cwd + '/photo_chem_data/AtChem2_APINENE/hiNOx/productionRates.output')

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


# open the PyCHAM tracking results
Pyfname = str(cwd + '/photo_chem_data/PyCHAM_APINENE/hiNOx/AtChem2_comp_hiNOx')

fname = str(Pyfname+'/HCHO_rate_of_change')
dydt = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header

fname = str(Pyfname+'/time')
Pyt = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header

PyCHAM_reac_num = dydt[0, :] # reaction numbers stored in PyCHAM results
Pyres = np.zeros((Atres.shape[0], Atres.shape[1]))

# for loss reactions
Pyres[:, 0] = Atres[:, 0] # times (s)
# reaction numbers, subtract one to account for Fortran indexing
Pyres[:, 1] = Atres[:, 1]-1

# PyCHAM indices for loss reactions, know that there are four loss reactions for HCHO
# in the MCM alpha-pinene ozonolysis scheme
Ati = 0 # count on AtChem  times

for i in range(len(Pyt)):
	
	if Pyt[i] == Atres[Ati*4, 0]:
		
		for ii in range(4):
			# -1 from Atres index to account for the Fortran indexing
			PyCHAM_ind = np.where(PyCHAM_reac_num==Atres[ii, 1]-1) 
			# record reaction rates
			Pyres[Ati*4+ii, 2] = dydt[i, PyCHAM_ind]
		Ati += 1
		
# repeat for production reactions
Pyres_prod = np.zeros((Atres_prod.shape[0], Atres_prod.shape[1]))		
Pyres_prod[:, 0] = Atres_prod[:, 0] # times (s)
# reaction numbers, subtract one to account for Fortran indexing
Pyres_prod[:, 1] = Atres_prod[:, 1]-1

# PyCHAM indices for production reactions, know that there are thirty five production 
# reactions for HCHO in the MCM alpha-pinene ozonolysis scheme
Ati = 0 # count on AtChem  times
for i in range(len(Pyt)):
	if Pyt[i] == Atres_prod[Ati*35, 0]:
		
		for ii in range(35):
			# -1 from Atres index to account for the Fortran indexing
			PyCHAM_ind = np.where(PyCHAM_reac_num==(Atres_prod[ii, 1]-1))

			# record reaction rates
			Pyres_prod[Ati*35+ii, 2] = dydt[i, PyCHAM_ind]
		Ati += 1

# ----------------------------------------------------------------------------------------
# comparative statistics

# number of loss reactions for this component
num_lossr = int((Atres[:, 0].shape[0])/(len(np.unique(Atres[:, 0]))))

# empty array for fractional deviation, first row for times, second for equation number
# and third for fractional deviation
frac_dev_loss = np.zeros((Atres[:, 0].shape[0], 1))
# loop through loss reactions
for i in range(num_lossr):
	# % difference, multiply AtChem result by -1 to represent loss
	frac_dev_loss[i::num_lossr, 0] = ((Pyres[i::num_lossr, 2]-(-1.0*Atres[i::num_lossr, 2]))/
									np.max(Atres[i::num_lossr, 2]))*100.0
									
# number of production reactions for this component
num_prodr = int((Atres_prod[:, 0].shape[0])/(len(np.unique(Atres_prod[:, 0]))))

# empty array for fractional deviation, first row for times, second for equation number
# and third for fractional deviation
frac_dev_prod = np.zeros((Atres_prod[:, 0].shape[0], 1))
# loop through loss reactions
for i in range(num_prodr):
	# % difference
	frac_dev_prod[i::num_prodr, 0] = ((Pyres_prod[i::num_prodr, 2]-(Atres_prod[i::num_prodr, 2]))/
									np.max(Atres_prod[i::num_prodr, 2]))*100.0

# identify the two loss and two production reactions with greatest deviation
max_dev1 = 0.0
max_dev0 = 0.0
loss_reac_i = np.zeros((2))
for i in range(num_lossr): # loop through loss reactions
	max_dev2 = np.max(np.abs(frac_dev_loss[i::num_lossr, 0]))
	if max_dev2>max_dev1 and max_dev2>max_dev0: # highest maximum seen so far
		max_dev1 = max_dev2
		loss_reac_i[1] = int(loss_reac_i[0]) # record second highest maximum
		loss_reac_i[0] = int(i) # record index of new highest maximum

# identify the two loss and two production reactions with greatest deviation
max_dev1 = 0.0
max_dev0 = 0.0
prod_reac_i = np.zeros((2))
for i in range(num_prodr): # loop through loss reactions
	max_dev2 = np.max(np.abs(frac_dev_prod[i::num_prodr, 0]))
	if max_dev2>max_dev1 and max_dev2>max_dev0: # highest maximum seen so far
		max_dev1 = max_dev2
		prod_reac_i[1] = int(prod_reac_i[0]) # record second highest maximum
		prod_reac_i[0] = int(i) # record index of new highest maximum


# plot results for hihg NOx case
fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10,6), sharex=True) # prepare plot
ax0.plot(Pyres[int(loss_reac_i[0])::num_lossr, 0]/3600.0, frac_dev_loss[int(loss_reac_i[0])::num_lossr, 0], linewidth=3, label=r'$\mathrm{OH+HCHO\rightarrow HO2+CO}$')
ax0.plot(Pyres[int(loss_reac_i[1])::num_lossr, 0]/3600.0, frac_dev_loss[int(loss_reac_i[1])::num_lossr, 0], linewidth=3, label=r'$\mathrm{NO3 + HCHO\rightarrow HNO3 + CO + HO2}$')
ax0.plot(Pyres_prod[int(prod_reac_i[0])::num_prodr, 0]/3600.0, frac_dev_prod[int(prod_reac_i[0])::num_prodr, 0], linewidth=3, label=r'$\mathrm{C621O\rightarrow HCHO + H1C23C4CHO + HO2}$')
ax0.plot(Pyres_prod[int(prod_reac_i[1])::num_prodr, 0]/3600.0, frac_dev_prod[int(prod_reac_i[1])::num_prodr, 0], linewidth=3, label=r'$\mathrm{C109O\rightarrow C89CO3 + HCHO}$')
ax0.text(x=-1.0, y=4.5, s='(a)', size=16)
ax0.set_ylabel(r'Deviation (%)', fontsize=16)
ax0.legend(fontsize=13, loc = [0.55, 0.5])
ax0.yaxis.set_tick_params(size=16)
ax0.xaxis.set_tick_params(size=16)
ax0.yaxis.set_tick_params(size=16)
ax0.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))

# ----------------------------------------------------------------------------------------
# Repeat for low NOx case

# open AtChem2 loss reaction results
Atfname = str(cwd + '/photo_chem_data/AtChem2_APINENE/loNOx/lossRates.output')

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
					
					
# open AtChem2 production reaction results
Atfname = str(cwd + '/photo_chem_data/AtChem2_APINENE/loNOx/productionRates.output')

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


# open the PyCHAM tracking results
Pyfname = str(cwd + '/photo_chem_data/PyCHAM_APINENE/loNOx/AtChem2_comp_loNOx')

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

# PyCHAM indices for loss reactions, know that there are four loss reactions for HCHO
# in the MCM alpha-pinene ozonolysis scheme
Ati = 0 # count on AtChem  times

for i in range(len(Pyt)):
	
	if Pyt[i] == Atres[Ati*4, 0]:
		
		for ii in range(4):
			# -1 from Atres index to account for the Fortran indexing
			PyCHAM_ind = np.where(PyCHAM_reac_num==Atres[ii, 1]-1) 
			# record reaction rates
			Pyres[Ati*4+ii, 2] = dydt[i, PyCHAM_ind]
		Ati += 1
		
# repeat for production reactions
Pyres_prod = np.zeros((Atres_prod.shape[0], Atres_prod.shape[1]))		
Pyres_prod[:, 0] = Atres_prod[:, 0] # times (s)
# reaction numbers, -1 to account for Fortran indexing
Pyres_prod[:, 1] = Atres_prod[:, 1]-1

# PyCHAM indices for production reactions, know that there are thirty five loss reactions
# for HCHO in the MCM alpha-pinene ozonolysis scheme
Ati = 0 # count on AtChem  times
for i in range(len(Pyt)):
	if Pyt[i] == Atres_prod[Ati*35, 0]:
		
		for ii in range(35):
			# -1 from Atres index to account for the Fortran indexing
			PyCHAM_ind = np.where(PyCHAM_reac_num==Atres_prod[ii, 1]-1) 
			# record reaction rates
			Pyres_prod[Ati*35+ii, 2] = dydt[i, PyCHAM_ind]
		Ati += 1

# ----------------------------------------------------------------------------------------
# comparative statistics

# number of loss reactions for this component
num_lossr = int((Atres[:, 0].shape[0])/(len(np.unique(Atres[:, 0]))))

# empty array for fractional deviation, first row for times, second for equation number
# and third for fractional deviation
frac_dev_loss = np.zeros((Atres[:, 0].shape[0], 1))
# loop through loss reactions
for i in range(num_lossr):
	
	# if maximum rate doesn't exceed 1 molecule/cc.s (air), then ignore
	if np.max(Atres[i::num_lossr, 2])<1.0:
		frac_dev_loss[i::num_lossr, 0] = 0.0
		continue
	
	# % difference, multiply AtChem result by -1 to represent loss
	frac_dev_loss[i::num_lossr, 0] = ((Pyres[i::num_lossr, 2]-(-1.0*Atres[i::num_lossr, 2]))/
									np.max(Atres[i::num_lossr, 2]))*100.0

							
# number of production reactions for this component
num_prodr = int((Atres_prod[:, 0].shape[0])/(len(np.unique(Atres_prod[:, 0]))))

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
max_dev1 = 0.0
max_dev0 = 0.0
loss_reac_i = np.zeros((2))
for i in range(num_lossr): # loop through loss reactions
	max_dev2 = np.max(np.abs(frac_dev_loss[i::num_lossr, 0]))
	if max_dev2>max_dev1 and max_dev2>max_dev0: # highest maximum seen so far
		max_dev1 = max_dev2
		loss_reac_i[1] = int(loss_reac_i[0]) # record second highest maximum
		loss_reac_i[0] = int(i) # record index of new highest maximum

# identify the two loss and two production reactions with greatest deviation
max_dev1 = 0.0
max_dev0 = 0.0
prod_reac_i = np.zeros((2))
for i in range(num_prodr): # loop through loss reactions
	max_dev2 = np.max(np.abs(frac_dev_prod[i::num_prodr, 0]))
	if max_dev2>max_dev1 and max_dev2>max_dev0: # highest maximum seen so far
		max_dev1 = max_dev2
		prod_reac_i[1] = int(prod_reac_i[0]) # record second highest maximum
		prod_reac_i[0] = int(i) # record index of new highest maximum

# print the reaction numbers where disagreement is greatest (Python indexing, starting at 
# 0)
# print('loNOx indices')
# print(Pyres[int(loss_reac_i[0]), 1])
# print(Pyres[int(loss_reac_i[1]), 1])
# print(Pyres_prod[int(prod_reac_i[0]), 1])
# print(Pyres_prod[int(prod_reac_i[1]), 1])

# plot results for hihg NOx case
ax1.plot(Pyres[int(loss_reac_i[0])::num_lossr, 0]/3600.0, frac_dev_loss[int(loss_reac_i[0])::num_lossr, 0], linewidth=3, label=r'$\mathrm{OH + HCHO\rightarrow HO2 + CO}$')
ax1.plot(Pyres[int(loss_reac_i[1])::num_lossr, 0]/3600.0, frac_dev_loss[int(loss_reac_i[1])::num_lossr, 0], linewidth=3, label=r'$\mathrm{HCHO\rightarrow H2 + CO}$')
ax1.plot(Pyres_prod[int(prod_reac_i[0])::num_prodr, 0]/3600.0, frac_dev_prod[int(prod_reac_i[0])::num_prodr, 0], linewidth=3, label=r'$\mathrm{HOCH2CO3\rightarrow HCHO + HO2}$')
ax1.plot(Pyres_prod[int(prod_reac_i[1])::num_prodr, 0]/3600.0, frac_dev_prod[int(prod_reac_i[1])::num_prodr, 0], linewidth=3, label=r'$\mathrm{CH3O2\rightarrow HCHO}$')


ax1.text(x=-1.0, y=0.8, s='(b)', size=16)
ax1.set_xlabel(r'Time of day (hours)', fontsize=16)
ax1.set_ylabel(r'Deviation (%)', fontsize=16)
ax1.legend(fontsize=13, loc = [0.70, -0.05])
ax1.yaxis.set_tick_params(size=16)
ax1.xaxis.set_tick_params(size=16)


fig.savefig('fig04.png')
plt.show()