'''code to plot results from AtChem2 and compare against PyCHAM results'''
# introduction
# aim is to calculate deviation of PyCHAM photochemistry output from AtChem2
# use the AtChem2_apinene_scheme.txt in the Results folder of the GMD_paper for chemical 
# scheme input
# use the xml file from the PyCHAM inputs folder
# use the fig03_mod_var_hiNOx.txt and fig03_mod_var_loNOx.txt of the 
# GMD_paper_ploting_scripts folder, for the high and low NOx simulations, respectively, where for 
# the low NOx case, initial concentration of NO2 is set to 0 ppb whilst for the high NOx 
# case it's set to 9.8 ppb
# important that exactly the same chemical scheme as used in the AtChem2 simulation is
# used, this saved as fig03_scheme.txt in the same folder as this file
# for the simulations to compare with AtChem2 use:
# update_step = 60.
# recording_time_step = 60.
# whereas, for the temporal resolution sensitivity simulations use the appropriate
# resolution for update_step and set recording_time_step to be the same

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

# can call from the PyCHAM home folder or the GMD paper Results folder

import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os

# get current working directory
cwd = os.getcwd()
# ensure modules can be seen 
sys.path.append(str(os.getcwd() + '/PyCHAM'))
import retr_out

# ----------------------------------------------------------------------------------------
# AtChem2 part
# open saved files
try: # in case calling from the PyCHAM folder
	Atfname = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig03_data/AtChem2_APINENE/hiNOx/speciesConcentrations.output')
	inputs = open(Atfname, mode='r') # open results
except: # in case calling from GMD paper/Results folder
	Atfname = str(cwd + '/fig03_data/AtChem2_APINENE/hiNOx/speciesConcentrations.output')
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
try: # in case calling from the PyCHAM folder
	Pyfname = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig03_data/PyCHAM_APINENE/hiNOx/PyCHAM_comp_hiNOx')
	# required outputs
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, xfm, t_array, PyCHAM_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(Pyfname)
except: # in case calling from the GMD paper results folder
	Pyfname = str(cwd + '/fig03_data/PyCHAM_APINENE/hiNOx/PyCHAM_comp_hiNOx')
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
elif (np.abs(sum(gconc[:, 0]-t_array)) > 1.e-3):
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


#-----------------------------------------------------------------------------------------
# AtChem2 part
# open saved files
try: # in case calling from the PyCHAM folder
	Atfname = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig03_data/AtChem2_APINENE/loNOx/speciesConcentrations.output')
	inputs = open(Atfname, mode='r') # open results
except: # in case calling from the GMD paper Results folder
	Atfname = str(cwd + '/fig03_data/AtChem2_APINENE/loNOx/speciesConcentrations.output')
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
try: # in case calling from the PyCHAM folder
	Pyfname = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig03_data/PyCHAM_APINENE/loNOx/PyCHAM_comp_loNOx')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array2, PyCHAM_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(Pyfname)
except: # in case calling from the GMD paper Results folder
	Pyfname = str(cwd + '/fig03_data/PyCHAM_APINENE/loNOx/PyCHAM_comp_loNOx')
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
