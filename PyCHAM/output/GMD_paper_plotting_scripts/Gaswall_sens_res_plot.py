'''Code to plot results for gas-wall partitioning sensitivity'''
# aim is to illustrate the effect of varying the mass transfer coefficient for gas-wall
# partitioning and the effective absorbing mass concentration of the wall on SOA mass
# concentration estimates

# all runs for plot (a) used the isoprene MCM sheme (isoprene_scheme_MCM.txt), whilst 
# plot (b) used Gaswall_simple_chem.txt
# PyCHAM inputs are saved as: Gaswall_sens_inputs.txt and Gaswall_sens_inputs_b.txt
#  

# import required modules
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as si
import os

# ----------------------------------------------------------------------------------------
#import result files

# first case is no gas-wall partitioning, where kgwt was set to 
# 1.0e-20 /s and Cw was set to 1.0e-10 g/m3 (air)

# get current working directory
cwd = os.getcwd()

# open saved files
output_by_sim = str(cwd + '/Gas_wall_partit_data/Gaswall_sens_Cw7e-20_kw1e-20')

# name of file where experiment constants saved (number of size bins and whether wall 
# included)
fname = str(output_by_sim+'/constants')

const_in = open(fname)
const = {} # prepare to create dictionary
for line in const_in.readlines():

	# convert to python list
	dlist = []
	for i in line.split(',')[1::]:
		if str(line.split(',')[0]) == 'num_sb':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'num_speci':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'mw' or  str(line.split(',')[0]) == 'mv':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			dlist.append(float(i))
		if str(line.split(',')[0]) == 'spec_namelist':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			i = i.strip('\'')
			dlist.append(str(i))
	const[str(line.split(',')[0])] = dlist

num_sb = (const['num_sb'])[0] # number of size bins
num_speci = int((const['num_speci'])[0]) # number of species
y_mw = const['mw']
y_MV = const['mv']
PyCHAM_names = const['spec_namelist']

# name of file where concentration (molecules/cc (air)) results saved
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

SOA0 = 0.0
for i in range(1, num_sb):
	# calculate SOA (*1.0E-12 to convert from g/cc (air) to ug/m3 (air))
	SOA0 += (((y[:,num_speci*i:num_speci*(i+1)-2]/si.N_A)*y_mw[0:-2]*1.0e12).sum(axis=1))
			
# ----------------------------------------------------------------------------------------
# second case is moderate gas-wall partitioning, where kgwt was set to 
# 1.0e-1 /s and Cw was set to 7.0e-5 g/m3 (air), as these values are comparable to those
# for particulate matter (for particulate matter, total concentration (core, water and
# SOA summed) was around 7.0e-4 g/m3 (air))

# open saved files
output_by_sim = str(cwd + '/Gas_wall_partit_data/Gaswall_sens_Cw7e-5_kw1e-1')

# name of file where experiment constants saved (number of size bins and whether wall 
# included)
fname = str(output_by_sim+'/constants')

const_in = open(fname)
const = {} # prepare to create dictionary
for line in const_in.readlines():

	# convert to python list
	dlist = []
	for i in line.split(',')[1::]:
		if str(line.split(',')[0]) == 'num_sb':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'num_speci':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'mw' or  str(line.split(',')[0]) == 'mv':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			dlist.append(float(i))
		if str(line.split(',')[0]) == 'spec_namelist':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			i = i.strip('\'')
			dlist.append(str(i))
	const[str(line.split(',')[0])] = dlist

num_sb = (const['num_sb'])[0] # number of size bins
num_speci = int((const['num_speci'])[0]) # number of species
y_mw = const['mw']
y_MV = const['mv']
PyCHAM_names = const['spec_namelist']

# name of file where concentration (molecules/cc (air)) results saved
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

SOA1 = 0.0
for i in range(1, num_sb):
	# calculate SOA (*1.0E-12 to convert from g/cc (air) to ug/m3 (air))
	SOA1 += (((y[:,num_speci*i:num_speci*(i+1)-2]/si.N_A)*y_mw[0:-2]*1.0e12).sum(axis=1))

# ----------------------------------------------------------------------------------------
# third case is high gas-wall partitioning, where kgwt was set to 
# 1.0e-1 /s and Cw was set to 7.0e-2 g/m3 (air)

# open saved files
output_by_sim = str(cwd + '/Gas_wall_partit_data/Gaswall_sens_Cw7e-2_kw1e-1')

# name of file where experiment constants saved (number of size bins and whether wall 
# included)
fname = str(output_by_sim+'/constants')

const_in = open(fname)
const = {} # prepare to create dictionary
for line in const_in.readlines():

	# convert to python list
	dlist = []
	for i in line.split(',')[1::]:
		if str(line.split(',')[0]) == 'num_sb':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'num_speci':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'mw' or  str(line.split(',')[0]) == 'mv':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			dlist.append(float(i))
		if str(line.split(',')[0]) == 'spec_namelist':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			i = i.strip('\'')
			dlist.append(str(i))
	const[str(line.split(',')[0])] = dlist

num_sb = (const['num_sb'])[0] # number of size bins
num_speci = int((const['num_speci'])[0]) # number of species
y_mw = const['mw']
y_MV = const['mv']
PyCHAM_names = const['spec_namelist']

# name of file where concentration (molecules/cc (air)) results saved
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

SOA2 = 0.0
for i in range(1, num_sb):
	# calculate SOA (*1.0E-12 to convert from g/cc (air) to ug/m3 (air))
	SOA2 += (((y[:,num_speci*i:num_speci*(i+1)-2]/si.N_A)*y_mw[0:-2]*1.0e12).sum(axis=1))

# ----------------------------------------------------------------------------------------
# fourth case is high gas-wall partitioning, where kgwt was set to 
# 1.0e2 /s and Cw was set to 7.0e-5 g/m3 (air)

# open saved files
output_by_sim = str(cwd + '/Gas_wall_partit_data/Gaswall_sens_Cw7e-5_kw1e2')

# name of file where experiment constants saved (number of size bins and whether wall 
# included)
fname = str(output_by_sim+'/constants')

const_in = open(fname)
const = {} # prepare to create dictionary
for line in const_in.readlines():

	# convert to python list
	dlist = []
	for i in line.split(',')[1::]:
		if str(line.split(',')[0]) == 'num_sb':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'num_speci':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'mw' or  str(line.split(',')[0]) == 'mv':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			dlist.append(float(i))
		if str(line.split(',')[0]) == 'spec_namelist':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			i = i.strip('\'')
			dlist.append(str(i))
	const[str(line.split(',')[0])] = dlist

num_sb = (const['num_sb'])[0] # number of size bins
num_speci = int((const['num_speci'])[0]) # number of species
y_mw = const['mw']
y_MV = const['mv']
PyCHAM_names = const['spec_namelist']

# name of file where concentration (molecules/cc (air)) results saved
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

SOA3 = 0.0
for i in range(1, num_sb):
	# calculate SOA (*1.0E-12 to convert from g/cc (air) to ug/m3 (air))
	SOA3 += (((y[:,num_speci*i:num_speci*(i+1)-2]/si.N_A)*y_mw[0:-2]*1.0e12).sum(axis=1))

# ----------------------------------------------------------------------------------------
# prepare plot
fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(8,5))

ax0.plot(t_array/3600.0, SOA0, linewidth=3, label='No wall')
ax0.plot(t_array/3600.0, SOA1, linewidth=3, label=r'$k_{w}=1x10^{-1}\, \mathrm{s^{-1}}$, $C_{w}=7x10^{1}\, \mathrm{\mu g \, m^{-3}}$')
ax0.plot(t_array/3600.0, SOA2, linewidth=3, label=r'$k_{w}=1x10^{-1}\, \mathrm{s^{-1}}$, $C_{w}=7x10^{4}\, \mathrm{\mu g \, m^{-3}}$')
ax0.plot(t_array/3600.0, SOA3, '--', linewidth=3, label=r'$k_{w}=1x10^{2}\, \mathrm{s^{-1}}$, $C_{w}=7x10^{1}\, \mathrm{\mu g \, m^{-3}}$')
ax0.set_xlabel(r'Time (hours of day)', fontsize=12)
ax0.set_ylabel(r'[SOA] ($\mathrm{\mu g\, m^{-3}}$)', fontsize=12)
ax0.yaxis.set_tick_params(size=12)
ax0.xaxis.set_tick_params(size=12)
ax0.legend(fontsize=10)
ax0.text(x=-1.1, y=31.7, s='(a)', size=12)


# ----------------------------------------------------------------------------------------
# for plot b, we show decay of individual components over time in a dark experiment with
# no lights or particles and no chemical interaction between the components of interest
# in the results below the Gaswall_simple_chem chemical scheme was used, and components 
# with an initial concentration of 50 ppb were (MCM names): MGLYOX (methylglyoxal) and 
# two-methylglyceric_acid, the latter had to be added manually to the xml file, but its
# production during isoprene ozonolysis is confirmed here: doi.org/10.1021/jp061734m
# we want to plot concentrations of all components against time to illustrate how to 
# constrain wall loss of vapour parameters 


# ----------------------------------------------------------------------------------------
# moderate Kw (1e-1 /s) and moderate Cw (7e-5 ug/m3)
# open saved files
output_by_sim = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v120/PyCHAM/output/GMD_paper/Results/Gas_wall_partit_data/Gaswall_sens_Cw7e-5_kw1e-1_2comp'

# name of file where experiment constants saved (number of size bins and whether wall 
# included)
fname = str(output_by_sim+'/constants')

const_in = open(fname)
const = {} # prepare to create dictionary
for line in const_in.readlines():

	# convert to python list
	dlist = []
	for i in line.split(',')[1::]:
		if str(line.split(',')[0]) == 'num_sb':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'num_speci':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'mw' or  str(line.split(',')[0]) == 'mv':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			dlist.append(float(i))
		if str(line.split(',')[0]) == 'spec_namelist':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			i = i.strip('\'')
			dlist.append(str(i))
		if str(line.split(',')[0]) == 'Cfactor':
			dlist.append(float(i))
	const[str(line.split(',')[0])] = dlist

num_sb = const['num_sb'] # number of size bins
num_speci = int((const['num_speci'])[0]) # number of species
y_mw = const['mw']
y_MV = const['mv']
PyCHAM_names = const['spec_namelist']
comp_index = [] # empty index list

# conversion factor to change gas-phase concentrations from molecules/cc 
# (air) into ppb
Cfactor = float((const['Cfactor'])[0])

# get index of three components of interest
for i in range(len(PyCHAM_names)):
	if PyCHAM_names[i] == 'MGLYOX':
		comp_index.append(i)
	if PyCHAM_names[i] == 'two-methylglyceric_acid':
		comp_index.append(i)
	
# name of file where concentration (molecules/cc (air)) results saved
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# ax1.plot(t_array[0:45], y[0:45, comp_index[0]]/Cfactor, '+r', label= r''+str('methylglyoxal')+' $k_{w}=1x10^{-1}\, \mathrm{s^{-1}}, C_{w}=7x10^{-5}\, \mathrm{g \, m^{-3}}$')
ax1.plot(t_array[0:45], y[0:45, comp_index[1]]/Cfactor, '+g', label= r''+' $k_{w}=1x10^{-1}\, \mathrm{s^{-1}}, C_{w}=7x10^{1}\, \mathrm{\mu g \, m^{-3}}$')

# ----------------------------------------------------------------------------------------
# moderate Kw (1e-1 /s) high Cw (7e-2 ug/m3)
# open saved files
output_by_sim = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v120/PyCHAM/output/GMD_paper/Results/Gas_wall_partit_data/Gaswall_sens_Cw7e-2_kw1e-1_2comp'

# name of file where experiment constants saved (number of size bins and whether wall 
# included)
fname = str(output_by_sim+'/constants')

const_in = open(fname)
const = {} # prepare to create dictionary
for line in const_in.readlines():

	# convert to python list
	dlist = []
	for i in line.split(',')[1::]:
		if str(line.split(',')[0]) == 'num_sb':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'num_speci':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'mw' or  str(line.split(',')[0]) == 'mv':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			dlist.append(float(i))
		if str(line.split(',')[0]) == 'spec_namelist':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			i = i.strip('\'')
			dlist.append(str(i))
		if str(line.split(',')[0]) == 'Cfactor':
			dlist.append(float(i))
	const[str(line.split(',')[0])] = dlist

num_sb = const['num_sb'] # number of size bins
num_speci = int((const['num_speci'])[0]) # number of species
y_mw = const['mw']
y_MV = const['mv']
PyCHAM_names = const['spec_namelist']
comp_index = [] # empty index list

# conversion factor to change gas-phase concentrations from molecules/cc 
# (air) into ppb
Cfactor = float((const['Cfactor'])[0])

# get index of three components of interest
for i in range(len(PyCHAM_names)):
	if PyCHAM_names[i] == 'MGLYOX':
		comp_index.append(i)
	if PyCHAM_names[i] == 'two-methylglyceric_acid':
		comp_index.append(i)
	
# name of file where concentration (molecules/cc (air)) results saved
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# ax1.plot(t_array[0:45], y[0:45, comp_index[0]]/Cfactor, '--r', label= r''+str('methylglyoxal')+' $k_{w}=1x10^{-1}\, \mathrm{s^{-1}}, C_{w}=7x10^{-2}\, \mathrm{g \, m^{-3}}$')
ax1.plot(t_array[0:45], y[0:45, comp_index[1]]/Cfactor, '--g', label= r''+' $k_{w}=1x10^{-1}\, \mathrm{s^{-1}}, C_{w}=7x10^{4}\, \mathrm{\mu g \, m^{-3}}$')


# ----------------------------------------------------------------------------------------
# fourth case is Kw (1e-2 /s) moderate Cw (7e-5 ug/m3)
# open saved files
output_by_sim = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v120/PyCHAM/output/GMD_paper/Results/Gas_wall_partit_data/Gaswall_sens_Cw7e-5_kw1e-2_2comp'

# name of file where experiment constants saved (number of size bins and whether wall 
# included)
fname = str(output_by_sim+'/constants')

const_in = open(fname)
const = {} # prepare to create dictionary
for line in const_in.readlines():

	# convert to python list
	dlist = []
	for i in line.split(',')[1::]:
		if str(line.split(',')[0]) == 'num_sb':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'num_speci':
			dlist.append(int(i))
		if str(line.split(',')[0]) == 'mw' or  str(line.split(',')[0]) == 'mv':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			dlist.append(float(i))
		if str(line.split(',')[0]) == 'spec_namelist':
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			i = i.strip('\'')
			dlist.append(str(i))
		if str(line.split(',')[0]) == 'Cfactor':
			dlist.append(float(i))
	const[str(line.split(',')[0])] = dlist

num_sb = const['num_sb'] # number of size bins
num_speci = int((const['num_speci'])[0]) # number of species
y_mw = const['mw']
y_MV = const['mv']
PyCHAM_names = const['spec_namelist']
comp_index = [] # empty index list

# conversion factor to change gas-phase concentrations from molecules/cc 
# (air) into ppb
Cfactor = float((const['Cfactor'])[0])

# get index of three components of interest
for i in range(len(PyCHAM_names)):
	if PyCHAM_names[i] == 'MGLYOX':
		comp_index.append(i)
	if PyCHAM_names[i] == 'two-methylglyceric_acid':
		comp_index.append(i)
	
# name of file where concentration (molecules/cc (air)) results saved
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# ax1.plot(t_array[0:45], y[0:45, comp_index[0]]/Cfactor, '.r', label= r''+str('methylglyoxal')+' $k_{w}=1x10^{-2}\, \mathrm{s^{-1}}, C_{w}=7x10^{-5}\, \mathrm{g \, m^{-3}}$')
ax1.plot(t_array[0:45], y[0:45, comp_index[1]]/Cfactor, '.g', label= r''+' $k_{w}=1x10^{-2}\, \mathrm{s^{-1}}, C_{w}=7x10^{1}\, \mathrm{\mu g \, m^{-3}}$')

ax1.set_xlabel(r'Time (seconds into simulation)', fontsize=12)
ax1.set_ylabel(r'[2-methylglyceric acid] (ppb)', fontsize=12)
ax1.yaxis.set_tick_params(size=12)
ax1.xaxis.set_tick_params(size=12)
ax1.legend(fontsize=10)
ax1.text(x=-4.0, y=53.0, s='(b)', size=12)
plt.subplots_adjust(hspace=0.34)
plt.show()
fig.savefig('fig06.png')
 