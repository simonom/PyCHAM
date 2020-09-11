'''code for plotting PyCHAM results'''
# this particular code used for analysing simulations investigating
# methane atmospheric lifetime in real-world troposphere

import numpy as np
import os
import matplotlib.pyplot as plt
from retr_out import retr_out as retr

# get current working directory

dir_path = os.getcwd() # current working directory
# obtain just part of the path up to PyCHAM home directory
for i in range(len(dir_path)):
	if dir_path[i:i+7] == 'PyCHAM':
		dir_path = dir_path[0:i+7]
		break
dir_path = str(dir_path+'/PyCHAM/output/scheme_simp/woCL')

# chamber condition ---------------------------------------------------------

Cfac_dict = {}

# retrieve results
(_, _, Cfac_dict['Cfac0'], _, _, _, _, _, _, _, _, _, _, _) = retr(dir_path)
Cfac = (Cfac_dict['Cfac0'])[0]

# concentrations ------------------------------------------------------------
# open results
res = open(str(dir_path+'/concentrations_all_components_all_times_gas_particle_wall'))

# separate into lines
resn = res.readlines()
res.close() # close file

# get component chemical scheme names and strip white space
comp = ((resn)[1]).split(',')
# remove excess punctuation
for i in range(len(comp)):
	comp[i] = comp[i].strip('#')
	comp[i] = comp[i].strip(' ')

# get concentrations (ppb for gas-phase), times in rows, components 
# in columns
y = np.zeros((len(resn)-2, len(comp)))
for i in range(2, len(resn)):
	y[i-2, :] = (resn[i]).split(',')

# tendencies to change ------------------------------------------------------

# open results
res = open(str(dir_path+'/CH4_rate_of_change'))

# separate into lines
resn = res.readlines()
res.close() # close file

# get reaction numbers (pythonic indexing) and strip white space
reacn = ((resn)[1]).split(',')
# remove excess punctuation
for i in range(len(reacn)):
	reacn[i] = reacn[i].strip('#')
	reacn[i] = reacn[i].strip(' ')

# get tendencies to change (molecules/cc/s), times in rows, reactions 
# in columns
dy = np.zeros((len(resn)-2, len(reacn)))
for i in range(2, len(resn)):
	dy[i-2, :] = (resn[i]).split(',')

# ----------------------------------------------------------------------------------------
# estimate atmospheric lifetime of methane
indx = comp.index('CH4_g') # get index
#Cfac = 24130931925.639053
tindx = 720
# average loss rate due to Cl
CL_loss = sum(dy[tindx::, 0])/len(dy[tindx::, 0])
# average loss rate due to OH
OH_loss = sum(dy[tindx::, 1])/len(dy[tindx::, 1])
print(CL_loss, OH_loss)
# lifetime in s
life = y[0, indx]*Cfac/(CL_loss+OH_loss)
# lifetime in years
life = life/(3600.*24.*365.25)
print('atmospheric lifetime: ', life, ' years')
# names of components we want to plot
names = ['CH4_g', 'O3_g', 'OH_g', 'APINENE_g', 'CL_g']

for i in names: # loop through components
	indx = comp.index(i) # get index
	plt.semilogy(y[:, indx], label = i)

plt.ylabel('[gas-phase] (ppb)')
plt.legend()
plt.show()
