'''module to test nuc.py module'''
print('module to test nuc.py')
print('importing required modules')
import os
import sys
dirpath = os.getcwd() # get current path
sys.path.append(os.path.split(dirpath)[0]) # add path to system path
import numpy as np
import scipy.constants as si
import ipdb
print('required modules imported okay, now importing nuc.py as in ode_gen.py')
from nuc import nuc
print('if no issue stated above, importing completed successfully')
print('now calling nuc.py as in ode_gen.py')

# preparing inputs
sumt = 400.0 # total time passed in simulation (s)
new_part_sum1 = 0.0 # new number of particles generated before this call 
N_perbin = np.array((10.0,12.0)) # number of particles per size bin (#/cc (air))
# number concentration of molecules (molecules/cc (air))
y = np.array((1.0e12, 1.0e12, 1.0e10, 1.0e16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
y_mw = (np.array((200.0, 48.0, 300.0, 18.0))).reshape(4,1) # molecular weight of components (g/mol)
y_dens = np.array((1.0e3, 1.0e3, 1.0e3, 1.0e3)) # liquid density of components (kg/m3)  
num_speci = 4 # number of components
x = (np.array((1.0e-2))).reshape(1) # current radius of particles in lowest size bin (um)
new_partr = 1.0e-7 # radius of two nucleating molecules together (cm)
t = 60.0 # time step (s)
MV = (np.array((200.0, 30.0, 300.0, 18.0))).reshape(4,1) # molar volume (cc/mol)
nucv1 = 31000.0
nucv2 = -55.0
nucv3 = 180.0
nuc_comp = 2 # index of the nucleating component
Varr = np.array((1.0e-7, 1.0e-5)) # volume of particles (um3)
[N_perbin, y, x[0], Varr[0], new_part_sum1] = nuc(sumt, new_part_sum1, 
							N_perbin, y, y_mw.reshape(-1, 1), np.squeeze(y_dens*1.0e-3),  
							num_speci, x[0], new_partr, t, MV, nucv1, nucv2, 
							nucv3, nuc_comp)

if int(N_perbin[0])!=89:
	print('issue with N_perbin')
if int(y[2]*1e-7)!=999 or int(y[6])!=672:
	print('issue with y')
if int(x[0]*1e5)!=96:
	print('issue with x')
if int(Varr[0]*1e11)!=372:
	print('issue with Varr')
if int(new_part_sum1)!=79:
	print('issue with new_part_sum1')

print('if no issues stated above, then code is working fine')
