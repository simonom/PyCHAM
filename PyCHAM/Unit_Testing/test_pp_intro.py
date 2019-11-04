# function to test pp_intro.py
print('function to test pp_intro.py and Size_distributions.py and init_water_partit.py')

import os
import sys
dirpath = os.getcwd() # get current path
sys.path.append(os.path.split(dirpath)[0]) # add path to system path
print('testing imports')
import numpy as np
import Size_distributions # custom library - see source code
from init_water_partit import init_water_partit
import scipy.constants as si
print('imports imported fine')

print('now calling pp_intro.py as in front.py')
from pp_intro import pp_intro

y = np.array((2.46e+11, 2.95e+11, 1.00e-40, 4.62e+17, 1.00e-40))
num_speci = 5
spec_list = ['CC1=CCC2CC1C2(C)C','[OH]','[O]OC(C)(C)C1CC=C(C)C(O)C1']
import pybel
Pybel_objects = []
for name_only in spec_list:
	Pybel_object = pybel.readstring('smi', name_only)
	# append to Pybel object list
	Pybel_objects.append(Pybel_object)
TEMP = 298.15
H2Oi=3
Psat_water = -1.5
mfp = (np.array((2.49e-8, 2.56e-7, 1.61e-8, 1.17e-7, 2.45e-8))).reshape(5,1)
y_mw = (np.array((130.0, 2.0, 200.0, 18.0, 132.14))).reshape(5,1)
surfT = 72.0
DStar_org = (np.array((7.40e-06, 1.20e-04, 5.56e-06, 2.77e-05, 7.32e-06))).reshape(5,1)
RH = 0.60
testf = 0
num_sb = 2
lowersize = 0.0
uppersize = 8.0e-2
pconc = 1.0e3
tmax = 60.0
nuc_comp = 2
voli = np.array((2, -1))
volP = np.array((1.0e-30, 1.0e-30))
testf = 2
std = 1.0
loc = 0.0
scale = 1.0e-2
kgwt = 332.0e-25
accom_coeff = np.ones((num_speci,1)) # particle accommodation coefficient
therm_sp = (np.array((0,0,0,218,0))).reshape(num_speci,1)
Cw = 301.0e13
y_dens = np.array((870, 1133, 1033, 1000, 1770)).reshape(num_speci,1)
Psat = (np.array((110*1e15, 526*1e20, 242/1e18, 778*1e15, 242/1e18))).reshape(num_speci,1)
core_diss = 3.0

[y, N_perbin, x, Varr, Vbou, 
							rad0, Vol0, rbou, 
							new_partr, MV, num_sb, nuc_comp] = pp_intro(y, 
							num_speci, spec_list, Pybel_objects, TEMP, H2Oi, 
							mfp, accom_coeff, y_mw, surfT, DStar_org, 
							RH, num_sb, lowersize, uppersize, pconc, tmax, nuc_comp, 
							voli, volP, testf, std, loc, scale,
							therm_sp, Cw, y_dens, Psat, core_diss, kgwt)
print('pp_intro called and returned fine, now checking returned values')

if int(y[5])!=0 or int(y[8]*1e-6)!=717 or int(y[9]*1e-6)!=252 or int(y[13]*1e-7)!=179 or int(y[14]*1e-6)!=476 or int(y[18]*1e-13)!=178:
	print('issue with y, possibly due to init_water_partit.py or Size_distributions.py')  
if int(N_perbin[0])!=934 or int(N_perbin[1])!=65:
	print('issue with N_perbin possibly due to Size_distributions')
if int(x[0]*1e4)!=237 or int(x[1]*1e4)!=744 :
	print('issue with x possibly due to Size_distributions')
if int(Varr[0]*1e7)!=564 or int(Varr[1]*1e5)!=172:
	print('issue with Varr possibly due to Size_distributions')
if int(Vbou[0])!=0 or int(Vbou[1]*1e6)!=268 or int(Vbou[2]*1e-2)!=214:
	print('issue with Vbou possibly due to Size_distributions')
if int(rad0[0]*1e2)!=2 or int(rad0[1]*1e2)!=6:
	print('issue with rad0 possibly due to Size_distributions')
if int(Vol0[0]*1e7)!=335 or int(Vol0[1]*1e6)!=904:
	print('issue with Vol0 possibly due to Size_distributions')
if int(rbou[0])!=0 or int(rbou[1]*1e4)!=400 or int(rbou[2])!=8:
	print('issue with rbou possibly due to Size_distributions')
if int(new_partr*1e10)!=535:
	print('issue with new_partr')
if int(MV[0])!=149 or int(MV[1]*1e2)!=176 or int(MV[2])!=193 or int(MV[3]*1e1)!=180 or int(MV[4]*1e1)!=746:   
	print('issue with MV')
if num_sb!=3:
	print('issue with num_sb')
print('testing finished, if no issues printed above, pp_intro.py working fine')