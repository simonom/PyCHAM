'''module to unit test kimt_prep.py, which is used in preparation of the gas-particle and gas-wall partitioning odes'''

# function to test kimt_prep.py
print('function to test kimt_prep.py')
import os
import sys
dirpath = os.getcwd() # get current path
sys.path.append(os.path.split(dirpath)[0]) # add path to system path

print('testing imports')

import numpy as np
import scipy.constants as si
from kimt_prep import kimt_prep
print('imports fine')
print('creating inputs to kimt_prep')

y_mw = np.array((130.0, 2.0, 200.0, 18.0, 132.14))
TEMP = 298.15
num_speci = 5
testf = 0
Cw = 1.0
kgwt = 1.0e5
print('calling kimt_prep')

[DStar_org, mfp, accom_coeff, therm_sp, surfT, Cw, kgwt] = kimt_prep(y_mw, TEMP, num_speci, 
																testf, Cw, kgwt)
print('checking outputs')
if int(DStar_org[0]*1e8)!=740 or int(DStar_org[1]*1e6)!=119 or int(DStar_org[2]*1e8)!=555 or int(DStar_org[3]*1e7)!=276 or int(DStar_org[4]*1e8)!=732:
	print('issue with Dstar_org')
if int(mfp[0]*1e10)!=249 or int(mfp[1]*1e9)!=256 or int(mfp[2]*1e10)!=161 or int(mfp[3]*1e9)!=117 or int(mfp[4]*1e10)!=245:
	print('issue with mfp')
if (accom_coeff[0])!=1. or (accom_coeff[2])!=1. or (accom_coeff[3])!=1. or (accom_coeff[4])!=1.:
	print('issue with accom_coeff')
if int(therm_sp[0])!=220 or int(therm_sp[1])!=1776 or int(therm_sp[2])!=177 or int(therm_sp[3])!=592 or int(therm_sp[4])!=218: 
	print('issue with therm_sp')
if surfT!=72.0:
	print('issue with surfT')
if int(Cw*1e-13)!= 301:
	print('issue with Cw')
if int(kgwt*1e13)!= 332:
	print('issue with kgwt')

print('if no issues printed above, code is fine')