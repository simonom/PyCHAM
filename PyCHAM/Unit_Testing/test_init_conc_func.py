# module to test init_conc_func, water_calc

print('tesint init_conc_func.py and eqn_parser.write_rate_file and water_calc.py and MCM_constants_auto.mcm_constants')
print('testing imports')

import os
import sys
dirpath = os.getcwd() # get current path
sys.path.append(os.path.split(dirpath)[0]) # add path for eqn_parser to system path

import numpy as np
from water_calc import water_calc
import MCM_const
import eqn_parser
import scipy.constants as si

print('init_conc_func imports imported fine, now importing init_conc_func')
from init_conc_func import init_conc_func

# test inputs for init_conc_func
num_speci = 3
init_SMIL = ['CC1=CCC2CC1C2(C)C','[OH]']
spec_list = ['CC1=CCC2CC1C2(C)C','[OH]','[O]OC(C)(C)C1CC=C(C)C(O)C1']
init_conc = np.array((10.0, 12.0))
TEMP = 298.15
RH = 0.60
NA = si.N_A
PInit = 101325.0
M_val = ((PInit/(8.3144598*TEMP)*6.0221409e+23)*1.0e-6)
N2_val = M_val*0.79
O2_val = M_val*0.2096
reac_coef = ['1.2e-11*numpy.exp(440/TEMP)*0.055']
fname = 'scheme_test'
start_sim_time = 11
lat = 12
lon = 13
import pybel
Pybel_objects = []
for name_only in spec_list:
	Pybel_object = pybel.readstring('smi', name_only)
	# append to Pybel object list
	Pybel_objects.append(Pybel_object)
testf = 2
pconc = 1000.0
print('submitting inputs to init_conc_func')

[y, H2Oi, Psat_water, y_mw, num_speci, 
	Cfactor, y_indx_plot, corei] = init_conc_func(num_speci, init_SMIL, spec_list, init_conc, 
							TEMP, RH, M_val, N2_val, O2_val, reac_coef, fname, PInit, 
							start_sim_time, lat, lon, Pybel_objects, testf, pconc)
							
print('checking init_conc_func outputs')				
if y[0]!=10.0*PInit*(NA/(8.3144598e6*TEMP))*1.0e-9:
	print('issue with y[0]')
if y[1]!=12.0*PInit*(NA/(8.3144598e6*TEMP))*1.0e-9:
	print('issue with y[1]')
if int(y[3]*1e-15)!=461:
	print('issue with y[3], possibly due to water_calc')
if H2Oi!=3:
	print('issue with H2Oi')
if int(Psat_water*100)!= -150:
	print('issue with Psat_water, possibly due to water_calc')
if int(y_mw[0])!=136 or int(y_mw[1])!=17 or int(y_mw[2])!=185 or int(y_mw[3])!=18 or int(y_mw[4])!=132:
	print('issue with y_mw')
if num_speci!=5:
	print(num_speci)
	print('issue with num_speci')
if Cfactor!=PInit*(NA/(8.3144598e6*TEMP))*1.0e-9:
	print('issue with Cfactor')
if y_indx_plot[0]!= 0 and y_indx_plot[1]!=1:
	print('issue with y_indx_plot')
if corei!=num_speci-1:
	print('issue with corei')

os.remove('Rate_coeffs.py') # remove rate coefficient file

print('if no issues stated above, the init_conc_func.py code and its called on functions (listed above) are fine')