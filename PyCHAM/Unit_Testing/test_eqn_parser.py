# module to test eqn_parser.py
import os
import sys
dirpath = os.getcwd() # get current path
sys.path.append(os.path.split(dirpath)[0]) # add path for eqn_parser to system path

print('checking imports for eqn_parser.py')

import re
import collections

import datetime
import formatting
import numpy as np
import pybel
import xmltodict # for opening and converting xml files to python dictionaries
import ipdb

print('imports okay')

import eqn_parser # module to test

# eqn_parser.extract_mechanism test inputs
fname = 'scheme_test.txt'
xmlname = 'xml_test.txt'
TEMP = 298.15
PInit = 101325
Comp0 = ['APINENE']
testf = 0
print('calling on eqn_parser.extract_mechanism, eqn_parser.write_RO2_indices and formatting with test inputs')
[rindx, pindx, rstoi, pstoi, reac_coef, spec_list, Pybel_objects, num_eqn, num_speci, 
		RO2_indices, nreac, nprod, prodn, 
		reacn, M_val, N2_val, O2_val, 
		init_SMIL] = eqn_parser.extract_mechanism(fname, xmlname, 
							TEMP, PInit, Comp0, testf)
print('successfully called, now checking outputs')

if rindx[0,0]!= 0 or rindx[0,1]!= 1:
	print('issue with rindx')
if pindx!=2:
	print('issue with pindx')
if rstoi[0,0]!=1.0 or rstoi[0,1]!=1.0:
	print('issue with rstoi')
if pstoi!=1:
	print('issue with pstoi')
if reac_coef[0]!='1.2e-11*numpy.exp(440/TEMP)*0.055':
	print('issue with reac_coef')
if spec_list[0]!='CC1=CCC2CC1C2(C)C':
	print('issue with spec_list')
if spec_list[1]!='[OH]':
	print('issue with spec_list')
if spec_list[2]!='[O]OC(C)(C)C1CC=C(C)C(O)C1':
	print('issue with spec_list')
if num_eqn!=1:
	print('issue with num_eqn')
if num_speci!=3:
	print('issue with num_speci')
if RO2_indices[0,0]!=4 or RO2_indices[0,1]!=2:
	print('issue with RO2_indices, possibly due to write_RO2_indices function inside eqn_parser')
if nreac[0]!=2:
	print('issue with nreac')
if nprod[0]!=1:
	print('issue with nprod')
if prodn!=1:
	print('issue with prodn')
if reacn!=2:
	print('issue with reacn')
if M_val!= ((PInit/(8.3144598*TEMP)*6.0221409e+23)*1.0e-6):
	print('issue with M_val')
if N2_val!= ((PInit/(8.3144598*TEMP)*6.0221409e+23)*1.0e-6)*0.79:
	print('issue with N2_val')
if O2_val!= ((PInit/(8.3144598*TEMP)*6.0221409e+23)*1.0e-6)*0.2096:
	print('issue with O2_val')
if init_SMIL[0]!='CC1=CCC2CC1C2(C)C':
	print('issue with init_SMIL')
	
print('if no issues stated above then eqn_parser.extract_mechanism, eqn_parser.write_RO2_indices and formatting functions working fine')