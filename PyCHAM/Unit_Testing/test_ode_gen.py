# function to test ode_gen.py is working fine
print('function to test ode_gen.py (and its called modules: kimt_calc, mov_cen_main, coag, wallloss and nuc) is working fine')

print('importing modules needed by ode_gen')
import os
import sys
dirpath = os.getcwd() # get current path
sys.path.append(os.path.split(dirpath)[0]) # add path to system path

import numpy as np
from assimulo.problem import Explicit_Problem
from assimulo.solvers import CVode
import numba
from numba import jit, f8
import matplotlib.pyplot as plt
import ipdb
from kimt_calc import kimt_calc
from recording import recording
from mov_cen_main import mov_cen_main as movcen # moving centre method for rebinning
from coag import coag
from wallloss import wallloss
from nuc import nuc
import scipy.constants as si

print('importing ode_gen.py as in front.py')
from ode_gen import ode_gen

print('if no issues stated above, imports are working fine')

tstep_len = 60.0
y = np.array((2.46e+11, 2.95e+11, 1.00e-40, 4.62e+17, 1.00e-40, 1.00e-40, 1.00e-40,
			 1.00e-40, 7.17e+08, 2.52e+08, 1.00e-40, 1.00e-40, 1.00e-40, 1.79e+09,
			 4.76e+08, 1.00e-40, 1.00e-40, 1.00e-40, 1.98e+20, 1.00e-40))
num_speci = 5
num_eqn = 1
rindx = np.array((0, 1)).reshape(1,2)
pindx = np.array((2)).reshape(1,1)
rstoi = np.array((1, 1)).reshape(1,2)
pstoi = np.array((1)).reshape(1,1)
H2Oi = 3
TEMP = 298.15
RO2_indices = np.array((4, 2)).reshape(1,2)
num_sb = 3
Psat = (np.array((110*1e15, 526*1e20, 242/1e18, 778*1e15, 242/1e18))).reshape(num_speci,1)
mfp = (np.array((2.49e-8, 2.56e-7, 1.61e-8, 1.17e-7, 2.45e-8))).reshape(5,1)
accom_coeff = np.ones((num_speci,1))
surfT = 72.0
y_dens = np.array((870, 1133, 1033, 1000, 1770)).reshape(num_speci,1)
N_perbin = np.array((934, 65))
DStar_org = (np.array((7.40e-06, 1.20e-04, 5.56e-06, 2.77e-05, 7.32e-06))).reshape(5,1)
y_mw = (np.array((130.0, 2.0, 200.0, 18.0, 132.14))).reshape(5,1)
x = np.array((237.0e-4, 744.0e-4))
core_diss = 3.0
Varr = np.array((564.0e-7, 172.0e-5))
Vbou = np.array((0.0, 268.0e-6, 214.0e2))
RH = 0.60
rad0 = np.array((2.0e-2, 6.0e-2))
Vol0 = np.array((335.0e-7, 904.0e-6))
end_sim_time = 60.0
pconc = 1000.0
save_step = 60.0
rbou = np.array((0.0, 400.0e-4, 8.0))
therm_sp = (np.array((0,0,0,218,0))).reshape(num_speci,1)
Kw = np.array((399.0, 675.0, 133.0, 430.0, 187.0)).reshape(num_speci,1)
kgwt = 332.0e-25
Cw = 301.0e13
light_time = np.array((0.0, 14400.0))
light_stat = np.array((0, 0))
nreac = np.array((1)).reshape(1,1)
nprod = np.array((1)).reshape(1,1)
prodn = 1
reacn = 2
new_partr = 535.0e-10
MV = (np.array((149.0, 176.0e-2, 193.0, 180.0e-1, 746.0e-1))).reshape(num_speci,1)
nucv1 = 31000.0
nucv2 = -55.0
nucv3 = 180.0
inflectDp = 2.0e-7
pwl_xpre = 0.0
pwl_xpro = -7.0e-5
inflectk = 9.5e-5
nuc_comp = -3
ChamR = 2.95
Rader = 0
PInit = 101325.0
testf = 2

print('now calling ode_gen.py with test inputs')

[t_out, y_mat, Nresult, x2] = ode_gen(tstep_len, y, num_speci, num_eqn, rindx, pindx, 
				rstoi, pstoi, H2Oi, TEMP, RO2_indices, 
				num_sb, Psat, mfp, accom_coeff, surfT, y_dens, N_perbin,
				DStar_org, y_mw, x, core_diss, Varr, Vbou, RH, rad0, Vol0,
				end_sim_time, pconc, save_step, 
				rbou, therm_sp, Cw, light_time, light_stat,
				nreac, nprod, prodn,
				reacn, new_partr, MV, nucv1, nucv2, nucv3, inflectDp, pwl_xpre, 
				pwl_xpro, inflectk, nuc_comp, ChamR, Rader, PInit, testf, kgwt)

print('now checking outputs')
if int(t_out[0])!=0 or int(t_out[1])!=65:
	print('issue with t_out')
if int(y_mat[0,0]*1e-9)!=246 or int(y_mat[0,1]*1e-9)!=295 or int(y_mat[1,0]*1e35)!=441 or int(y_mat[1,1]*1e-8)!=489 or int(y_mat[1,2]*1e-9)!=123 or int(y_mat[1,12]*1e-9)!=122:
	print('issue with y_mat, possibly due to ode solver or mov_cen_main')
if int(Nresult[0,0])!=934 or int(Nresult[0,1])!=65 or int(Nresult[1,0])!=0 or int(Nresult[1,1])!=994:
	print('issue with Nresult')
if int(x2[0,0]*1e4)!=237 or int(x2[0,1]*1e4)!=743 or int(x2[1,0]*1e2)!=2 or int(x2[1,1]*1e2)!=127:
	print('issue with x2')
print('if no issue stated above, then ode_gen working fine')