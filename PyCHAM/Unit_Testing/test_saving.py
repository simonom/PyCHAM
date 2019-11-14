# function to test saving module
print('function to test saving.py')
print('importing required modules')
print('importing modules needed by ode_gen')
import os
import sys
dirpath = os.getcwd() # get current path
sys.path.append(os.path.split(dirpath)[0]) # add path to system path
import pandas as pd
import os
import sys
import shutil
import numpy as np
print('import saving.py as in front.py')
from saving import saving
print('if no issue stated above, importing completed successfully')
print('now calling saving as in front')
fname = 'test'
y_mat = np.array(((2.46e+11, 2.95e+11, 1.00e-40, 4.62e+17, 1.00e-40, 1.00e-40, 1.00e-40,
			 1.00e-40, 7.17e+08, 2.52e+08, 1.00e-40, 1.00e-40, 1.00e-40, 1.79e+09,
			 4.76e+08, 1.00e-40, 1.00e-40, 1.00e-40, 1.98e+20, 1.00e-40),
			 (2.46e+11, 2.95e+11, 1.00e-40, 4.62e+17, 1.00e-40, 1.00e-40, 1.00e-40,
			 1.00e-40, 7.17e+08, 2.52e+08, 1.00e-40, 1.00e-40, 1.00e-40, 1.79e+09,
			 4.76e+08, 1.00e-40, 1.00e-40, 1.00e-40, 1.98e+20, 1.00e-40)))
t_out = np.array((0.0, 65.0))
Nresult = np.array(((934.0,65.0),(1.0e-40, 993.0)))
x2 = np.array(((237.0e-4,743.0e-4),(2.0e-2, 210.0e-3)))
num_sb=3
num_speci = 5
y_mw = (np.array((136.0, 17.0, 185.0, 18.0, 132.0))).reshape(num_speci,1)
resfname = 'test_out'

# remove temporary output folder
def handleRemoveReadonly(func, path, exc):
		excvalue = exc[1]
		if not os.access(path, os.W_OK):
			# Is the error an access error ?
			os.chmod(path, stat.S_IWUSR)
			func(path)
		else:
			raise

# remove if directory already exists
if os.path.exists(dirpath +'/output/'+fname+'/'+resfname)==1:
	shutil.rmtree(dirpath +'/output', ignore_errors=False, onerror=handleRemoveReadonly) # remove existing folder, onerror will change permission of directory if needed.
	
rbou = np.array((0.0, 400.0e-4, 8.0))
import scipy.constants as si
Cfactor = 101325.0*(si.N_A/(8.3144598e6*298.15))*1.0e-9
MV = (np.array((149.0, 176.0e-2, 193.0, 180.0e-1, 746.0e-1))).reshape(num_speci,1)
testf = 2
output_by_sim = saving(fname, y_mat, t_out, Nresult, x2, num_sb, y_mw, num_speci, 
							resfname, rbou, Cfactor, MV, testf)
print('now checking output')
if os.path.exists(dirpath +'/output/'+fname+'/'+resfname)==0:
	print('results directory not created')
else:
	print('test results directory created fine, now removing')


shutil.rmtree(dirpath +'/output', ignore_errors=False, onerror=handleRemoveReadonly) # remove existing folder, onerror will change permission of directory if needed. 


print('if no issues stated above then saving.py is working fine')