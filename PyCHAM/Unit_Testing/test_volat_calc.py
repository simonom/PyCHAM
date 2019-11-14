# function to test volat_calc.py which calls UManSysProp
print('function to test volat_calc.py')
import os
import sys
dirpath = os.getcwd() # get current path
sys.path.append(os.path.split(dirpath)[0]) # add path to system path

print('importing required modules, note for volat_calc this requires internet connection to download the latest UManSysProp repository')

import numpy as np
import sys
import os
from git import Repo
import shutil
import scipy.constants as si

# download latest version of umansysprop
cwd = os.getcwd() # address of current working directory

if os.path.isdir(cwd + '/umansysprop'): # check if there is an existing umansysprop folder
	def handleRemoveReadonly(func, path, exc):
		excvalue = exc[1]
		if not os.access(path, os.W_OK):
			# Is the error an access error ?
			os.chmod(path, stat.S_IWUSR)
			func(path)
		else:
			raise
	shutil.rmtree(cwd + '/umansysprop', ignore_errors=False, onerror=handleRemoveReadonly) # remove existing folder, onerror will change permission of directory if needed. 

sys.path.insert(1, (cwd + '/umansysprop')) # address for updated version
git_url = 'https://github.com/loftytopping/UManSysProp_public.git'
Repo.clone_from(git_url, (cwd + '/umansysprop'))

from umansysprop import boiling_points
from umansysprop import vapour_pressures
from umansysprop import liquid_densities

print('if no issues stated above, imported fine')
print('now calling volat_calc as in front.py')
from volat_calc import volat_calc
spec_list = ['CC1=CCC2CC1C2(C)C','[OH]','[O]OC(C)(C)C1CC=C(C)C(O)C1']
import pybel
Pybel_objects = []
for name_only in spec_list:
	Pybel_object = pybel.readstring('smi', name_only)
	# append to Pybel object list
	Pybel_objects.append(Pybel_object)
TEMP = 298.15
H2Oi = 3
num_speci = 5
Psat_water = -1.5
voli = np.array((2,4))
volP = np.array((1.0e-30, 1.0e-30))
testf = 0
corei = 4
[Psat, y_dens, Psat_Pa] = volat_calc(spec_list, Pybel_objects, TEMP, H2Oi, num_speci,  
								Psat_water, voli, volP, testf, corei)

print('now asserting returned values are as expected')
if int(Psat[0]*1e-15)!=110 or int(Psat[1]*1e-20)!=526 or int(Psat[2]*1e18)!=242 or int(Psat[3]*1e-15)!=778 or int(Psat[4]*1e18)!=242:
	print('issue with y_dens, possibly due to UManSysProp')
if int(Psat_Pa[0]*1e0)!=456 or int(Psat_Pa[1]*1e-6)!=216 or int(Psat_Pa[2]*1e32)!=100 or int(Psat_Pa[3]*1e-1)!=320 or int(Psat_Pa[4]*1e32)!=100:
	print('issue with y_dens, possibly due to UManSysProp')
if int(y_dens[0])!=870 or int(y_dens[1])!=1133 or int(y_dens[2])!=1033 or int(y_dens[3])!=1000 or int(y_dens[4])!=1770:
	print('issue with y_dens, possibly due to UManSysProp')
	
print('now removing temporary UManSysProp directory')

if os.path.isdir(cwd + '/umansysprop'): # check if there is an existing umansysprop folder
	def handleRemoveReadonly(func, path, exc):
		excvalue = exc[1]
		if not os.access(path, os.W_OK):
			# Is the error an access error ?
			os.chmod(path, stat.S_IWUSR)
			func(path)
		else:
			raise
	shutil.rmtree(cwd + '/umansysprop', ignore_errors=False, onerror=handleRemoveReadonly) # remove existing folder, onerror will change permission of directory if needed. 


print('if no issues stated above, volat_calc is working fine, test complete')