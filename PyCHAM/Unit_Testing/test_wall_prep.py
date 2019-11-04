# function to test wall_prep.py
print('function to test wall_prep.py')
import os
import sys
dirpath = os.getcwd() # get current path
sys.path.append(os.path.split(dirpath)[0]) # add path to system path

print('importing modules required for wall_prep')
import numpy as np
import scipy.constants as si

print('if no issues stated above, importing successful')
print('now calling wall_prep.py as in front.py')
from wall_prep import wall_prep
print('wall_prep imported okay')
TEMP = 298.15
num_speci = 5
y_mw = (np.array((136.0, 17.0, 185.0, 18.0, 132.0))).reshape(num_speci,1)
Psat_Pa = (np.array((456, 216, 100, 320, 100))).reshape(num_speci,1)
wall_accom = 1.0
testf = 0
[wall_accom, Kw] = wall_prep(TEMP, Psat_Pa, y_mw, num_speci, wall_accom, testf)
print('now asserting returned values correct')
if (wall_accom!=1.0).sum()>0:
	print('issue with wall_accom')
if (int(Kw[0]*1e4)!=399) or (int(Kw[1]*1e3)!=675) or (int(Kw[2]*1e3)!=133) or (int(Kw[3]*1e3)!=430) or (int(Kw[4]*1e3)!=187):
	print('issue with Kw')
print('if no issues stated above, then wall_prep is working fine and test is complete')