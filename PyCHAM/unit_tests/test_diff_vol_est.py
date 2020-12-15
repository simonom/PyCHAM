'''unit test for diff_vol_est'''
# the module to be tested - diff_vol_est is responsible for estimating the 
# diffusion volumes of components
# assumes calling from the PyCHAM home folder
print('unit test for estimating gas-phase diffusion volumes of components, assumed calling from the PyCHAM home folder')

import os
import sys
dir_path = os.getcwd() # current working directory
# temporarily add the PyCHAM folder to path
sys.path.append(str(dir_path+'/PyCHAM'))

import pybel
import diff_vol_est



# the SMILE strings of components to be tested
name_SMILE = []
name_SMILE.append('O=C=O') # CO2
name_SMILE.append('C') # CH4
name_SMILE.append('CC') # propane
name_SMILE.append('CCC') # ethane
name_SMILE.append('c1ccccc1') # benzene
name_SMILE.append('CC1=CCC2CC1C2(C)C') # alpha-pinene
name_SMILE.append('[H]N([H])[H]') # NH3

# pybel objects of components
Pybel_object = []

# generate pybel objects from SMILES
for i in name_SMILE:	
	Pybel_object.append(pybel.readstring('smi', i))


# get diffusion volumes
diff_vol = diff_vol_est.diff_vol_est(Pybel_object)

# compare output against Table 4.1 of  the Taylor (1993) textbook 
# Multicomponent Mass Transfer, ISBN: 0-471-57417-1
print(diff_vol)