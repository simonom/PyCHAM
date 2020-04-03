'''module to test the res_plot_super.py module'''
print('module to the res_plot_super.py module, note that test_res_plot_super.py must be called from the PyCHAM home directory')

# import necessary packages
import os
import sys

# import the required module
print('importing res_plot_super as done from PyCHAM.py')


cwd = os.getcwd() # address of current working directory
sys.path.insert(1, str(cwd+'/PyCHAM')) # ensure code can see res_plot_super
# setting test flag input
testf = 2

import res_plot_super as plotter
plotter.run(testf)