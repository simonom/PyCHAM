'''allows continuous integration testing using Travis CI without invoking the GUI (which doesn't work in the Travis application)'''

# address for Travis CI: https://travis-ci.org/account/repositories

import sys
import pickle # for storing inputs
import threading # for calling on further PyCHAM functions
import os
import numpy as np


dirpath = os.getcwd() # get current path
		
fname = dirpath+'/PyCHAM/inputs/Example_Run.txt' # chemical scheme file

with open(dirpath+'/fname.txt','w') as f:
	f.write(fname)
f.close()

xmlname = dirpath+'/PyCHAM/inputs/Example_Run_xml.xml' # chemical scheme xml file
		
with open(dirpath+'/xmlname.txt','w') as f:
	f.write(xmlname)
f.close()

# open the file
inname = dirpath+'/PyCHAM/inputs/Example_Run_inputs_TravisC1_test.txt' # name of model inputs file
inputs = open(inname, mode='r')

# read the file and store everything into a list
in_list = inputs.readlines()
inputs.close()

if len(in_list) != 36:
	print('Error: The number of model inputs is incorrect, should be 36, but is ' + str(len(in_list)) )
	sys.exit()
for i in range(len(in_list)):
	key, value = in_list[i].split('=')
	key = key.strip() # a string
	if key == 'Res_file_name':
		resfname = str(value.strip())
	if key == 'Total_model_time':
		end_sim_time = float(value.strip())
	if key == 'Time_step':
		tstep_len = float(value.strip())
	if key == 'Recording_time_step':
		save_step = float(value.strip())
	if key == 'Number_size_bins':
		num_sb = int(value.strip())
	if key == 'lower_part_size':
		lowersize = float(value.strip())
	if key == 'upper_part_size':
		uppersize = float(value.strip())
	if key == 'mass_trans_coeff':
		kgwt = float(value.strip())
	if key == 'eff_abs_wall_massC':
		Cw = float(value.strip())
	if key == 'Temperature':
		TEMP = float(value.strip())
	if key == 'PInit':
		PInit = float(value.strip())
	if key == 'RH':
		RH = float(value.strip())
	if key == 'lat':
		lat = float(value.strip())
	if key == 'lon':
		lon = float(value.strip())
	if key == 'daytime_start':
		dt_start = float(value.strip())
	if key == 'ChamSA':
		ChamSA = float(value.strip())
	if key == 'nucv1':
		nucv1 = float(value.strip())
	if key == 'nucv2':
		nucv2 = float(value.strip())
	if key == 'nucv3':
		nucv3 = float(value.strip())
	if key == 'nuc_comp':
		nuc_comp = int(value.strip())
	if key == 'inflectDp':
		inflectDp = float(value.strip())
	if key == 'Grad_pre_inflect':
		pwl_xpre = float(value.strip())
	if key == 'Grad_post_inflect':
		pwl_xpro = float(value.strip())
	if key == 'Kern_at_inflect':
		inflectk = float(value.strip())
	if key == 'Rader_flag':
		Rader = int(value.strip())
	if key == 'C0':
		C0 = [float(i) for i in (value.split(','))]
	if key == 'Comp0': # strip removes white space
		Comp0 = [str(i).strip() for i in (value.split(','))]
	if key == 'voli':
		
		if value.split(',')==['\n']:
			voli = np.empty(0)
		else:
			voli = [int(i) for i in (value.split(','))]	
	if key == 'volP':
		if value.split(',')==['\n']:
			volP = np.empty(0)
		else:
			volP = [float(i) for i in (value.split(','))]
	if key == 'pconc':
		pconc = float(value.strip())
	if key == 'std':
		std = float(value.strip())
	if key == 'loc':
		loc = float(value.strip())
	if key == 'scale':
		scale = float(value.strip())
	if key == 'core_diss':
		core_diss = float(value.strip())
	if key == 'light_time':
		light_time = [float(i) for i in (value.split(','))]
	if key == 'light_stat':
		light_stat = [int(i) for i in (value.split(','))]

# can't have a recording time step (s) less than the ode solver time step (s),
# so force to be equal and tell user
if save_step<tstep_len:
	print('Recording time step cannot be less than the ode solver time step, so increasing recording time step to be equal to input ode solver time step')
	save_step = tstep_len
# get names of chemical scheme and xml files
# set path prefix based on whether this is a test or not (i.e. whether test flag 
# file exists)
dirpath = os.getcwd() # get current path
	
f = open(dirpath+'/fname.txt','r')
content = f.readlines()
f.close()
fname = str(content[0])
f = open(dirpath+'/xmlname.txt','r')
content = f.readlines()
f.close()
xmlname = str(content[0])
# remove temporary files
os.remove(dirpath+'/fname.txt')
os.remove(dirpath+'/xmlname.txt')

# write variable values to pickle file
list_vars = [fname, num_sb, lowersize, uppersize, end_sim_time, 
resfname, tstep_len, TEMP, PInit, RH, lat, lon, dt_start, Cw,  
save_step, ChamSA, nucv1, nucv2, nucv3, nuc_comp,   
inflectDp, pwl_xpre, pwl_xpro, inflectk, Rader, xmlname, C0, Comp0, 
voli, volP, pconc, std, loc, scale, core_diss, light_stat, light_time,
kgwt]
	
if os.path.isfile(dirpath+'/testf.txt'):
	print('Model input buttons work successfully')
	with open('test_var_store.pkl','wb') as f:
		pickle.dump(list_vars,f)
		print('Pickle file dumped successfully')
		os.remove('test_var_store.pkl')
	f.close()
else:
	with open('PyCHAM/var_store.pkl','wb') as f:
		pickle.dump(list_vars,f)
	f.close()
	
# now run model automatically
import front as model
dirpath = os.getcwd() # get current path
if os.path.isfile(dirpath+'/testf.txt'):
	testf=2	
else:
	testf=0
# call on model to run
t = threading.Thread(target=model.run(testf))
t.daemon = False
t.start()