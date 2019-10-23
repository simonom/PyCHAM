# function to test user_input.py module
# set path to the PyCHAM folder
import os
import sys
import numpy as np
dirpath = os.getcwd() # get current path
sys.path.append(os.path.split(dirpath)[0]) # add path to system path

print('testing pickling of files as done in PyCHAM.py')
# dump dummy variables in a pickle file in the same way that PyCHAM module does
import pickle # for storing inputs
list_vars = ('test',1,2,3,4,'test2',6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
			25,26,27,28,29,30,31,32,33)
with open('var_store.pkl','wb') as f:
		pickle.dump(list_vars,f)
print('pickled variables fine')
print('now retrieving variables from user_input.py as done in front.py')
# call user_input
import user_input as ui
# module to ask, receive and return required inputs as in front
[fname, num_sb, lowersize, uppersize, end_sim_time, resfname, tstep_len, 
tmax, TEMP, PInit, RH, lat, lon, start_sim_time, save_step, lightt, Cw, Ke, 
cham_dim, ChamR, wall_accom, nucv1, nucv2, nucv3, nuc_comp, inflectDp, pwl_xpre,  
pwl_xpro, inflectk, xmlname, init_conc, Comp0, Rader, voli, volP, pconc, 
	std] = ui.run(0, 2)

if fname!='test':
	print('issue with fname variable')
if num_sb!=2:
	print('issue with num_sb variable')
if lowersize!=2:
	print('issue with lowersize variable')
if uppersize!=3:
	print('issue with uppersize variable')
if end_sim_time!=4:
	print('issue with end_sim_time variable')
if resfname!='test2':
	print('issue with resfname variable')
if tstep_len!=6:
	print('issue with tstep_len')
if tmax!=6:
	print('issue with tmax variable')
if TEMP!=7:
	print('issue with TEMP variable')
if PInit!=8:
	print('issue with PInit variable')
if RH!=9:
	print('issue with RH variable')
if lat!=10:
	print('issue with lat variable')
if lon!=11:
	print('issue with lon variable')
if start_sim_time!=12:
	print('issue with start_sim_time variable')
if save_step!=14:
	print('issue with save_step variable')
if lightt!=16:
	print('issue with lightt variable')
if Cw!=13:
	print('issue with Cw variable')
if Ke!=18:
	print('issue with Ke variable')
if cham_dim!=15:
	print('issue with cham_dim variable')
if ChamR!=((16/(4.0*np.pi))**0.5):
	print('issue with ChamR variable')
if wall_accom!=17:
	print('issue with wall_accom variable')
if nucv1!=18:
	print('issue with nucv1 variable')
if nucv2!=19:
	print('issue with nucv2 variable')
if nucv3!=20:
	print('issue with nucv3 variable')
if nuc_comp!=21:
	print('issue with nuc_comp variable')
if inflectDp!=22:
	print('issue with inflectDp variable')
if pwl_xpre!=23:
	print('issue with pwl_xpre variable')
if pwl_xpro!=24:
	print('issue with pwl_xpro variable')
if inflectk!=25:
	print('issue with inflectk variable')
if xmlname!=27:
	print('issue with xmlname variable')
if init_conc!=28:
	print('issue with init_conc variable')
if Comp0!=29:
	print('issue with Comp0 variable')
if Rader!=26:
	print('issue with Rader variable')
if voli!=30:
	print('issue with voli variable')
if volP!=31:
	print('issue with volP variable')
if pconc!=32:
	print('issue with pconc variable')
if std!=33:
	print('issue with std variable')

	
print('testing when called by front.py finished, if no issues printed above, then call and return is fine')

print('now testing user_input in the manner of pickling files by front.py and loading files by res_plot_super.py')
# dump dummy variables in pickle file as front.py does
import pickle
list_vars = ['test1',1,2,3,4,'test2',5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
			 23,'test3',24,25,26]

with open('var_store.pkl','wb') as f:
	pickle.dump(list_vars,f)
print('pickled variables fine, now calling user_input.py to retrieve variables')

# open pickle file variables as in res_plot_super.py
[fname, num_sb, lowersize, uppersize, end_sim_time, resfname, tstep_len, 
	tmax, TEMP, PInit, RH, lat, lon, start_sim_time, save_step, lightt, Cw, Ke, 
	cham_dim, wall_accom, nucv1, nucv2, nucv3, nuc_comp, inflectDp, pwl_xpre, pwl_xpro, 
	inflectk, xmlname, Comp0, y_indx_plot] = ui.run(1,2)

if fname!='test1':
	print('issue with fname variable')
if num_sb!=2:
	print('issue with num_sb')
if lowersize!=2:
	print('issue with lowersize variable')
if uppersize!=3:
	print('issue with uppersize variable')
if end_sim_time!=4:
	print('issue with end_sim_time variable')
if resfname!= 'test2':
	print('issue with resfname variable')
if tstep_len!=5:
	print('issue with tstep_len variable')
if tmax!=5:
	print('issue with tmax variable')
if TEMP!=6:
	print('issue with TEMP variable')
if PInit!=7:
	print('issue with PInit variable')
if RH!=8:
	print('issue with the RH variable')
if lat!=9:
	print('issue with the lat variable')
if lon!=10:
	print('issue with the lon variable')
if start_sim_time!=11:
	print('issue with the start_sim_time variable')
if save_step!=13:
	print('issue with the save_step variable')
if lightt!=16:
	print('issue with the lightt variable')
if Cw!=12:
	print('issue with the Cw variable')
if Ke!=18:
	print('issue with Ke variable')
if cham_dim!=14:
	print('issue with cham_dim variable')
if wall_accom!=15:
	print('issue with wall_accom variable')
if nucv1!=16:
	print('issue with nucv1 variable')
if nucv2!=17:
	print('issue with nucv2 variable')
if nucv3!=18:
	print('issue with nucv3 variable')
if nuc_comp!=19:
	print('issue with nuc_comp variable')
if inflectDp!=20:
	print('issue with inflectDp variable')
if pwl_xpre!=21:
	print('issue with pwl_xpre variable')
if pwl_xpro!=22:
	print('issue with pwl_xpro variable')
if inflectk!=23:
	print('issue with inflectk variable')
if xmlname!='test3':
	print('issue with xmlname variable')
if Comp0!=25:
	print('issue with Comp0 variable')
if y_indx_plot!=26:
	print('issue with y_indx_plot variable')

print('testing when called by res_plot_super.py complete, if no issues printed above then code is fine')