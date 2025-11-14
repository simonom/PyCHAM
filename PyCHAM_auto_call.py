'''code to automatically select PyCHAM inputs from INGENIOUS 
data, and call PyCHAM. Note that model variable files are 
created by mod_var_setup.py in OneDrive/SOAPRA/EMEP'''

def autorun(): # define function

	# import dependencies
	import os
	import platform
	import sys
	from PyQt6.QtWidgets import QApplication, QWidget	

	# get operating system
	if 'Darwin' in platform.system() or 'Linux' in platform.system():
		pj = '/'
		PyCHAM_path = '/Users/user/Documents/GitHub/PyCHAM/PyCHAM'
	if 'Win' in platform.system():
		pj = '\\'
		PyCHAM_path = 'C:\\Users\\Psymo\\Desktop\\PyCHAM\\PyCHAM\\PyCHAM'
		# the directory with inputs
		directory = 'C:\\Users\\Psymo\\OneDrive - The University of Manchester\\SOAPRA\\EMEP\\PyCHAM_inputs'

	# create dictionary to pass to PyCHAM
	param_const = {}
	# filler
	param_const['param_ranges'] = []		
	param_const['sim_num'] = 1	

	# simulation type
	param_const['sim_type'] = 'standard_call'

	# loop through model runs
	for runi in os.listdir(directory):
	
		if ('mod_var' in (runi)):

			# get the name of folder to save results to
			mod_var_path = str(directory + pj + runi)
			inputs = open(mod_var_path, mode= 'r' ) # open model variables file
			in_list = inputs.readlines() # read file and store everything into a list
			inputs.close() # close file
			# get name of folder to save results to
			for listi in in_list: # loop through list items
				if ('res_file_name' in listi):
					# get directory to save results to
					rfn = listi[listi.index('=')+1::]
					break # break out of list loop

			# continue onto next simulation if results 
			# already exist
			if (os.path.isdir(rfn)):
				break 			

			# get the model variables file
			param_const['mod_var_name'] = mod_var_path
			print(mod_var_path)
			# call PyCHAM
			# ensure gui module can be seen
			sys.path.append(PyCHAM_path)
			
			import gui
			app = QApplication(sys.argv)

			ex = gui.PyCHAM(param_const)
		
			sys.exit() # fully disconnect from QApplication
			print(str('completed ' + str(runi)))

	return() # end function

autorun() # call function