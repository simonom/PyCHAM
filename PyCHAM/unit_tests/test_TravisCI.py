'''module to test on the Travis CI platform'''
# Travis Continuous Integration tests the software to ensure it
# runs on every update to the repository

def TravisCI(): # define function

	import sys
	import os
	# ensure modules can be seen 
	# (assumes calling function is in the home folder)
	sys.path.append(str(os.getcwd() + '/PyCHAM'))
	# the default model variables
	import def_mod_var
	def_mod_var.def_mod_var(1)
	import middle # the main call to program
	middle.middle()
	
	return()

TravisCI() # call the function
