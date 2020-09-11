'''determines number of photochemical reactions'''
# based on the photochemistry file path, the number of reactions is
# inferred

import os

# define function
def photo_num(photo_file):

	# inputs: -------------------------------------------------------------
	# photo_file - path to file containing abosrption cross-sections and
	# quantum yields
	# ---------------------------------------------------------------------
	# number of photolysis reactions, if this relevant
	cwd = os.getcwd() # address of current working directory
	if photo_file == str(cwd + '/PyCHAM/photofiles/MCMv3.2'):
		Jlen = 62 # for MCM (default name of photolysis parameters)
	else: # need to find out number of photolysis reactions
		# use Fortran indexing to be consistent with MCM photochemical reaction numbers
		Jlen = 1 
		# open file to read
		f = open(str(photo_file), 'r')
		for line in f: # loop through line
			if line.strip() == str('J_'+str(Jlen) + '_axs'):
				Jlen += 1
	return(Jlen)
