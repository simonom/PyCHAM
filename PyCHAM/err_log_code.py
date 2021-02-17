'''logs any problems reported by numpy or scipy'''
# issues found by numpy or scipy are logged here rather than printed to terminal

import os

def write(msg): # define function

	err_log = str(os.getcwd() + '/PyCHAM/err_log.txt')
	
	# open file for logging errors
	with open(err_log, 'a') as el:
		el.write(msg)
		el.close()
	
	return()