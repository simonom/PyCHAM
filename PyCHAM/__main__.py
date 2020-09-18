'''module called to run the PyCHAM model, produces graphical user interface'''
from PyQt5.QtWidgets import QApplication, QWidget
import sys
import os
# to push to Github repo:
# switch to wanted branch: git checkout -b branch_name
# if not already linked to github: git init
# followed by: git remote add origin https://github.com/simonom/PyCHAM.git
# then use command above to switch to wanted branch
# git add .
# git commit -m "commit comment"
# git push origin master

def main(): 

	import PyCHAM
	import signal
	
	# what to do if user breaks in terminal using ctrl+c and the 
	# programme is running, specifically when it is not waiting for
	# a button on the GUI to be pressed
	def sig_han(signum, frame):
		sys.exit()
	
	app = QApplication(sys.argv)
	signal.signal(signal.SIGINT, sig_han) # in case of ctrl+c	
	ex = PyCHAM.PyCHAM()
	app.exec_()
	
	return()
	

if __name__ == '__main__':
	main() # call function
	
sys.exit()
