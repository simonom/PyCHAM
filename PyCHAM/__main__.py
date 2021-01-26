'''module called to run the PyCHAM programme via a graphical user interface'''
from PyQt5.QtWidgets import QApplication, QWidget
import sys
import os

def main(): 

	import gui
	import signal
	
	# what to do if user breaks in terminal using ctrl+c and the 
	# programme is running, specifically when it is not waiting for
	# a button on the GUI to be pressed
	def sig_han(signum, frame):
		sys.exit()
	
	app = QApplication(sys.argv)
	signal.signal(signal.SIGINT, sig_han) # in case of ctrl+c	
	ex = gui.PyCHAM()
	app.exec_()
	
	return()

if __name__ == '__main__':
	main() # call function
	
sys.exit()
