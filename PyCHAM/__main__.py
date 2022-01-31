##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2022 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
#                                                                                       											 #
#    All Rights Reserved.                                                                									 #
#    This file is part of PyCHAM                                                         									 #
#                                                                                        											 #
#    PyCHAM is free software: you can redistribute it and/or modify it under              						 #
#    the terms of the GNU General Public License as published by the Free Software       					 #
#    Foundation, either version 3 of the License, or (at your option) any later          						 #
#    version.                                                                            										 #
#                                                                                        											 #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT                						 #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       			 #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              				 #
#    details.                                                                            										 #
#                                                                                        											 #
#    You should have received a copy of the GNU General Public License along with        					 #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                 							 #
#                                                                                        											 #
##########################################################################################
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
