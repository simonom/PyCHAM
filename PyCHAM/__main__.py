'''module called to run the PyCHAM model, produces graphical user interface'''
from PyQt5.QtWidgets import QApplication
import sys
# to push to Github repo:
# git add .
# git commit -m "commit comment"
# git push origin master

def main(): 
	
	from PyCHAM import PyCHAM
	app = QApplication(sys.argv)
	ex = PyCHAM()
	app.exec_()
main()