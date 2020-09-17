'''module called to run the PyCHAM model, produces graphical user interface'''
from PyQt5.QtWidgets import QApplication, QWidget
import sys
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
	app = QApplication(sys.argv)
	ex = PyCHAM.PyCHAM()
	app.exec_()
	
	
main() # call function
sys.exit()
