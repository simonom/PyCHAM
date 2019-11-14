# function to test the front.py module
# set path to the PyCHAM folder
import os
import sys
import threading # for calling on further functions

dirpath = os.getcwd() # get current path
sys.path.append(os.path.split(dirpath)[0]) # add path to system path
print('importing front.py as in PyCHAM.py')
import front as model # import module to be tested
print('front.py imported fine, now calling as in PyCHAM.py to test')
# flag to tell module this is a test
testf = 1
# call on module as in PyCHAM.py
t = threading.Thread(target=model.run(testf))
t.daemon = False
t.start()

os.remove('test_var_store.pkl') # remove pickle file

print('front.py has passed testing')