'''unit test for estimating natural light intensity'''
# will estimate and plot intensity as a function of time
# assumes calling from the PyCHAM home folder
print('unit test that will plot the photochemical rate as a function of day time for the MCM equation with greatest photosensitivity to natural sunlight, when sunlight intensity change during a day is at a maximum for Earth')

import os
import sys
dir_path = os.getcwd() # current working directory
# temporarily add the PyCHAM folder to path
sys.path.append(str(dir_path+'/PyCHAM'))

import zenith
import matplotlib.pyplot as plt
import numpy as np

# define function
def test_zenith():

	# for a conservative estimate of how natural light intensity varies
	# select the latitude as the tropic of cancer (23.5 degrees north)
	# at the summer solstice, which is when the maximum angle of 
	# incidence is 90 degrees and the minimum is 0 degrees, and 
	# therefore change to photon intensity is greatest

	# time of day (s)
	timer = np.arange(3.6e3*24)
	lat = 23.5
	lon = 0.
	DayOfYear = 172
	secx = np.empty((len(timer)))
	cosx = np.empty((len(timer)))

	cnt = 0
	for time in timer:
		(secx[cnt], cosx[cnt]) = zenith.zenith(time, lat, lon, DayOfYear)
		cnt += 1
	# estimate photochemistry reaction rate coefficients for the MCM
	# chemical scheme
	J = np.empty((35, len(cosx)))
	J[0, :]=6.073E-05*cosx**(1.743)*np.exp(-1.0*0.474*secx)
	J[1, :]=4.775E-04*cosx**(0.298)*np.exp(-1.0*0.080*secx)
	J[2, :]=1.041E-05*cosx**(0.723)*np.exp(-1.0*0.279*secx)
	J[3, :]=1.165E-02*cosx**(0.244)*np.exp(-1.0*0.267*secx)
	J[4, :]=2.485E-02*cosx**(0.168)*np.exp(-1.0*0.108*secx)
	J[5, :]=1.747E-01*cosx**(0.155)*np.exp(-1.0*0.125*secx)
	J[6, :]=2.644E-03*cosx**(0.261)*np.exp(-1.0*0.288*secx)
	J[7, :]=9.312E-07*cosx**(1.230)*np.exp(-1.0*0.307*secx)
	J[8, :]=4.642E-05*cosx**(0.762)*np.exp(-1.0*0.353*secx)
	J[9, :]=6.853E-05*cosx**(0.477)*np.exp(-1.0*0.323*secx)
	J[10, :]=7.344E-06*cosx**(1.202)*np.exp(-1.0*0.417*secx)
	J[11, :]=2.879E-05*cosx**(1.067)*np.exp(-1.0*0.358*secx)
	J[12, :]=2.792E-05*cosx**(0.805)*np.exp(-1.0*0.338*secx)
	J[13, :]=1.675E-05*cosx**(0.805)*np.exp(-1.0*0.338*secx)
	J[14, :]=7.914E-05*cosx**(0.764)*np.exp(-1.0*0.364*secx)
	J[15, :]=1.140E-05*cosx**(0.396)*np.exp(-1.0*0.298*secx)
	J[16, :]=1.140E-05*cosx**(0.396)*np.exp(-1.0*0.298*secx)
	J[17, :]=7.992E-07*cosx**(1.578)*np.exp(-1.0*0.271*secx)
	J[18, :]=5.804E-06*cosx**(1.092)*np.exp(-1.0*0.377*secx)
	J[19, :]=1.836E-05*cosx**(0.395)*np.exp(-1.0*0.296*secx)
	J[20, :]=1.836E-05*cosx**(0.395)*np.exp(-1.0*0.296*secx)
	J[21, :]=6.845E-05*cosx**(0.130)*np.exp(-1.0*0.201*secx)
	J[22, :]=1.032E-05*cosx**(0.130)*np.exp(-1.0*0.201*secx)
	J[23, :]=3.802E-05*cosx**(0.644)*np.exp(-1.0*0.312*secx)
	J[24, :]=1.537E-04*cosx**(0.170)*np.exp(-1.0*0.208*secx)
	J[25, :]=3.326E-04*cosx**(0.148)*np.exp(-1.0*0.215*secx)
	J[26, :]=7.649E-06*cosx**(0.682)*np.exp(-1.0*0.279*secx)
	J[27, :]=1.588E-06*cosx**(1.154)*np.exp(-1.0*0.318*secx)
	J[28, :]=1.907E-06*cosx**(1.244)*np.exp(-1.0*0.335*secx)
	J[29, :]=2.485E-06*cosx**(1.196)*np.exp(-1.0*0.328*secx)
	J[30, :]=4.095E-06*cosx**(1.111)*np.exp(-1.0*0.316*secx)
	J[31, :]=1.135E-05*cosx**(0.974)*np.exp(-1.0*0.309*secx)
	J[32, :]=7.549E-06*cosx**(1.015)*np.exp(-1.0*0.324*secx)
	J[33, :]=3.363E-06*cosx**(1.296)*np.exp(-1.0*0.322*secx)
	J[34, :]=7.537E-04*cosx**(0.499)*np.exp(-1.0*0.266*secx)
	
	maxJ = np.empty((J.shape[0]))
	# the maximum photochemistry rate
	for i in range(0, 35):
		maxJ[i] = max(J[i, :])
	maxi = np.where(maxJ == max(maxJ))
	plt.semilogy(timer, J[maxi[0], :][0, :], label = str('Eq. # ' + str(maxi[0][0])))
	plt.xlabel('Time (s)')
	plt.ylabel('Photochemistry rate (/s)')
	plt.legend()
	plt.show()

	return()

# call function
test_zenith()
