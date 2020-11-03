'''module to produce examples of sensitivity of Gompertz nucleation function to user-defined parameters'''
# self-sufficient - does not need to import any data, can be called from the PyCHAM home folder or 
# the GMD paper Results folder

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import matplotlib.gridspec as gspec # setting subplots on a grid

def nuc_sens(t, nucv1, nucv2, nucv3):

	Pnew = np.zeros((len(t), len(nucv1)))
	for i in range(len(nucv1)): # loop through combinations
		Pnew[:, i] = nucv1[i]*np.exp(nucv2[i]*(np.exp(-t/nucv3[i])))

	fig, (ax0) = plt.subplots(1, 1, figsize=(6,5))
	ax0.plot(t, Pnew[:, 0],label=r'$nuc_{v1}=50, nuc_{v2}=-10, nuc_{v3}=40$')
	ax0.plot(t, Pnew[:, 1],label=r'$nuc_{v1}=70, nuc_{v2}=-10, nuc_{v3}=40$')
	ax0.plot(t, Pnew[:, 2],label=r'$nuc_{v1}=50, nuc_{v2}=-20, nuc_{v3}=40$')
	ax0.plot(t, Pnew[:, 3],label=r'$nuc_{v1}=50, nuc_{v2}=-10, nuc_{v3}=80$')
	ax0.set_xlim([0.0, 600.0])
	ax0.set_xlabel(r'Time (s)', fontsize=10)
	ax0.set_ylabel(r'Concentration of new particles (#/cc (air))', fontsize=10)
	ax0.legend()
	fig.savefig('fig10.png')
	plt.show()			

# ----------------------------------------------------------------------------------------
# inputs
t = np.linspace(0, 3600, 100)# time (s)
nucv1 = np.array((50.0, 70.0, 50.0, 50.0))
nucv2 = np.array((-10.0, -10.0, -20.0, -10.0))
nucv3 = np.array((40.0, 40.0, 40.0, 80.0))



nuc_sens(t, nucv1, nucv2, nucv3) # call function with inputs
