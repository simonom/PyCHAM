'''module to produce examples of sensitivity of Gompertz nucleation function to user-defined parameters'''

import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd # for reading in .xls spreadsheets
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
# 	
# 	# ------------------------------------------------------------------------------------
# 	# second part is to show how nucleation parameters can be fitted to measurements
# 	
# 	# first, open MAC nucleation experiment measurements
# 	fnamepre = '/Users/Simon_OMeara/Documents/Manchester/postdoc_stuff/box-model/PyCHAM_Gitw/PyCHAM/output/'
# 	# file name for number size distribution contours
# 	fname = str(fnamepre+'MACobs/dndlogD.xlsx')
# 	sht_name = str('dNdlogD') # sheet name
# 	dNdlogDm = pd.read_excel(fname, sheet_name=sht_name)
# 	dNdlogDm = (dNdlogDm.values)[1::, 1::] # convert to numpy array from pandas data frame
# 	
# 	# file name for number size distribution times (hours)
# 	fname = str(fnamepre+'MACobs/dndlogDt.xlsx')
# 	sht_name = str('dndlogDt') # sheet name
# 	dNdlogDmt = pd.read_excel(fname, sheet_name=sht_name)
# 	dNdlogDmt = (dNdlogDmt.values)[1::, 1] # convert to numpy array from pandas data frame
# 	
# 	# file name for number size distribution diameters (nm)
# 	fname = str(fnamepre+'MACobs/Dp.xlsx')
# 	sht_name = str('Dp') # sheet name
# 	dNdlogDmDp = pd.read_excel(fname, sheet_name=sht_name)
# 	dNdlogDmDp = (dNdlogDmDp.values)[1::, 1] # convert to numpy array from pandas data frame
# 	
# 	# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
# 	colors = [(0.60, 0.0, 0.70), (0, 0, 1), (0, 1.0, 1.0), (0, 1.0, 0.0), (1.0, 1.0, 0.0), (1.0, 0.0, 0.0)]  # R -> G -> B
# 	n_bin = 100  # Discretizes the colormap interpolation into bins
# 	cmap_name = 'my_list'
# 	# Create the colormap
# 	cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
# 	
# 	# set contour levels
# 	levels = (MaxNLocator(nbins = 100).tick_values(np.min(dNdlogDm[~np.isnan(dNdlogDm)]), 
# 			np.max(dNdlogDm[~np.isnan(dNdlogDm)])))
# 
# 	# associate colours and contour levels
# 	norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)
# 	
# 	fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10,5), sharey=True, 
#         gridspec_kw={'width_ratios':[1,1,1.2]})
# 	CS1 = ax1.pcolormesh(dNdlogDmt, dNdlogDmDp, np.transpose(dNdlogDm[:,0:-1]), 
# 			cmap=cm, norm=norm1)
# 	ax1.set_ylim([0.0, 180.0])
# 	ax1.set_xlim([0.0, 4.0])
# 	ax1.set_xlabel('Time (hrs)', size=18)
# 	ax1.set_ylabel('Diameter (nm)', size=18)
# 	ax1.set_title('$\mathrm{a)\; nuc_{v1,v2,v3}=15x10^3,0,60}$', size=10)
# 	
# 	
# 	# ------------------------------------------------------------------------------------
# 	# second iteration plot
# 	CS2 = ax2.pcolormesh(dNdlogDmt, dNdlogDmDp, np.transpose(dNdlogDm[:,0:-1]), 
# 						cmap=cm, norm=norm1)
# 	ax2.set_ylim([0.0, 180.0])
# 	ax2.set_xlim([0.0, 4.0])
# 	ax2.set_xlabel('Time (hrs)', size=18)
# 	ax2.set_title('$\mathrm{b)\; nuc_{v1,v2,v3}=15x10^3,-70,60}$', size=10)
# 	
# 	# ------------------------------------------------------------------------------------
# 	# third iteration plot
# 	CS3 = ax3.pcolormesh(dNdlogDmt, dNdlogDmDp, np.transpose(dNdlogDm[:,0:-1]), 
# 						cmap=cm, norm=norm1)
# 	ax3.set_ylim([0.0, 180.0])
# 	ax3.set_xlim([0.0, 4.0])
# 	ax3.set_xlabel('Time (hrs)', size=18)
# 	ax3.set_title('$\mathrm{c)\; nuc_{v1,v2,v3}=25x10^3,-70,120}$', size=10)
# 
# 	# ------------------------------------------------------------------------------------	
# 	# colour bar
# 	# function for doing colorbar tick labels in standard notation
# 	def fmt(x, pos):
# 		a, b = '{:.1e}'.format(x).split('e')
# 		b = int(b)
# 		return r'${} \times 10^{{{}}}$'.format(a, b)
# 	
# 	cb = fig.colorbar(CS3, format=ticker.FuncFormatter(fmt))
# 	# colour bar label
# 	cb.set_label('dN/d$\mathrm{log_{10}}$($\mathrm{D_{p}}$) $\mathrm{(cm^{-3})}$', 
# 					size=18, rotation=270, labelpad=20)
# 	
# 	# ------------------------------------------------------------------------------------
# 	# simulation results
# 	
# 	# first iteration plot
# 	# open
# 	output_by_sim = '/Users/Simon_OMeara/Documents/Manchester/postdoc_stuff/box-model/PyCHAM_Gitw/PyCHAM/output/nuc_sens_chem/iter0'
# 	fname = str(output_by_sim+'/N') # particle number concentration (# particles/cc (air))
# 	N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
# 	
# 	# size bin radius bounds (um) for calculating dN/dlog10D (/cc (air))
# 	# name of file where particle size results saved
# 	fname = str(output_by_sim+'/sbb')
# 	sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
# 	# convert from um to nm and from radius to diameter
# 	sbb = sbb*2.0e3
# 
# 	# withdraw times
# 	fname = str(output_by_sim+'/t') # (s)
# 	t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
# 	# convert seconds to hours
# 	t_array = t_array/3600.0
# 	
# 	# diameters at particle size bin centres (um)
# 	fname = str(output_by_sim+'/x')
# 	x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
# 	# convert um to nm and radius to diameter
# 	x = x*2.0e3
# 	
# 	# calculate dlog10(Dp)
# 	# two-point moving sum of particle number concentrations (# particles/cc(air))
# 	N2p = (N[:,0:-1]+N[:,1::])
# 	# change first element in sbb from zero to very low value
# 	if sbb[0]==0.0:
# 		sbb[0] = 1.0e-40 # very low filler, effectively zero
# 	# log10 diameter bin bounds
# 	logDp = np.log10(sbb)
# 	# note we use 2 in the index as using two-point moving sum for N2p
# 	dlogDp = (logDp[2::]-logDp[0:-2])
# 	# quotient of particle number concentration per size bin and bin width
# 	dNdlogDp = N2p/dlogDp
# 	
# 	# add contour lines to measurement plot
# 	CS1sim = ax1.contour(t_array, (x[0,0:-1]+x[0,1::])/2.0, np.transpose(dNdlogDp), cmap=cm, norm=norm1)
# 	# ------------------------------------------------------------------------------------
# 	
# 	# second iteration plot
# 	# open
# 	output_by_sim = '/Users/Simon_OMeara/Documents/Manchester/postdoc_stuff/box-model/PyCHAM_Gitw/PyCHAM/output/nuc_sens_chem/iter1'
# 	fname = str(output_by_sim+'/N') # particle number concentration (# particles/cc (air))
# 	N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
# 	
# 	# size bin radius bounds (um) for calculating dN/dlog10D (/cc (air))
# 	# name of file where particle size results saved
# 	fname = str(output_by_sim+'/sbb')
# 	sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
# 	# convert from um to nm and from radius to diameter
# 	sbb = sbb*2.0e3
# 
# 	# withdraw times
# 	fname = str(output_by_sim+'/t') # (s)
# 	t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
# 	# convert seconds to hours
# 	t_array = t_array/3600.0
# 	
# 	# diameters at particle size bin centres (um)
# 	fname = str(output_by_sim+'/x')
# 	x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
# 	# convert um to nm and radius to diameter
# 	x = x*2.0e3
# 	
# 	# calculate dlog10(Dp)
# 	# two-point moving sum of particle number concentrations (# particles/cc(air))
# 	N2p = (N[:,0:-1]+N[:,1::])
# 	# change first element in sbb from zero to very low value
# 	if sbb[0]==0.0:
# 		sbb[0] = 1.0e-40 # very low filler, effectively zero
# 	# log10 diameter bin bounds
# 	logDp = np.log10(sbb)
# 	# note we use 2 in the index as using two-point moving sum for N2p
# 	dlogDp = (logDp[2::]-logDp[0:-2])
# 	# quotient of particle number concentration per size bin and bin width
# 	dNdlogDp = N2p/dlogDp
# 	
# 	# add contour lines to measurement plot
# 	CS2sim = ax2.contour(t_array, (x[0,0:-1]+x[0,1::])/2.0, np.transpose(dNdlogDp), 
# 							cmap=cm, norm=norm1)
# 	# ------------------------------------------------------------------------------------
# 	
# 	# third iteration plot
# 	# open
# 	output_by_sim = '/Users/Simon_OMeara/Documents/Manchester/postdoc_stuff/box-model/PyCHAM_Gitw/PyCHAM/output/nuc_sens_chem/iter2'
# 	fname = str(output_by_sim+'/N') # particle number concentration (# particles/cc (air))
# 	N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
# 	
# 	# size bin radius bounds (um) for calculating dN/dlog10D (/cc (air))
# 	# name of file where particle size results saved
# 	fname = str(output_by_sim+'/sbb')
# 	sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
# 	# convert from um to nm and from radius to diameter
# 	sbb = sbb*2.0e3
# 
# 	# withdraw times
# 	fname = str(output_by_sim+'/t') # (s)
# 	t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
# 	# convert seconds to hours
# 	t_array = t_array/3600.0
# 	
# 	# diameters at particle size bin centres (um)
# 	fname = str(output_by_sim+'/x')
# 	x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
# 	# convert um to nm and radius to diameter
# 	x = x*2.0e3
# 	
# 	# calculate dlog10(Dp)
# 	# two-point moving sum of particle number concentrations (# particles/cc(air))
# 	N2p = (N[:,0:-1]+N[:,1::])
# 	# change first element in sbb from zero to very low value
# 	if sbb[0]==0.0:
# 		sbb[0] = 1.0e-40 # very low filler, effectively zero
# 	# log10 diameter bin bounds
# 	logDp = np.log10(sbb)
# 	# note we use 2 in the index as using two-point moving sum for N2p
# 	dlogDp = (logDp[2::]-logDp[0:-2])
# 	# quotient of particle number concentration per size bin and bin width
# 	dNdlogDp = N2p/dlogDp
# 	
# 	# add contour lines to measurement plot
# 	CS3sim = ax3.contour(t_array, (x[0,0:-1]+x[0,1::])/2.0, np.transpose(dNdlogDp), 
# 							cmap=cm, norm=norm1)
# 
# 	fig.savefig('nuc_sens_vsobs.png')
# 	plt.show() # showing contour plot of number size distribution

# ----------------------------------------------------------------------------------------
# inputs
t = np.linspace(0, 3600, 100)# time (s)
nucv1 = np.array((50.0, 70.0, 50.0, 50.0))
nucv2 = np.array((-10.0, -10.0, -20.0, -10.0))
nucv3 = np.array((40.0, 40.0, 40.0, 80.0))

nuc_sens(t, nucv1, nucv2, nucv3) # call function with inputs