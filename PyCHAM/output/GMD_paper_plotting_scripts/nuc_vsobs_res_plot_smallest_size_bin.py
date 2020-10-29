'''script to plot results for the comparison of modelled nucleation against observations'''
# using the output produced by nuc_vsobs_inputs.txt and limonene_simp_scheme.txt, 
# assumed calling from the PyCHAM home folder - note that 24 size bin used for
# moving-centre and 12 size bins used for full-moving, with the save file name
# and size strucure inputs changed accordingly also

# import required modules
import  matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import numpy as np
import scipy.constants as si
import sys
import os
# ensure modules can be seen 
# (assumes calling from the home folder)
sys.path.append(str(os.getcwd() + '/PyCHAM'))
import retr_out

# ----------------------------------------------------------------------------------------
# function for doing colorbar tick labels in standard notation
def fmt(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	return r'${} \times 10^{{{}}}$'.format(a, b)

# ----------------------------------------------------------------------------------------
# create figure to plot results
#fig, (ax0, ax00, ax1, ax11) = plt.subplots(4, 1, figsize=(12, 6.5), sharex='all')
fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(12, 6.5), sharex='all')
fig.subplots_adjust(hspace = 0.05)

#par0 = ax0.twinx() # first parasite axis uses right-hand vertical axis
#par1 = par0.twiny() # first parasite axis uses upper horizontal axis

# ------------------------------------------------------------------------------
# full-moving simulation results
cwd = os.getcwd() # get current working directory
output_by_sim = str(cwd + '/PyCHAM/output/limonene_MCM_PRAM/nuc_vsobs_output_fm_12sb')

# required outputs from full-moving
(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, xfm, time_array, comp_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)


# isolate only particle number concentrations with a diameter between 10 nm and 46.9 nm,
# note xfm in um
Nsmall = np.zeros((len(time_array), 1))
lr = 22.e-3 # smallest radius of interest (um)
hr = 23.45e-3 # largest radius of interest (um)
rdiff = hr-lr
for it in range(len(time_array)):
	for isb in range(num_sb-wall_on):
		# cross-over into relevant radius space
		if rbou_rec[it, isb+1]<=lr or rbou_rec[it, isb]>=hr:
			continue
		# complete encompassing by model bounds
		if rbou_rec[it, isb]<=lr and rbou_rec[it, isb+1]>=hr: 
			rfrac = rdiff/(rbou_rec[it, isb+1]-rbou_rec[it, isb])
		# complete encompassing by reference bounds
		if rbou_rec[it, isb]>lr and rbou_rec[it, isb+1]<hr:
			rfrac = 1.
		# lower bound within required bounds
		if rbou_rec[it, isb]>=lr and rbou_rec[it, isb]<hr:
			rfrac = (hr-rbou_rec[it, isb])/(rbou_rec[it, isb+1]-rbou_rec[it, isb])
		# upper bound within required bounds
		if rbou_rec[it, isb+1]>=lr and rbou_rec[it, isb+1]<hr:
			rfrac = (rbou_rec[it, isb+1]-lr)/(rbou_rec[it, isb+1]-rbou_rec[it, isb])
		
		Nsmall[it, 0] += Ndry[it, isb]*rfrac

ax0.plot(time_array[time_array<1.], Nsmall[time_array<1.], '--b')

# observation part ---------------------------------------------------------------
import xlrd # for opening xlsx file

wb = xlrd.open_workbook(str(cwd+'/PyCHAM/output/GMD_paper_plotting_scripts/nuc_vsobs_obs.xlsx'))
wb = wb.sheet_by_index(0)

sr = 7#14 # starting row for desired time since O3 injection

obst = np.zeros((46-sr, 1)) # empty array for observation times (s)
obsN = np.zeros((46-sr, 23)) # empty array for particle number concentration (#/cc (air)) 
obsDp= np.zeros((1, 23)) # empty array for diameters (nm)

# get data
rc = 0 # row count
for r in range(0, wb.nrows): # row loop
	if (r > sr): # desired times (s)
		
		obstn = wb.cell_value(r, 1)
		if (isinstance(obstn, str)): # convert string type times to fraction of day
			
			splt = obstn.split(':')
			obstn = float(splt[0])*(1/24)+float(splt[1])*(1/(24*60))+float(splt[2])*(1/(24*60*60))
			if ((splt[0][0] == '-') and obstn > 0): # check for negative times
				obstn = -1*obstn
		obst[rc, 0] = obstn*(3600.*24.) # convert to seconds

	cc = 0 # column count
	for c in range(2, wb.ncols): # column loop
		if (r == 0):
			obsDp[0, cc] = wb.cell_value(r, c)
		if (r > sr):
			obsN[rc, cc] = wb.cell_value(r, c)
		cc += 1 # column count

	if (r > sr):
		rc += 1 # row count

# number concentration in smallest size bin - note values in excel already normalised by d(log10(Dp))
Nsmall_obs = obsN[:, 0]*(np.log10(46.9)-np.log10(30.))

ax0.plot(obst[obst[:, 0]<3600., 0]/3600., Nsmall_obs[obst[:, 0]<3600.], ':k')
plt.show()
import ipdb; ipdb.set_trace()

# ----------------------------------------------------------------------------------------
# moving-centre simulation results
cwd = os.getcwd() # get current working directory
output_by_sim = str(cwd + '/PyCHAM/output/limonene_simp_scheme/nuc_vsobs_output_mc')

# required outputs from full-moving
(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, xfm, time_array, comp_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

# in case we want to check on any gas-phase concentrations
#indx = comp_names.index('LIMONENE')
#ax0.plot(time_array, yrec[:, indx]*Cfac, '-r')

dlog10D = np.log10(rbou_rec[:, 1::]*2.0)-np.log10(rbou_rec[:, 0:-1]*2.0)
dNdD = Ndry/dlog10D # normalised number size distribution (#/cc 9air))
# two-point moving average
dNdD = (dNdD[:, 1::]+dNdD[:, 0:-1])/2.
rbou_rec = (rbou_rec[:, 1::]+rbou_rec[:, 0:-1])/2. 

# note - using same colour scheme as for full-moving above
for ti in range(len(time_array)-1): # loop through times
	p1 = ax1.pcolormesh(time_array[ti:ti+2], (rbou_rec[ti, :]*2*1e3), dNdD[ti, :].reshape(-1, 1), cmap=cm, norm=norm1)

ax1.set_yscale("log")
ax1.set_ylabel('$D_p$ (nm)', size = 14)
ax1.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
ax1.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
ax1.set_xlabel(r'Time since O3 injected (hours)', fontsize = 14)		
# observation part ---------------------------------------------------------------
# contour line plot of observations
p2 = ax1.contour(obst[:, 0]/3600., obsDp[0, :], dNdDp.transpose(), cmap=cm, norm=norm1)

ax1.text(txtl, 4.e2, '(b)', size = 14) # plot label
ax1.set_ylim(1.e1, 4.e2) # set vertical axis limits
# repeat for the close up --------------------------------------------------------
# note - using same colour scheme as for full-moving above
#for ti in range(len(time_array)-1): # loop through times
#	p1 = ax11.pcolormesh(time_array[ti:ti+2], (rbou_rec[ti, :]*2*1e3), dNdD[ti, :].reshape(-1, 1), cmap=cm, norm=norm1)

#ax11.set_yscale("log")
#ax11.set_ylabel('$D_p$ (nm)', size = 14)
#ax11.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
#ax11.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
#ax11.set_xlabel(r'Time since O3 injected (hours)', fontsize = 14)

#p2 = ax11.contour(obst[:, 0]/3600., obsDp[0, :], dNdDp.transpose(), cmap=cm, norm=norm1)

#ax11.set_ylim(7.e1, 3.e2) # set vertical axis limits

#ax11.text(txtl, 3.e2, '(d)', size = 14) # plot label

fig.savefig('fig11.png')


plt.show()
