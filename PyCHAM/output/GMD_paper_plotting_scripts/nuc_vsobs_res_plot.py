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
output_by_sim = str(cwd + '/PyCHAM/output/limonene_simp_scheme/nuc_vsobs_output_mc_n1_2e4_n2_-4e_n3_1e2')

# required outputs from full-moving
(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, xfm, time_array, comp_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

# in case we want to check on any gas-phase concentrations
#indx = comp_names.index('LIMONENE')
#ax0.plot(time_array, yrec[:, indx]*Cfac, '-r')

dlog10D = np.log10(rbou_rec[:, 1::]*2.0)-np.log10(rbou_rec[:, 0:-1]*2.0)
dNdD = Ndry/dlog10D # normalised number size distribution (#/cc 9air))
# moving average
#dNdD = (dNdD[:, 3::]+dNdD[:, 2:-1]+dNdD[:, 1:-2]+dNdD[:, 0:-3])/3.
#rbou_rec = (rbou_rec[:, 3::]+rbou_rec[:, 2:-1]+rbou_rec[:, 1:-2]+rbou_rec[:, 0:-3])/3.

# two-point moving average
#dNdD = (dNdD[:, 1::]+dNdD[:, 0:-1])/2.
#rbou_rec = (rbou_rec[:, 1::]+rbou_rec[:, 0:-1])/2. 

# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
colors = [(0.60, 0.0, 0.70), (0, 0, 1), (0, 1.0, 1.0), (0, 1.0, 0.0), (1.0, 1.0, 0.0), (1.0, 0.0, 0.0)]  # R -> G -> B
n_bin = 100  # discretizes the colormap interpolation into bins
cmap_name = 'my_list'
# Create the colormap
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N = n_bin)
	
# set contour levels
levels = (MaxNLocator(nbins = 100).tick_values(np.min(dNdD), np.max(1.8e5)))
	
# associate colours and contour levels
norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)

for ti in range(len(time_array)-1): # loop through times
	p0 = ax0.pcolormesh(time_array[ti:ti+2], (rbou_rec[ti, :]*2*1e3), dNdD[ti, :].reshape(-1, 1), cmap=cm, norm=norm1)

ax0.set_yscale("log")
ax0.set_ylabel('$D_p$ (nm)', size = 14)
ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
ax0.set_yticks([1.e2])
ax0.set_yticklabels(['$10^2$'])
ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
# colorbar
cax = plt.axes([0.30, 0.92, 0.40, 0.01])
cb = plt.colorbar(p0, cax = cax, ticks=[0., 6.0e4, 1.2e5, 1.8e5], format = ticker.FuncFormatter(fmt), orientation = 'horizontal')
cb.ax.tick_params(labelsize = 12)
# colour bar label
cb.set_label('dN $\mathrm{(cm^{-3})}$/dlog10($D_p$)', size = 12, rotation = 0, labelpad = -44.)

# observation part ---------------------------------------------------------------
import xlrd # for opening xlsx file

wb = xlrd.open_workbook(str(cwd+'/PyCHAM/output/GMD_paper_plotting_scripts/nuc_vsobs_obs.xlsx'))
wb = wb.sheet_by_index(0)

sr = 14 # starting row for desired time since O3 injection

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

# normalised number concentration - note values in excel already normalised by d(log10(Dp))
dNdDp = obsN

# contour line plot of observations
p2 = ax0.contour(obst[:, 0]/3600., obsDp[0, :], dNdDp.transpose(), cmap=cm, norm=norm1)

txtl = -0.2 # horizontal positioning of plot labels
ax0.text(txtl, 4.e2, '(a)', size = 14) # plot label
ax0.set_ylim(1.e1, 4.e2) # set vertical axis limits
# repeat for the close up ---------------------------------------------------------------
#for ti in range(len(time_array)-1): # loop through times
#	p0 = ax00.pcolormesh(time_array[ti:ti+2], (rbou_rec[ti, :]*2*1e3), dNdD[ti, :].reshape(-1, 1), cmap=cm, norm=norm1)

#ax00.set_yscale("log")
#ax00.set_ylabel('$D_p$ (nm)', size = 14)
#ax00.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
#ax00.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')

# contour line plot of observations
#p2 = ax00.contour(obst[:, 0]/3600., obsDp[0, :], dNdDp.transpose(), cmap=cm, norm=norm1)

#ax00.set_ylim(7.e1, 3.e2) # set vertical axis limits
#ax00.text(txtl, 3.e2, '(b)', size = 14) # plot label

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
