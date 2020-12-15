'''script to plot results for the comparison of modelled nucleation against observations'''
# using the output produced by fig11_mod_var.txt and fig11_simp_scheme.txt, 
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
# (assumes calling from the home folder, but also works from the GMD paper Results 
# folder)
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

fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(12, 6.5), sharex='all')
fig.subplots_adjust(hspace = 0.05)


# ------------------------------------------------------------------------------
# results - note that all results for table 1 saved in fig11_data
cwd = os.getcwd() # get current working directory
try: # if calling from the PyCHAM home directory
	output_by_sim = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig11_data/nuc_vsobs_output_fm_n1_2e4_n2_-4e2_n3_1e2')
	# required outputs
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, xfm, time_array, comp_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

except: # if calling from the GMD paper Results folder
	output_by_sim = str(cwd + '/fig11_data/nuc_vsobs_output_fm_n1_2e4_n2_-4e2_n3_1e2')
	# required outputs
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, xfm, time_array, comp_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

	
dlog10D = np.log10(rbou_rec[:, 1::]*2.0)-np.log10(rbou_rec[:, 0:-1]*2.0)
dNdD = Ndry/dlog10D # normalised number size distribution (#/cc (air))

# observation part ---------------------------------------------------------------
import xlrd # for opening xlsx file

try: # if calling from the PyCHAM home folder
	wb = xlrd.open_workbook(str(cwd+'/PyCHAM/output/GMD_paper_plotting_scripts/fig11_data/obs.xlsx'))
	wb = wb.sheet_by_index(0)
except: # if calling from the GMD paper Results folder
	wb = xlrd.open_workbook(str(cwd+'/fig11_data/obs.xlsx'))
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

# consider only the observations from the first hour of the experiment
t_indx_max = int(sum((obst/3600.) <= 1.))

obst = obst[0:t_indx_max, :]
dNdDp = dNdDp[0:t_indx_max, :]

# make simulation and observation results comparable by mapping simulation results onto 
# observation times and sizes (the latter by interpolation)
dNdD_interp = np.zeros((dNdDp.shape[0], dNdDp.shape[1]))

# observation diameters (um)
obss = obsDp*1.e-3

for it in range(len(obst)-1): # loop through observation time steps
	# observation time (hours)
	obstn = obst[it]/3600.

	# index of simulation result at this time
	simi = np.where(time_array==((time_array[time_array>=obstn])[0]))[0]
	
	# fill model matrix - convert simulation radii to diameter
	dNdD_interp[it, :] = np.interp(obss[0, :], np.squeeze(xfm[simi, :]*2.), np.squeeze(dNdD[simi, :]))


# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
colors = [(0.60, 0.0, 0.70), (0, 0, 1), (0, 1.0, 1.0), (0, 1.0, 0.0), (1.0, 1.0, 0.0), (1.0, 0.0, 0.0)]  # R -> G -> B
n_bin = 100  # discretizes the colormap interpolation into bins
cmap_name = 'my_list'
# Create the colormap
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N = n_bin)
	
# set contour levels
levels = (MaxNLocator(nbins = 100).tick_values(np.min(dNdDp), np.max(dNdDp)))
	
# associate colours and contour levels
norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)
p1 = ax0.contour(obst[:, 0]/3600., obsDp[0, :], dNdD_interp.transpose(), cmap=cm, norm=norm1)
ax0.set_yscale("log")
ax0.set_ylim(1.e1, 4.e2) # set vertical axis limits
p2 = ax1.contour(obst[:, 0]/3600., obsDp[0, :], dNdDp.transpose(), cmap=cm, norm=norm1)
ax1.set_yscale("log")
ax1.set_ylim(1.e1, 4.e2) # set vertical axis limits

# colorbar
cax = plt.axes([0.30, 0.92, 0.40, 0.01])
cb = plt.colorbar(p1, cax = cax, ticks=[0., 6.0e4, 1.2e5, 1.8e5], format = ticker.FuncFormatter(fmt), orientation = 'horizontal')

plt.show()

res = np.zeros((obst.shape[0], obsDp.shape[1]))
# loop through times and size bins to get residuals
for it in range(obst.shape[0]):
	for isb in range(obsDp.shape[1]):
		if ((dNdD_interp[it, isb] >= 0.) and (dNdDp[it, isb] >= 0.)):
			res[it, isb] = (np.abs(dNdD_interp[it, isb]-dNdDp[it, isb]))

# denominator is the sum of number concentrations from all size bins at all times for each simulation
# results to use
dNdD_interp_sp = dNdD_interp[0:obst.shape[0], :]
den_cnt = (sum(dNdD_interp_sp[dNdD_interp_sp>0.]))
# results to use
dNdDp_sp = dNdDp[0:obst.shape[0], :]
den_cnt += (sum(dNdDp_sp[dNdDp_sp>0.]))

# sum over times and size bins to get percentage error
resT = ((sum(sum(res)))/(den_cnt))*100.

print(resT)