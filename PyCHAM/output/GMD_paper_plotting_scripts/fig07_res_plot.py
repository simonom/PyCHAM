'''Code to plot results for gas-wall partitioning sensitivity'''
# aim is to illustrate the effect of varying the mass transfer coefficient for gas-wall
# partitioning and the effective absorbing mass concentration of the wall on SOA mass
# concentration estimates

# all runs for first and second plot used the isoprene MCM sheme (fig07ab_scheme.txt), whilst 
# the final plot uses fig07c_scheme.txt
# PyCHAM inputs are saved as: fig07ab_mod_var.txt and fig07c_mod_var.txt
# note that in order to produce plot b, umansysprop, xml_interr.py and 
# example_xml.xml need to be available 

# import required modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import scipy.constants as si
import os
import sys
# ensure modules can be seen 
# (assumes calling from the home folder, but will also work when in the GMD paper Results folder)
sys.path.append(str(os.getcwd() + '/PyCHAM'))
import retr_out


# ----------------------------------------------------------------------------------------
# function for doing colorbar tick labels in standard notation
def fmt(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	return r'${} \times 10^{{{}}}$'.format(a, b)

# ----------------------------------------------------------------------------------------
# prepare plot
fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(10,7))
fig.subplots_adjust(hspace = 0.7)


# ----------------------------------------------------------------------------------------
#import result files

# first case is no gas-wall partitioning, where kw set to 
# 1.e2 /s and Cw was set to 0. g/m3 (air)

cwd = os.getcwd() # get current working directory
try: # when calling from PyCHAM home directory
	output_by_sim = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig07_data/Gaswall_sens_Cw0_kw1e2')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

except: # when calling from GMD paper Results folder
	output_by_sim = str(cwd + '/fig07_data/Gaswall_sens_Cw0_kw1e2')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

SOA0 = 0.
for i in range(1, num_sb):
	# calculate SOA (*1.0E-12 to convert from g/cc (air) to ug/m3 (air))
	SOA0 += (((y[:, num_comp*i:num_comp*(i+1)-2]/si.N_A)*y_mw[0:-2]*1.0e12).sum(axis=1))
			
# ----------------------------------------------------------------------------------------
# second case is moderate gas-wall partitioning, where kw was set to 
# 1.0e-3 /s and Cw set to 1.e-4 g/m3 (air), as these values are comparable to those
# for particulate matter (for particulate matter, total concentration (core, water and
# SOA summed) was around 7.85e-6 g/m3 (air))

cwd = os.getcwd() # get current working directory

try: # when calling from PyCHAM home directory
	output_by_sim = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig07_data/Gaswall_sens_Cw1e-4_kw1e-3')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

except: # when calling from GMD paper Results folder
	output_by_sim = str(cwd + '/fig07_data/Gaswall_sens_Cw1e-4_kw1e-3')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)


# required outputs from full-moving
(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

SOA1 = 0.
for i in range(1, num_sb):
	# calculate SOA (*1.0E-12 to convert from g/cc (air) to ug/m3 (air))
	SOA1 += (((y[:,num_comp*i:num_comp*(i+1)-2]/si.N_A)*y_mw[0:-2]*1.0e12).sum(axis=1))

# ----------------------------------------------------------------------------------------
# third case is high gas-wall partitioning due to high Cw, where kw set to 
# 1.0e-3 /s and Cw set to 1.e-2 g/m3 (air)
cwd = os.getcwd() # get current working directory

try: # when calling from PyCHAM home directory
	output_by_sim = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig07_data/Gaswall_sens_Cw1e-2_kw1e-3')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

except: # when calling from GMD paper Results folder
	output_by_sim = str(cwd + '/fig07_data/Gaswall_sens_Cw1e-2_kw1e-3')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)


# required outputs from full-moving
(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

SOA2 = 0.
for i in range(1, num_sb):
	# calculate SOA (*1.0E-12 to convert from g/cc (air) to ug/m3 (air))
	SOA2 += (((y[:,num_comp*i:num_comp*(i+1)-2]/si.N_A)*y_mw[0:-2]*1.0e12).sum(axis=1))

# ----------------------------------------------------------------------------------------
# fourth case is high gas-wall partitioning due to high kw, where kw set to 
# 1.0e2 /s and Cw was set to 1.e-2 g/m3 (air)

cwd = os.getcwd() # get current working directory

try: # when calling from PyCHAM home directory
	output_by_sim = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig07_data/Gaswall_sens_Cw1e-4_kw1e2')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

except: # when calling from GMD paper Results folder
	output_by_sim = str(cwd + '/fig07_data/Gaswall_sens_Cw1e-4_kw1e2')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)


# required outputs from full-moving
(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

SOA3 = 0.
for i in range(1, num_sb):
	# calculate SOA (*1.0E-12 to convert from g/cc (air) to ug/m3 (air))
	SOA3 += (((y[:,num_comp*i:num_comp*(i+1)-2]/si.N_A)*y_mw[0:-2]*1.0e12).sum(axis=1))

# ----------------------------------------------------------------------------------------

ax0.plot(t_array, SOA0, linewidth=3, label='No wall')
ax0.plot(t_array, SOA1, linewidth=3, label=r'$k_{w}=1\rm{x}10^{-3}\, \mathrm{s^{-1}}$, $C_{w}=1\rm{x}10^{-4}\, \mathrm{\mu g \, m^{-3}}$')
ax0.plot(t_array, SOA2, linewidth=3, label=r'$k_{w}=1\rm{x}10^{-3}\, \mathrm{s^{-1}}$, $C_{w}=1\rm{x}10^{-2}\, \mathrm{\mu g \, m^{-3}}$')
ax0.plot(t_array, SOA3, '--', linewidth=3, label=r'$k_{w}=1\rm{x}10^{2}\, \mathrm{s^{-1}}$, $C_{w}=1\rm{x}10^{-4}\, \mathrm{\mu g \, m^{-3}}$')
ax0.set_xlabel(r'Time (hours of day)', fontsize=12)
ax0.set_ylabel(r'[SPM] ($\mathrm{\mu g\, m^{-3}}$)', fontsize=12)
ax0.yaxis.set_tick_params(direction = 'in')
ax0.xaxis.set_tick_params(direction = 'in')
ax0.legend(fontsize=10)
ax0.text(x=-0.8, y=9.2, s='(a)', size=12)

# ----------------------------------------------------------------------------------------
# for the next plot we exemplify displaying the volatility basis set interpretation
# of particle-phase concentrations
import pybel
import xml_interr

# use the no gas-wall partitioning case, where kw set to 
# 1.e2 /s and Cw was set to 0. g/m3 (air)

cwd = os.getcwd() # get current working directory

try: # when calling from PyCHAM home directory
	output_by_sim = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig07_data/Gaswall_sens_Cw0_kw1e2')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

except: # when calling from GMD paper Results folder
	output_by_sim = str(cwd + '/fig07_data/Gaswall_sens_Cw0_kw1e2')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)


# required outputs from full-moving
(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

# get all SMILE strings from the xml file
# interrogate xml to list component names and SMILES
try: # when calling from the PyCHAM home directory
	[comp_smil, comp_name] = xml_interr.xml_interr(str(cwd + '/PyCHAM/input/example_xml.xml'))
except: # when calling from the GMD paper Results folder
	[comp_smil, comp_name] = xml_interr.xml_interr(str(cwd + '/example_xml.xml'))

# convert chemical scheme component names into SMILEs
comp_smiles = [] # holder
for name in comp_names[0:-2]: # omit H20 and core at end of comp_names
	comp_smiles.append(comp_smil[comp_name.index(name)])

SOA0 = 0.
for i in range(1, num_sb):
	# calculate SOA (*1.0E-12 to convert from g/cc (air) to ug/m3 (air))
	SOA0 += (((y[:, num_comp*i:num_comp*(i+1)-2]/si.N_A)*y_mw[0:-2]*1.0e12).sum(axis=1))

Pybel_objects = [] # holder for Pybel object names
for i in range(num_comp-2): # component loop
	# generate pybel object
	Pybel_objects.append(pybel.readstring('smi', comp_smiles[i]))

# point to umansysprop folder
sys.path.insert(1, (cwd + '/umansysprop')) # address for updated version
	
from umansysprop import boiling_points
from umansysprop import vapour_pressures
from umansysprop import liquid_densities

NA = si.Avogadro # Avogadro's number (molecules/mol)
# vapour pressures of components, excluding water and core at end
Psat = np.zeros((1, num_comp-2))
TEMP = 298.15 # temperature (K) 


for i in range(num_comp-2): # component loop
	# vapour pressure (log10(atm)) (# eq. 6 of Nannoolal et al. (2008), with dB of 
	# that equation given by eq. 7 of same reference)
	Psat[0, i] = ((vapour_pressures.nannoolal(Pybel_objects[i], TEMP, 
			boiling_points.nannoolal(Pybel_objects[i]))))

# convert from list to array
y_mw = np.array((y_mw))

# convert vapour pressures in log10(atm) to saturation concentrations in ug/m3
# using eq. 1 of O'Meara et al. 2014
Psat = 1.e6*y_mw[0:-2].reshape(1, -1)*(10**(Psat))/(8.2057e-5*TEMP)

# convert particle-phase concentrations from molecules/cc to ug/m3 to 
# be comparable with total secondary particulate matter concentration
# isolate just particle-phase concentrations (molecules/cc)
# note just one size bin in this simulation and want to exclude water 
# and core
SPMi = y[:, num_comp:(num_comp*(num_sb)-2)]
# convert to mol/cc
SPMi = SPMi/si.N_A
# convert to ug/m3
SPMi = SPMi*(y_mw[0:-2].reshape(1, -1)*1.e12)

# the saturation concentrations to consider (log10(C* (ug/m3)))
# note these will be values at the centre of the volatility size bins
sc = np.arange(10)-2.5

# empty array for normalised mass contributions
nmc = np.zeros((len(sc), len(t_array)))
# loop through saturation concentrations and find normalised mass contributions
# to secondary particulate matter
for i in range(len(sc)):
	if (i == 0):
		indx = Psat < 10**(sc[i]+0.5)
	if (i > 0 and i < len(sc)-2):
		indx = (Psat >= 10**(sc[i-1]+0.5))*(Psat <10**(sc[i]+0.5))
	if (i == len(sc)-1):
		indx = Psat >= 10**(sc[i-1]+0.5)
	
	for it in range(len(t_array)): # loop over times
		if (SOA0[it] > 0.):
			nmc[i, it] = (SPMi[it, indx[0, :]].sum())/SOA0[it]

# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
colors = [(0.60, 0.0, 0.70), (0, 0, 1), (0, 1.0, 1.0), (0, 1.0, 0.0), (1.0, 1.0, 0.0), (1.0, 0.0, 0.0)]  # R -> G -> B
n_bin = 100  # discretizes the colormap interpolation into bins
cmap_name = 'my_list'
# Create the colormap
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N = n_bin)
	
# set contour levels
levels = (MaxNLocator(nbins = 100).tick_values(np.min(nmc), np.max(nmc)))
	
# associate colours and contour levels
norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)


ptindx = (SOA0 > 0) # indices of times where secondary material present
p0 = ax1.pcolormesh(t_array[ptindx], sc, nmc[:, ptindx], cmap=cm, norm=norm1)

cax = plt.axes([0.875, 0.40, 0.02, 0.18]) # specify colour bar position
cb = plt.colorbar(p0, cax = cax, ticks=[0.00, 0.25, 0.50, 0.75, 1.00], orientation = 'vertical')
cb.ax.tick_params(labelsize = 12)
cb.set_label('mass fraction', size = 12, rotation = 270, labelpad = 10.)

ax1.set_xlabel(r'Time (hours of day)', fontsize=12)
ax1.set_ylabel(r'$\rm{log_{10}(}$$C*_{\mathrm{298.15 K}}$$\rm{\, (\mu g\, m^{-3}))}$', fontsize=12, labelpad = 10.)
ax1.yaxis.set_tick_params(direction = 'in', which = 'both')
ax1.xaxis.set_tick_params(direction = 'in', which = 'both')
ax1.text(x=-0.8, y=7.5, s='(b)', size=12)
ax1.set_xlim(-0.6, 12.6)
ax1.set_yticks([-2, 0, 2, 4, 6])
ax1.set_yticklabels(['<-2', '0', '2', '4', '$\geq$6'])

# ----------------------------------------------------------------------------------------
# for the next plot, we show decay of individual components over time in a dark experiment with
# no lights or particles and no chemical interaction between the components of interest
# in the results below the fig07c_scheme chemical scheme was used, and components 
# with an initial concentration of 50 ppb were (MCM names): MGLYOX (methylglyoxal) and 
# two-methylglyceric_acid, the latter had to be added manually to the xml file, but its
# production during isoprene ozonolysis is confirmed here: doi.org/10.1021/jp061734m
# we want to plot concentrations of all components against time to illustrate how to 
# constrain wall loss of vapour parameters 


# ----------------------------------------------------------------------------------------
# moderate Kw (1e-1 /s) and moderate Cw (7e-5 ug/m3)
# open saved files
cwd = os.getcwd() # get current working directory

try: # when calling from PyCHAM home directory
	output_by_sim = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig07_data/Gaswall_sens_Cw7e-5_kw1e-1_2comp')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

except: # when calling from GMD paper Results folder
	output_by_sim = str(cwd + '/fig07_data/Gaswall_sens_Cw7e-5_kw1e-1_2comp')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)


# required outputs from full-moving
(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

comp_index = [] # holder for component indices
# get index of three components of interest
for i in range(len(comp_names)):
	if comp_names[i] == 'MGLYOX':
		comp_index.append(i)
	if comp_names[i] == 'two-methylglyceric_acid':
		comp_index.append(i)
	

ax2.plot(t_array[0:45]*3600., y[0:45, comp_index[1]], '+g', label= r'$k_{w}=1\rm{x}10^{-1}\, \mathrm{s^{-1}}$, $C_{w}=7\rm{x}10^{1}\, \mathrm{\mu g \, m^{-3}}$')

# ----------------------------------------------------------------------------------------
# moderate Kw (1e-1 /s) high Cw (7e-2 ug/m3)
# open saved files
cwd = os.getcwd() # get current working directory

try: # when calling from PyCHAM home directory
	output_by_sim = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig07_data/Gaswall_sens_Cw7e-2_kw1e-1_2comp')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

except: # when calling from GMD paper Results folder
	output_by_sim = str(cwd + '/fig07_data/Gaswall_sens_Cw7e-2_kw1e-1_2comp')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)



# required outputs from full-moving
(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

comp_index = [] # holder for component indices
# get index of three components of interest
for i in range(len(comp_names)):
	if comp_names[i] == 'MGLYOX':
		comp_index.append(i)
	if comp_names[i] == 'two-methylglyceric_acid':
		comp_index.append(i)


ax2.plot(t_array[0:45]*3600., y[0:45, comp_index[1]], '--g', label= r'$k_{w}=1\rm{x}10^{-1}\, \mathrm{s^{-1}}$, $C_{w}=7\rm{x}10^{4}\, \mathrm{\mu g \, m^{-3}}$')


# ----------------------------------------------------------------------------------------
# fourth case is Kw (1e-2 /s) moderate Cw (7e-5 ug/m3)
# open saved files
cwd = os.getcwd() # get current working directory

try: # when calling from PyCHAM home directory
	output_by_sim = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig07_data/Gaswall_sens_Cw7e-5_kw1e-2_2comp')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

except: # when calling from GMD paper Results folder
	output_by_sim = str(cwd + '/fig07_data/Gaswall_sens_Cw7e-5_kw1e-2_2comp')
	# required outputs
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

comp_index = [] # holder for component indices
# get index of three components of interest
for i in range(len(comp_names)):
	if comp_names[i] == 'MGLYOX':
		comp_index.append(i)
	if comp_names[i] == 'two-methylglyceric_acid':
		comp_index.append(i)
	
ax2.plot(t_array[0:45]*3600., y[0:45, comp_index[1]], '.g', label= r'$k_{w}=1\rm{x}10^{-2}\, \rm{s^{-1}}$, $C_{w}=7\rm{x}10^{1}\, \rm{\mu g \, m^{-3}}$')

ax2.set_xlabel(r'Time (seconds into simulation)', fontsize=12)
ax2.set_ylabel(r'[2-methylglyceric acid] (ppb)', fontsize=12)
ax2.yaxis.set_tick_params(direction = 'in', which = 'both')
ax2.xaxis.set_tick_params(direction = 'in', which = 'both')
ax2.legend(fontsize=10)
ax2.text(x=-3.0, y=57.0, s='(c)', size=12)
plt.subplots_adjust(hspace=0.34)

# show to screen and save
plt.show()
fig.savefig('fig07.png')
