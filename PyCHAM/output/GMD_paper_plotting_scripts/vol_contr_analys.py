'''script for plotting volatility basis set mass fraction of particle phase with time and tabulating component volatilities and particle-phase concentrations'''
# aids interpretation of gas-particle partitioning results, 
# assumes calling from the PyCHAM home directory. Note that path to results
# has to be manually set below

# import required modules
import os
import sys
# ensure modules can be seen 
# (assumes calling from the home folder)
sys.path.append(str(os.getcwd() + '/PyCHAM'))
import numpy as np
import pybel
import xml_interr
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import scipy.constants as si
import retr_out

# ----------------------------------------------------------------------------------------
# prepare plot
fig, (ax1) = plt.subplots(1, 1, figsize=(10,7))
fig.subplots_adjust(hspace = 0.7)

# ----------------------------------------------------------------------------------------
# Prepare the volatility basis set interpretation
# of particle-phase concentrations

cwd = os.getcwd() # get current working directory
output_by_sim = str(cwd + '/PyCHAM/output/isoprene_scheme_MCM/Gaswall_sens_Cw0_kw1e2')

# required outputs from full-moving
(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, comp_names, 
		y_mw, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)

print(((63.4*Cfac[0]*1.e6)/si.N_A)*(68.12*1e6))
import ipdb; ipdb.set_trace()
# number of particle size bins without wall
num_asb = (num_sb-wall_on)

# get all SMILE strings from the xml file
# interrogate xml to list component names and SMILES
[comp_smil, comp_name] = xml_interr.xml_interr(str(cwd + '/PyCHAM/input/example_xml.xml'))

# convert chemical scheme component names into SMILEs
comp_smiles = [] # holder
for name in comp_names[0:-2]: # omit H20 and core at end of comp_names
	comp_smiles.append(comp_smil[comp_name.index(name)])

SOA0 = 0.
for i in range(1, (num_asb+1)):
	# calculate SOA (*1.0E-12 to convert from g/cc (air) to ug/m3 (air))
	# excluding water and core
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

print(comp_names.index('C536OOH'))
import ipdb; ipdb.set_trace()

# convert from list to array
y_mw = np.array((y_mw))

# convert vapour pressures in log10(atm) to saturation concentrations in ug/m3
# using eq. 1 of O'Meara et al. 2014
Psat_Cst = (1.e6*y_mw[0:-2].reshape(1, -1))*(10**(Psat))/(8.2057e-5*TEMP)

# saturation vapour pressure in Pa
Psat_Pa = (10**(Psat))*101325.

# convert particle-phase concentrations from molecules/cc to ug/m3 to 
# be comparable with total secondary particulate matter concentration
# isolate just particle-phase concentrations (molecules/cc)
# note just one size bin in this simulation and want to exclude water 
# and core
SPMi = y[:, num_comp:(num_comp*(num_asb+1)-2)]
# convert to mol/cc
SPMi = SPMi/si.N_A
# convert to ug/m3
SPMi = SPMi*(y_mw[0:-2].reshape(1, -1)*1.e12)

# the saturation concentrations to consider (log10(C* (ug/m3)))
# note these will be values at the centre of the volatility size bins
sc = np.arange(-2.5, 7.5, 1.)

# empty array for normalised mass contributions
nmc = np.zeros((len(sc), len(t_array)))

for it in range(len(t_array)): # loop through times
	# loop through saturation concentrations and find normalised mass contributions
	# to secondary particulate matter
	for i in range(len(sc)):

		if (i == 0):
			indx = Psat_Cst < 10**(sc[i]+0.5)
		if (i > 0 and i < len(sc)-1):
			indx = (Psat_Cst >= 10**(sc[i-1]+0.5))*(Psat_Cst < 10**(sc[i]+0.5))
		if (i == len(sc)-1):
			indx = (Psat_Cst >= 10**(sc[i-1]+0.5))

		# repeat over size bins
		indx = np.repeat(indx, num_asb, axis = 1)
		
	
		if (SOA0[it] > 0.):
			nmc[i, it] = (SPMi[it, indx[0, :]].sum())/SOA0[it]
			#print(i, it, SPMi[it, indx[0, :]].sum(), SOA0[it], nmc[i, it])
			#import ipdb;ipdb.set_trace()

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

# show figure
plt.show()

# ----------------------------------------------------------------------------------------
# next create table with components in columns along with their saturation vapour pressure
# and the particle-phase concentration as a function of time in rows

# particle-phase concentrations (molecules/cc)
y_mat = y[:, num_comp:(num_comp*(num_asb+1))]
# convert to ug/m3
y_mat = ((y_mat/si.N_A)*1.e6)*y_mw

# prepare header for concentrations with time file
# note header element begins with here, as
# first column used for times
y_header = str('time (hours)/component name')
# empty matrix for holding times and concentrations, excluding water and core
y_matn = np.zeros((y_mat.shape[0], 1))

for i in range(num_asb):

	end = str('_p' + str(i+1))

	# exclude water and core
	y_matn = np.append(y_matn, y_mat[:, i*(num_comp):(i+1)*(num_comp-2)], axis = 1)	

	for ii in range(num_comp-2): # exclude water and core
		
		start = ', ' 
		y_header = str(y_header+str(start+comp_names[ii])+end)

# insert times (hours)
y_matn[:, 0] = t_array
# add row to hold saturation vapour pressures per component (Pa)
y_matn = np.append(np.zeros((1, y_matn.shape[1])), y_matn, axis = 0)
y_matn[0, 1::] = np.repeat(Psat_Pa.reshape(1, -1), num_asb, axis = 1) 

# saving both gas- and particle-phase concentrations of species
np.savetxt(os.path.join(output_by_sim, 'particle_phase_concentrations_and_volatility_for_all_components_at_all_times'), y_matn, delimiter=',', header=str('time (given in first column (hours)) changes with rows, components in columns, with _pi representing particle phase where i is the size bin number (starting at 1) (ug/m3 (air)); saturation vapour pressures (Pa) at 298.15 K given in first row\n'+y_header))
