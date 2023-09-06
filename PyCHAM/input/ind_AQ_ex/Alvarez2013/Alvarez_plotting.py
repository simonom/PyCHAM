'''code for generating plots from the PyCHAM/PyCHAM/inputs/ind_AQ_ex/Alvarez2013 simulation,
in the style of Figure 3 of Alvarez et al. 2013: 
https://dx.doi.org/10.1073/pnas.1308310110. The purpose is to show that PyCHAM is able to 
reproduce observations of indoor air. Note that the user should adjust the file path to suit
their directory'''
# Simon O'Meara (University of Manchester, National Centre for Atmospheric Science) 2023

# import dependencies
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as si

# define class
class self:

	# path where retr_out module saved (note that by default this is saved in 'PyCHAM home directory/PyCHAM' folder
	ret_path = '/Users/user/Documents/GitHub/PyCHAM/PyCHAM'
	
	# allow access to retr_out module
	sys.path.insert(0, ret_path)

	def plotting(self): # define function
		
		# retrieving outputs ---------------------------
		# retrieve all outputs and include in self
		from retr_out import retr_out
		
		# path where gas-phase chemistry results saved	
		self.dir_path = '/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/PyCHAM/inputs/EAC23_poster_input_output/Alvarez2013/output/Alvarez2013_output'

		for prog in retr_out(self): # call on modules to get results
			
			if (isinstance(prog, str)): # check if it's a message
				mess = prog
					
				
		# finished retrieving outputs -------------------
		
		# just the relevant results
		comp_names = self.ro_obj.names_of_comp
		temper = self.ro_obj.env_cond[:, 0]
		timehr = self.ro_obj.thr
		Cfac = np.array((self.ro_obj.cfac))
		y_MW = self.ro_obj.comp_MW
		yrec = self.ro_obj.yrec[:, :]
		num_comp = self.ro_obj.nc

		# get indices of all results between hours 37 and 48
		tindx = (timehr>=37.0)*(timehr<48.)

		# number of time points
		tnum = sum(tindx)

		# prepare for holding J(HONO)[HONO] (# molecules/cm3/s) 
		# concentrations of OH (# molecules/cm3)
		res = np.zeros((tnum, 2))
		
		# get indices of OH and HONO in gas phase
		OHi = comp_names.index('OH')

		# get relevant gas-phase concentrations (ppb)
		res[:, 0] = yrec[tindx, OHi]

		# multiply by Cfac to convert concentrations from 
		# ppb to # molecules/cm3
		res[:, 0] = res[:, 0]*(Cfac[tindx].reshape(-1))

		# get the photolysis rate of HONO at these times
		# (# molecules/cm3/s)
		fname = str(self.dir_path + '/' +'HONO_rate_of_change')
		dydt = np.loadtxt(fname, delimiter = ',', skiprows = 1) # skiprows = 1 omits header
		# multiply by -1 to make rate positive, whereas results 
		# report a negative, to represent a loss term for HONO			
		res[:, 1] = dydt[tindx, dydt[0, :] == 41]*-1

		# state the units
		munit = str('(molecules' + ' cm' + u'\u207B' + u'\u00B3' + ')')
		runit = str('(molecules' + ' cm' + u'\u207B' + u'\u00B3' + ' s' + u'\u207B' + u'\u00B9' + ')')
		
				
		# setup figure
		fig = plt.figure(figsize=(14, 7))

		ax0 = fig.add_subplot(1, 1, 1)
		
		xobs = [0.0, 8.e6, 1.6e7, 2.4e7]
		yobs = [3.e5, 6.e5, 9.e5, 1.2e6]


		p3, = ax0.plot(res[:, 1], res[:, 0], '.k', label = 'sim.')
		p4, = ax0.plot(xobs, yobs, '-k', label = 'obs.')		

		ax0.set_ylabel(r'$\mathrm{[OH]}$ ' + munit, fontsize = 18)
		ax0.set_xlabel(r'$J$(HONO)[HONO] ' + munit, fontsize = 18)
		ax0.yaxis.set_tick_params(labelsize = 18, direction = 'in', which = 'both')
		ax0.xaxis.set_tick_params(labelsize = 18, direction = 'in', which = 'both')
		
		# figure title
		fig.suptitle('Figure 3 of Alvarez et al. 2013 (doi.org/10.1073/pnas.1308310110)')

		fig.tight_layout() # space out subplots
		plt.legend(fontsize=18, handles=[p3, p4] , loc=0, fancybox=True, framealpha=0.5)
		plt.show() # show figure
		
		return() # end function

# instantiate
plotter = self()
plotter.plotting() # call function