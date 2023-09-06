'''code for generating plots from the PyCHAM/PyCHAM/inputs/ind_AQ_ex/Collins2018 simulation,
in the style of Figure 4 and Figure 5 of Lunderberg et al. 2018: 
https://dx.doi.org/10.1021/acs.est.8b04512. The purpose is to show that PyCHAM is able to 
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

	def collins_plotting(self): # define function
		
		# retrieving outputs ---------------------------
		# retrieve all outputs and include in self
		from retr_out import retr_out
		
		# path where gas-phase chemistry results saved	
		self.dir_path = '/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/PyCHAM/inputs/EAC23_poster_input_output/Collins2018/output/Collins2018_allprocesses_gasphaseHONOprod'

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

		# get indices of all results between hours 2 and 3.12
		tindx = (timehr>=2.0)*(timehr<=3.12)
		

		# number of time points
		tnum = sum(tindx)

		# prepare for holding concentrations of NO2 (for gas-phase 
		# chemistry case only) and HONO (for all cases) (ppb)
		res = np.zeros((tnum, 6))
		
		# get indices of NO2 and HONO in gas phase
		NO2i = comp_names.index('NO2')
		HONOi = comp_names.index('HONO')

		# get relevant gas-phase concentrations (ppb) for gas-phase 
		# chemistry case
		res[:, 0] = yrec[tindx, NO2i]
		res[:, 1] = yrec[tindx, HONOi]

		# path where surface-phase chemistry results saved	
		self.dir_path = '/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/PyCHAM/inputs/EAC23_poster_input_output/Collins2018/output/Collins2018_allprocesses_surfacephaseHONOprod'

		for prog in retr_out(self): # call on module to get results
			
			if (isinstance(prog, str)): # check if it's a message
				mess = prog
		yrec = self.ro_obj.yrec[:, :]
		timehr1 = self.ro_obj.thr
		
		# because the other runs use the 2 hour point of the gas-phase 
		# reaction runs as the starting point, get first 0 to 1.12 hour 
		# indices for these runs
		tindx1 = (timehr1>=0.0)*(timehr1<=1.12)
	
		res[:, 2] = yrec[tindx1, HONOi]

		# path where no chemistry results saved	
		self.dir_path = '/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/PyCHAM/inputs/EAC23_poster_input_output/Collins2018/output/Collins2018_allprocesses_nochem'

		for prog in retr_out(self): # call on module to get results
			
			if (isinstance(prog, str)): # check if it's a message
				mess = prog
		yrec = self.ro_obj.yrec[:, :]	
		res[:, 3] = yrec[tindx1, HONOi]

		# path where no equilibrium results saved	
		self.dir_path = '/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/PyCHAM/inputs/EAC23_poster_input_output/Collins2018/output/Collins2018_allprocesses_noequilib'

		for prog in retr_out(self): # call on module to get results
			
			if (isinstance(prog, str)): # check if it's a message
				mess = prog
		yrec = self.ro_obj.yrec[:, :]	
		res[:, 4] = yrec[tindx1, HONOi]

		# path where air exchange only results saved	
		self.dir_path = '/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/PyCHAM/inputs/EAC23_poster_input_output/Collins2018/output/Collins2018_airexchangeonly'

		for prog in retr_out(self): # call on module to get results
			
			if (isinstance(prog, str)): # check if it's a message
				mess = prog
		yrec = self.ro_obj.yrec[:, :]	
		res[:, 5] = yrec[tindx1, HONOi]

		munit = '(ppb)'

				
		# setup figure
		fig = plt.figure(figsize=(14, 7))

		ax0 = fig.add_subplot(1, 1, 1)
		par1 = ax0.twinx() # first parasite axis


		xobs = [0., 1.e3, 2.e3, 3.e3, 4.e3]
		yobs = [35., 15., 10., 7., 5.]

		p4, = ax0.plot(xobs, yobs, marker = 's', color = 'orange', label = 'obs. NO2')



		p3, = ax0.plot((timehr[tindx]-2.)*3.6e3, res[:, 0], '-k', label = 'NO2 gas-phase chem.')
		
		xobs = [0., 1.e3, 2.e3, 3.e3, 4.e3]
		yobs = [31., 13., 10., 9., 8.]

		p10, = par1.plot(xobs, yobs, 'g', marker = 'o', label = 'obs. HONO')	


		p5, = par1.plot((timehr[tindx]-2.)*3.6e3, res[:, 1], '-k', linewidth=3.0, label = 'HONO gas-phase chem.')		
		p6, = par1.plot((timehr1[tindx1])*3.6e3, res[:, 2], '-m', linewidth=3.0, label = 'HONO surface-phase chem.')
		p7, = par1.plot((timehr1[tindx1])*3.6e3, res[:, 3], '--k', linewidth=1.0, label = 'HONO no chem.')			
		p8, = par1.plot((timehr1[tindx1])*3.6e3, res[:, 4], '-c', linewidth=1.0, label = 'HONO no equilib.')
		p9, = par1.plot((timehr1[tindx1])*3.6e3, res[:, 5], '--c', linewidth=1.0, label = 'HONO air exchange only')

		ax0.set_ylabel(r'$\mathrm{NO_{2}}$ ' + munit, fontsize = 18)
		ax0.set_xlabel('Elapsed Time (s)', fontsize = 18)
		ax0.yaxis.set_tick_params(labelsize = 18, direction = 'out', which = 'both')
		ax0.xaxis.set_tick_params(labelsize = 18, direction = 'out', which = 'both')
		par1.set_ylabel(r'$\mathrm{HONO}$ ' + munit, rotation = 270, fontsize = 18, labelpad=20)		
		par1.yaxis.set_tick_params(labelsize = 18, direction = 'out', which = 'both')

		# figure title
		fig.suptitle('Figure 4 of Collins et al. 2018 (doi.org/10.1021/acs.est.8b04512)')

		fig.tight_layout() # space out subplots
		plt.legend(fontsize=18, handles=[p4, p3, p10, p5, p6, p7, p8, p9] , loc=1, fancybox=True, framealpha=0.5)
		plt.show() # show figure
		
		return() # end function

# instantiate
plotter = self()
plotter.collins_plotting() # call function