'''code for generating plots from the PyCHAM/PyCHAM/inputs/ind_AQ_ex/Tran2017 simulation,
in the style of Figure 3 of Tran et al. 2017: 
https://doi.org/10.1177/1420326X15610798. The purpose is to show that PyCHAM is able to 
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
		self.dir_path = '/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/PyCHAM/inputs/EAC23_poster_input_output/Tran2017/output/Tran2017_output'

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
		num_sb = self.ro_obj.nsb
		rbou_rec = np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
		rbou_rec[:, :] = self.ro_obj.rad[:, :]

		# number of actual particle size bins
		num_asb = (num_sb-self.ro_obj.wf)

		import scipy.constants as si
		# particle-phase concentrations of all components (# molecules/cm3)
		if (self.ro_obj.wf > 0): # wall on
			ppc = yrec[:, num_comp:-num_comp*self.ro_obj.wf]
		if (self.ro_obj.wf == 0): # wall off
			ppc = yrec[:, num_comp::]

		# tile molar weights over size bins and times
		y_mwt = np.tile(np.array((self.ro_obj.comp_MW)).reshape(1, -1), (1, num_asb))
		y_mwt = np.tile(y_mwt, (ppc.shape[0], 1))
		# convert from # molecules/cm3 to ug/m3
		ppc = (ppc/si.N_A)*y_mwt*1.e12

		# get index of size bins with diameters between 2 um and 10 um
		coarse_index = (rbou_rec[:, 1::]*2.>2.)*(rbou_rec[:, 1::]*2.<=10.)

		# repeat index over components
		coarse_index = np.repeat(coarse_index, num_comp, axis=1)

		# zero components outside the coarse size interval
		ppc = ppc*coarse_index

		# sum over components and size bins (ug/m3)
		ppc = np.sum(ppc, axis=1)

		# get indices of all results between hours 0 and 10
		tindx = (timehr>=0.0)*(timehr<=10.)

		# state the units
		munit = str('(' + u'\u03BC' + 'g/m' + u'\u00B3' + ')')
		
				
		# setup figure
		fig = plt.figure(figsize=(14, 7))

		ax0 = fig.add_subplot(1, 1, 1)

		xobs = [120., 220., 320., 420., 520., 580.]
		yobs = [168., 45., 20., 10., 5., 4.]
		
		p4, = ax0.plot(xobs, yobs, 'sk', label = '$\mathrm{Obs. Indoor\, PM_{2-10}}$')

		p3, = ax0.plot(timehr[tindx]*60., ppc[tindx], '.k', label = '$\mathrm{Sim. Indoor\, PM_{2-10}}$')
		p5, = ax0.plot(timehr[tindx]*60., np.ones(51)*8., '-k', label = '$\mathrm{Obs. Outdoor\, PM_{2-10}}$')
		
		ax0.set_ylabel(r'$\mathrm{PM_{2-10} concentrations}$ ' + munit, fontsize = 18)
		ax0.set_xlabel(r'Time (minute)', fontsize = 18)
		ax0.yaxis.set_tick_params(labelsize = 18, direction = 'in', which = 'both')
		ax0.xaxis.set_tick_params(labelsize = 18, direction = 'in', which = 'both')
		
		# figure title
		fig.suptitle('Figure 3 of Tran et al. 2013 (doi.org/10.1177/1420326X15610798)')

		fig.tight_layout() # space out subplots
		plt.legend(fontsize=18, handles=[p4, p3, p5] , loc=1, fancybox=True, framealpha=0.5)
		plt.show() # show figure
		
		return() # end function

# instantiate
plotter = self()
plotter.plotting() # call function