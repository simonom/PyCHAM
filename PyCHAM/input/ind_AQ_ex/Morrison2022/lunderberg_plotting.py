'''code for generating plots from the PyCHAM/PyCHAM/inputs/ind_AQ_ex/Lunderberg2020 simulation,
in the style of Figure 2 and Figure 3 of Lunderberg et al. 2020: 
https://dx.doi.org/10.1021/acs.est.0c00966. The purpose is to show that PyCHAM is able to 
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

	# path where results saved	
	dir_path = 'C:\\Users\\Psymo\\Desktop\\PyCHAM\\PyCHAM\\PyCHAM\\output\\Lunderberg_scheme\\Lundenberg2020'

	# path where retr_out module saved
	ret_path = 'C:\\Users\\Psymo\\Desktop\\PyCHAM\\PyCHAM\\PyCHAM'
	
	# allow access to retr_out module
	sys.path.insert(0, ret_path)

	def lunderberg_plotting(self): # define function
		
		# retrieving outputs ---------------------------
		# retrieve all outputs and include in self
		from retr_out import retr_out
			
		for prog in retr_out(self): # call on modules to solve problem
			
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

		# prepare for holding concentrations (ug/m3)
		res = np.zeros((10, 3))
		
		# prepare for holding gradients (ug/m3/K)
		m = np.zeros((10))

		tempers = [288., 290., 292.] # temperatues to consider (K)

		# get gas+particle mass concentrations of C13-C22 ready for abundance vs. temperature plot

		compcnt = 0 # count on components

		# just gas-phase concentrations first (ppb)
		for Cnum in ['C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22']:
			indx = comp_names.index(Cnum) # get index
			
			tcnt = 0 # count on temperatures
			for tempern in tempers: # loop through temperatures
				last_time_of_this_temp = timehr[temper == tempern][-1]
		
				# concentration (# molecules/cm3)
				res[compcnt, tcnt] = yrec[timehr == last_time_of_this_temp, indx]*Cfac[timehr == last_time_of_this_temp]

				# add particle-phase concentration (# molecules/cm3)
				res[compcnt, tcnt] += yrec[timehr == last_time_of_this_temp, indx+num_comp]

				# convert to ug/m3
				res[compcnt, tcnt] = (res[compcnt, tcnt]/si.N_A)*y_MW[indx]*1.e12

				tcnt += 1 # count on temperatures

			# estimate gradient (ug/m3/K)
			m[compcnt] = (res[compcnt, -1] - res[compcnt, 0])/(tempers[-1] - tempers[0]) 

			compcnt += 1 # count on components

		# unicode for units
		munit = str('(' + u'\u03BC' + 'g m' + u'\u207B' + u'\u00B3' + ' K' + u'\u207B' + u'\u00B9' + ')')
		cunit = str('(' + u'\u03BC' + 'g m' + u'\u207B' + u'\u00B3' + ')')
		
		# setup figure
		fig = plt.figure(figsize=(14, 7))

		ax0 = fig.add_subplot(5, 3, 1)
		ax0.plot(tempers, res[0, :], 'x-')
		ax0.set_title('C13 bin')
		ax0.text(288., 12.5, str('m = ' + str(round(m[0], 1))))

		ax1 = fig.add_subplot(5, 3, 2)
		ax1.plot(tempers, res[1, :], 'x-')
		ax1.set_title('C14 bin')
		ax1.text(288., 9., str('m = ' + str(round(m[1], 1))))

		ax3 = fig.add_subplot(5, 3, 4)
		ax3.plot(tempers, res[2, :], 'x-')
		ax3.set_title('C15 bin')
		ax3.text(288., 6.5, str('m = ' + str(round(m[2], 1))))

		ax4 = fig.add_subplot(5, 3, 5)
		ax4.plot(tempers, res[3, :], 'x-')
		ax4.set_title('C16 bin')
		ax4.text(288., 5., str('m = ' + str(round(m[3], 1))))

		ax2 = fig.add_subplot(1, 3, 3)
		ax2.plot(range(13, 23), m, 'x-')
		ax2.set_ylabel(r'Fit Slope ' + munit, fontsize = 14)
		ax2.set_xlabel('Alkane-Equivalent Volatility Bin', fontsize = 14)
		ax2.yaxis.set_tick_params(labelsize = 14, direction = 'out', which = 'both')
		ax2.xaxis.set_tick_params(labelsize = 14, direction = 'out', which = 'both')
		ax2.text(15., 1., str('m = ' + str(round(((m[-1]-m[0])/(23.-13.)), 1))))

		ax5 = fig.add_subplot(5, 3, 7)
		ax5.plot(tempers, res[4, :], 'x-')
		ax5.set_title('C17 bin')
		ax5.set_ylabel(r'SVOC Concentration ' + cunit, fontsize = 14)
		ax5.text(288., 3., str('m = ' + str(round(m[4], 1))))

		ax6 = fig.add_subplot(5, 3, 8)
		ax6.plot(tempers, res[5, :], 'x-')
		ax6.set_title('C18 bin')
		ax6.text(288., 3., str('m = ' + str(round(m[5], 1))))

		ax7 = fig.add_subplot(5, 3, 10)
		ax7.plot(tempers, res[6, :], 'x-')
		ax7.set_title('C19 bin')
		ax7.text(288., 2.5, str('m = ' + str(round(m[6], 1))))

		ax8 = fig.add_subplot(5, 3, 11)
		ax8.plot(tempers, res[7, :], 'x-')
		ax8.set_title('C20 bin')
		ax8.text(288., 1.5, str('m = ' + str(round(m[7], 1))))

		ax9 = fig.add_subplot(5, 3, 13)
		ax9.plot(tempers, res[8, :], 'x-')
		ax9.set_title('C21 bin')
		ax9.set_xlabel('Temperature (K)', fontsize = 14)
		ax9.text(288., 0.75, str('m = ' + str(round(m[8], 1))))

		ax10 = fig.add_subplot(5, 3, 14)
		ax10.plot(tempers, res[9, :], 'x-')
		ax10.set_title('C22 bin')
		ax10.set_xlabel('Temperature (K)', fontsize = 14)
		ax10.text(288., 0.4, str('m = ' + str(round(m[9], 1))))		

		# figure title
		fig.suptitle('Figure 2 of Lunderberg et al. 2020 (doi.org/10.1021/acs.est.0c00966)')

		fig.tight_layout() # space out subplots
		plt.show() # show figure
		

		# prepare for holding concentrations (ug/m3)
		res = np.zeros((10, 3))
		
		# prepare for holding gradients (ug/m3/K)
		m = np.zeros((10))
		# get gas+particle mass concentrations of C24-C31 ready for abundance vs. [PM2.5] plot

		compcnt = 0 # count on components

		# components of interest
		for Cnum in ['C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30', 'C31']:
			indx = comp_names.index(Cnum) # get index
			
			tcnt = 0 # count on [PM2.5]s

			PM2p5 = np.zeros((3)) # prepare to hold particle mass concentrations
			# loop through PM2.5 concentrations at 292 K - note that the 
			# PM2.5 concentration varies every 5 hours in the final 15 hours of simulation
			for last_time_of_this_temp in [15., 20., 25.]:
		
				# get the PM2.5 concentration (ug/m3)
				PM2p5[tcnt] = np.sum(((yrec[timehr == last_time_of_this_temp, num_comp:num_comp*2]/si.N_A)*y_MW)*1.e12)

				# gas-phase concentration (# molecules/cm3)
				res[compcnt, tcnt] = yrec[timehr == last_time_of_this_temp, indx]*Cfac[timehr == last_time_of_this_temp]

				# add particle-phase concentration (# molecules/cm3)
				res[compcnt, tcnt] += yrec[timehr == last_time_of_this_temp, indx+num_comp]

				# convert to ug/m3
				res[compcnt, tcnt] = (res[compcnt, tcnt]/si.N_A)*y_MW[indx]*1.e12

				tcnt += 1 # count on [PM2.5]s

			# estimate gradient (ug/m3/K)
			m[compcnt] = (res[compcnt, -1] - res[compcnt, 0])/(PM2p5[-1] - PM2p5[0]) 

			compcnt += 1 # count on components

		# setup figure
		fig2 = plt.figure(figsize = (14, 7))

		ax0 = fig2.add_subplot(2, 4, 1)
		ax0.plot(PM2p5, res[0, :], 'x-')
		ax0.set_title('C24 bin')
		ax0.set_ylim((0., 0.15))
		ax0.text(1., 0.13, str('m = ' + str(round(m[0], 3))))

		ax1 = fig2.add_subplot(2, 4, 2)
		ax1.plot(PM2p5, res[1, :], 'x-')
		ax1.set_title('C25 bin')
		ax1.set_ylim((0., 0.15))
		ax1.text(1., 0.13, str('m = ' + str(round(m[1], 3))))

		ax2 = fig2.add_subplot(2, 4, 3)
		ax2.plot(PM2p5, res[2, :], 'x-')
		ax2.set_title('C26 bin')
		ax2.set_ylim((0., 0.15))
		ax2.text(1., 0.13, str('m = ' + str(round(m[2], 3))))

		ax3 = fig2.add_subplot(2, 4, 4)
		ax3.plot(PM2p5, res[3, :], 'x-')
		ax3.set_title('C27 bin')
		ax3.set_ylim((0., 0.15))
		ax3.text(1., 0.13, str('m = ' + str(round(m[3], 3))))

		ax4 = fig2.add_subplot(2, 4, 5)
		ax4.plot(PM2p5, res[4, :], 'x-')
		ax4.set_title('C28 bin')
		ax4.set_ylim((0., 0.15))
		ax4.text(1., 0.13, str('m = ' + str(round(m[4], 3))))
		ax4.set_ylabel(r'SVOC Concentration ' + cunit, fontsize = 14)		

		ax5 = fig2.add_subplot(2, 4, 6)
		ax5.plot(PM2p5, res[5, :], 'x-')
		ax5.set_title('C29 bin')
		ax5.set_ylim((0., 0.15))
		ax5.text(1., 0.13, str('m = ' + str(round(m[5], 3))))		

		ax6 = fig2.add_subplot(2, 4, 7)
		ax6.plot(PM2p5, res[6, :], 'x-')
		ax6.set_title('C30 bin')
		ax6.set_ylim((0., 0.15))
		ax6.text(1., 0.13, str('m = ' + str(round(m[6], 3))))
		ax6.set_xlabel(r'PM2.5 Concentration ' + cunit, fontsize = 14)
		

		ax7 = fig2.add_subplot(2, 4, 8)
		ax7.plot(PM2p5, res[7, :], 'x-')
		ax7.set_title('C31 bin')
		ax7.set_ylim((0., 0.15))
		ax7.text(1., 0.13, str('m = ' + str(round(m[7], 3))))

		# figure title
		fig2.suptitle('Figure 3 of Lunderberg et al. 2020 (doi.org/10.1021/acs.est.0c00966)')

		fig2.tight_layout() # space out subplots
		plt.show() # show figure

		return() # end function

# instantiate
plotter = self()
plotter.lunderberg_plotting() # call function