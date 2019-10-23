# function called to run the PyCHAM model
# produces graphical user interface

import os
import sys
from tkinter import*
from tkinter.filedialog import askopenfilename
import pickle # for storing inputs
import threading # for calling on further functions
import numpy as np


if __name__ == '__main__':
        
	master = Tk() # intiate gui window
	heightv=550
	# to produce a vertical scrollbar need to embed button and entry widgets into a canvas
	# which is embedded into a frame
	frame = Frame(master, width=550, height=heightv) # establish frame
	frame.grid(row=0,column=0) # state frame position
	canvas = Canvas(frame, width=550, height=heightv, scrollregion=(0,0,2000,heightv*1.5))
	vbar = Scrollbar(frame, orient=VERTICAL)
	vbar.pack(side=RIGHT,fill=Y)
	vbar.config(command=canvas.yview)
	canvas.config(width=550, height=heightv)
	canvas.config(yscrollcommand=vbar.set)
	canvas.pack(side=LEFT,expand=True,fill=BOTH)
	master.title("PyCHAM") # goes in title bar


	def chemscm(): # function called when button pressed for chemical scheme
		fname = askopenfilename()
		
		def xml():
			xmlname = askopenfilename() # function called when button pressed for xml file
			
			def File():
				inname = askopenfilename()
				
				# open the file
				inputs = open(inname, mode='r')
				
				# read the file and store everything into a list
				in_list = inputs.readlines()
				inputs.close()
				
				if len(in_list) != 38:
					print('The number of model inputs is incorrect, should be 38, but is ' + str(len(in_list)) )
				for i in range(len(in_list)):
					key, value = in_list[i].split('=')
					key = key.strip() # a string
					if key == 'Res_file_name':
						resfname = str(value.strip())
					if key == 'Total_model_time':
						end_sim_time = float(value.strip())
					if key == 'Time_step':
						tstep_len = float(value.strip())
					if key == 'Recording_time_step':
						save_step = float(value.strip())
					if key == 'Number_size_bins':
						num_sb = int(value.strip())
					if key == 'lower_part_size':
						lowersize = float(value.strip())
					if key == 'upper_part_size':
						uppersize = float(value.strip())
					if key == 'eff_abs_wall_massC':
						Cw = float(value.strip())
					if key == 'Temperature':
						TEMP = float(value.strip())
					if key == 'PInit':
						PInit = float(value.strip())
					if key == 'RH':
						RH = float(value.strip())
					if key == 'lat':
						lat = float(value.strip())
					if key == 'lon':
						lon = float(value.strip())
					if key == 'daytime_start':
						dt_start = float(value.strip())
					if key == 'Cham_SA_to_Vol_ratio':
						cham_dim = float(value.strip())
					if key == 'ChamSA':
						ChamSA = float(value.strip())
					if key == 'Wall_accomm_coeff':
						wall_accom = float(value.strip())
					if key == 'nucv1':
						nucv1 = float(value.strip())
					if key == 'nucv2':
						nucv2 = float(value.strip())
					if key == 'nucv3':
						nucv3 = float(value.strip())
					if key == 'nuc_comp':
						nuc_comp = int(value.strip())
					if key == 'inflectDp':
						inflectDp = float(value.strip())
					if key == 'Grad_pre_inflect':
						pwl_xpre = float(value.strip())
					if key == 'Grad_post_inflect':
						pwl_xpro = float(value.strip())
					if key == 'Kern_at_inflect':
						inflectk = float(value.strip())
					if key == 'Rader_flag':
						Rader = int(value.strip())
					if key == 'C0':
						C0 = [float(i) for i in (value.split(','))]
					if key == 'Comp0': # strip removes white space
						Comp0 = [str(i).strip() for i in (value.split(','))]
					if key == 'voli':
						voli = [int(i) for i in (value.split(','))]
					if key == 'volP':
						volP = [float(i) for i in (value.split(','))]
					if key == 'test':
						test = int(value.strip())
					if key == 'pconc':
						pconc = float(value.strip())
					if key == 'std':
						std = float(value.strip())
					if key == 'loc':
						loc = float(value.strip())
					if key == 'scale':
						scale = float(value.strip())
					if key == 'core_diss':
						core_diss = float(value.strip())
					if key == 'light_time':
						light_time = [float(i) for i in (value.split(','))]
					if key == 'light_stat':
						light_stat = [int(i) for i in (value.split(','))]
				# function to write namelist to pickle file, which is read in by model
				def clicked():
					
					# write variable values to pickle file
					list_vars = [fname, num_sb, lowersize, uppersize, end_sim_time, 
					resfname, tstep_len, TEMP, PInit, RH, lat, lon, dt_start, Cw,  
					save_step, cham_dim, ChamSA, wall_accom, nucv1, nucv2, nucv3, nuc_comp,   
					inflectDp, pwl_xpre, pwl_xpro, inflectk, Rader, xmlname, C0, Comp0, 
					voli, volP, pconc, std, loc, scale, core_diss, light_stat, light_time]
					
					if test==1:
						print('Run Model button works successfully')
						with open('var_store.pkl','wb') as f:
							pickle.dump(list_vars,f)
							print('pickle file dumped successfully')
							print('now press the "Plot Results" button to ensure this works')
								
							
						# plots one figure, a superimposition of the:
						# particle number concentration, SOA mass, and number size
						# distribution
						def plotting():
							print('Plotting button works successfully')
							print('PyCHAM.py test complete and successful')
							exit()
							
						# vertical coordinate for plot results button
						b6 = Button(canvas, text='Plot Results', command=plotting)
						b6.pack()
						canvas.create_window(45, ystart+60, window=b6)
					
					if test!=1:
						with open('var_store.pkl','wb') as f:
							pickle.dump(list_vars,f)
							
						import front as model
						# flag to tell front to operate normally
						testf = 0
						# call on model to run
						t = threading.Thread(target=model.run(testf))
						t.daemon = False
						t.start()
						
						# plots one figure, a superimposition of the:
						# particle number concentration, SOA mass, and number size
						# distribution
						def plotting():
							
							import res_plot_super as plotter
							# pass the name of the folder where results are saved
							t = threading.Thread(target=plotter.run())
							t.daemon = False
							t.start()
							
	
						# vertical coordinate for plot results button
						b6 = Button(canvas, text='Plot Results', command=plotting)
						b6.pack()
						canvas.create_window(45, ystart+60, window=b6)
				
				# vertical coordinate for manual input button
				# button that calls PyCHAM model
				b5 = Button(canvas, text='Run Model', command=clicked)
				b5.pack()
				canvas.create_window(40, ystart+30, window=b5)				
				
			# button to choose text file with variable parameters stated
			b4 = Button(canvas, text='Model Variables .txt File', command=File)
			b4.pack()
			canvas.create_window(84, ystart-10, window=b4)
		
		# set button to select chemical scheme from file
		b2 = Button(canvas, text='Chemical Scheme .xml File', command=xml)
		b2.pack()
		canvas.create_window(92, ystart-40, window=b2)
	
	# vertical coordinate for manual input button
	ystart = 150
	
	# set button to select chemical scheme from file
	b1 = Button(canvas, text='Chemical Scheme .txt File', command=chemscm)
	b1.pack()
	canvas.create_window(89, ystart-70, window=b1)
	
	l21 = Label(canvas, text="PyCHAM Inputs",font=("Helvetica",30))
	l21.pack()
	canvas.create_window(280, ystart-130, window=l21)
	
	
	l20 = Label(canvas, text="Please see README file for help and Example Run for an example")
	l20.pack()
	canvas.create_window(234, ystart-100, window=l20)
	
	mainloop( )