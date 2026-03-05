Welcome to the input and analysis of PyCHAM simulations for the Zhang et al. paper on the chemical characteristics and environmental drivers of nitrogen-containing organic aerosol formation in coastal and inland urban atmospheres.

To run a simulation:
i) Find the model variable files for the six simulations used in the paper inside the  folder of PyCHAM/output/Myanmar_NOCs_paper_Zhang_etal/Myanmar_input. Model variable file names are something like Mandalay_summer_40RH_mv.txt
ii) For the simulation of interest, make sure its model variable file is the only model variable file with mod_var in the file name (whilst the other model variable files do not contain the mo_var string). For example, for Mandalay winter at 80 % set the file name to Mandalay_winter_80RH_mod_var.txt. Doing this means PyCHAM will automatically select this model variable file.
iii) Inside the model variable file of of interest, provide the path to the wanted chemical scheme, and continuous influx rate file, for example, for Mandalay winter at 80 %:
	chem_sch_name = mcm_export.fac
	cont_infl = cont_infl_Mandalay_winter_80RH.xlsx
If these files are in the same folder as the model variable file, their full path does not need stating, just their file names. Similarly, PyCHAM will automatically detect the xml file in the same folder if it has the file name xml.
iv) Open the PyCHAM GUI, select the PyCHAM/output/Myanmar_NOCs_paper_Zhang_etal/Myanmar_input folder and check that the correct files have been selected by checking the dialogue boxes on the right hand side of the PyCHAM GUI Simulate tab. Then run PyCHAM with the prepared simulation.
v) New results will (unless the res_file_name variable in the model variable file is changed) be stored in /PyCHAM/output/mcm_export
vi) To plot results for relevant gas-phase or particle-phase components averaged over the final 24 hours of simulation, use the python script stored at PyCHAM/output/Myanmar_NOCs_paper_Zhang_etal/Myanmar_results/ind_component_time_integrated_conc.py. Set the names of wanted components using the plot_name variable and set the phase variable to 'g' for gas-phase and 'p' for particle-phase.
vi) To plot time-series use the relevant options in the PyCHAM GUI Plot tab.