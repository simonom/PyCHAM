# in this example, the user has set the photolysis 
# rates (J) explicitly, i.e. they do not need to be 
# calculated from actinic fluxes combined with
# absorption cross-sections and quantum yields,
# nor from the Hayman approach for solar
# radiation provided in photolysisRates.py
# note that to set this simulaion up, 
# light_status has been set to the path of the
# file containing the photolysis rates. These 
# photolysis rates are from spectral analysis of
# Manchester Aerosol Chamber lights, combined
# with the Mainz database. The photolysis rates
# file contains a change of photolysis with time,
# starting with the full intensity of lights,
# followed by a change to half intensity after
# 1800 seconds (half an hour). The change can 
# be demonstrated in PyCHAM by first selecting
# this example simulation and then asking to
# display photolysis rates using the 
# 'Check Values' button of the 'Simulate' tab
# and selecting 'Photolysis Rates' from the drop-
# down options. In addition, one can see the
# effect on the component O3, which is
# photolysed in this simulation, by running the 
# simulation and then selecting the 'Plot' tab
# followed by the 'Flux' tab followed by writing
# O3 into the text bar and selecting 'Plot change 
# tendency due to each chemical reaction'
res_file_name = set_J_directly_example_output
total_model_time = 3600.
update_step = 900.
recording_time_step = 900.
size_structure = 0
number_size_bins = 16
lower_part_size = 0.
upper_part_size = 5.e-2
space_mode = log
mass_trans_coeff = 1.e-4
eff_abs_wall_massC = 1.e2
temperature = 288.
tempt = 0.
p_init = 101300
rh = 0.65
rht = 0
light_status = J_values.xlsx
trans_fac = 1.0
DayOfYear = 183
daytime_start = 3.24e4
ChamSA = 109.85
coag_on = 1
nucv1 = 31000.
nucv2 = -55.
nucv3 = 180.
inflectDp = 2.0e-7
Grad_pre_inflect = 0.0
Grad_post_inflect = -7.0e-5
Rate_at_inflect = 9.5e-5
McMurry_flag = 0
C0 = 25., 30., 0.
Comp0 = O3, APINENE, ELVOC_O3
vol_Comp = ELVOC_O3
volP = 1e-30
umansysprop_update = 0
tracked_comp = O3
chem_scheme_markers = {, RO2, +, C(ind_, ), , &, , , :, }, ;,
