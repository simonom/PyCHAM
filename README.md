# PyCHAM: CHemistry with Aerosol Microphysics in Python Box Model

Welcome to the PyCHAM software for modelling of aerosol chambers.  Funding has been provided by the [EUROCHAMP-2020 research project](http://www.eurochamp.org) and the National Centre for Atmospheric Science ([NCAS](https://www.ncas.ac.uk/en/)).  Please open an issue above or contact Simon O'Meara (simon.omeara@manchester.ac.uk) with any issues, comments or suggestions.

PyCHAM is an open-access computer code (written in Python) for simulating aerosol chambers.  It is supplied under the GNU General Public License v3.0.

# Table of Content
1. [Documentation](#Documentation)
2. [Installation](#Installation)
3. [Running](#Running)
4. [Testing](#Testing)
5. [Inputs](#Inputs)

## Documentation

The README file you are now viewing serves as the PyCHAM manual, explaining how to setup the software and use it.

The [article](https://doi.org/10.21105/joss.01918) published in the Journal for Open Source Software explains the underlying mechanisms of PyCHAM and its purpose.  This article was reviewed using v0.2.4 of PyCHAM, which is stored in the releases of the PyCHAM github repository and [here](https://doi.org/10.5281/zenodo.3752677).

Version numbers of PyCHAM try to adhere to the semantics described by [semver](https://semver.org).

## Installation

There are two options for installing, via conda and via pip.  The pip method takes longer as the openbabel package has to be installed separately.  The instructions below for the pip method currently apply only to linux and macOS, whilst the conda instructions apply to windows, linux and macOS.


## Install via conda

1. Download the PyCHAM repository from github.com/simonom/PyCHAM

2. Download and install the package manager Anaconda using the following address and selecting the appropriate operating system version: https://www.anaconda.com/distribution/#download-section

3. Ensure conda is operating correctly, the method varies between operating systems and is explained [here](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html)

The following steps are at the command line:

4. cd into the directory where the PyCHAM package is stored.

5. Use the following command to install: conda env create -f PyCHAM_OSenv.yml -n PyCHAM, where OS is replaced by your operating system name (win (for Windows), lin (for Linux), mac (for macOS)).

6. Activate the environment with the command: conda activate PyCHAM

Install is complete, to run PyCHAM please see [Running](#Running).

## Install via pip

1) Ensure that swig is installed on your system.  For example, for macOS, a command at the command line like: brew install swig , and for linux a command at the command line like: sudo apt install swig

2) Ensure that eigen is available on your system.  For example, for macOS, a command at the command line like: brew install eigen, and for linux, download the latest stable eigen release [here](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download), then unzip and move the unzipped folder into the /usr/local/include folder using a command at the command line like: sudo mv eigen-3.3.9 /usr/local/include

3) Create a virtual environment in a suitable location running with at least python3.6.  For example, if your command line recognises python3.6 as python3.6, the command to make a virtual environment called 3env on macOS is: python3.6 -m venv 3env.  Note that python3.6 in this example should be replaced with the appropriate command for recognising python on your machine (often just python).

4) Activate this environment.  For example, for a virtual environment called 3env on macOS and linux the command at the command line is: source 3env/bin/activate.  For example, for a virtual environment called 3env on Windows the command at the command line is: .\3env\Scripts\activate

5) openbabel must be installed separately to PyCHAM, begin by downloading and unzipping the tar file (file name containing .tar.) for the latest version of openbabel on [github](https://github.com/openbabel/openbabel/releases).  Note that the unzipped version can be stored in the Downloads folder for the installation process.

The following steps are at the command line:

6) At the command line cd into the unzipped openbabel folder, for example: cd openbabel-3.1.1

7) Create a build directory: mkdir build

8) Change into the build directory: cd build

9) Using cmake, prepare the openbabel build files.  This requires several specifications, the -DCMAKE_INSTALL_PREFIX specification should be the path to the site-packages folder of the virtual environment created above: cmake -DRUN_SWIG=ON -DPYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=~/3env/lib/python3.6/site-packages ..

9a) If the above command causes a message at the command prompt that eigen cannot be found this can be fixed by stating its location, for example: 
cmake .. -DRUN_SWIG=ON -DPYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=~/3env/lib/python3.6/site-packages -DEIGEN3_INCLUDE_DIR=/usr/local/Cellar/eigen/3.3.9/include/eigen3

10) Complete installation with: make install

11) Test that openbabel is functioning with: python

11a) Then inside the python interpreter use: import openbabel

11b) If this works fine (no error message) continue to step 12 If this returns: ModuleNotFoundError: No module named 'openbabel', then quit the python interpreter: quit()

11c) If error seen in step above, add the openbabel path to the python path, for example: export PYTHONPATH=~/3env/lib/python3.6/site-packages/lib/python3.6/site-packages/openbabel:$PYTHONPATH

12) Test that pybel is functioning with: python

12a) Then inside the python interpreter use: import pybel

12b) If no error message seen continue to next step.  During testing this threw a relative import error which was corrected by changing the relevant line (given in the error message) in the pybel.py file (file location given in the error message) from "from . import openbabel as ob" to "import openbabel as ob"

13) Ensure pip and wheel up to date: pip install --upgrade pip wheel

14) Install PyCHAM and its dependencies in the virtual environment: python -m pip install --upgrade PyCHAM

Install is complete, to run PyCHAM please see [Running](#Running).

## Running

1. For model inputs, ensure you have: a .txt file chemical reaction scheme, a .xml file for converting the component names used in the chemical reaction scheme file to SMILE strings and a .txt file stating values of model variables (e.g. temperature) - see details about these three files below and note that example files are available in PyCHAM/input

2. Once [Installation](#Installation) is complete and the appropriate environment has been activated (see [Installation](#Installation)), use the command line to change into the top level directory PyCHAM (the directory above the PyCHAM __main__ file).

3. Run the model from the command line: python PyCHAM

4. The first three buttons of the GUI allow identification of the input files.  These can be used in any order and if no file is selected for any or all of the buttons, default files will be used (defaults: Chemical Scheme .txt File - example_scheme.txt, Chemical Scheme .xml File - example_xml.xml, Model Variables .txt File - example_model_var.txt) (see below for details on the chemical scheme, xml and model variables input files).

5. The 'Run Model' button starts the simulation - results will be saved in the output folder.  The time through the simulation will be displayed in the console along with updates of programme progress.

6. The 'Plot Results' button produces (and saves in the output folder) two plots: one with the particle number distribution, secondary aerosol mass, and particle number concentration against time, and another plot that shows the gas-phase concentrations of specified components with time (the specified components are those with initial concentrations given in the model variables file).  This button will not operate correctly until the simulation is complete.  The simulation is complete when the console says so.

7. The 'Quit' button will stop the programme.  If the programme is running and Quit does not work, the ctrl+z key combination in the console window can cease operations safely, though without results being saved and therefore making the 'Plot Results' button redundant.

## Testing

Unit tests for PyCHAM modules can be found in the PyCHAM/unit_tests folder.  Call these tests from the home folder for PyCHAM, with: python test_module.py with module replaced by the name of the PyCHAM module to be tested.  For some unit tests example inputs are required, the chemical scheme files for these are stored in unit_tests/input with file names beginning with test_ ..., therefore we recommend users do not use chemical schemes with the same naming convention to prevent confusion.  Where required, model variables for unit tests either use the default values or those given in the unit test script and use the xml file provided in PyCHAM/input.

Continuous integration testing can be completed using the '.travis.yml' (home folder) and 'test_TravisCI.py' (unit_tests folder) files at the [Travis CI website](https://travis-ci.com).

Example run output is saved in the PyCHAM/output/example_scheme folder.  To reproduce this, select from PyCHAM/input example_scheme.txt for the chemical scheme, example_xml.xml for the xml file and example_model_var.txt for the model variables.  Note that the example output may vary between releases so please check correspondence.

## Inputs

## Chemical Scheme file

The chemical scheme file states the reactions and their rate cofficients in the gas- and aqueous-phases.

An example chemical scheme file is given in the PyCHAM/input folder, called 'example_scheme.txt', which has been obtained from the [Master Chemical Mechanism (MCM) website](http://mcm.leeds.ac.uk/MCM/) (KPP version) and modified.  

Results are automatically saved in PyCHAM/output/name_of_chemical_scheme_file/name_given_in_model_variables_input_file_for_saving.  

The unit tests described above save results with the prefix 'test_', therefore we recommend using a different convention for chemical scheme names to prevent confusion.

Markers are required to recognise different sections of the chemical scheme.  The default markers are for the MCM KPP format, however, others can be specified using the chem_scheme_markers input in the model variables input file.  A guide to chem_scheme_markers is given in the Model Variables .txt File section below.  This includes how to distinguish between gas- and aqueous-phase reactions.

Reaction rate coefficients for chemical reactions and generic rate coefficients must adhere to the following rules:
The expression for the rate coefficient can use Fortran type scientific notation or python type; acceptable math functions: EXP, exp, dsqrt, dlog, LOG, dabs, LOG10, numpy.exp, numpy.sqrt, numpy.log, numpy.abs, numpy.log10; rate coefficients may be functions of TEMP, RH, M, N2, O2 where TEMP is temperature (K), RH is relative humidity (0-1), M, N2 and O2 are the concentrations of third body, nitrogen and oxygen, respectively (# molecules/cc (air)).

Inside the chemical scheme file, the expression for the reaction rate coefficient of a chemical reaction and the reaction itself must be contained on the same line of the file, with some marker (described above with chem_scheme_markers) separating them.

## Chemical Scheme .xml file

An example is given in the inputs folder (of the Github repository), called 'examples_xml.xml'.  It has a two line header, the first states that the mechanism is beginning (`<mechanism>`) and the second states that the species definition is beginning (`<species_defs>`).  The end of the species list must be marked (`</species_defs>`) and finally, the end of the mechanism must be marked (`</mechanism>`). 

Beneath this, every component included in the reactions of the chemical scheme must have its SMILES string given.  To add new components, use this three line example per new component:

`<species species_number="s6058" species_name="O2">`

`<smiles>O=O</smiles>`

`</species>`

Here the first line states that the species definition is beginning and gives a unique code, s6058 in this case, and its chemical scheme name, in this case O2.  The second line provides the SMILES string, in this case O=O.  The third line states that the definition is finished.  For information on SMILES please see: [SMILES website](https://daylight.com/smiles/index.html).


## Model Variables .txt File

An example is provided in the inputs folder (of the Github repository), called 
'example_model_var.txt' , this can include the following variables separated by a 
return (so one line per variable), 
note that if a variable is irrelevant for your simulation, it can be omitted and will be replaced by the default.

| Input Name | Description|
| ---------- | ---------- |
| res_file_name = | Name of folder to save results to |
| total_model_time = | Total experiment time to be simulated (s) |
| update_step =  | Time (s) interval for updating integration constants (specifically natural light intensity (if applicable) and particle number concentration due to its change during any of: coagulation, particle loss to wall and/or nucleation).  Defaults to 1 s.  Can be set to more than the total_model_time variable above to prevent updates. |
| recording_time_step =  | Time interval for recording results (s).  Must be at least the value of update_step if particles are present (number_size_bins variable below greater than zero).  Defaults to 60 s.|
| size_structure = | The size structure for the sectional approach to particles of varying size.  Set to 0 for moving-centre (default) and 1 for full-moving |
| number_size_bins = | Number of size bins (excluding wall); to turn off particle considerations set to 0 (which is also the default), likewise set pconc and seed_name variables below off.  Must be integer (e.g. 1) not float (e.g. 1.0) |
| lower_part_size = | Radius of smallest size bin boundary (um) |
| upper_part_size = | Radius of largest size bin boundary (um) |
| space_mode = | lin for linear spacing of size bins in radius space, or log for logarithmic spacing of size bins in radius space, if empty defaults to linear spacing|
| wall_on = | 1 to consider wall for gas-wall partitioning and particle deposition to wall, 0 to neglect these processes. |
| mass_trans_coeff = | Mass transfer coefficient of vapour-wall partitioning (/s), if left empty defaults to zero |
| eff_abs_wall_massC = | Effective absorbing wall mass concentration (g/m3 (air)), if left empty defaults to zero |
| temperature = | Air temperature inside the chamber (K).  At least one value must be given for the experiment start (times corresponding to temperatures given in tempt variable below).  If multiple values, representing changes in temperature at different times, then separate by a comma.  For example, if the temperature at experiment start is 290.0 K and this increases to 300.0 K after 3600.0 s of the experiment, inputs are: temperature = 290.0, 300.0, tempt = 0.0, 3600.0 |
| tempt = | Times since start of experiment (s) at which the temperature(s) set by the temperature variable above, are reached.  Defaults to 0.0 if left empty as at least the temperature at experiment start needs to be known.  If multiple values, representing changes in temperature at different times, then separate by a comma.  For example, if the temperature at experiment start is 290.0 K and this increases to 300.0 K after 3600.0 s of the experiment, inputs are: temperature = 290.0, 300.0; tempt = 0.0, 3600.0 |
| p_init =  | Pressure of air inside the chamber (Pa) |
| rh = | Relative Humidity (fraction, 0-1) |
| lat = | latitude (degrees) for natural light intensity (if applicable, leave empty if not (if experiment is dark set light_status below to 0 for all times)) |
| lon = | longitude (degrees) for natural light intensity (if applicable, leave empty if not (if experiment is dark set light_status below to 0 for all times)) |	
| DayOfYear = | day of the year for natural light intensity (if applicable, leave empty if not (if experiment is dark set light_status below to 0 for all times)), must be integer between 1 and 365|
| daytime_start = | Time of the day for natural light intensity (if applicable, leave empty if not (if experiment is dark set light_status below to 0 for all times)) (s since midnight) |
| act_flux_file = | Name of csv file stored in PyCHAM/photofiles containing the actinic flux values; use only if artificial lights inside chamber are used during experiment.  The file should have a line for each wavelength, with the first number in each line representing the wavelength in nm, and the second number separated from the first by a comma stating the flux (Photons /cm2/nm/s) at that wavelength.  No headers should be present in this file.  Example of file given by /PyCHAM/photofiles/Example_act_flux and example of the act_flux_path variable is: act_flux_path = Example_act_flux.csv.  Note, please include the .csv in the variable name if this is part of the file name.  Defaults to null file |
| photo_par_file = | Name of txt file stored in PyCHAM/photofiles containing the wavelength-dependent absorption cross-sections and quantum yields for photochemistry.  If left empty defaults to MCMv3.2, and is only used if act_flux_path variable above is stated.  File must be of .txt format with the formatting: <br> J_n_axs <br> wv_m, axs_m <br> J_n_qy <br> wv_M, qy_m <br> J_end <br> where n is the photochemical reaction number, axs represents the absorption cross-section (cm2/molecule), wv is wavelength (nm), _m is the wavelength number, and qy represents quantum yield (fraction).  J_end marks the end of the photolysis file.  An example is provided in PyCHAM/photofiles/example_inputs.txt.  Note, please include the .txt in the file name. |
| ChamSA = | Chamber surface area (m2), used if the Rader and McMurry wall loss of particles option (Rader_flag) is set to 1 (on) below|
| coag_on = | set to 1 (default if left empty) for coagulation to be modelled, or set to zero to omit coagulation|
| nucv1 = | Nucleation parameterisation value 1 to control the total number of newly formed particles|
| nucv2 = | Nucleation parameterisation value 2 to control the start time of nucleation|
| nucv3 = | Nucleation parameterisation value 3 to control the duration of nucleation|
| nuc_comp = | Name of component contributing to nucleation (only one allowed), must correspond to a name in the chemical scheme file.  Deafults to empty.  If empty, the nucleation module (nuc.py) will not be called.|
| new_partr = | Radius of newly nucleated particles (cm), if empty defaults to 2.0e-7 cm. |
| inflectDp = | Particle diameter wall deposition rate at inflection point (m). |
| Grad_pre_inflect = | Negative log10 of the gradient of particle wall deposition rate against the log10 of particle diameter before inflection (/s).  For example, for the rate to decrease by an order of magnitude every order of magnitude increase in particle diameter, set to 1.|
| Grad_post_inflect = | log10 of the gradient of particle wall deposition rate against the log10 of particle diameter after inflection (/s).  For example, for the rate to increase by an order of magnitude every order of magnitude increase in particle diameter, set to 1.|
| Rate_at_inflect = | Particle deposition rate to wall at inflection (/s). |
| part_charge_num = | Average number of charges per particle, only required if the McMurry and Rader (1985) model for particle deposition to walls is selected.|
| elec_field = | Average electric field inside the chamber (g.m/A.s3), only required if the McMurry and Rader (1985) model for particle deposition to walls is selected |
| McMurry_flag = | 0 to use the particle wall loss parameter values given above or 1 to the McMurry and Rader (1985, doi: 10.1080/02786828508959054) method for particle wall loss, which uses the chamber surface area given by ChamSA above, average number of charges per particle (part_charge_num above) and average electric field inside chamber (elec_field above), defaults to no particle wall loss if empty, similarly -1 turns off particle wall loss |
| C0 = | Initial concentrations of any trace gases input at the experiment's start (ppb), must correspond to component names in Comp0 variable below |
| Comp0 = | Names of trace gases present at start (in the order corresponding to their concentrations in C0).  Note, this is case sensitive, with the case matching that in the chemical scheme file|
| Ct = | Concentrations of component achieved when injected at some time after experiment start (ppb), if multiple values (representing injection at multiple times), separate with commas.  If multiple components are  injected after the start time, then this input should comprise the injected concentrations of components with times separated by commas and components separated by semicolons.  E.g., if k ppb of component A injected after m seconds and j ppb of component B injected after n (n>m) seconds, then injectt should be m, n and Compt should be A, B and Ct should be k,0;0,j.  The value here is the increase in concentration from the moment before the injection to the moment after (ppb).  Note this is for components with concentrations allowed to change, see const_comp for those with invariable concentrations |
| Compt = | Name of component injected at some time after experiment start.  Note, this is case sensitive, with the case matching that in the chemical scheme file - note this for components with concentrations allowed to change, see const_comp for those with invariable concentrations|
| injectt = | Time(s) at which injections occur (seconds), which correspond to the concentrations in Ct, if multiple values (representing injection at multiple times), please separate with commas.  If multiple components are  injected after the start time, then this input should still consist of just one series of times as these will apply to all components.  E.g., if k ppb of component A injected after m seconds and j ppb of component B injected after n (n>m) seconds, then this input should be m, n and Compt should be A, B and the Ct should be k,0;0,j Note this is for components with concentrations allowed to change, see const_comp for those with invariable concentrations|
| const_comp = | Name of component with continuous gas-phase concentration inside chamber.  Note, this is case sensitive, with the case matching that in the chemical file.  Defaults to nothing if left empty.  To specifically account for constant influx, see const_infl variable below.|
| const_infl = | Name of component(s) with continuous gas-phase influx to chamber. Note, this is case sensitive, with the case matching that in the chemical file. Defaults to nothing if left empty. For constant gas-phase concentration see const_comp variable above. Should be one dimensional array covering all components. For example, if component A has constant influx of K ppb/s from 0 s to 10 s and component B has constant influx of J ppb/s from 5 s to 20 s, the input is: const_infl = A, B Cinfl = K, K, 0, 0; 0, J, J, 0 const_infl_t = 0, 5, 10, 20 therefore, the semicolon in Cinfl is used to distinguish the influxes of different components |
| const_infl_t = | Times during which constant influx of each component given in the const_infl variable occurs, with the rate of their influx given in the Cinfl variable.  Should be one dimensional array covering all components.  For example, if component A has constant influx of K ppb/s from 0 s to 10 s and component B has constant influx of J ppb/s from 5 s to 20 s, the input is: const_infl = A, B Cinfl = K, K, 0, 0; 0, J, J, 0 const_infl_t = 0, 5, 10, 20 therefore, the semicolon in Cinfl is used to distinguish the influxes of different components |
| Cinfl = | Rate of gas-phase influx of components with constant influx (stated in the const_infl variable above).  In units of ppb/s.  Defaults to zero if left empty.  If multiple components affected, their influx rate should be separated by a semicolon, with a rate given for all times presented in const_infl_t (even if this is constant from the previous time step for a given component).  For example, if component A has constant influx of K ppb/s from 0 s to 10 s and component B has constant influx of J ppb/s from 5 s to 20 s, the input is: const_infl = A, B Cinfl = K, K, 0, 0; 0, J, J, 0 const_infl_t = 0, 5, 10, 20 therefore, the semicolon in Cinfl is used to distinguish the influxes of different components |
| dens_Comp = | Chemical scheme names of components with a specified density, if more than one name then separate with comma. The number of names must match the number of densities provided in the dens input.  Default is to estimate density based on the SMILE string of each component and the Girolami method contained in UManSysProp. |
| dens = | The density of components specified in the dens_Comp input above (g/cc), if more than one density then separate with a comma. The number of densities must match the number of names provided in the dens_Comp input.  Default is to estimate density based on the SMILE string of each component and the Girolami method contained in UManSysProp. |
| vol_Comp = | Names of components with vapour pressures to be manually assigned from volP, names must correspond to those in the chemical scheme file and if more than one, separated by commas.  Can be left empty, which is the default. |
| volP = | Vapour pressures (Pa) of components with names given in vol_Comp variable above, where one vapour pressure must be stated for each component named in vol_Comp and multiple values should be separated by a comma.  Acceptable for inputs to use e for standard notation, such as 1.0e-2 for 0.01 Pa. |
| act_comp = | Names of components (names given in the chemical scheme) with activity coefficients stated in act_user variable below (if multiple names, separate with a comma).  Must have same length as act_user.|
| act_user = | Activity coefficients of components with names given in act_comp variable above, if multiple values then separate with a comma.  Must have same length as act_comp. |
| accom_coeff_comp = | Names of components (corresponding to names in chemical scheme file) with accommodation coefficients set by the user in the accom_coeff_user variable below, therefore length must equal that of accom_coeff_user.  Multiple names must be separated by a comma.  For any components not mentioned in accom_coeff_comp, accommodation coefficient defaults to 1.0 |
| accom_coeff_user = | Accommodation coefficients (dimensionless) of the components with names given in the variable accom_coeff_comp variable, therefore number of accommodation coefficients must equal number of names, with multiple coefficients separated by a comma.  Can be a function of radius (m), in which case use the variable name radius, e.g: for NO2 and N2O5 with accommodation coefficients set to 1.0 and 6.09e-08/Rp, respectively, where Rp is radius of particle at a given time (m), the inputs are: accom_coeff_comp = NO2, N2O5 accom_coeff_user = 1., 6.09e-08/radius.  For any components not mentioned in accom_coeff_comp, accommodation coefficient defaults to 1.. |
| partit_cutoff = | The product of component saturation vapour pressure and activity coefficient (Pa) above which gas-particle partitioning assumed to be zero.  Defaults to empty list meaning that all components allowed to gas-particle partition.  Will only accept one value, which is applied to all components. |
| pconct = | Times (seconds) at which seed particles of number concentration given in pconc are introduced to the chamber.  If introduced at multiple times, separate times by a semicolon.  For example, for a two size bin simulation with 10 and 5 particles/cc in the first and second size bin respectively introduced at time 0 s, and later at time 120 s seed particles of concentration 6 an 0 particles/cc in the first and second size bin respectively are introduced, the pconc input is: pconc = 10, 5; 6, 0 and the pconct input is: pconct = 0; 120 and the number_size_bins input is: number_size_bins = 1 |
| pconc = | Either total particle concentration per mode (modes separated by a colon), or particle concentration per size bin, in which case length should equal number of particle size bins and  values should be separated by a comma (# particles/cc (air)).  If total particle concentration per mode, particles will be spread across size bins based on the values in the std and mean_rad inputs.  If seed aerosol introduced at mutiple times during the simulation, separate times using a semicolon, however maintain consistency between times as to whether number size distributions are being expressed in terms of modes or explicitly via the concentration per size bin.  For example, for a two size bin simulation with 10 and 5 particles/cc in the first and second size bin respectively introduced at time 0 s, and later at time 120 s seed particles of concentration 6 an 0 particles/cc in the first and second size bin respectively are introduced, the pconc input is: pconc = 10, 5; 6, 0 and the pconct input is: pconct = 0; 120 and the number_size_bins input is: number_size_bins = 1 | 
| seed_name = | Name of component(s) comprising seed particles, can either be core for a component not present in the chemical scheme, or a name from the chemical scheme, note no quotation marks needed.  If more than one component then separate names with commas.  PyCHAM will error if a name given here is not included in the chemical scheme (and isn't core).  It can be included in the chemical scheme with a zero reaction rate coefficient, for example, for AMM_SUL (ammonium sulphate) the text inside the quotation marks: "% 0 : AMM_SUL = AMM_SUL ;".  IMPORTANT: if the component(s) given in seed_name are volatile, the PyCHAM ODE solver may error since particles will evaporate, to prevent this, set the vapour pressures of seed components using the vol_Comp and volP model variables explained above.|
| seed_mw = | Molecular weight of seed component (g/mol), if empty defaults to that of ammonium sulphate (132.14 g/mol).  This only needs to be specified if seed_name input contains core.  If seed_name is a component(s) from the chemical scheme, then its molecular weight is estimated by Pybel (based on the component SMILE strings). |
| seed_dens = | Density of seed material (g/cc), defaults to 1.0 g/cc if left empty.  This only needs to be specified if the seed_name contains core.  If seed_name is a component(s) from the chemical scheme, then density should be specified using the dens_comp and dens inputs, otherwise density is estimated by UManSysProp (based on the component SMILE string). |
| seedVr = | Volume ratio of component(s) in seed particles, must match length of seed_name with the ratio of different components separated by a comma.  The same ratio will be applied to all size bins for all seed injections.  Defaults to equal volume contributions from each component given in seed_name.  For example, for two components with a 1:4 ratio please use 1, 4 (or equivalent, e.g., 25, 100). |
| mean_rad = | Mean radius of particles (um).  mean_rad should represent the mean radius of the lognormal size distribution per mode (modes separated by a colon).  Defaults to mean radius of the particle size bin radius bounds given by lower_part_size and upper_part_size inputs.  If seed particles are introduced at more than one time, then mean_rad for the different times should be separated by a semicolon.  For example, if seed particle with a mean_rad of 1.0e-2 um introduced at start and with mean_rad of 1.0e-1 um introduced after 120 s, the mean_rad input is: mean_rad = 1.0e-2; 1.0e-1 and the pconct input is pconct = 0; 120 |	
| std = | Geometric mean standard deviation of seed particle number concentration (dimensionless) when total particle number concentration(s) per mode provided in pconc variable.  If more than one mode, separate modes with a colon.  Role explained online in scipy.stats.lognorm page, under pdf method: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html.  If left empty defaults to 1.2.  If seed particles introduced after the experiment start, then separate std for different times using a semicolon.  For example, if seed particle with a standard deviation of 1.2 introduced at start and with standard deviation of 1.3 introduced after 120 s, the std input is: std = 1.2; 1.3 and the pconct input is: pconct = 0; 120 |
| seed_diss = | Dissociation constant(s) for seed component(s) (dimensionless), if empty defaults to one.  If more than one component comprising seed particles please separate their dissociation constants with a comma and ensure that number of constants matches number of components named in seed_name input. |
| light_time = | Times (s) for light status, corresponding to the elements of light_status (below), if empty defaults to lights off for whole experiment.  Use this setting regardless of whether light is natural or artificial (chamber lamps).  For example, for a 4 hour experiment, with lights on for first half and lights off for second, use: light_time = 0.0, 7200.0 light_status = 1, 0. If light_time doesn't include the experiment start (0.0 s), default is lights off at experiment start. |
| light_status = | 1 for lights on and 0 for lights off, with times given in light_time (above), if empty defaults to lights off for whole experiment.  Setting to off (0) means that even if variables defining light intensity above, the simulation will be dark.  Use this setting regardless of whether light is natural or artificial (chamber lamps).  The setting for a particular time is recognised when the time step will surpass the time given in light_time.  For example, for a 4 hour experiment, with lights on for first half and lights off for second, use: light_time = 0.0, 7200.0 light_status = 1, 0.  If status not given for the experiment start (0.0 s), default is lights off at experiment start. |
| tracked_comp = | Name of component(s) to track rate of concentration change (molecules/cc.s); must match name given in chemical scheme, and if multiple components given they must be separated by a comma.  Can be left empty and then defaults to tracking no components. |
| umansysprop_update = | Flag to update the UManSysProp module via internet connection: set to 1 to update and 0 to not update.  If empty defaults to no update.  In the case of no update, the module PyCHAM checks whether an existing UManSysProp module is available and if not tries to update via the internet.  If update requested and either no internet or UManSysProp repository page is down, code stops with an error. |
| chem_scheme_markers = | markers denoting various sections of the user's chemical scheme.  If left empty defaults to Kinetic Pre-Processor (KPP) formatting.  If filled, must have following elements separated with commas (brackets at start of description give pythonic index): (0) marker for start of gas-phase reaction lines (just the first element), note this must be different to that for aqueous-phase reaction, (1) marker for peroxy radical list starting, note that this should occur at the start of the peroxy radical list in the chemical scheme file, (2) marker between peroxy radical names, (3) prefix to peroxy radical name, (4) string after peroxy radical name, (5) marker for end of peroxy radical list (if no marker, then leave empty), (6) marker for RO2 list continuation onto next line, note this may be the same as marker between peroxy radical names, (7) marker at the end of each line containing generic rate coefficients, (8) marker for start of aqueous-phase reaction lines (just the first element), note this must be different to that for gas-phase reaction, (9) marker for start of reaction rate coefficient section of an equation line (note this must be the same for gas- and aqueous-phase reactions), (10) marker for start of equation section of an equation line (note this must be the same for gas- and aqueous-phase reactions), (11) final element of an equation line (should be constant for all phases of reactions).  For example, for the MCM KPP format (which only includes gas-phase reactions): chem_scheme_markers = {, RO2, +, C(ind_, ), , &, , , :, }, ; |
| int_tol = | Integration tolerances, with absolute tolerance first followed by relative tolerance, if left empty defaults to the maximum required during testing for stable solution: 1.0e-3 for absolute and 1.0e-4 for relative. |
| dil_fac = |Volume fraction per second chamber is diluted by, should be just a single number.  Defaults to zero if left empty.|
		
 
This project has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement No 730997.  Simon O'Meara has received funding from National Centre for Atmospheric Science (NCAS).
