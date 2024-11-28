# PyCHAM: CHemistry with Aerosol Microphysics in Python Box Model

Welcome to the PyCHAM software for modelling of aerosol chambers.  Funding has been provided by the [EUROCHAMP-2020 research project](http://www.eurochamp.org) and the National Centre for Atmospheric Science ([NCAS](https://www.ncas.ac.uk/en/)).  Please open an issue on the GitHub repository or contact Simon O'Meara (simon.omeara@manchester.ac.uk) with any issues, comments or suggestions.

PyCHAM is an open-source computer code (written in Python) for simulating aerosol chambers.  It is supplied under the GNU General Public License v3.0.  The license document is provided with the software (LICENSE) and contains information around modification and conveyancing.

# Table of Content
1. [Documentation](#Documentation)
2. [Installation](#Installation)
3. [Running](#Running)
4. [Testing](#Testing)
5. [Inputs](#Inputs)
6. [Outputs](#Outputs)
7. [Photochemistry](#Photochemistry)
8. [Gas-particle Partitioning](#Gas-particle-Partitioning)
9. [Numerical Considerations](#Numerical-Considerations)
10. [Quick Plotting Tab](#Quick-Plotting-Tab)
11. [Flow Mode](#Flow-Mode)
13. [Indoor Air Quality Modelling](#Indoor-Air-Quality-Modelling)
14. [Frequently Asked Questions](#Frequently-Asked-Questions)
15. [Acknowledgements](#Acknowledgements)

## Documentation

The README file you are now viewing serves as the PyCHAM manual, explaining how to setup the software and use it.  As an additional resource, we also provide an [introductory video](https://www.youtube.com/watch?v=W8NbcU8WHeg&t=506s).

The [article](https://doi.org/10.21105/joss.01918) published in the Journal for Open Source Software explains the underlying mechanisms of PyCHAM and its purpose.  This article was reviewed using v0.2.4 of PyCHAM.  Additionally, the [article](https://doi.org/10.5194/gmd-14-675-2021) published in Geophysical Model Development provides a detailed introduction of PyCHAM and its use.  This article was reviewed using v2.1.1 of PyCHAM    The DOI for all PyCHAM releases is: [10.5281/zenodo.3752676](https://www.doi.org/10.5281/zenodo.3752676).



Version numbers of PyCHAM try to adhere to the semantics described by [semver](https://semver.org).

## Installation

There are two options for installing, via conda and via pip.  The pip method takes longer as the openbabel package has to be installed separately.  The instructions below for the pip method currently apply only to linux and macOS, whilst the conda instructions apply to windows, linux and macOS.


## Install via conda

1. Download the PyCHAM repository from github.com/simonom/PyCHAM

2. Download and install the package manager Miniconda (Anaconda is also suitable but takes more memory and takes longer to install) using the following address and selecting the appropriate operating system version: https://docs.conda.io/en/latest/miniconda.html.

3. Ensure conda is operating correctly, the method varies between operating systems and is explained [here](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html)

The following steps are at the command line:

4. cd into the directory where the PyCHAM package is stored.

5. Use the following command to install: conda env create -f PyCHAM_OSenv.yml -n PyCHAM, where OS is replaced by your operating system name (win (for Windows), lin (for Linux), mac (for macOS)).

6. Activate the environment with the command: conda activate PyCHAM

6a. For Windows (please skip this step if on Linux or Mac) it has been noted that when following the above steps openbabel does not always link to its data directory.  You can test this by looking for open babel warning messages when running the PyCHAM examples (to run PyCHAM please see [Running](#Running)).  A message such as 'the aromatic.txt file could not be found' indicates the data library is not set.  To set, search for Control Panel in the Windows search bar then User Accounts and User Accounts again.  Then select 'Change my environment variables' for User variables select New then create the new variable with Variable name 'BABEL_DATADIR' and with Variable value set as the path to the relevant open babel data directory.  To identify your open babel data directory you should look for the open babel folder containing aromatic.txt at a location such as: C:\Users\your_account_name\\.conda\pkgs\openbabel-2.4.1-py37_5\share\openbabel.  Once the path is suppled to the Variable value box select OK to create the variable then open a new anaconda prompt to restart PyCHAM.

6b. On Mac (please skip this step if on Linux or Windows) it has been noted that when following the above steps and then running PyCHAM (to run PyCHAM please see [Running](#Running)), an error is given stating that Developer Tools have not been installed. If this error displays please follow the on-screen instructions to install Developer Tools.

7. Install is complete for Linux, Mac and Windows; to run PyCHAM please see [Running](#Running).  

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

3. There are two choices for starting up the programme.  If you are new to PyCHAM, and/or have simulations that can be practically selected manually, then begin the programme from the command line to use PyCHAM via the graphical user interface: python PyCHAM.  Alternatively, if you wish to automate simulation setup and run, then edit the module 'automated_setup_and_call.py' to your needs and call this from the command line: python automated_setup_and_call

4. The PyCHAM graphical user interface (GUI) should now display on your screen.  Using the 'Simulate' tab, one can select the folder containing all input files using the 'Select Folder Containing Input Files' button.  This will search the selected folder for the input files (chemical reaction scheme, xml and model variables).  For the chemical scheme, files with filenames including 'chem' will be identified.  For the xml, files with filenames including 'xml' will be identified.  For the model variables, files with filenames including 'var' will be identified.  Any identified files will then be displayed in the GUI (see below for details on the contents of the chemical scheme, xml and model variables input files).

5. To select any of the input files individually, one can use the corresponding GUI button.

6. The GUI will display the found inputs provided in the selected model variables file.  For inputs not stated in this file, the displayed variables are default.

7. Problems with the input files will be displayed in the GUI - this functionality is under development, meaning that not all problems are currently captured.

8. Once the first simulation is ready (through selection of the desired combination of correct input files described above), the user chooses between a single simulation or adding to batch, with the latter allowing multiple simulations to be queued.

9a.  If the user chooses a single simulation to run, a progress bar will show, which represents the time through the experiment as a fraction of the total experiment time.

9b. If the user chooses to add to batch, then further simulations can be chosen by repeating steps 4-7 above.  When ready, the batch can be run with the start series of simulations button.  The progress bar then represents individual experiments and the current simulation is shown in the GUI.  Note that when adding to batch input files should be located in different folders (rather than changing the inputs inside a folder already selected for batch between adding to batch).

10. The 'Plot' tab allows multiple plotting options.  The Standard Results Plot produces two sub-plots in one figure: one with the particle number distribution, secondary aerosol mass, and particle number concentration against time, and another plot that shows the gas-phase concentrations of specified components with time (the specified components are those with initial concentrations given in the model variables file).

11. The 'Quit' button will stop the programme.  If it does not work, the ctrl+z key combination in the console window can cease operations safely.  In both cases Python will release all memory associated with the simulation.

## Testing

Unit tests for PyCHAM modules can be found in the PyCHAM/unit_tests folder.  Call these tests from the home folder for PyCHAM, with: python test_module.py with module replaced by the name of the PyCHAM module to be tested.  For some unit tests example inputs are required, the chemical scheme files for these are stored in unit_tests/input with file names beginning with test_ ..., therefore we recommend users do not use chemical schemes with the same naming convention to prevent confusion.  Where required, model variables for unit tests either use the default values or those given in the unit test script and use the xml file provided in PyCHAM/input.

Continuous integration testing can be completed using the '.travis.yml' (home folder) and 'test_TravisCI.py' (unit_tests folder) files at the [Travis CI website](https://travis-ci.com).

Example run output is saved in the PyCHAM/output/example_scheme folder.  To reproduce this, select from PyCHAM/input example_scheme.txt for the chemical scheme, example_xml.xml for the xml file and example_model_var.txt for the model variables.  Note that the example output may vary between releases so please check correspondence.

## Inputs

## Chemical Scheme file

The chemical scheme file states the reactions and their rate coefficients in the gas- and aqueous-phases.

An example chemical scheme file is given in the PyCHAM/input folder, called 'example_scheme.txt', which has been obtained from the [Master Chemical Mechanism (MCM) website](http://mcm.leeds.ac.uk/MCM/) (KPP version) and modified.  

Results are automatically saved in PyCHAM/output/name_of_chemical_scheme_file/name_given_in_model_variables_input_file_for_saving.  

The unit tests described above save results with the prefix 'test_', therefore we recommend using a different convention for chemical scheme names to prevent confusion.

Markers are required to recognise different sections of the chemical scheme.  The default markers are for the MCM KPP format, however, others can be specified using the chem_scheme_markers input in the model variables input file.  A guide to chem_scheme_markers is given in the Model Variables .txt file section below.  This includes how to distinguish between gas- and aqueous-phase reactions.

Reaction rate coefficients for chemical reactions and generic rate coefficients must adhere to the following rules:
The expression for the rate coefficient can use Fortran type scientific notation or python type; acceptable math functions: EXP, exp, dsqrt, dlog, LOG, dabs, LOG10, numpy.exp, numpy.sqrt, numpy.log, numpy.abs, numpy.log10; rate coefficients may be functions of TEMP, RH, M, N2, O2 where TEMP is temperature (K), RH is relative humidity (0-1), M, N2 and O2 are the concentrations of third body, nitrogen and oxygen, respectively (# molecules/cm3 (air)).

Inside the chemical scheme file, the expression for the reaction rate coefficient of a chemical reaction and the reaction itself must be contained on the same line of the file, with some delimiter (described above with chem_scheme_markers) separating them.

## Chemical Scheme .xml file

An example is given in the inputs folder (of the Github repository), called 'examples_xml.xml'.  It has a two line header, the first states that the mechanism is beginning (`<mechanism>`) and the second states that the species definition is beginning (`<species_defs>`).  The end of the species list must be marked (`</species_defs>`) and finally, the end of the mechanism must be marked (`</mechanism>`). 

Beneath this, every component included in the reactions of the chemical scheme must have its SMILES string given.  To add new components, use this three line example per new component:

`<species species_number="s6058" species_name="O2">`

`<smiles>O=O</smiles>`

`</species>`

Here the first line states that the species definition is beginning and gives a unique code, s6058 in this case, and its chemical scheme name, in this case O2.  The second line provides the SMILES string, in this case O=O.  The third line states that the definition is finished.  For information on SMILES please see: [SMILES website](https://daylight.com/smiles/index.html).

For Master Chemical Mechanism chemical schemes the associated xml file can be acquired from the MCM website.


## Model Variables .txt file

An example is provided in the inputs folder (of the Github repository), called 
'example_model_var.txt' , this can include the following variables separated by a 
return (so one line per variable), 
note that if a variable is irrelevant for your simulation, it can be omitted and will be replaced by the default.
In addition, you can find more information around photochemistry and flow mode instruments in the corresponding section of this README file.

| Input Name | Description|
| ---------- | ---------- |
| res_file_name = | Name of folder to save results to |
| total_model_time = | Total experiment time to be simulated (s) |
| update_step =  | Time (s) interval for updating integration constants (specifically natural light intensity (if applicable) and particle number concentration due to its change during any of: coagulation, particle loss to wall and/or nucleation).  Defaults to 1 s.  Can be set to more than the total_model_time variable above to prevent updates. |
| recording_time_step =  | Time interval for recording results (s).  Must be at least the value of update_step if particles are present (number_size_bins variable below greater than zero).  Defaults to 60 s.  Note that recorded values represent values at the recording interval.  Therefore, instantaneous injections at the recording interval are registered. |
| size_structure = | The size structure for the sectional approach to particles of varying size.  Set to 0 for moving-centre (default) and 1 for full-moving. |
| number_size_bins = | Number of size bins (excluding wall); to turn off particle considerations set to 0 (which is also the default), likewise set pconc and seed_name variables below off.  Must be integer (e.g. 1) not float (e.g. 1.0) |
| lower_part_size = | Radius of smallest size bin boundary (um), defaults to 0.0 um |
| upper_part_size = | Radius of largest size bin boundary (um), defaults to 0.5 um |
| space_mode = | lin for linear spacing of size bins in radius space, or log for logarithmic spacing of size bins in radius space, if empty defaults to linear spacing. |
| wall_on = | 1 to consider wall for gas-wall partitioning and particle deposition to wall, 0 to neglect these processes.  Defaults to 1.  However, in addition, please note the described defaults below for model variables relevant to wall processes. |
| number_wall_bins = | The number of wall bins to use, defaults to one.  Note that if wall_on set to zero, then this variable will be ignored and no wall will be simulated. |
| mass_trans_coeff = | Mass transfer coefficients for vapour-wall partitioning (/s), if left empty defaults to zero (which implies no partitioning with wall).  If using multiple wall bins, then separate values by a comma and ensure alignment with multiple values in the eff_abs_wall_massC variable (described below). To specify mass transfer coefficients for vapour-wall partitioning (/s) for specific components on specific wall, use the syntax: ComponentNameInChemicalScheme_walln_mtc, where ComponentNameInChemicalScheme is the name of the component as used in the chemical scheme (case sensitive), n is the wall number (beginning at 1), and mtc is the vapour-wall partitioning mass transfer coefficient (/s) for that chemical on that surface. To separate different components use a semi-colon. For example, to specify NO2 and HONO transfer coefficients onto the first surface, whilst all other components have the same transfer coefficient to that surface of 1x10-4 /s, and all components (including NO2 and HONO) have a transfer coefficient to a second surface of 1x10-5 /s: mass_trans_coeff = 1e-4; NO2_wall1_2e-4; HONO_wall1_3e-4, 1e-5|
| eff_abs_wall_massC = | Effective absorbing wall mass concentrations (g/m3 (air)), if left empty defaults to zero (which implies no partitioning with wall). If using multiple wall bins, then separate values by a comma and ensure alignment with multiple values in the mass_trans_coeff variable (described above). |
| temperature = | Air temperature inside the chamber (K).  At least one value must be given for the experiment start (times corresponding to temperatures given in tempt variable below).  If multiple values, representing changes in temperature at different times, then separate by a comma.  For example, if the temperature at experiment start is 290.0 K and this increases to 300.0 K after 3600.0 s of the experiment, inputs are: temperature = 290.0, 300.0, tempt = 0.0, 3600.0.  A change in temperature during the simulation will automatically cause relative humidity, chamber pressure, component volatilities and gas-phase diffusivities to change accordingly. |
| tempt = | Times since start of experiment (s) at which the temperature(s) set by the temperature variable above, are reached.  Defaults to 0.0 if left empty as at least the temperature at experiment start needs to be known.  If multiple values, representing changes in temperature at different times, then separate by a comma.  For example, if the temperature at experiment start is 290.0 K and this increases to 300.0 K after 3600.0 s of the experiment, inputs are: temperature = 290.0, 300.0; tempt = 0.0, 3600.0 |
| p_init =  | Pressure of air inside the chamber (Pa) |
| rh = | Relative Humidity (fraction, 0-1), if this changes during the simulation, values at different times should be separated by a comma, with the corresponding times provided in the rht model variable. Defaults to 0.65. If this model variable is used, a relative humidity at experiment start must be provided.  For example, for an experiment starting at relative humidity 0.8 and increasing to 0.9 after 30 minutes, inputs would be: rh = 0.8, 0.9 and rht = 0., 1800..  Note that relative humidity cannot be changed through Compt and associated model variables for instantaneous injection of gas-phase components to prevent conflicts with this rh model variable.  Note also that relative humidity will change in response to changing temperature (temperature model variable). It is possible to specify a change in relative humidity that is too great for PyCHAM to maintain stability in the ODE solver, see the [Numerical Considerations](#Numerical-Considerations) section below for more information on this. |
| rht = | Times (s) through simulation at which the relative humidities stated in the rh model variable are reached.  Defaults to 0, which implies a constant relative humidity.  If times provided, a time of 0 (experiment start) must also be provided along with a corresponding relative humidity in the rh model variable.  For example, for an experiment starting at relative humidity 0.8 and increasing to 0.9 after 30 minutes, inputs would be: rh = 0.8, 0.9 and rht = 0., 1800.|
| lat = | Latitude (degrees) for natural light intensity (if applicable, leave empty if not (if experiment is dark set light_status below to 0 for all times)). |
| lon = | Longitude (degrees) for natural light intensity (if applicable, leave empty if not (if experiment is dark set light_status below to 0 for all times)). |	
| DayOfYear = | Day of the year for natural light intensity (if applicable, leave empty if not (if experiment is dark set light_status below to 0 for all times)), must be integer between 1 and 365. |
| daytime_start = | Time of day experiment starts, for natural light intensity (if applicable, leave empty if not (if experiment is dark set light_status below to 0 for all times)) (Greenwich Mean Time (GMT)/Coordinated Universal Time (UTC) in seconds (not hours:minutes:seconds)). |
| act_flux_path = | Path to the csv file containing the actinic flux values; use only if you wish to specify actinic fluxes.  The file should have a line for each wavelength, with the first number in each line representing the wavelength in nm, and the second number separated from the first by a comma stating the flux (Photons/cm2/nm/s) at that wavelength.  No headers should be present in this file.  Example of file given by /PyCHAM/photofiles/Example_act_flux.csv and an example of the act_flux_path variable is: act_flux_path = /PyCHAM/photofiles/Example_act_flux.csv.  Note, please include the .csv in the variable name if this is part of the file name.  If the chamber light status is set to illuminated and a Master Chemical Mechanism chemical scheme is used, PyCHAM defaults to estimating the MCM photolysis reactions based on natural solar radiation using the parameterisation of Hayman (1997) which is described in [Saunders et al. (2003)](https://doi.org/10.5194/acp-3-161-2003) and which requires estimation of the solar zenith angle as described by the textbook chapter "The Atmosphere and UV-B Radiation at Ground Level" by S. Madronich (in 'Environmental UV Photobiology' textbook, 1993).  It is not necessary to state the full path to the actinic flux file - if the file is saved in the photofiles folder of PyCHAM, it is only necessary to state the name of the file. |
| photo_par_file = | Name of txt file stored in PyCHAM/photofiles containing the wavelength-dependent absorption cross-sections and quantum yields for photochemistry.  If left empty defaults to MCMv3.2 recommended values (http://mcm.leeds.ac.uk/MCMv3.3.1/parameters/photolysis.htt), which come as part of PyCHAM.  File must be of .txt format with the formatting: <br> J_n_axs <br> wv_m, axs_m <br> J_n_qy <br> wv_M, qy_m <br> J_end <br> where n is the photochemical reaction number, axs represents the absorption cross-section (cm2/molecule), wv is wavelength (nm), _m is the wavelength number, and qy represents quantum yield (fraction).  J_end marks the end of the photolysis file.  An example is provided in PyCHAM/photofiles/example_inputs.txt.  Note, please include the .txt in the file name. |
| ChamSA = | Chamber surface area (m2), used if the Rader and McMurry wall loss of particles option (McMurry_flag) is set to 1 (on) below.  Note that the model will convert this to a chamber radius by assuming the chamber is a sphere.|
| coag_on = | set to 1 (default if left empty) for coagulation to be modelled, or set to zero to omit coagulation|
| nucv1 = | Nucleation parameterisation value 1 to control the total number of newly formed particles|
| nucv2 = | Nucleation parameterisation value 2 to control the start time of nucleation|
| nucv3 = | Nucleation parameterisation value 3 to control the duration of nucleation|
| nuc_comp = | Name of component contributing to nucleation (only one allowed), must correspond to a name in the chemical scheme file, or 'core' for a generic zero vapour pressure component.  Defaults to 'core'.|
| new_partr = | Radius of newly nucleated particles (cm), if empty defaults to 2.0e-7 cm. |
| inflectDp = | Particle diameter (m) at which the particle deposition to wall rate function has an inflection point. Defaults to 1e-6 m (== 1 um). |
| Grad_pre_inflect = | Gradient of the logarithm of particle wall deposition rate against the logarithm of particle diameter before inflection.  For example, for the rate to decrease by an order of magnitude every order of magnitude increase in particle diameter, set to 1. Defaults to 0.|
| Grad_post_inflect = | Gradient of the logarithm of particle wall deposition rate against the logarithm of particle diameter after inflection .  For example, for the rate to increase by an order of magnitude every order of magnitude increase in particle diameter, set to 1. Defaults to 0.|
| Rate_at_inflect = | Particle deposition rate to wall at inflection (/s). Defaults to 0.|
| part_charge_num = | Average number of charges per particle, only required if the McMurry and Rader (1985) model for particle deposition to walls is selected.|
| elec_field = | Average electric field inside the chamber (g.m/A.s3), only required if the McMurry and Rader (1985) model for particle deposition to walls is selected |
| McMurry_flag = | 0 to use a user-defined particle to wall deposition rate as a function of particle size.  With the user defining this function through the model variables (inflectDp, Grad_pre_inflect, Grad_post_inflect, Rate_at_inflect, which are described above).  1 to use the McMurry and Rader (1985, doi: 10.1080/02786828508959054) method for particle wall loss, which uses the chamber surface area given by ChamSA above, average number of charges per particle (part_charge_num above) and average electric field inside chamber (elec_field above), defaults to no particle wall loss if empty, similarly -1 turns off particle wall loss.  Note, that if using the McMurry and Rader approach, it will be assumed that the surface is spherical.|
| C0 = | Initial concentrations of components present at the experiment start (ppb), must correspond to component names in Comp0 variable below. Can affect the gas-phase and/or surface (e.g. wall), see the Comp0 dexcription below for how to distinguish between reservoirs. |
| Comp0 = | Names of components present at experiment start (in the order corresponding to their concentrations in C0).  Note, this is case sensitive, with the case matching that in the chemical scheme file.  For components in the gas-phase, only their name (chemical scheme name) is needed. For surfaces, postfix chemical scheme names with '_walln', where n is the surface (e.g. wall) number starting from 1.|
| Ct = | The concentrations of components following instantaneous injection at some time after experiment start (ppb).  Seperate injections at different times with commas.  Seperate different components with a semicolon.  E.g., if k ppb of component A injected after m seconds and j ppb of component B injected after n (n>m) seconds, then injectt should be m, n and Compt should be A, B and Ct should be k,0;0,j.  Note, this is for components with concentrations allowed to change, see const_comp for those with invariable concentrations. |
| Compt = | Chemical scheme name of component injected instantaneously at some time after experiment start.  Note, this is case sensitive, with the case matching that in the chemical scheme file - note this for components with concentrations allowed to change, see const_comp for those with invariable concentrations.  Also note that water should not be stated here, rather, for varying relative humidity, use the rh and rht model variables. Separate components with a comma. |
| injectt = | Time(s) at which instantaneous injections occur (seconds), which correspond to the concentrations in Ct.  Separate multiple values (representing injection at multiple times) with commas.  If multiple components are  injected after the start time, then this input should still consist of just one series of times as these will apply to all components.  E.g., if k ppb of component A injected after m seconds and j ppb of component B injected after n (n>m) seconds, then this input should be m, n and Compt should be A, B and the Ct should be k,0;0,j Note this is for components with concentrations allowed to change, see const_comp for those with invariable concentrations. |
| const_comp = | Name of component with continuous gas-phase concentration inside chamber.  Note, this is case sensitive, with the case matching that in the chemical file.  Defaults to nothing if left empty.  To specifically account for constant influx, see const_infl variable below.|
| obs_file = | Name of xlsx file containing concentrations (molecules/cm3) of components and times (s) through experiment.  Component chemical scheme names must be in the first row, time must vary with later rows and the first column must contain times (s), whilst later columns contain the concentrations (molecules/cm3) of the components given in the first row.  If this variable is provided PyCHAM will automatically fix the component concentrations to those provided in this file.  The name of the file must also include the path from the PyCHAM home directory (and the file must be contained within the PyCHAM home directory or its sub-folders).  The sheet name to extract information from must be called PyCHAMobs.|
| const_infl = | Name of component(s) with continuous gas-phase influx to chamber. Note, this is case sensitive, with the case matching that in the chemical file. Defaults to nothing if left empty. For constant gas-phase concentration see const_comp variable above. Should be one dimensional array covering all components. For example, if component A has constant influx of K ppb/s from 0 s to 10 s and component B has constant influx of J ppb/s from 5 s to 20 s, the input is: const_infl = A, B Cinfl = K, K, 0, 0; 0, J, J, 0 const_infl_t = 0, 5, 10, 20 therefore, the semicolon in Cinfl is used to distinguish the influxes of different components.  If information on influxing component names, times and concentrations is more neatly held in an excel spreadsheet, then state the path to that spreadsheet as the value for const_infl. Do not contain the path inside quotation marks.  An example of stating the path inside a model variables file can be seen in PyCHAM/input/ind_AQ_ex/model_var_test.txt.  An example of such a spreadsheet is available in PyCHAM/input/ind_AQ_ex/const_infl.xlsx. Note that inside this example spreadsheet, the following formatting rules are exemplified: i) in the first row, first column cell, state the abundance unit of the emission rate: ppb for parts per billion, molecules/cm3 for # molecules/cm3, ii) in the first row, column 2 onwards, state the times (s) through simulation that emissions relate to, iii) in the first column, row 2 onwards |
| const_infl_t = | Times during which constant influx of each component given in the const_infl variable occurs, with the rate of their influx given in the Cinfl variable.  Should be one dimensional array covering all components.  For example, if component A has constant influx of K ppb/s from 0 s to 10 s and component B has constant influx of J ppb/s from 5 s to 20 s, the input is: const_infl = A, B Cinfl = K, K, 0, 0; 0, J, J, 0 const_infl_t = 0, 5, 10, 20 therefore, the semicolon in Cinfl is used to distinguish the influxes of different components |
| Cinfl = | Rate of gas-phase influx of components with constant influx (stated in the const_infl variable above).  In units of ppb/s.  Defaults to zero if left empty.  If multiple components affected, their influx rate should be separated by a semicolon, with a rate given for all times presented in const_infl_t (even if this is constant from the previous time step for a given component).  For example, if component A has constant influx of K ppb/s from 0 s to 10 s and component B has constant influx of J ppb/s from 5 s to 20 s, the input is: const_infl = A, B Cinfl = K, K, 0, 0; 0, J, J, 0 const_infl_t = 0, 5, 10, 20 therefore, the semicolon in Cinfl is used to distinguish the influxes of different components.  Cannot be an expression, e.g. 1.e-1, must be number, e.g. 0.1 instead of 1.e-1.  |
| dens_Comp = | Chemical scheme names of components with a specified density, if more than one name then separate with comma. The number of names must match the number of densities provided in the dens input.  Default is to estimate density based on the SMILE string of each component and the Girolami method contained in UManSysProp. |
| dens = | The density of components specified in the dens_Comp input above (g/cm3), if more than one density then separate with a comma. The number of densities must match the number of names provided in the dens_Comp input.  Default is to estimate density based on the SMILE string of each component and the Girolami method contained in UManSysProp. |
| vol_Comp = | Names of components with vapour pressures to be manually assigned from volP, names must correspond to those in the chemical scheme file and if more than one, separated by commas.  Can be left empty, which is the default (in which case vapour pressures are estimated from a vapour pressure estimation method by UManSysProp).  To specify a group of components based on their estimated vapour pressures (Pa), use inequalities, e.g., for all components with vapour pressures less than 1.e-2 Pa: all_<1.e-2.  And to specify for a certain wall, postfix with \_walln, where n is the wall number (starting at 1 for the first wall), e.g. for all components with estimated vapour pressure less than 1e-2 Pa and for the second wall: all_<1.e-2_wall2. |
| volP = | Vapour pressures (Pa) of components with names given in vol_Comp variable above, where one vapour pressure must be stated for each component named in vol_Comp and multiple values should be separated by a comma.  Acceptable to use e for standard notation, such as 1.e-2 for 0.01 Pa.  To specify the vapour pressure associated with a particular wall (wall order given by the order for the mass_trans_coeff and mass_trans_coeff variables described above), inside the vol_Comp variable described above, use _walln to postfix the relevant component/vapour pressure category name, where n is the wall number (starting at 1 for the first wall). |
| act_comp = | Names of components (names given in the chemical scheme) with activity coefficients stated in act_user variable below (if multiple names, separate with a comma).  Must have same length as act_user.|
| act_user = | Activity coefficients of components with names given in act_comp variable above, if multiple values then separate with a comma.  Must have same length as act_comp. |
| accom_coeff_comp = | Names of components (corresponding to names in chemical scheme file) with accommodation coefficients set by the user in the accom_coeff_user variable below, therefore length must equal that of accom_coeff_user.  Multiple names must be separated by a comma.  For any components not mentioned in accom_coeff_comp, accommodation coefficient defaults to 1.0.  For an introduction to accommodation coefficients, the recommended reading is page 525 of [Seinfeld and Pandis 2016](https://www.wiley.com/en-us/Atmospheric+Chemistry+and+Physics%3A+From+Air+Pollution+to+Climate+Change%2C+3rd+Edition-p-9781118947401).  |
| accom_coeff_user = | Accommodation coefficients (dimensionless) of the components with names given in the variable accom_coeff_comp variable, therefore number of accommodation coefficients must equal number of names, with multiple coefficients separated by a comma.  Can be a function of radius (m), in which case use the variable name radius, e.g: for NO2 and N2O5 with accommodation coefficients set to 1.0 and 6.09e-08/Rp, respectively, where Rp is radius of particle at a given time (m), you would use the inputs: accom_coeff_comp = NO2, N2O5 accom_coeff_user = 1., 6.09e-08/radius.  For any components not mentioned in accom_coeff_comp, accommodation coefficient defaults to 1.  See the description for the accom_coeff_comp variable for recommended reading on accommodation coefficients. |
| pconct = | Times (seconds) at which seed particles of number concentration given in pconc are introduced to the chamber (by default this assumed to be instantaneous injection but a continuous injection can be specified using the pcont variable).  If introduced at multiple times, separate times by a semicolon.  For example, for a two size bin simulation with 10 and 5 particles/cm3 in the first and second size bin respectively introduced at time 0 s, and later at time 120 s seed particles of concentration 6 and 0 particles/cm3 in the first and second size bin respectively are introduced, the pconc input is: pconc = 10, 5; 6, 0 and the pconct input is: pconct = 0; 120 and the number_size_bins input is: number_size_bins = 2.  Only one initial (pconct = 0.0 s) number size distribution is allowed.  If you wish to have injection of particles after experiment start but close to time = 0 s please use a pconct value that is greater than 0.0 s but small compared to the recording time step. |
| pconc = | Either total particle concentration per mode (modes separated by a colon), or particle concentration per size bin, in which case length should equal number of particle size bins and  values should be separated by a comma (# particles/cm3 (air)).  If total particle concentration per mode, the particles per mode will be spread across size bins, with the degree and location of spread based on the values in the std and mean_rad inputs, respectively.  If seed aerosol introduced at multiple times during the simulation, separate times using a semicolon, however maintain consistency between times as to whether number size distributions are being expressed in terms of modes or explicitly via the concentration per size bin.  For example, for a two size bin simulation with 10 and 5 particles/cm3 in the first and second size bin respectively introduced at time 0 s, and later at time 120 s seed particles of concentration 6 an 0 particles/cm3 in the first and second size bin respectively are introduced, the pconc input is: pconc = 10, 5; 6, 0 and the pconct input is: pconct = 0; 120 and the number_size_bins input is: number_size_bins = 2.  When transferring measured number concentration from a particle counter the unit of #particles/cm3 used for the pconc input refers to total number concentration (dN) rather than the density spectral (dN/dLogDp). |
| pcont = | Flag for whether the injection of particles given by pconct, pconc and associated inputs is continuous or instantaneous.  Defaults to instantaneous (flag = 0), in which case units of pconc are # particles/cm3 or can be set to 1 for continuous, in which case units of pconc are # particles/cm3.s. E.g., to change the example given in pconc description from the default two instantaneous injections to a continuous injection followed by an instantaneous injection: pcont = 1; 0. |
| seed_name = | Name of components comprising seed particles, can either be core for a component not present in the chemical scheme, or a name from the chemical scheme.  If more than one component then separate names with a comma.  Defaults to core.  Do not include water in seed_name, instead use the chamber relative humidity setting along with the Vwat_inc and seed_eq_wat model variables to determine water's contribution and timing of contribution to seed particles.  PyCHAM will error if a name given here is not included in the chemical scheme (and isn't core).  It can be included in the chemical scheme with a zero reaction rate coefficient.  IMPORTANT: if the components given in seed_name are volatile, the PyCHAM ODE solver may error since particles will evaporate.  To prevent this, set the vapour pressures of seed components using the vol_Comp and volP model variables explained above.  Note that seed components without a manually assigned vapour pressure (using the vol_Comp and volP model variables) will be automatically assigned a low vapour pressure to stop evaporation, but this can be changed through manual assignment using the vol_Comp and volP model variables. |
| seed_mw = | Molecular weight of seed component (g/mol), if empty defaults to that of ammonium sulphate (132.14 g/mol).  This only needs to be specified if seed_name input contains core.  If seed_name is a component(s) from the chemical scheme, then its molecular weight is estimated by Pybel (based on the component SMILE strings). |
| seed_dens = | Density of seed material (g/cm3), defaults to 1.0 g/cm3 if left empty.  This only needs to be specified if the seed_name contains core.  If seed_name is a component(s) from the chemical scheme, then density should be specified using the dens_comp and dens inputs, otherwise density is estimated by UManSysProp (based on the component SMILE string). |
| seedx = | Mole fraction of non-water components in dry (no water) seed particles.  Defaults to equal mole fractions for all non-water components stated in the seed_name model variable.  Must match length of seed_name with the mole fractions of different components separated by a semicolon.  If just one mole fraction provided for each component it will be assumed that this applies to all size bins.  If mole fractions are to be provided for each size bin, then separate values with a comma and the number of values per component must equal the number_size_bins variable above (the first value representing the smallest size bin and so on).  For example, for two components (seed_name = AMM_SUL, SOOT) and two size bins with mole fractions for the first component (AMM_SUL) of 0.1 and 0.2 in size bin 1 and size bin 2 and mole fractions for the second component (SOOT) of 0.9 and 0.8 in size bin 1 and size bin 2: seedx = 0.1, 0.2; 0.9, 0.8, where the component order (first, second, etc.), is given by the order stated in seed_name.  Please note that PyCHAM does not automatically estimate the effect of particle composition on activity coefficients (which affect gas-particle partitioning), nor reactions at the particle surface, nor reactions in the particle bulk.  Please see the [Gas-particle Partitioning](#Gas-particle-Partitioning) section for more related information. |
| Vwat_inc = | Flag to say whether (set to 1) or not (set to 0) water volume is accounted for in the seed particle number size distribution.  Default is 1. |
| mean_rad = | Mean radius of particles (um).  If in number size distributions given in modal mode, then mean_rad should represent the mean radius of the lognormal size distribution per mode (modes separated by a colon). Whereas if particle number concentrations given per size bin, mean_rad .  Defaults to mean radius of the particle size bin radius bounds given by lower_part_size and upper_part_size inputs.  If seed particles are introduced at more than one time, then mean_rad for the different times should be separated by a semicolon.  For example, if seed particle with two modes of mean radii of 1.e-2 and 1.e-1 um introduced at start and with mean radii of 2.0e-2 and 2.e-1 um introduced after 120 s, the mean_rad input is: mean_rad = 1.e-2 : 1.e-1 ; 2.e-2 : 2.e-1 and the pconct input is pconct = 0. ; 120. |
| seed_eq_wat = | If water is not included in the provided seed particle number size distribution (determined by the Vwat_inc model variable), this variable determines whether (1) or not (0) to allow water vapour in chamber (determined by the relative humidity) to equilibrate with seed particles prior to the experiment starting.  Defaults to 1 which allows equilibrium. |
| std = | Geometric mean standard deviation of seed particle number concentration (dimensionless) when total particle number concentrations per mode provided in pconc variable.  If more than one mode, separate modes with a colon.  Role explained online in scipy.stats.lognorm page, under pdf method: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html.  If left empty defaults to 1.2.  If seed particles introduced at multiple times, then separate std for different times using a semicolon.  For example, if seed particle with modes of standard deviation  1.2 and 1.3 introduced at start and with standard deviations of 1.4 and 1.5 introduced after 120 s, the std input is: std = 1.2 : 1.3 ; 1.4 : 1.5 and the pconct input is: pconct = 0. ; 120..  Note that std must be greater than 1.0, with wider distributions given by greater values.  Will error if a number less than or equal to 1.0 supplied.|
| seed_diss = | Dissociation constant(s) for seed component(s) (dimensionless), if empty defaults to one.  If more than one component comprising seed particles please separate their dissociation constants with a comma and ensure that number of constants matches number of components named in seed_name input. |
| z_prt_coeff = | Fraction of the total gas-particle partitioning coefficient below which partitioning of components (including water) to a size bin is treated as zero.  Defaults to one billionth (1.0e-9).  This setting is necessary for ODE solver stability when some particle size bins have relatively very small surface area.  To the best of our knowledge, the default value of one billionth has no significant effect on model results.  Please note that this is a variable required for numerical practicality and not for improved representation of simulated processes.  See the [Numerical Considerations](#Numerical-Considerations) section for further information. |
| light_time = | Times (s) for light status, corresponding to the elements of light_status (below), defaults to 0.0s (start of experiment).  Use this setting regardless of whether light is natural or artificial (chamber lamps).  For example, for a 4 hour experiment, with lights on for first half and lights off for second, use: light_time = 0.0, 7200.0 light_status = 1, 0. light_time must include 0 (experiment start) and a corresponding light status in the light_status model variable. |
| light_status = | 1 for lights on and 0 for lights off, with times given in light_time (above), if empty defaults to lights off for whole experiment.  Setting to off (0) means that even if variables defining light intensity above, the simulation will be dark.  Use this variable for both natural and artificial (chamber lamps) light.  The setting for a particular time is recognised when the time through experiment reached the time given in light_time.  For example, for a 4 hour experiment, with lights on for first half and lights off for second, use: light_time = 0.0, 7200.0 light_status = 1, 0.  If status not given for the experiment start (0.0 s), default is lights off at experiment start. |
| trans_fac = | Transmission factor (0-1) for light.  Use 1 for unhindered light and lower values to represent resistances to sunlight transmission, such as clouds or windows.  Note, if using the Hayman approach to estimating MCM photolysis rates (for natural light), a single value should be supplied, which will be applied to all MCM photolysis reactions (i.e. not wavelength-dependent).  For wavelength-dependent transmission factors (which are available when the user supplies their own actinic flux, see act_flux_path model variable), starting with the longest wavelength region, state the wavelength (nm) followed by the transmission factor (0-1 fraction) with an underscore separating the two.  For example, for 49.2 \% for wavelengths above UV and 15 \% for wavelengths in and below UV range: 400_0.492, 0_0.15   |
| tf_UVC = | Fraction (0-1) of UVC light (100-280 nm) (where relevant) stated in the provided actinic flux file (specified in the act_flux_file model variable) allowed into chamber.  E.g. when a UV-C lamp has variable input. |
| tf_UVCt = | Times (s) through experiment when values for the tf_UVC model variable are valid.  Defaults to 0.0 s (start of experiment), provide values in the same manner as described for the light_time model variable. |
| tracked_comp = | Name of component(s) to track rate of concentration change (molecules/cm3/s); must match name given in chemical scheme (description of how to track multiple components with a group name given later in this section), and if multiple components given they must be separated by a comma.  Can be left empty and then defaults to tracking no components.  Use RO2_ind and RO_ind to track all individual alkyl peroxy radicals and alkoxy radicals, respectively. |
| umansysprop_update = | Flag to update the UManSysProp module via internet connection: set to 1 to update and 0 to not update.  If empty defaults to no update.  In the case of no update, the module PyCHAM checks whether an existing UManSysProp module is available and if not tries to update via the internet.  If update requested and either no internet or UManSysProp repository page is down, code stops with an error. |
| chem_scheme_markers = | Markers denoting various sections of the user's chemical scheme.  If left empty defaults to Kinetic Pre-Processor (KPP) formatting.  If provided, must have following elements separated with commas (brackets at start of description give pythonic index): (0) marker for start of gas-phase reaction lines (just the first element), note this must be different to that for aqueous-phase reaction, (1) marker for peroxy radical list starting, note that this should occur at the start of the peroxy radical list in the chemical scheme file, (2) marker between peroxy radical names, (3) prefix to peroxy radical name, (4) string after peroxy radical name, (5) marker for end of peroxy radical list (if no marker, then leave empty), (6) marker for RO2 list continuation onto next line, note this may be the same as marker between peroxy radical names, (7) marker at the end of each line containing generic rate coefficients, (8) marker for start of aqueous-phase reaction lines (just the first element), note this must be different to that for gas-phase reaction, (9) marker for start of reaction rate coefficient section of an equation line (note this must be the same for gas- and aqueous-phase reactions), (10) marker for start of equation section of an equation line (note this must be the same for gas- and aqueous-phase reactions), (11) final element of an equation line (should be constant for all phases of reactions), (12) marker for start of reaction lines corresponing to surface (e.g. wall) (just the first element), note this must be different to that for gas-phase and aqueous-phase reaction.  For example, for the MCM KPP format (which only includes gas-phase reactions): chem_scheme_markers = {, RO2, +, C(ind_, ), , &, , , :, }, ;,  |
| int_tol = | Integration tolerances, with absolute tolerance first followed by relative tolerance. Separate absolute and relative tolerance with a comma, for example: 1.e-6, 1.e-7.  Defaults to the maximum required during testing for stable solution: 1.e-3 for absolute and 1.e-4 for relative. |
| dil_fac = | Volume fraction per second chamber is diluted by, should be just a single number.  Defaults to zero if left empty.|
| H2O_hist = | Flag for particle-phase history with respect to water partitioning: 0 for dry and therefore on the deliquescence curve, 1 for wet and therefore on the efflorescence curve.  Defaults to 1 if left empty.|
| drh_ft = | Expression for deliquescence relative humidity (fraction between 0-1) as a function of temperature, where the usual python math symbols should be used for mathematical functions and TEMP should be used to represent temperature which has units K.  E.g. for a deliquescence relative humidity at 298.15 K of 0.5 and an increase/decrease of 0.001 for every unit decrease/increase in temperature: drh_ft = 0.5-(1.e-3*(TEMP-298.15)).  Defaults to a deliquescence relative humidity of 0.0 at all temperatures if left empty (which combined with the default H2O_hist model variable of 1 would result in the assumption of no crystallisation and therefore particle-phase always treated as a solution).  For general information on deliquescence and efflorescence, please see page 410 of [Seinfeld and Pandis 2016](https://www.wiley.com/en-us/Atmospheric+Chemistry+and+Physics%3A+From+Air+Pollution+to+Climate+Change%2C+3rd+Edition-p-9781118947401) and references therein.  Research into the gas-particle partitioning of water for various types of particles is ongoing, and a literature search is recommended for a particular PyCHAM simulation setup.  Please see the [Gas-particle Partitioning](#Gas-particle-Partitioning) section for information on how PyCHAM uses this model variable.|
| erh_ft = | Expression for efflorescence relative humidity (fraction between 0-1) as a function of temperature, where the usual python math symbols should be used for mathematical functions and TEMP should be used to represent temperature which has units K.  E.g. for an efflorescence relative humidity at 298.15 K of 0.5 and an increase/decrease of 0.001 for every unit decrease/increase in temperature: erh_ft = 0.5-(1.e-3*(TEMP-298.15)).  Defaults to an efflorescence relative humidity of 0.0 at all temperatures if left empty (which combined with the default H2O_hist model variable of 1 would result in the assumption of no crystallisation and therefore particle-phase always treated as a solution).  For general information on deliquescence and efflorescence, please see page 410 of [Seinfeld and Pandis 2016](https://www.wiley.com/en-us/Atmospheric+Chemistry+and+Physics%3A+From+Air+Pollution+to+Climate+Change%2C+3rd+Edition-p-9781118947401) and references therein.  Research into the gas-particle partitioning of water for various types of particles is ongoing, and a literature search is recommended for a particular PyCHAM simulation setup.  Please see the [Gas-particle Partitioning](#Gas-particle-Partitioning) section for information on how PyCHAM uses this model variable.|
| ser_H2O = | Integer value for whether to separate the integration of the water partitioning between vapour and particle problem from integration of other processes.  Set to 0 to turn off separation and set to 1 to turn on (1 is default).  See [Numerical Considerations](#Numerical-Considerations) for more information. |

## Outputs

Model results are saved to the folder specified by the user in the model variables file, using the model variable res_file_name (described above).  PyCHAM automatically places this folder inside the PyCHAM/output/name of chemical scheme/ folder.  Several files are stored in this output folder.  The concentrations of components is stored in the file with name beginning 'concentrations_all_components_all_times'.  This is a comma separated value (csv) file.  If you would like suggestions for code to open and view results, please see below, along with the plotter_gp.py file and the retr_out.py files.

A minimum working example (for plotting the time profile of the gas-phase concentration of a given component) is (note that you may need to activate the PyCHAM environment in order to have the necessary packages available (numpy and matplotlib)):

```
# state path to output folder (for Windows Operating System use \\ to separate folders rather than /))
output_by_sim = 'path to your output folder'

# combine folder path with specific file name for component names
fname = str(output_by_sim + '/model_and_component_constants')

# open file containing component name
const_in = open(fname)
# create empty dictionary to hold component names
const = {}

# loop through lines of file containing component names
for line in const_in.readlines():
		
	dlist = [] # empty list to hold values
	for i in line.split(',')[1::]:
			
		if (str(line.split(',')[0]) == 'chem_scheme_names') or (str(line.split(',')[0]) == 'SMILES') or (str(line.split(',')[0]) == 'space_mode'):
			i = i.strip('\n')
			i = i.strip('[')
			i = i.strip(']')
			i = i.strip(' ')
			i = i.strip('\'')
			dlist.append(str(i))
			
	const[str(line.split(',')[0])] = dlist

# close file with component names
const_in.close()

# isolate component names from dictionary
comp_names = const['chem_scheme_names']

# withdraw times (s) -----------------
fname = str(output_by_sim+'/time')
# import numpy package
import numpy as np
# load times
t_array = np.loadtxt(fname, delimiter=',', skiprows=1)
timehr = t_array/3600.0 # convert from s to hr
# ------------------------------------------

# combine folder path with specific file name for component concentrations
fname = str(output_by_sim + '/concentrations_all_components_all_times_gas_particle_wall')
# load file, omitting headers
y = np.loadtxt(fname, delimiter=',', skiprows=1)

# state name of component you want to plot
comp_names_to_plot = 'APINENE'

# get index of this component
indx_plot = [comp_names.index(comp_names_to_plot.strip())]

# import packages for plotting
import matplotlib.pyplot as plt

# make plot
plt.plot(timehr, y[:, indx_plot])

# show plot
plt.show()

```

## Photochemistry
Chemical schemes may include photochemical reactions where the rate of reaction is dependent on light intensity.  Several of the model variables described here in the Model Variables .txt file section are relevant to correct modelling of photochemistry and these will be further detailed here.  

The input variables light_status and light_time determine when the chamber is illuminated or dark.  

The input variable act_flux_file states the actinic flux (photon/cm2/nm/s) as a function of wavelength (nm).  For chambers with artificial light (lamps) it is necessary to supply this file so that PyCHAM knows the light intensity spectrum.  PyCHAM will automatically interpolate the wavelengths and corresponding actinic fluxes given in act_flux_file to unit wavelength resolution (every 1 nm) to ensure correct integration of photolysis rate across the spectrum.  Inside the photofiles folder are examples of act_flux_file (e.g. Example_act_flux.csv), including the file for Manchester Aerosol Chamber (MAC) (MAC_Actinic_Flux_Spectrum.csv).  The required format is a comma separated value file with wavelength (nm) in the first column and the corresponding actinic flux (photon/cm2/nm/s) in the second column.  No headers are allowed.

For chambers with natural light (open roof), users may also supply an act_flux_file representing the relevant solar light intensity spectrum.  However, if natural light is present and the chemical scheme is derived from the Master Chemical Mechanism then PyCHAM will use the parameterisation of Hayman (1997), described in [Saunders et al. (2003)](https://doi.org/10.5194/acp-3-161-2003), to estimate the photolysis rates of the Master Chemical Mechanism.  In this model setup users may also supply the day number of the year (# days) (DayOfYear model variabled) time of day (Greenwich Mean Time (GMT)/Coordinated Universal Time (UTC) in seconds (not hours:minutes:seconds)) that the experiment starts (daytime_start model variable), the latitude (lat model variable) (degrees) and longitude (lon model variable) (degrees).  These inputs allow the solar zenith angle to be calculated according to the first chapter of Environmental UV Photobiology (1993): "The Atmosphere and UV-B Radiation at Ground Level" by S. Madronich (Environmental UV Photobiology, 1993).  This setting of solar radiation and deriving the photolysis rates for MCM is the default setting when light_status is set to illuminated.

Photolysis rate is the product of actinic flux, component absorption cross-sections (wavelength dependent) and quantum yield (wavelength dependent) integrated over the relevant range of the light spectrum.  By default PyCHAM assumes the Master Chemical Mechanism photolysis reaction rate coefficients require estimation.  For this reason, the PyCHAM software comes with the component absorption cross-sections and quantum yields as recommended by the Master Chemical Mechanism v3.3.1 website: http://mcm.leeds.ac.uk/MCMv3.3.1/parameters/photolysis.htt.

For photolysis reactions in chemical schemes other than Master Chemical Mechanism, the user can supply their own file for component absorption cross-section and quantum yield (the values are to be contained in the same file).  The file name should be stated in the photo_par_file model variable and the file should be stored in PyCHAM/photofiles.  A short example is given in the photofiles folder, called example_inputs.txt.  File must be of .txt format with the formatting: <br> J_n_axs <br> wv_m, axs_m <br> J_n_qy <br> wv_M, qy_m <br> J_end <br> where n is the photochemical reaction number, axs represents the absorption cross-section (cm2/molecule), wv is wavelength (nm), _m is the wavelength number, and qy represents quantum yield (fraction).  J_end marks the end of the photolysis file.

If a lamp power is given in terms of watts (J (kg m^2/s^2)/s), and a spectrum per unit wavelength is provided then it is possible to convert to actinic flux in units of photons/s.  To do this use the photon energy formula: E (J (kg m^2/s^2)) = h (J (kg m^2/s^2) s)*c (m/s)/lamda (m), where h is Planck's constant (6.626e-34 Js), c is the speed of light (3e8 m/s) and lambda is the wavelength (m), to get the energy of one photon.  Then, divide the provided Watts (J/s) per unit wavelength value by the result to get the actinic flux (photons/s) (note that by definition the actinic flux corresponds to a unit wavelength and a unit area).

## Gas-particle Partitioning
Whilst the numerical treatment of gas-particle partitioning is described in the [Geophysical Model Development (GMD) paper](https://doi.org/10.5194/gmd-14-675-2021), here the link between model variables and gas-particle partitioning is further described.

Setting the seed particle composition through the model variables seed_name and seedx has no automatic effect on the activity coefficient of components (and therefore on gas-particle partitioning via activity), since PyCHAM does not automatically estimate activity coefficients.

To estimate gas-particle partitiong, PyCHAM uses the difference in concentration of a component between the gas and particle phase, which depends on the mole fraction of that component in the particle phase.  Particle-phase mole fractions are a function of the concentrations of all components present in a particle and, for seed particle components, depends on the product of their particle-phase concentration and their dissociation constants (seed_diss model variable).  The dissociation constant is the number of ions an inorganic component dissociates into at infinite dilution per molecule of that component, e.g. for ammonium sulphate the constant is 3 and for sodium chloride it is 2.  Gas-particle partitioning is also a function of the activity coefficient of a component with respect to particles.  The activity coefficient is not estimated automatically in PyCHAM, but maybe set by the user through the act_user and act_comp model variables.  This is currently limited to one value per component for all size bins throughout a whole simulation.  More information on activity coefficients can be found on page 407 of [Seinfeld and Pandis 2016](https://www.wiley.com/en-us/Atmospheric+Chemistry+and+Physics%3A+From+Air+Pollution+to+Climate+Change%2C+3rd+Edition-p-9781118947401).

The default model treatment for gas-particle partitioning of water is the same as for other components, including the ability of the user to set the activity coefficient for water.  The user may additionally (for water only) use the drh_ft, erh_ft and H2O_hist  model variables to set whether water is able to partition to particles or not at a given relative humidity and water partitioning history.  If partitioning is not allowed, particles are modelled as completely dry (zero mole fraction of water), whilst if partitioning is allowed, water is allowed to partition as described above in this section and in the [Geophysical Model Development (GMD) paper] (https://doi.org/10.5194/gmd-14-675-2021).  PyCHAM cannot automatically determine deliquescence or efflorescence relative humidities, nor activity coefficients of water.   

## Numerical Considerations
In this section aspects of PyCHAM affecting numerical stability and computation time speed-up are discussed.

Speed-up of computation time can be achieved through solving gas-particle partitioning of water separately to other processes (the other processes are gas-wall partitioning of water, partitioning of non-water components between gas-particle and gas-wall, and chemical reactions).  For a system with ~1000 chemical reactions and 32 particle size bins a speed-up of a factor ~500 was seen when this separation was introduced.  Separation is done by default but can be turned off by setting the ser_H2O model variable to 0.  To the best of our knowledge, solving water gas-particle partitioning separately has negligible effect on integration estimates.

The Ordinary Differential Equation (ODE) solver package is [solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) from Scipy.  For the integration of the vapour-particle partitioning of water problem the Radau integration method is used as testing indicates this gives least computation time.  For integration of other processes (vapour-particle partitioning of non-water components and chemistry) problems, the backward differentiation formula (BDF) method is used as it is well suited to stiff problems.

The user can supply their own integration tolerances (int_tol) in the model variables file.  By default PyCHAM uses tolerances that were found to suit the problems presented in the GMD software decription paper (cited above).  However, non-stiff problems can be solved with less computation using higher integration tolerances, whilst stiffer problems may become unstable unless lower tolerances are used.

The ODE solver can fail to find a stable solution.  Sometimes this is due to a problem being too stiff.  If this occurs a message is printed to the PyCHAM GUI.  The initial message describes that instability has been noted and a smaller integration time step is being tried to improve stability.  Below a certain small time step the message changes and the simulation stops.  It has been decided that a stable solution cannot be found.  In this situation, to help the user identify the cause of instability, a file called ODE_solver_break_relevant_fluxes is produced which contains the fluxes of components with negative concentrations output by the ODE solver (where negative concentrations are physically unrealistic and therefore indicative of causing/resulting from integration instability).  The user could identify small or large fluxes (relative to other fluxes in this file) to gain indication of the processes causing instability and therefore determine whether a change of inputs is possible or whether PyCHAM is unsuitable for the problem in question.

Stability of the ODE solver may be compromised by a change in gas-particle partitioning of water that is too great for the solver. We recommend that instantaneous changes to relative humidity therefore do not surpass 0.01. In the example input titled rh_change_example is an example of relativly large changes to relative humidity that PyCHAM can cope with, however much greater changes would cause instability.

If the user-supplied model variables includes a seed component that is volatile then it may evaporate.  This can cause instability in the ODE solver.  Because this scenario is (in our experimental experience) unlikely in reality, PyCHAM automatically makes seed components non-volatile.  However, this can be over-ridden by specifying a pure component saturation vapour pressure for seed components in the user-supplied model variables. 

Computation time speed-up and ODE solver stability may be gained by increasing the value of the user-supplied model variable z_prt_coeff.  This variable determines for which size bins gas-particle partitioning of components (including water) is allowed.  If the gas-particle partitioning coefficient of a particle size bin is sufficiently low that it comprises less than this fraction of the total gas-particle partitioning coefficient (sum of coefficients across all size bins), then for this size bin partitioning is treated as zero.  If a size bin represents a relatively very small fraction of total partitioning the partitioning problem can become too stiff and prevent the ODE solver from concluding.  To the best of our knowledge the default value of 1.e-9 for z_prt_coeff has negligible effect on model estimates.  However, if users state a higher number, we recommend a sensitivity test for their simulation between this user-supplied number and one an order of magnitude higher to check whether a significant effect is seen - if not then the user-supplied number is shown to be conservative as it does not affect results - if it is then the user must decide whether the loss of accuracy is acceptable.

A change in boundary conditions, such as lighting, can cause concentrations to tend toward zero.  In turn, this can cause negative concentrations of components.  Therefore, limited negativity is allowed when checking for instability in the solution.  This can be seen inside the ode_updater module below the call to the ode_solv module.  An example input folder called neg_conc_example can invoke the allowable negative concentrations.  Gas-phase (and possibly particle-phase) output for the components with MCM names: O, O1D and PROPALO show a tendency toward zero after 50 minutes of experiment, with fluctuations above and below zero.

## Quick Plotting Tab

Particle mass concentration (whether total of all components or excluding certain components) is calculated based on the concentration of components and their molecular weight.  At a given time, for a given component in a given particle size bin (note that the concentration of a component represents the sum over all particles present in a size bin) the formula is: ((# molecules/cm3)/(Avogadro's Number (# molecules/mol)))*(g/mol)*1.e12 = ug/m3.  This equation can be extended over certain size bins and certain groupings of components to attain the desired particle-phase mass concentration.

## Flow Mode

When preparing model variable inputs for an instrument in flow mode (e.g. a flow tube or a chamber in flow reactor mode), the dilution factor and continuous influx of component model variables can be used.  To simulate removal of a constant fraction of the chamber's volume per second, set the dil_fac model variable accordingly.  For example, if the residence time in the instrument is 10 seconds, then 0.1 of the volume is removed per second, so use dil_fac = 0.1.  If components of interest (including potentially water) are injected to the instrument to replace the components lost through chamber air being extracted, then use the model variables: const_infl, const_infl_t and Cinfl to describe their continuous influx.

## Indoor Air Quality Modelling

For simulating atmospheres inside a building, an example set of input files is provided at PyCHAM/input/ind_AQ_ex.  When simulating buildings for first time, it is recommended to first read the notes above regarding Flow Mode.  This is because ventilation of the building will extract a certain fraction of air per second, which will be replaced by outdoor air.

Other features that are relevant to indoor air simulation include: the wavelength-dependent attenuation of solar radiation (see the model variable trans_fac above); the difference rates of deposition of gases and vapours to surfaces (possibly due to differing reaction rates on surfaces) (see the vol_Comp model variable above); the emission of gases/vapours into the building either from outdoors, or from indoor surfaces can be dealt with through the const_infl model variable (described above), which has the option of specifying a path to a spreadsheet containing emissions into the gas phase (e.g. if the number of components or/and the number of time points is sufficiently great that it is more neatly contained in a spreadsheet).

PyCHAM has been tested against observations from homes with the purpose of verifying that it can include the key aspects of atmospheric science that are particularly important for the indoor environment: surface-gas-particle parititioning of organics (verified against [Lunderberg et al. (2020)] (https://dx.doi.org/10.1021/acs.est.0c00966)), ozone reaction on surfaces to generate oxidised organics (verified against [Morrison et al. (2022)] (https://pubs.rsc.org/en/content/articlelanding/2022/EM/D2EM00307D)), the effect of ventilation; the effect of indoor emissions on indoor aerosol (gas and particle), infiltration of outdoor pollution indoors, infiltration of indoor pollution outdoors. 

The [Lunderberg et al. (2020)] (https://dx.doi.org/10.1021/acs.est.0c00966) study observes that components with volatilities comparable to C13-C22 alkanes show increased mixing ratios in the gas-phase+particle-phase with increased temperature, with decreased dependency as volatility decreases. They conclude that increased temperature decreases the fraction of these components on surfaces. Furthermore, components with volatilities comparable to C23-C31 alkanes show increased mixing ratios in the gas-phase+particle-phase with increased mass concentrations of particulate matter with diameter below 2.5 um. They conclude that a reservoir for these components on surfaces acts to modulate the particle-phase concentration. The PyCHAM inputs to generate the Lunderberg et al. (2020) results and the code to generate the PyCHAM results below are stored in PyCHAM/input/ind_AQ_ex/Lunderberg2020. Note that to run this code for your system, you will need to set the dir_path and ret_path variables to the path for your system. The plots below show results from PyCHAM simulations and reproduce the Lunderberg et al. (2020) observations.

![Lunderberg Fig. 2](https://github.com/simonom/PyCHAM/blob/master/Images_for_readme/Lunderberg2020Fig2.png "Figure 2")
![Lunderberg Fig. 3](https://github.com/simonom/PyCHAM/blob/master/Images_for_readme/Lunderberg2020Fig3.png "Figure 3")

## Frequently Asked Questions

**Why does PyCHAM crash without an error message?**
This has been observed using the conda install on a Windows 10 operating system with a Intel(R) Core(TM) i7-8500y CPU @ 1.50Ghz processor with processor speed 1400 MHz.  Checking the Event Viewer application on Windows, under the Windows Logs/Application tab, showed that libblas was crashing.  To correct, the PyCHAM virtual environment was activated in the command line, then from the command line conda was used to uninstall libblas and its dependents, then conda was used to install libblas again followed by its dependents.  This solved the issue.

**How are seed component concentrations initialised?**
The concentration of seed particles is based on the seed particle properties supplied by the user in the model variable file.  The molar volume (cm3/mol) of the seed component is calculated using: (g/mol)/(g/cm3) or (molar mass)/(component density), where density is estimated from the Girolami method of UManSysProp.  The volume of seed particles per size bin is calculated from the number size distribution stated in the model variable file.  Seed particle volumes per size bin are then divided by the molar volume of seed components to estimate the concentration of the seed components per size bin: # molecules/cm3 (seed components per size bin) = (cm3/cm3 (seed particle per size bin)/(cm3/mol) (seed components))*Avogadro's constant.

**What if I need an automated method for generating the xml file?**
A suggested code is contained in SMILES_generator.py - which allows xml file generation for the PRAM chemical mechanism, as well as accretion product formation from two organic peroxy radicals.

## Acknowledgements
This project has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement No 730997.  Simon O'Meara received funding support from the Natural Environment Research Council through the National Centre for Atmospheric Science.
