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

The README file you are now reading serves as a manual explaining how to setup the software and use it.  

The [article](https://doi.org/10.21105/joss.01918) published in the Journal for Open Source Software explains the underlying mechanisms of PyCHAM and its purpose.  This article was reviewed using v0.2.4 of PyCHAM, which is stored in the releases of the PyCHAM github repository and [here](https://doi.org/10.5281/zenodo.3752677).

Version numbers of PyCHAM try to adhere to the semantics described by [semver](https://semver.org).

## Installation

There are two options for installing, via conda and from source.  Experience indicates that the conda install is more straightforward than the pip, therefore we recommend this.


## Install from conda

1. Download the PyCHAM repository from github.com/simonom/PyCHAM

2. Download and install the package manager Anaconda using the following address and selecting the appropriate operating system version: https://www.anaconda.com/distribution/#download-section

The following steps need to be at the command line:

3. cd into the directory where the PyCHAM package is stored.

4. Use the following command to install: conda env create -f PyCHAM_OSenv.yml -n PyCHAM, where OS is replaced by your operating system name (win (for Windows), lin (for Linux), mac (for Mac)).

5. Now the environment is setup you can activate it with the command: conda activate PyCHAM

Install is complete

## Install from source

1) open your terminal/command prompt 

2) cd to directory where you want the PyCHAM environment stored (here we will use the example Documents)

3) create an environment called myenv: python3 -m venv myenv

4) activate the environment: source myenv/bin/activate

5) cd to the environment’s site packages: cd lib/python3.x/site-packages

6) install PyCHAM: python3 -m pip install --upgrade PyCHAM

7) make directory to contain sundials build and install (inside site-packages): mkdir sundials

8) make directory to build sundials: mkdir sundials/builddir

9) make directory to install sundials: mkdir sundials/installdir

10) download .tar file for sundials-3.2.1 from: https://github.com/LLNL/sundials/releases/tag/v3.2.1

11) unzip and move to the site-packages folder in the environment you have created above

12) in the terminal/command prompt, from inside the site-packages directory: cd sundials/builddir

13) this next step requires that cmake is installed on your system (https://cmake.org/install/), it allows you to configure sundials: ccmake /Documents/myenv/lib/python3.x/site-packages/sundials-3.2.1

14) press c to view install options

15) using the i key, set the CMAKE_INSTALL_PREFIX and EXAMPLES_INSTALL_PATH to your installdir path, e.g. Documents/myenv/lib/python3.x/site-packages/sundials/installdir

16) press ‘c’ (causes configuration) then ‘g’ (generation)

17) back in the terminal/command window: make

18) finally, to complete installation of sundials: make install

19) download the .tgz file for BLAS from: http://www.netlib.org/blas/

20) unzip the BLAS download and move to the site-packages folder

21) cd into the BLAS folder

22) into terminal type: make

23) download the .tar.gz file for LAPACK from: 	http://www.netlib.org/lapack/

24) unzip and move to the site-packages folder

25) copy the blas_LINUX.a from the BLAS folder to the LAPACK folder

26) inside LAPACK folder copy the make.inc.template (or make.inc.example) file and rename make.inc and state address of the blas_LINUX.a beside the BLASLIB variable, e.g. (note the one space between $ and the path): BLASLIB = $ Documents/myenv/lib/python3.x/site-packages/lapack3.9.0/blas_LINUX.a

27) in terminal, inside the LAPACK folder type: make

28) still inside the LAPACK folder copy BLAS .a file to the system folder: sudo cp blas_LINUX.a /usr/local/lib/ 

29) still inside the LAPACK folder copy lapack .a. file to the system fodler: sudo cp liblapack.a /usr/local/lib/
 
30) install Cython: pip3 install Cython

31) in a text editor open Cython/Compiler/main.py and find language_level =, set to: language_level = 3

32) save and close Main.py and cd back to the site-packages directory

33) download the .tar file for Assimulo-3.0: https://github.com/modelon/Assimulo/releases

34) unzip this and move to the site-packages folder

35) cd to the new Assimulo folder

36) install assimulo, stating the path to sundials, blas and lapack: e.g.: python setup.py install --sundials-home=/Users/Simon_OMeara/Documents/Manchester/postdoc_stuff/box-model/PyCHAM/myenv/lib/python3.6/site-packages/sundials/installdir --blas-home=/Users/Simon_OMeara/Documents/Manchester/postdoc_stuff/box-model/PyCHAM/myenv/lib/python3.6/site-packages/BLAS-3.8.0 --lapack-home=/Users/Simon_OMeara/Documents/Manchester/postdoc_stuff/box-model/PyCHAM/myenv/lib/python3.6/site-packages/lapack-3.9.0

37) set the environment variable so that assimulo can link to the sundials library: export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Users/Simon_OMeara/Documents/Manchester/postdoc_stuff/box-model/PyCHAM/myenv/lib/python3.6/site-packages/sundials/installdir/lib/

38) cd out of the assimulo folder: cd ..

39) install openbabel: pip3 install openbabel

Install is complete

## Running

1. For model inputs, ensure you have: a .txt file chemical reaction scheme, a .xml file for converting species names to SMILE strings and a .txt file stating values of model variables (e.g. temperature) - see details about these three files below and note that example files are available in PyCHAM/inputs

2. Now you are ready to run the model: python PyCHAM

3. Follow the gui directions (see below for details on the chemical scheme, xml and model input files)

4. The 'run model' button starts the simulation - results will be saved in the output folder in your PyCHAM directory

5. The 'plot results' button produces (and saves in the output folder) two plots: the particle number distribution, SOA mass and particle number concentration against time, and another that shows the gas-phase concentrations of specified components with time.

## Testing

Unit tests for PyCHAM modules can be found in the PyCHAM/Unit_Testing folder in the Github repository.  To use, cd to this folder and use python test_module.py with module replaced by the name of the module to be tested.

Integration testing can be completed using the '.travis.yml' and 'PyCHAM_CI_test.py' files at the [Travis CI website](https://travis-ci.org).

Example run output is saved in the PyCHAM/output/Example_Run folder.  To reproduce this, select from PyCHAM/inputs Example_Run for the chemical scheme, Example_Run_xml for the xml file and Example_Run_inputs for the model variables.  Note that the example output may vary between releases so check correspondence.

## Inputs

## Chemical Scheme .txt file

An example chemical scheme .txt file is given in the inputs folder (of the Github repository), called 'Example_Run.txt', which has been obtained
from the [Master Chemical Mechanism (MCM) website](http://mcm.leeds.ac.uk/MCM/) (FACSIMILE version) and modified.

Identifiers are required to recognise different sections of the chemical scheme (to be stated in the model variables input file in the chem_scheme_markers input).

Reaction rate coefficients for chemical reactions and generic rate coefficients must adhere to the following rules:
The expression for the rate coefficient can use Fortran type scientific notation or python type; acceptable math functions: EXP, exp, dsqrt, dlog, LOG, dabs, LOG10, numpy.exp, numpy.sqrt, numpy.log, numpy.abs, numpy.log10; rate coefficients may be functions of TEMP, RH, M, N2, O2 where TEMP is temperature (K), RH is relative humidity (0-1), M, N2 and O2 are the concentrations of third body, nitrogen and oxygen (# molecules/cc (air)).

Inside the chemical scheme file, the expression for the reaction rate coefficient of a chemical reaction and the reaction itself must be contained on the same line of the file, with some identifier separating them.

## Chemical Scheme .xml file

An example is given in the inputs folder (of the Github repository), called 'Examples_Run_xml.xml'.  It has a two line header, the first states that the mechanism is beginning (`<mechanism>`) and the second states that the species definition is beginning (`<species_defs>`).  The end of the species list must be marked (`</species_defs>`) and finally, the end of the mechanism must be marked (`</mechanism>`). 

Beneath this, every species included in the reactions of the chemical scheme must have its SMILES string given.


## Model Variables .txt File

An example is provided in the inputs folder (of the Github repository), called 
'Example_Run_inputs.txt' , this must include the following variables separated by a 
return (so one line per variable, an error message will show if a variable is missing), 
note that if a variable is irrelevant for your simulation, it can be left empty 
(e.g. vol_Comp = ):

| Input Name | Description|
| ---------- | ---------- |
| res_file_name = | Name of folder to save results to |
| total_model_time = | Total experiment time to be simulated (s) |
| update_step =  | Time (s) interval for updating integration constants (specifically natural light intensity (if applicable) and particle number concentration due to its change during any of: coagulation, particle loss to wall and/or nucleation).  Default to 60 s.  Can be set to more than the total_model_time variable above to prevent updates. |
| recording_time_step =  | Time interval for recording results (s).  Must be at least the value of update_step if particles are present (number_size_bins variable below greater than zero).  Defaults to 60 s.|
| number_size_bins = | Number of size bins (excluding wall); to turn off particle considerations set to 0 (which is also the default), likewise set pconc and seed_name variables below off.  Must be integer (e.g. 1) not float (e.g. 1.0) |
| lower_part_size = | Radius of smallest size bin boundary (um) |
| upper_part_size = | Radius of largest size bin boundary (um) |
| space_mode = | lin for linear spacing of size bins in radius space, or log for logarithmic spacing of size bins in radius space, if empty defaults to linear spacing|
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
| Ct = | Concentrations of component achieved when injected at some time after experiment start (ppb), if multiple values (representing injection at multiple times), please separate with commas.  If multiple components are  injected after the start time, then this input should comprise the injected concentrations of components with times separated by commas and components separated by semicolons.  E.g., if k ppb of component A injected after m seconds and j ppb of component B injected after n (n>m) seconds, then injectt should be m, n and Compt should be A, B and Ct should be k,0;0,j.  The value here is the increase in concentration from the moment before the injection to the moment after (ppb).  Note this is for components with concentrations allowed to change, see const_comp for those with invariable concentrations |
| Compt = | Name of component injected at some time after experiment start.  Note, this is case sensitive, with the case matching that in the chemical scheme file - note this for components with concentrations allowed to change, see const_comp for those with invariable concentrations|
| injectt = | Time(s) at which injections occur (seconds), which correspond to the concentrations in Ct, if multiple values (representing injection at multiple times), please separate with commas.  If multiple components are  injected after the start time, then this input should still consist of just one series of times as these will apply to all components.  E.g., if k ppb of component A injected after m seconds and j ppb of component B injected after n (n>m) seconds, then this input should be m, n and Compt should be A, B and the Ct should be k,0;0,j Note this is for components with concentrations allowed to change, see const_comp for those with invariable concentrations|
| const_comp = | Name of component with continuous gas-phase concentration inside chamber.  Note, this is case sensitive, with the case matching that in the chemical file.  Defaults to nothing if left empty.  To specifically account for constant influx, see const_infl variable below.|
| const_infl = | Name of component(s) with continuous gas-phase influx to chamber. Note, this is case sensitive, with the case matching that in the chemical file. Defaults to nothing if left empty. For constant gas-phase concentration see const_comp variable above. Should be one dimensional array covering all components. For example, if component A has constant influx of K ppb/s from 0 s to 10 s and component B has constant influx of J ppb/s from 5 s to 20 s, the input is: const_infl = A, B Cinfl = K, K, 0, 0; 0, J, J, 0 const_infl_t = 0, 5, 10, 20 therefore, the semicolon in Cinfl is used to distinguish the influxes of different components |
| const_infl_t = | Times during which constant influx of each component given in the const_infl variable occurs, with the rate of their influx given in the Cinfl variable.  Should be one dimensional array covering all components.  For example, if component A has constant influx of K ppb/s from 0 s to 10 s and component B has constant influx of J ppb/s from 5 s to 20 s, the input is: const_infl = A, B Cinfl = K, K, 0, 0; 0, J, J, 0 const_infl_t = 0, 5, 10, 20 therefore, the semicolon in Cinfl is used to distinguish the influxes of different components |
| Cinfl = | Rate of gas-phase influx of components with constant influx (stated in the const_infl variable above).  In units of ppb/s.  Defaults to zero if left empty.  If multiple components affected, their influx rate should be separated by a semicolon, with a rate given for all times presented in const_infl_t (even if this is constant from the previous time step for a given component).  For example, if component A has constant influx of K ppb/s from 0 s to 10 s and component B has constant influx of J ppb/s from 5 s to 20 s, the input is: const_infl = A, B Cinfl = K, K, 0, 0; 0, J, J, 0 const_infl_t = 0, 5, 10, 20 therefore, the semicolon in Cinfl is used to distinguish the influxes of different components |
| vol_Comp = | Names of components with vapour pressures to be manually assigned from volP, names must correspond to those in the chemical scheme file and if more than one, separated by commas.  Can be left empty, which is the default. |			
| volP = | Vapour pressures (Pa) of components with names given in vol_Comp variable above, where one vapour pressure must be stated for each component named in vol_Comp and multiple values should be separated by a comma.  Acceptable for inputs to use e for standard notation, such as 1.0e-2 for 0.01 Pa |
| act_comp = | Names of components (names given in the chemical scheme) with activity coefficients stated in act_user variable below (if multiple names, separate with a comma).  Must have same length as act_user.|
| act_user = | Activity coefficients of components with names given in act_comp variable above, if multiple values then separate with a comma.  Must have same length as act_comp. |
| accom_coeff_comp = | Names of components (corresponding to names in chemical scheme file) with accommodation coefficients set by the user in the accom_coeff_user variable below, therefore length must equal that of accom_coeff_user.  Multiple names must be separated by a comma.  For any components not mentioned in accom_coeff_comp, accommodation coefficient defaults to 1.0 |
| accom_coeff_user = | Accommodation coefficients (dimensionless) of the components with names given in the variable accom_coeff_comp variable, therefore number of accommodation coefficients must equal number of names, with multiple coefficients separated by a comma.  Can be a function of radius (m), in which case use the variable name radius, e.g: for NO2 and N2O5 with accommodation coefficients set to 1.0 and 6.09e-08/Rp, respectively, where Rp is radius of particle at a given time (m), the inputs are: accom_coeff_comp = NO2, N2O5 accom_coeff_user = 1.0, 6.09e-08/radius.  For any components not mentioned in accom_coeff_comp, accommodation coefficient defaults to 1.0 | 
| pconct = | Times (seconds) at which seed particles of number concentration given in pconc are introduced to the chamber.  If introduced at multiple times, separate times by a semicolon.  For example, for a two size bin simulation with 10 and 5 particles/cc in the first and second size bin respectively introduced at time 0 s, and later at time 120 s seed particles of concentration 6 an 0 particles/cc in the first and second size bin respectively are introduced, the pconc input is: pconc = 10, 5; 6, 0 and the pconct input is: pconct = 0; 120 and the number_size_bins input is: number_size_bins = 1 |
| pconc = | Either total particle concentration, in which case should be a scalar, or particle concentration per size bin, in which case length should equal number of particle size bins (# particles/cc (air)).  If an array of numbers, then separate numbers by a comma.  If a scalar, the particles will be spread across size bins based on the values in the std and mean_rad inputs.  To turn off particle considerations leave empty.  If seed aerosol introduced at mutiple times during the simulation, separate times using a semicolon.  For example, for a two size bin simulation with 10 and 5 particles/cc in the first and second size bin respectively introduced at time 0 s, and later at time 120 s seed particles of concentration 6 an 0 particles/cc in the first and second size bin respectively are introduced, the pconc input is: pconc = 10, 5; 6, 0 and the pconct input is: pconct = 0; 120 and the number_size_bins input is: number_size_bins = 1 | 
| seed_name = | name of component comprising the seed particles, can either be core for a component not present in the equation file, or a name from the equation list or H2O for water, note no quotation marks needed |		
| seed_mw = | molecular weight of seed component (g/mol), if empty defaults to that of ammonium sulphate - 132.14 g/mol |
| seed_dens = | Density of seed material (g/cc), defaults to 1.0 g/cc if left empty. |			
| mean_rad = | Mean radius of particles (um), defaults to flag that tells software to estimate mean radius from the particle size bin radius bounds given by lower_part_size and upper_part_size inputs above.  If more than one size bin the default is the mid-point of each.  If the lognormal size distribution is being found (using the std input below), mean_rad should be a scalar representing the mean radius of the lognormal size distribution.  If seed particles are introduced at more than one time, then mean_rad for the different times should be separated by a semicolon.  For example, if seed particle with a mean_rad of 1.0e-2 um introduced at start and with mean_rad of 1.0e-1 um introduced after 120 s, the mean_rad input is: mean_rad = 1.0e-2; 1.0e-1 and the pconct input is pconct = 0; 120 |	
| std = | Geometric mean standard deviation of seed particle number concentration (dimensionless) when scalar provided in pconc variable above, role explained online in scipy.stats.lognorm page, under pdf method: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html.  If left empty defaults to 1.1.  If seed particles introduced after the experiment start, then separate std for different times using a semicolon.  For example, if seed particle with a standard deviation of 1.2 introduced at start and with standard deviation of 1.3 introduced after 120 s, the std input is: std = 1.2; 1.3 and the pconct input is: pconct = 0; 120 |
| core_diss = | Core dissociation constant (for seed component) (dimensionless) (1), if empty defaults to one|
| light_time = | Times (s) for light status, corresponding to the elements of light_status (below), if empty defaults to lights off for whole experiment.  Use this setting regardless of whether light is natural or artificial (chamber lamps).  For example, for a 4 hour experiment, with lights on for first half and lights off for second, use: light_time = 0.0, 7200.0 light_status = 1, 0. If light_time doesn't include the experiment start (0.0 s), default is lights off at experiment start. |
| light_status = | 1 for lights on and 0 for lights off, with times given in light_time (above), if empty defaults to lights off for whole experiment.  Setting to off (0) means that even if variables defining light intensity above, the simulation will be dark.  Use this setting regardless of whether light is natural or artificial (chamber lamps).  The setting for a particular time is recognised when the time step will surpass the time given in light_time.  For example, for a 4 hour experiment, with lights on for first half and lights off for second, use: light_time = 0.0, 7200.0 light_status = 1, 0.  If status not given for the experiment start (0.0 s), default is lights off at experiment start. |
| tracked_comp = | Name of component(s) to track rate of concentration change (molecules/cc.s); must match name given in chemical scheme, and if multiple components given they must be separated by a comma.  Can be left empty and then defaults to tracking no components. |				 
| umansysprop_update = | Flag to update the UManSysProp module via internet connection: set to 1 to update and 0 to not update.  If empty defaults to no update.  In the case of no update, the module PyCHAM checks whether an existing UManSysProp module is available and if not tries to update via the internet.  If update requested and either no internet or UManSysProp repository page is down, code stops with an error. |
| chem_scheme_markers = | markers denoting various sections of the user's chemical scheme.  If left empty defaults to Kinetic Pre-Processor (KPP) formatting.  If filled, must have following elements separated with commas (brackets at start of description give pythonic index): (0) marker for punctuation at start of gas-phase reaction lines (just the first element), (1) marker for peroxy radical list starting, (2) punctuation between peroxy radical names, (3) prefix to peroxy radical name, (4) string after peroxy radical name, (5) marker for end of peroxy radical list (if no marker, then use the marker for RO2 list continuation onto next line), (6) marker for RO2 list continuation onto next line, (7) marker at the end of each line containing generic rate coefficients, (8) first element of aqueous-phase reaction lines (just the first element), (9) marker for start of reaction rate coefficient section of an equation line, (10) marker for start of equation section of an equation line, (11) final element of an equation line (should be constant for all phases of reactions).  For example, for the MCM KPP format: chem_scheme_markers = {, RO2, +, C(ind_, ), , &, , , :, }, ; |
| int_tol = | Integration tolerances, with absolute tolerance first followed by relative tolerance, if left empty defaults to the maximum required during testing for stable solution: 1.0e-3 for absolute and 1.0e-4 for relative. |
| dil_fac = |Volume fraction per second chamber is diluted by, should be just a single number.  Defaults to zero if left empty.|
		
 
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 730997.  Simon O'Meara has received funding from National Centre for Atmospheric Science (NCAS).