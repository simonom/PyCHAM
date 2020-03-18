# PyCHAM: Python CHemistry with Aerosol Microphysics Box Model

Welcome to the PyCHAM.  Funding has been provided from the [EUROCHAMP-2020 research project](http://www.eurochamp.org).  Please contact Simon O'Meara (simon.omeara@manchester.ac.uk) with any issues, comments or suggestions.

PyCHAM is an open-access computer code (written in Python) for simulating aerosol chambers.  It is supplied under the GNU General Public License v3.0.


# Table of Content
1. [Installation](#Installation)
2. [Running](#Running)
3. [Testing](#Testing)
4. [Inputs](#Inputs)

## Installation

There are two options for installing, via conda and from source.  Experience indicates that the conda install is more straightforward than the pip, therefore we recommend this.


## Install from conda

1. Download the PyCHAM repository from github.com/simonom/PyCHAM or data.eurochamp.org/modelling-tools/

2. Download and install the package manager Anaconda using the following address and selecting the appropriate operating system version: https://www.anaconda.com/distribution/#download-section

3. To set-up PyCHAM, use the terminal/command prompt, cd to the directory where the PyCHAM package is stored (likely called PyCHAM-master), then use the following command to install: conda env create -f PyCHAM_OSenv.yml -n PyCHAM, where OS is replaced by your operating system name (win (for Windows), lin (for Linux), mac (for Mac))

4. This will install all PyCHAM dependencies, no need for further installations.

5. Now the environment is set up you can activate it by typing into terminal: conda activate PyCHAM

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

5. The 'plot results' button produces (and saves in the output folder) two figures, the first called contours.png that shows the particle number distribution, SOA mass and particle number concentration against time, and the second called gas_pbb.png that shows the gas-phase concentrations of initial components with time  

## Testing

Unit tests for PyCHAM modules can be found in the PyCHAM/Unit_Testing folder in the Github repository.  To use, cd to this folder and use python test_module.py with module replaced by the name of the module to be tested.

Integration testing can be completed using the '.travis.yml' file at the website: travis-ci.org.

The example run output is saved in output/Example_Run/Example_output in the Github repository.  To reproduce this (e.g. for testing), use the inputs/Example_Run.txt for the chemical scheme, inputs/Example_Run_xml.xml for the xml file and inputs/Example_Run_inputs.txt for the model inputs.  Note that the example output was produced using version 0.2.1 of PyCHAM.

## Inputs

## Chemical Scheme .txt file

An example chemical scheme .txt file is given in the inputs folder (of the Github repository), called 'Example_Run.txt', which has been obtained
from the [Master Chemical Mechanism website](http://mcm.leeds.ac.uk/MCM/) and modified.

// can be used for preceding comments, as in the example.

Equations should adhere to the syntax in the following standard expression for an equation: % rate coefficient equation : reactant1 + reactant2 = product1 + product2 ;
Each equation should be represented on a new line.


The rate coefficients must adhere to the following rules:
The expression for the rate coefficient can use Fortran type scientific notation or python type; acceptable math functions: EXP, dsqrt, dlog, LOG, dabs, LOG10, numpy.exp, numpy.sqrt, numpy.log, numpy.abs, numpy.log10; MCM (master chemical mechanism) rate constants are accepted, e.g. KMT01, as are MCM photolysis rates, e.g. J<1> ;rate coefficients may be functions of TEMP, RH, M, N2, O2 where TEMP is temperature (K), RH is relative humidity (0-1), M, N2 and O2 are the concentrations of third body, nitrogen and oxygen (# molecules/cc (air)) - these concentrations are calculated automatically as a function of temperature and pressure inside eqn_parser.py; rate coefficients may also be a function of MCM names such as RO2.


## Chemical Scheme .xml file

An example is given in the inputs folder (of the Github repository), called 'Examples_Run_xml.xml'.  It has a two line header, the first states that the mechanism is beginning (`<mechanism>`) and the second states that the species definition is beginning (`<species_defs>`).  The end of the species list must be marked (`</species_defs>`) and finally, the end of the mechanism must be marked (`</mechanism>`). 

Beneath this, every species included in the reactions of the chemical scheme must have its SMILES string given.


## Model Variables .txt File

An example is provided in the inputs folder (of the Github repository), called 'Example_Run_inputs.txt' , this must include the values for the following variables separated by a return (so one line per variable), note that if a variable is irrelevant for your simulation, it can be left empty (e.g. voli = ):

Res_file_name = Name of folder to save results to

Total_model_time = Simulation time (s)

Time_step = Maximum time interval for ode (s)

Recording_time_step = Time interval for recording results (s)

Number_size_bins = Number of size bins (excluding wall), please use zero if no particles to be modelled

lower_part_size = Radius of smallest size bin boundary (um)

upper_part_size = Radius of largest size bin boundary (um)

kgwt = mass transfer coefficient of vapour-wall partitioning (m3/g.s)

eff_abs_wall_massC = effective absorbing wall mass concentration (g/m3 (air))

Temperature = Temperature (K)

PInit = the chamber pressure (Pa)

RH = Relative Humidity (fraction, 0-1)

lat = representative latitude (degrees) for light intensity (if lights used)

lon = representative longitude (degrees) for light intensity (if lights used)

daytime_start = representative time of the day for light intensity (s since midnight)

ChamSA = Chamber surface area (m2)

nucv1 = Nucleation parameterisation value 1

nucv2 = Nucleation parameterisation value 2

nucv3 = Nucleation parameterisation value 3

nuc_comp = Index of component contributing to nucleation (only one index allowed, can be absolute or relative)

inflectDp = Particle diameter wall loss kernel inflection point (m)

Grad_pre_inflect = Gradient of particle wall loss before inflection (/s)

Grad_post_inflect = Gradient of particle wall loss after inflection (/s)

Kern_at_inflect = Wall loss kernel value at inflection (/s)

Rader_flag = 0 to use the particle wall loss parameter values given above or
			 1 to the Rader and McMurry (1983) method for particle wall loss, which
			 uses the chamber surface area given by ChamSA above

C0 = Initial concentrations of any trace gases input at the experiment's start (ppb)

Comp0 = Names of trace gases present at start (in the order corresponding to their 
		concentrations in C0).  Note, this is case sensitive, with the case matching that 
		in the xml file

voli = index of components with vapour pressures stated in volP (multiple indices allowed, can be absolute or relative)

volP = vapour pressures (Pa) of components with indices given in voli

pconc = total particle concentration at start of experiment (seed particles) 
		(# particles/cc (air))

std = geometric mean standard deviation of seed particle number concentration 
		(dimensionless)

loc = shift of lognormal probability distribution function for seed particles 
	number-size distribution (um)

scale = scaling factor of lognormal probability distribution function for seed 
	particles (dimensionless)

core_diss = core dissociation constant (for seed component) (dimensionless) (1)

light_time = times (s) for light status, corresponding to the elements of light_status
				(below), minimum requirement is start (0.0) and total time of the 
				experiment 

light_status = 1 for lights on and 0 for lights off, with times given in light_time 
				(above), minimum requirement is stating the status at the start and end
				time of the experiment
	
## Acknowledgements

This project has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 730997.
