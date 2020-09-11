---
title: 'PyCHAM: CHemistry with Aerosol Microphysics in Python'
tags:
  - Python
  - atmospheric science
  - aerosol
  - aerosol chambers
  - air quality
  - climate change
  - EUROCHAMP
authors:
  - name: Simon O’Meara
    orcid: 0000-0001-7197-0651
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Shuxuan Xu
    affiliation: 1
  - name: David Topping
    orcid: 0000-0002-2104-4577
    affiliation: 1
  - name: Gerard Capes
    affiliation: 3
  - name: Douglas Lowe
    affiliation: 3
  - name: M. Rami Alfarra
    orcid: 0000-0002-3925-3780
    affiliation: "1, 2"
  - name: Gordon McFiggans
    orcid: 0000-0002-3423-7896
    affiliation: "1, 2"
affiliations:
 - name: Department of Earth and Environmental Sciences, University of Manchester
   index: 1
 - name: National Centre for Atmospheric Science
   index: 2
 - name: Research Computing Services, University of Manchester
   index: 3
date: 15 November 2019
bibliography: paper.bib
---

# Summary

PyCHAM (CHemistry with Aerosol Microphysics in Python) is an open-access 0-D box model for aerosol chamber studies.  Aerosol chambers provide a means for atmospheric scientists to interrogate aerosol processes [@Finlayson:2000; @Schwantes:2017; @Hidy:2019].  Allying PyCHAM with chamber measurements allows quantification of unknown parameters, such as branching ratios for oxidation schemes [@Chen:2005].  Although several box models have already been published [@Pierce:2008; @Roldin:2014; @Sunol:2018; @Roldin:2019], PyCHAM is novel in its ease of accessibility and utility.  To the best of our knowledge, PyBOX [@Topping:2018], which formed the basis of PyCHAM, and AtChem [@AtCHEM:2020], are the only other open source box models with graphical user interfaces.  However, PyCHAM provides a considerable upgrade to the functionality of both prior models, as AtChem only simulates gas-phase photochemistry and PyBOX has limited particle effects, whilst PyCHAM models significantly more particle microphysics (detailed below) and wall effects, in addition to gas-phase photochemistry.  The intention therefore, is that it will be readily employed by research groups undertaking aerosol chamber measurements.

With air quality and climate models increasingly important to guiding sustainable societies, the accuracy of simulations must suffice [@Tong:2019].  However, research shows that the simulated aerosol effects in these models provides a relatively high amount of uncertainty [@Johnson:2018].  The combination of box models like PyCHAM with aerosol chamber measurements to better constrain aerosol processes is therefore necessary to ultimately improve societal sustainability.

At the time of writing, PyCHAM is being used to investigate the autoxidation of organic vapours in the atmosphere.  The autoxidation process has recently been discovered to play a significant role in the formation of airbourne particulates [@Ehn:2014], however its exact chemical mechanism is yet to be elucidated.  Through comparison of chamber measurements with PyCHAM outputs using various mechanism possibilities, a constrained autoxidation chemical scheme is being generated.

The model employs non-equilibrium equations to simulate the known processes occurring in aerosol chambers.  At its core is integration of ordinary differential equations (ODEs) for gas-phase photochemistry and gas partitioning to particles and walls.  Here, the CVode function of the Assimulo package [@Andersson:2015] for ODE solvers is called on, using the backward differentiation formula of the Sundials solvers [@hindmarsh:2005sundials], which studies have shown is most reliable for solution of these equations [@Jacobson:2005].  The general equation for chemical reactions is [@Jacobson:2005]: 

$$\frac{d[i_{g}]}{dt} = \pm i_{s}k_{n}[a_{g}]^{a_{s}}[b_{g}]^{b_{s}},
$$

where square brackets represent concentrations, with $_{g}$ representing the gas phase, $i$ is the affected component, $t$ is time, $k$ is the reaction rate coefficient for reaction $n$, $a$ and $b$ are example reactants, with stoichiometries $_{s}$.  The equation is positive for products of reactions and negative for reactants.

The gas-particle partitioning equation is [@Jacobson:2005]:

$$\frac{d[i_{g}]}{dt} = -k_{i}([i_{g}]-p^{0}_{i}x_{i}K_{v}),
$$

where $k$ is the mass transfer coefficient, $p^{0}$ is the liquid saturation vapour pressure, $x$ is particle-phase mole fraction and $K_{v}$ is the Kelvin effect.

Gas partitioning to walls is an area of ongoing research [@Zhang:2015], therefore we provide an equation analogous to gas-particle partitioning given above:

$$\frac{d[i_{g}]}{dt} = -k_{gwt}C_{w}([i_{g}]-p^{0}_{i}\frac{[i_{w}]}{C_{w}}),
$$

where $_{w}$ represents the wall, $k_{gwt}$ is the gas-wall mass transfer coefficient, $C_{w}$ is the effective absorbing concentration of the wall.  Users set the values of $k_{gwt}$ and $C_{w}$ as wall effects vary significantly between chambers.

Outside the ODE solver, particle loss to walls, coagulation and nucleation are also solved, with equations for the former two given in @Jacobson:2005 and the parameters inside the nucleation expression (for relevent experiments) tuned by the user.  PyCHAM takes a sectional approach to particulates, dividing particles into a number of size bins and treating their changing size using the moving-centre approach [@Jacobson:2005].  It builds upon PyBOX [@Topping:2018], which did not include coagulation, nucleation or a sectional method.

Several variables change between different aerosol chambers and different experiments; therefore, the software is designed to allow the user to set these with ease.

# Acknowledgements

Funding has been provided by the EUROCHAMP-2020 research project [@EUROCHAMP:2020].  
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 730997 and the National Centre for Atmospheric Science.

# References
