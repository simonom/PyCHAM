---
# title: 'PyCHAM: CHemistry with Aerosol Microphysics in Python'
tags:
  - Python
  - atmospheric science
  - aerosol
  - environmental chambers
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

 
#affiliations:
 - name: Department of Earth and Environmental Sciences, University of Manchester
   index: 1
 - name: National Centre for Atmospheric Science
   index: 2
 - name: Research Computing Services, University of Manchester}
   index: 3

date: 5 November 2019
bibliography: paper_refs.bib
---

# Summary

PyCHAM (CHemistry with Aerosol Microphysics in Python) is an open-access O-D box model for environmental chamber studies.  Environmental chambers provide a means for atmospheric scientists to interrogate aerosol processes [@Finlayson:2000, @Schwantes:2017, @Hidy:2019].  Allying PyCHAM with chamber measurements allows quantification of unknown parameters, such as branching ratios for oxidation schemes [@Chen:2005].  Although several box models have already been published [@Pierce:2008, @Roldin2014, @Sunol2018, @Roldin2019], PyCHAM is novel in its ease of accessibility and utility.  The intention therefore, is that it will be readily employed by research groups undertaking environmental chamber measurements.

With air quality and climate models increasingly important to guiding sustainable societies, the accuracy of simulations must suffice [e.g. @Tong:2019].  However, research shows that the simulated aerosol effects in these models provides a relatively high amount of uncertainty [@Johnson:2018].  The combination of box models like PyCHAM with environmental chamber measurements to better constrain aerosol processes is therefore necessary to ultimately improve societal sustainability.

Funding for model development has been provided by the EUROCHAMP-2020 research project [@EUROCHAMP2020].  At the time of writing, PyCHAM is being used to investigate the autoxidation of organic vapours in the atmosphere and the partitioning of semi-volatile vapours to soot aerosol.  Here is the repository: github.com/simonom/PyCHAM.

The model employs non-equilibrium equations to simulate the known processes occurring in environmental chambers: gas partitioning to particles and walls, particle loss to walls, gas-phase photochemistry, coagulation and nucleation.  It takes a sectional approach to particulates, dividing particles into a number of size bins and treating their changing size using the moving-centre approach [@Jacobson:2005].  It builds upon PyBOX [@Topping2018] which did not include coagulation, nucleation or a sectional method.

Several variables change between different environmental chambers and different experiments; therefore, the software is designed to allow the user to set these with ease.

Central to the model is solution of the ordinary differential equations for the rate of gas-phase chemistry, gas-particle partitioning and gas-wall partitioning.  Here, the CVode function of the Assimulo package for ODE solvers is called on [@Andersson2015], using the backward differentiation formula, which studies have shown is most reliable for solution of these equations [@Jacobson2005].

# Acknowledgements

This project has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 730997 and the National Centre for Atmospheric Science.

# References
