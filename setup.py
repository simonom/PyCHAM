'''module to state PyCHAM package metadata'''
from setuptools import setup
from setuptools import find_packages


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="PyCHAM",
    version="2.1.9",
    author="Simon O'Meara",
    author_email="simon.omeara@manchester.ac.uk",
    description="PyCHAM: CHemistry with Aerosol Microphysics in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords=['aerosol', 'chamber', 'model'],
    url="https://github.com/simonom/PyCHAM",
    packages=find_packages(include=['PyCHAM', 'PyCHAM.*']),
    include_package_data=True,
    package_data={'PyCHAM': ['input/*', 'output/example_scheme/Example_Run_output/*', 'output/example_scheme_inc_aq/Example_Run_output_inc_aq/*', 'output/GMD_paper_plotting_scripts/*', 'output/plotting_scripts/*', 'photofiles/*', 'photofiles/MCMv3.2/*', 'photofiles/MCMv3.2/(CH3)3CONO2/*', 'photofiles/MCMv3.2/BIACET/*', 'photofiles/MCMv3.2/C2H5CHO/*', 'photofiles/MCMv3.2/C5HPALD1/*', 'photofiles/MCMv3.2/CH3CH2ONO2/*', 'photofiles/MCMv3.2/CH3CHO/*', 'photofiles/MCMv3.2/CH3COCH3/*', 'photofiles/MCMv3.2/(CH3COCHO/*', 'photofiles/MCMv3.2/CH3ONO2/*', 'photofiles/MCMv3.2/CH3OOH/*', 'photofiles/MCMv3.2/CHOCHO/*', 'photofiles/MCMv3.2/H2O2/*', 'photofiles/MCMv3.2/HCHO/*', 'photofiles/MCMv3.2/HNO3/*', 'photofiles/MCMv3.2/HONO/*', 'photofiles/MCMv3.2/i_C3H7CHO/*', 'photofiles/MCMv3.2/MACR/*', 'photofiles/MCMv3.2/MEK/*', 'photofiles/MCMv3.2/MVK/*', 'photofiles/MCMv3.2/n_C3H7CHO/*', 'photofiles/MCMv3.2/NO2/*', 'photofiles/MCMv3.2/NO3/*', 'photofiles/MCMv3.2/NOA/*', 'photofiles/MCMv3.2/O3/*', 'unit_tests/*', 'unit_tests/input/*']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    
# state dependencies for PyCHAM, note that openbabel is installed separately, as described in README
install_requires=["PyQt5","pickleshare","numpy","backcall","certifi","cycler","decorator","gitpython","ipdb","matplotlib","scipy","xmltodict","requests"]
)
