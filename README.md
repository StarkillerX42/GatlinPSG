# GatlinPSG
This is a project analyzing exoplanet transit spectra using Eric T. Wolf's climate models and NASA's Planetary Spectrum Generator (PSG). The generator, found at https://psg.gsfc.nasa.gov/ was made cheifly by Geronimo L. Villanueva. This project was created as a wrapper for the PSG for 3D cliamte models (as NETCDFs). This work was used to produce my honors thesis on the subject.

### Permissions
Currently, psginput files are available on GitHub, which can be used to send to the PSG manually, but Eric's climate models are only available upon request. The PSG can still be run using the psginput.txt files, although some early steps in the pipeline will need to be run differently, see below.

## Running The PSG:
The PSG is a flexible and scientifically very useful tool, but it has undergone many changes since it has been created. The changes are seldom backwards compatible, and for that reason, a number of subdirectories here will not run as is, but may run with a little modification. The correct way to run the PSG is via the following code.

```
# Initializes the object planet
planet = PSG.PSG(planet_name)
# Downloads NASA Exoplanet archive to populate variables
planet.fetch_archive(is_earth)
# Read the netCDF file atmosphere model
planet.from_cdf(cdf_file, phase)
# Add necessary values and add the upper layers and produce psginput.txt
planet.calculate(atmosphere_ceiling, n_uplayers)
# Write a file to send to the NASA PSG
planet.write(scope, exp_time, n_exposures, rad_units)
# Send to the PSG
planet.send(keep_files, run)
# Prepare plots
planet.plot_setup()
planet.depth_plot()
```

## Integrating with Jupyter Notebook
In jupyter notebook, you have some added flexibility to interact with the results. PSG.py has been optimized for this purpose. In order to interact with the data, some python dictionaries are created to store meaningful values. `planet_data` and `star_data`. There is also an astropy table of the atmosphere profile called `atmosphere`. Those are the most useful values, but there is also `n_layers`, `n_downlayers`, and `n_uplayers`. Each of which describes the layers in each section of the atmosphere. The downlayers are given by the input files, the uplayers are isothermal layers in the upper atmosphere. A very useful technique is to run calculate, modify something about the atmosphere (like remove water), then run write, which will use the modified values.

After running plot_setup(), you will be able to access any results obtained by the PSG, so variables like Wavelengths, Transit, Stellar, Planet, and Noise will all be accessible. Attributes without a label in front of them came from the radiance file PSG returns, which includes all observables. Attributes with a t in front came from the transmittance file, which described the behavior of each particle in the atmosphere. Attributes with an n in front of it came from the noise file, which breaks down the noise into components.

## Directories

This project has a number of sub-directories, which are variably old. Some of them haven't been used since early version of the pipeline, and I make no promises about how they will behave. The most up to date directory is Thesis, which contians all the code used in my honors thesis.

### Thesis
The most polished directory with a number of scientific results discussed in my thesis. This is likely the most interesting as it has the most examples and the most meaningful results.

### Standard Transits
This is the main pipeline with nothing fancy added to it, focusing on the 1barN2 0.4barCO2 model using MIRI-MRS. It produces the most reliable results. This has not been used in a while since the creation of Thesis.

### Raw Reduction
This is an attempt at more closely matching an actual observation of TRAPPIST-1 by only using the total spectra received, instead of the cleaned data that PSG returns. Currently, it proves to be almost identical to Standard Transits

### Dry TRAPPIST-1 e
This file set has the abundances of any waters artificially set to zero, so the transit is waterless. This is an attempt of determining how measurable water is in the atmosphere.

### NIRCamBadhan
This was an attempt at reproducing results in Badhan 2018. Some meaningful results agree, most particularly their output spectra. Unfortunately this was abandoned due to other commitments at the time it was created

### PSGExample Earth
This represents an Earth-like input to the PSG

### ALMA
This was an attempt at computing thermal phase curves using ALMA. Interesting work, but incomplete and old

### NoiseAnalysis
This is my first introduction to analyzing PSG noise, and it is far from functional

### Phase Curves
This is previous work on phase curves and out of date relative to Thesis
