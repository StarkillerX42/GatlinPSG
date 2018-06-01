# Wolf
This is a project analyzing exoplanet transit spectra using Eric T. Wolf's climate models and NASA's Planetary Spectrum Generator. The generator, found at https://psg.gsfc.nasa.gov/ was made cheifly by Geronimo L. Villanueva. Currently, this project has 3 sections on github, each with their own section explaining them here. All of them use PSG.py, which needs to be added to the $PYTHONPATH when running any versions of the PSG_Shells, where all the code is run. It is recommended to use jupyter notebooks for running the code as they offer some flexibility in analysis.

### Permissions
Currently, psginput files are available on github, which can be used to send to the PSG manually, but Eric's terminator profiles are not available, so the code will not run properly as-is.

## Running The PSG via PSG.py:
The PSG is a cool tool, but there are many steps which need to be completed before you put things in as inputs to the PSG, and there's a lot of work to do before you get reasonable results back, which is what this repository is for. In order to run the PSG, run the line.

```
planet = PSG.PSG(planet_name, file_name, is_earth, astmosphere_ceiling, n_uplayers, phase)
# Initializes the PSG
planet.calculate(skprows)  # Computes necessary values
planet.write(scope, exposure_time, exposure_count, rad_units)
# Writes a file to be sent
planet.send(run)  # Sends the file to the PSG
planet.plot_setup()  # Prepares for plotting
planet.<plot_function>  # Each creates a plot and writes it to a file
```

## Running  PSG via jupyter notebook
In jupyter notebook, you have some added flexibility to interact with the results. Fortunately, PSG.py has been optimized for this purpose. In order to interact with the data, some python dictionaries are created to store meaningful values. `planet_data` and `star_data`. There is also an astropy table of the atmosphere profile called `atmosphere`. Those are the most useful values, but there is also `n_layers`, `n_downlayers`, and `n_uplayers`. Each of which describes the layers in each section of the atmosphere. The downlayers are given by the input files, the uplayers are isothermal layers in the upper atmosphere. A very useful technique is to run calculate, modify something about the atmosphere (like remove water), then run write, which will use the modified values.

After running plot_setup(), you will be able to access any results obtained by the PSG, so variables like Wavelengths, Transit, Stellar, Planet, and Noise will all be accessible. Attributes without a label in front of them came from the radiance file PSG returns, which includes all observables. Attributes with a t in front came from the transmittance file, which described the behavior of each particle in the atmosphere. Attributes with an n in front of it came from the noise file, which breaks down the noise into components.
## Standard Transits
This is the main pipeline with nothing fancy added to it, focusing on the 1barN2 0.4barCO2 model using MIRI-MRS. It produces the most reliable results.

## Raw Reduction
This is an attempt at more closely matching an actual observation of TRAPPIST-1 by only using the total spectra received, instead of the cleaned data that PSG returns. Currently, it proves to be almost identical to Standard Transits

## Dry TRAPPIST-1 e
This file set has the abundances of any waters artificially set to zero, so the transit is waterless. This is an attempt of determining how measurable water is in the atmosphere.
