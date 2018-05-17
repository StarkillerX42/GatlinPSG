# Wolf
This is a project analyzing exoplanet transit spectra using Eric T. Wolf's climate models and NASA's Planetary Spectrum Generator. The generator, found at https://psg.gsfc.nasa.gov/ was made cheifly by Geronimo L. Villanueva. Currently, this project has 3 sections on github, each with their own section explaining them here. All of them use PSG.py, which needs to be added to the $PYTHONPATH when running any versions of the PSG_Shells, where all the code is run. It is recommended to use jupyter notebooks for running the code as they offer some flexibility in analysis.

### Permissions
Currently, psginput files are available on github, which can be used to send to the PSG manually, but Eric's terminator profiles are not available, so the code will not run properly as-is.

## Running The PSG via PSG.py:
The PSG is a cool tool, but there are many steps which need to be completed before you put things in as inputs to the PSG, and there's a lot of work to do before you get reasonable results back, which is what this repository is for. In order to run the PSG, run the line
```foo = PSG.PSG(*args)
foo.calculate() # This line creates things like a planet_data dictionary, a star_data dictionary
# an astropy table of the atmosphere called atmosphere, and much more
foo.write() # This line creates a psginput.txt file which will be sent to the PSG
foo.send(*args) # This sends the file to the PSG using the curl API and deletes unwanted files
foo.plot_setup() # This prepares you for running the plot commands by unpacking the PSG output files
```

## Standard Transits
This is the main pipeline with nothing fancy added to it, focusing on the 1barN2 0.4barCO2 model using MIRI-MRS. It produces the most reliable results.

## Raw Reduction
This is an attempt at more closely matching an actual observation of TRAPPIST-1 by only using the total spectra received, instead of the cleaned data that PSG returns. Currently, it proves to be almost identical to Standard Transits

## Dry TRAPPIST-1 e
This file set has the abundances of any waters artificially set to zero, so the transit is waterless. This is an attempt of determining how measurable water is in the atmosphere.
