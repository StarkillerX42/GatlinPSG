import csv
import glob
import os
import time
import matplotlib.pyplot as plt
import numpy as np
import astropy.table

__author__ = "Dylan_Gatlin"
__version__ = 3.6


class PSG(object):
    """Parent class for all PSG Operations, including config file writing,
    sending files to the PSG, and plotting. To properly use this, run the
    commands in the following order: x=PSG.PSG(*args,**kwargs) x.calculate()
    # This function defines everything we need for PSG x.write()  # This
    function writes what PSG needs to a PSG-friendly file x.Run(**kwargs)  #
    This sends the data to PSG and gets a response x.plot_setup()  # This
    does the background unpacking to plot things x.<PlotFunction>  # Each of
    these makes a single plot where you replace <PlotFunction> with
    depth_plot, emission, absorption, raw, star, signal_and_noise,
    signal_noise_ratio, trn_total,trn_co2, trn_n2, trn_h2o, trn_ice,
    trn_cloud, trn_cia, trn_species,trn_aero,trn_all, noise_components. You
    can run all of them sequentially if intended, but not before plot_setup()

    Required Files, Arguments, and Keyword Arguments: PSG.PSG.write():

        This function will take a planet name and return a config file in the
        formats necessary to run in the NASA Goddard Planetary Spectrum
        Generator. The planet name must be in the same format as NASA's
        Exoplanet API archive. This means "TRAPPIST-1 e" is valid,
        but "trappist1e" is not.

        Required Files:

        In the folder you"re running this, one file must be accessible.
        file_name, a text file following Eric Wolf's formatting file_name must
        be the terminator profile ending in _terminator.txt or another similar
        formatted file. If the observation phase is different than 180, the
        corresponding atmosphere profile must be given.

        Arguments:

        planet: (str) The name of the planet according to NASA"s API.
        For a list of known planets, run this function with a dummy planet name
        and read through exoplanets.txt, which should appear in your directory.

        file_name: (str) The name of the file which you would like to read in,
        described above

        Keyword Arguments:

        scope: (str) The name of the scope you are observing with. "MIRI-LRS",
        "MIRI-MRS","NIRISS-SOSS","NIRCam-Grism","NIRSpec-1000",
        "NIRSpec-2700", "NIRSpec-Prism", "Hubble" and "Spitzer-Short-High"are
        currently supported. "MIRI-MRS" is the default

        is_earth: (bool)Whether or not this planet is a fake exoplanet with
        the same mass and radius as Earth If is_earth is True, the planet
        will be imagined as if it was around the star of the variable "planet"

        atmosphere_ceiling: (float)The pressure where the atmosphere ends
        The PSG will only produce useful results if there are layers at
        extremely low pressure_profile. If the atmosphere layers don't go to low
        pressure_profile, imaginary isothermal layers can be added. Default is
        zero, which will use only the given profile

        n_uplayers: (int) The number of isothermal layers between the top given
        layer and the top of the atmosphere

        exposure_time: (float) The exposure time of an image

        exposure_count: (integer) The number of exposures
        
        phase: (float) The angle of the observer, with 180 being a perfect
        transit or alternatively (str) "file", if the angle is in the file_name

        skprow: (int)The number of rows to skip in input file

        Returns: A file ending in psginput.txt in your directory.

    Required Files, Arguments, and Keyword Arguments: PSG.PSG.Run():
    This function send the the file file_name to the NASA GSFC PSG for analysis.
    It will not return anything, but will write files in the curent directory.
    If you want it to return a filetype it doesn"t, hashtag out one of the
    lower lines so that the necessary file isn"t deleted.

        Required Files:

        A single file _psginput.txt or any other config file if you specify
        self._psginput_name to be something else
        Arguments:
        None

        pln,cfg,lay,noi,rad,str,trn,all,err (bool) True if you would like
        that type of file as output, False if you do not. The most useful
        files default to True, the unusual ones default to False. noi, rad,
        and trn are all True pln, cfg,lay,str, all are all False. lay if
        probably the most useful of the ones not left True by default. err is
        only referenced if there is an error. It should be left True always.
        If there is a err.txt file, you should read it

    Required Files, Arguments, and Keyword Arguments: PSG.PSG.plot_setup():
        None, If you have run the previous commands in order, it will run

    Requires LaTeX installed, siunitx, cmbright, type1cm, l3kernel, l3packages,
    beamer, zhmetrics to run all the plots

    """

    def __init__(self, planet_name: str,
                 file_name: str,
                 scope: str = "MIRI-MRS",
                 is_earth: bool = False,
                 atmosphere_ceiling=0,
                 n_uplayers: int = 0,
                 exposure_time=1000,
                 exposure_count: int = 110,
                 phase=180):
        super(PSG, self).__init__()
        self.planet = planet_name
        self.file_name = file_name
        self.scope = scope
        self.is_earth = is_earth
        self.atmosphere_ceiling = atmosphere_ceiling
        self.n_uplayers = n_uplayers
        self.exposure_time = exposure_time
        self.exposure_count = exposure_count
        self.phase = phase

        # Planet Variables
        self.planet_data = {}
        self.star_data = {}
        self.atmosphere = None
        self.n_downlayers = None
        self.n_layers = None

        # PSG Variables
        self.is_transit = None
        self.returned_files = None
        # Plot variables
        self._psginput_name = None
        self._plot_range = None
        self._file_stem = None
        self._title_stem = None
        self.Wavelengths = None
        self.Total = None
        self.Noise = None
        self.Stellar = None
        self.Planet = None
        self.Transit = None
        self.Thermal = None
        self.tWavelengthSpec = None
        self.tTotalSpec = None
        self.tN2Spec = None
        self.tCO2Spec = None
        self.tH2OSpec = None
        self.tIceSpec = None
        self.tCloudSpec = None
        self.tCIASpec = None
        self.nWavelengthSpec = None
        self.nTotal = None
        self.nSource = None
        self.nDetector = None
        self.nTelescope = None
        self.nBackground = None

    def calculate(self, skprow: int = 11):

        """See PSG parent class docstring for details"""
        print("Calculating Planet Data")
        # See if we need a new/updated exoplanet"s list
        if not os.path.isfile("exoplanets.txt"):
            need_file = True
        else:
            st = os.stat("exoplanets.txt")
            age = time.time() - st.st_mtime
            if age > (3600 * 24 * 2.0):
                need_file = True
            else:
                need_file = False
        if need_file:
            print("Retrieving planet variables from NASA's Exoplanet Archive")
            import urllib3
            import certifi
            http = urllib3.PoolManager(cert_reqs="CERT_REQUIRED",
                                       ca_certs=certifi.where())
            r = http.request("GET",
                             "https://exoplanetarchive.ipac.caltech.edu/cgi"
                             "-bin/nstedAPI/nph-nstedAPI?table=exoplanets"
                             "&select=pl_name,pl_masse,pl_rade,pl_orbsmax,"
                             "pl_orbincl,pl_trandep,st_teff,st_rad,st_radv,"
                             "st_dist,st_optmag,pl_insol&format=csv")
            lines = str(r.data)[2:].split("\\n")
            with open("exoplanets.txt", "w") as fil:
                for line in lines:
                    fil.write(line + "\n")
        # Extracts necessary details about the exoplanet from NASA"s API
        (planet_name, planet_emass, planet_erad,
         sma, inclination, transit_depth,
         star_temp, star_srad, star_velocity,
         star_distance, star_magnitude,
         insolation) = [None] * 12
        with open("exoplanets.txt", "r") as file:
            planet_list = csv.reader(file, delimiter=",")
            for line in planet_list:
                if line[0] == self.planet:
                    (planet_name, planet_emass, planet_erad,
                     sma, inclination, transit_depth,
                     star_temp, star_srad, star_velocity,
                     star_distance, star_magnitude,
                     insolation) = line
        if planet_name == "":
            print("    Planet not found")
            exit()
        if star_velocity == "":
            star_velocity = 0.
        if star_magnitude == "nan":
            star_magnitude = 0.
        if self.is_earth:
            planet_name = "ExoEarth like " + planet_name
            planet_emass = 1.
            planet_erad = 1.
            transit_depth = str(float(planet_erad)
                                / float(star_srad) * 0.009154)
        # Converts NASA"s values to proper units
        self.planet_data["Name"] = planet_name
        self.planet_data["Mass"] = float(planet_emass) * 5.9736e24
        self.planet_data["ERadius"] = float(planet_erad)
        self.planet_data["Radius"] = float(planet_erad) * 6371.0
        self.planet_data["Diameter"] = self.planet_data["Radius"] * 2
        self.planet_data["SemiMajorAxis"] = float(sma)
        self.planet_data["Inclination"] = float(inclination) * 1361.
        self.planet_data["TransitDepth"] = float(transit_depth) / 100
        self.planet_data["Insolation"] = float(insolation)
        self.star_data["SRadius"] = float(star_srad)
        self.star_data["Radius"] = self.star_data["SRadius"] * 1.39e6
        self.star_data["Velocity"] = float(star_velocity)
        self.star_data["Temperature"] = float(star_temp)
        self.star_data["Distance"] = float(star_distance)
        self.star_data["Magnitude"] = float(star_magnitude)
        self.planet_data["DepthEff"] = (self.planet_data["Radius"]
                                        / self.star_data["Radius"]) ** 2
        self.planet_data["Gravity"] = (6.67e-11
                                       * self.planet_data["Mass"]
                                       / (self.planet_data["Radius"]*1000) ** 2)
        if 2500 <= self.star_data["Temperature"] <= 3800:
            self.star_data["Type"] = "M"
        elif 3800 <= self.star_data["Temperature"] <= 5240:
            self.star_data["Type"] = "K"
        elif 5240 <= self.star_data["Temperature"] <= 5930:
            self.star_data["Type"] = "G"
        elif 5930 <= self.star_data["Temperature"] <= 7020:
            self.star_data["Type"] = "F"
        else:
            print("    star temp: {}".format(self.star_data["Temperature"]))
            self.star_data["Type"] = input("    Star type not recognized,"
                                           " please specify star type:")
        # Relevant if there are phase curves, this method isn't very good
        # unless you're using the terminator
        if "lon" not in self.file_name:
            if 178 < self.phase <= 182:
                self.file_name = self.file_name
                print("    Using terminator profile")
            if 135 < self.phase <= 178 or 182 < self.phase <= 225:
                self.file_name = self.file_name.split("_term")[0] \
                                 + "_antistellar.txt "
                print("    Using antistellar profile")
            if 45 < self.phase <= 135 or 225 < self.phase <= 315:
                self.file_name = self.file_name.split("_term")[0] \
                                 + "_globalmean.txt"
                print("    Using globalmean profile")
            if 0 <= self.phase <= 45 or 315 < self.phase <= 360:
                self.file_name = self.file_name.split("_term")[0] \
                                 + "_substellar.txt"
                print("    Using substellar profile")
        elif "lon" in self.file_name:
            print("    Using provided file, phase is {}".format(self.phase))

        # GCM Inputs
        with open(self.file_name) as fp:
            for i, line in enumerate(fp):
                if i == skprow - 2:
                    b = line
                    c = b.split()
                    self.planet_data["SurfacePressure"] = float(c[0]) / 1000.0
                    self.planet_data["SurfaceTemperature"] = float(c[1])
                    self.planet_data["Albedo"] = float(c[2])
                    try:
                        self.planet_data["MWeightDry"] = float(c[3])
                    except IndexError:
                        print("    No dry molecular weight")
                        exit()
        self.planet_data["EffectiveTemp"] = (self.planet_data["Insolation"]
                                             * (1 - self.planet_data["Albedo"])
                                             / (4 * 5.67e-8)) ** 0.25
        self.planet_data["ScaleHeight"] = (8.3144598
                                           * self.planet_data[
                                               "SurfaceTemperature"]
                                           / (self.planet_data["MWeightDry"]
                                              * self.planet_data["Gravity"]))

        a = np.loadtxt(self.file_name, skiprows=skprow, unpack=True)
        layers_profile = a[0] + 1
        self.n_downlayers = len(layers_profile)
        heights_profile = np.flipud(a[1])  # m
        pressure_profile = np.flipud(a[2] / 1000)  # bar
        temperature_profile = np.flipud(a[3])  # K
        n2_profile = np.flipud(a[4])  # Kg/Kg
        co2_profile = np.flipud(a[5])  # Kg/Kg
        h2o_profile = np.flipud(a[6])  # Kg/Kg
        liquid_cloud_profile = np.flipud(a[7])  # Kg/Kg
        ice_cloud_profile = np.flipud(a[8])  # Kg/Kg
        liquid_cloud_size_profile = np.flipud(a[9])  # um
        ice_cloud_size_profile = np.flipud(a[10])  # um
        self.atmosphere = astropy.table.Table([layers_profile,
                                               heights_profile,
                                               pressure_profile,
                                               temperature_profile,
                                               n2_profile,
                                               co2_profile,
                                               h2o_profile,
                                               liquid_cloud_profile,
                                               ice_cloud_profile,
                                               liquid_cloud_size_profile,
                                               ice_cloud_size_profile],
                                              names=("Layer", "Height",
                                                     "Pressure",
                                                     "Temperature",
                                                     "N2", "CO2", "H2O",
                                                     "LiquidCloud",
                                                     "IceCloud",
                                                     "LiquidCloudSize",
                                                     "IceCloudSize"))

        h0 = 0.0
        p0 = self.planet_data["SurfacePressure"]
        t0 = self.planet_data["SurfaceTemperature"]
        ntot = 0.0
        mtot = 0.0
        nlyrs = []
        for i, lvl in enumerate(self.atmosphere):
            # Height of the layer [m] is the height of the middle
            # of the layer
            dh = lvl["Height"] - h0
            # Layer molecular density [mol/m3]
            mol = 1e5 * (p0 + lvl["Pressure"]) / (
                    8.3144598 * (t0 + lvl["Temperature"]))
            nlyr = mol * dh  # Layer integrated density [mol/m2]
            # Layer molecular weight [g/mol]
            mlyr = (lvl["N2"] * 28.0134
                    + lvl["CO2"] * 44.01
                    + lvl["H2O"] * 18.01528
                    + lvl["LiquidCloud"] * 18.01528
                    + lvl["IceCloud"] * 18.01528)
            # N2 abundance  [mol/mol]
            self.atmosphere["N2"][i] = lvl["N2"] * (mlyr / 28.0134)
            # CO2 abundance [mol/mol]
            self.atmosphere["CO2"][i] = lvl["CO2"] * (mlyr / 44.01)
            # H2O abundance [mol/mol]
            self.atmosphere["H2O"][i] = lvl["H2O"] * (mlyr / 18.01528)
            # Liq abundance [mol/mol]
            self.atmosphere["LiquidCloud"][i] = (lvl["LiquidCloud"]
                                                 * (mlyr / 18.01528))
            # ice abundance [mol/mol]
            self.atmosphere["IceCloud"][i] = (lvl["IceCloud"] * (
                    mlyr / 18.01528))
            h0 += dh
            p0 = self.atmosphere["Pressure"][i]
            t0 = self.atmosphere["Temperature"][i]
            ntot += nlyr
            mtot += mlyr * nlyr
            nlyrs.append(nlyr)
        self.atmosphere["Height"] /= 1e3

        # Compute representative quantities
        self.planet_data["MWeightTotal"] = mtot / ntot
        # This needs to exist in case aerosols are absent,
        # ATMOSPHERE-AABUN must be 0
        self.planet_data["LiquidCloudAbundance"] = int(
            any(i > 0. for i in self.atmosphere["LiquidCloud"]))
        self.planet_data["IceCloudAbundance"] = int(
            any(i > 0. for i in self.atmosphere["IceCloud"]))
        self.planet_data["MeanLiquidCloudSize"] = np.average(
            self.atmosphere["LiquidCloudSize"],
            weights=self.atmosphere["LiquidCloudSize"])
        self.planet_data["MeanIceCloudSize"] = np.average(
            self.atmosphere["IceCloudSize"],
            weights=self.atmosphere["IceCloudSize"])

        # Adding layers
        if self.atmosphere_ceiling != 0:
            # This handles the imaginary upper atmosphere levels from
            # the top layer to the specified pressure atmosphere_ceiling
            # Defines upper atmosphere abundances
            ind = self.n_downlayers - 1
            (lyr, upheight, uppres, uptemp, upN2, upCO2, upH2O, upLiq,
             upIce, uplsize, upisize) = self.atmosphere[ind]
            # Planet radius at top of atmosphere
            # Hypsometric equation to find top of atmosphere
            uptop = (-np.log(self.atmosphere_ceiling / uppres)
                     * self.planet_data["ScaleHeight"]
                     + upheight)

            uplayers = np.geomspace(upheight+2, uptop, self.n_uplayers)
            uppressure = uppres * np.exp(-(uplayers - upheight)
                                         / self.planet_data["ScaleHeight"])
            i = 1
            for h, p in zip(uplayers, uppressure):
                self.atmosphere.add_row([self.n_downlayers+i, h, p, uptemp,
                                         upN2, upCO2, upH2O, upLiq, upIce,
                                         uplsize, upisize])
                i += 1
            self.n_layers = self.n_downlayers + self.n_uplayers

    def write(self):
        # This concludes the calculations, the rest is psg input writing
        self._psginput_name = ""
        if "aqua" in self.file_name:
            self._psginput_name += self.file_name.split("_aqua")[0]
            end = "_psginput.txt"
        else:
            print("    Assuming the file is a dry file")
            self._psginput_name = self.file_name.split("_dry")[0]
            end = "_dry_psginput.txt"
        if self.phase != 180:
            self._psginput_name += "_" + str(self.phase)
        if self.scope != "MIRI-MRS":
            self._psginput_name += "_" + self.scope
        self._psginput_name += end
        with open(self._psginput_name, "w") as results:

            # OBJECT
            results.write("<OBJECT-SUMMARY>Exoplanet transit simulation\n")
            results.write("<OBJECT>Exoplanet\n")
            results.write("<OBJECT-TYPE>Exoplanet\n")
            results.write("<OBJECT-NAME>{}\n".format(self.planet_data["Name"]))
            results.write("<OBJECT-DATE>2017/05/25 15:22\n")
            results.write("<OBJECT-DIAMETER>{}\n".format(
                self.planet_data["Diameter"]))
            results.write("<OBJECT-GRAVITY>{}\n".format(
                self.planet_data["Mass"]))
            results.write("<OBJECT-GRAVITY-UNIT>kg\n")
            results.write("<OBJECT-STAR-DISTANCE>{}\n".format(
                self.planet_data["SemiMajorAxis"]))
            results.write("<OBJECT-STAR-VELOCITY>0\n")
            results.write("<OBJECT-SOLAR-LONGITUDE>0\n")
            results.write("<OBJECT-SOLAR-LATITUDE>0\n")
            results.write("<OBJECT-SEASON>{}\n".format(self.phase))
            results.write("<OBJECT-SEASON-YEAR>2016\n")
            results.write("<OBJECT-STAR-TYPE>{}\n".format(
                self.star_data["Type"]))
            results.write("<OBJECT-STAR-TEMPERATURE>{}\n".format(
                self.star_data["Temperature"]))
            results.write("<OBJECT-STAR-RADIUS>{}\n".format(
                self.star_data["SRadius"]))
            results.write("<OBJECT-OBS-VELOCITY>{}\n".format(
                self.star_data["Velocity"]))
            results.write("<OBJECT-OBS-LONGITUDE>0.0\n")
            results.write("<OBJECT-OBS-LATITUDE>{}\n".format(
                self.planet_data["Inclination"]))

            # GEOMETRY
            results.write("<GEOMETRY-SUMMARY>Variables taken from 'Exoplanet "
                          "from JWST' config file\n")
            results.write("<GEOMETRY>Observatory\n")
            results.write("<GEOMETRY-OFFSET-UNIT>arcsec\n")
            results.write("<GEOMETRY-OBS-ALTITUDE>{}\n".format(
                self.star_data["Distance"]))
            results.write("<GEOMETRY-ALTITUDE-UNIT>pc\n")
            results.write("<GEOMETRY-DIAMETER-FOV>0.4\n")
            results.write("<GEOMETRY-FOV-UNIT>arcsec\n")
            results.write("<GEOMETRY-OBS-ANGLE>0.0\n")
            results.write("<GEOMETRY-SOLAR-ANGLE>180.0\n")
            results.write("<GEOMETRY-PLANET-FRACTION>1\n")
            results.write("<GEOMETRY-STAR-DISTANCE>0.0\n")
            results.write("<GEOMETRY-STAR-FRACTION>{}\n".format(
                self.planet_data["TransitDepth"]))
            results.write("<GEOMETRY-SOLAR-ANGLE>180\n")
            results.write("<GEOMETRY-PHASE>180\n")

            # GENERATOR
            if self.scope == "MIRI-LRS":
                results.write(
                    "<GENERATOR-INSTRUMENT>JWST_MIRI-LRS: JWST/MIRI is the "
                    "Mid-Infrared Instrument onboard the James Webb Space "
                    "Telescope. The low resolution spectroscopy (LRS) "
                    "configuration samples the 5 to 12 um spectral region "
                    "with a resolving power of 100 (at 7.5um), "
                    "and observations can be performed with a 0.5um slit or "
                    "via slitless spectroscopy. Fringing and other systematic "
                    "sources of noise have been identified in this "
                    "instrument, yet these effects are not included in the "
                    "noise simulator.\n")
                results.write("<GENERATOR-RANGE1>5.00\n")
                results.write("<GENERATOR-RANGE2>12.00\n")
                results.write("<GENERATOR-RANGEUNIT>um\n")
                results.write("<GENERATOR-RESOLUTION>100\n")
                results.write("<GENERATOR-RESOLUTIONUNIT>RP\n")
                results.write("<GENERATOR-TELESCOPE>SINGLE\n")
                results.write("<GENERATOR-DIAMTELE>6.5\n")
                results.write("<GENERATOR-BEAM>1.0\n")
                results.write("<GENERATOR-BEAM-UNIT>diffrac\n")
                results.write("<GENERATOR-TELESCOPE1>1\n")
                results.write("<GENERATOR-TELESCOPE2>0.0\n")
                results.write("<GENERATOR-TELESCOPE3>1.0\n")
                results.write("<GENERATOR-NOISE>CCD\n")
                results.write("<GENERATOR-NOISE1>32.6\n")
                results.write("<GENERATOR-NOISE2>0.005\n")
                results.write("<GENERATOR-NOISEOTEMP>50\n")
                results.write("<GENERATOR-NOISEOEFF>0.3\n")
                results.write("<GENERATOR-NOISEOEMIS>0.10\n")
                results.write("<GENERATOR-NOISETIME>{}\n"
                              .format(self.exposure_time))
                results.write("<GENERATOR-NOISEFRAMES>{}\n"
                              .format(self.exposure_count))
                results.write("<GENERATOR-NOISEPIXELS>8\n")
                results.write("<GENERATOR-TRANS-APPLY>N\n")
                results.write("<GENERATOR-TRANS-SHOW>N\n")
                results.write("<GENERATOR-LOGRAD>Y\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>Wsrm2um\n")
                results.write(
                    "<GENERATOR-SUMMARY>Wavelengths range 5.00-12.00 um with a "
                    "resolution of 100 RP. Molecular radiative-transfer "
                    "enabled; Continuum flux module enabled;\n")
                results.write("<GENERATOR-TRANS>02-01\n")

            elif self.scope == "MIRI-MRS":
                results.write(
                    "<GENERATOR-INSTRUMENT>JWST_MIRI-MRS: JWST/MIRI is the "
                    "Mid-Infrared Instrument onboard the James Webb Space "
                    "Telescope. The medium resolution spectroscopy (MRS) "
                    "configuration covers wavelengths from 5 to 28.3 um, "
                    "enabled by four Integral Field Units (IFUs) with FOVs of "
                    "3.9"" to 7.7"" (#1: 5 to 7.7 um, #2: 7.7 to 11.9 um, "
                    "#3: 11.9 to 18.4 um, #4: 18.4 to 28.3 um). This "
                    "configuration has an effective throughput of 0.3 and an "
                    "average resolving power of 2400.\n")
                results.write("<GENERATOR-RANGE1>5.00\n")
                results.write("<GENERATOR-RANGE2>28.30\n")
                results.write("<GENERATOR-RANGEUNIT>um\n")
                results.write("<GENERATOR-RESOLUTION>2400\n")
                results.write("<GENERATOR-RESOLUTIONUNIT>RP\n")
                results.write("<GENERATOR-TELESCOPE>SINGLE\n")
                results.write("<GENERATOR-DIAMTELE>6.5\n")
                results.write("<GENERATOR-BEAM>1.0\n")
                results.write("<GENERATOR-BEAM-UNIT>diffrac\n")
                results.write("<GENERATOR-TELESCOPE1>1\n")
                results.write("<GENERATOR-TELESCOPE2>0.0\n")
                results.write("<GENERATOR-TELESCOPE3>1.0\n")
                results.write("<GENERATOR-NOISE>CCD\n")
                results.write("<GENERATOR-NOISE1>32.6\n")
                results.write("<GENERATOR-NOISE2>0.005\n")
                results.write("<GENERATOR-NOISEOTEMP>50\n")
                results.write("<GENERATOR-NOISEOEFF>0.3\n")
                results.write("<GENERATOR-NOISEOEMIS>0.10\n")
                results.write("<GENERATOR-NOISETIME>{}\n".format(
                    self.exposure_time))
                results.write("<GENERATOR-NOISEFRAMES>{}\n".format(
                    self.exposure_count))
                results.write("<GENERATOR-NOISEPIXELS>8\n")
                results.write("<GENERATOR-TRANS-APPLY>N\n")
                results.write("<GENERATOR-TRANS-SHOW>N\n")
                results.write("<GENERATOR-LOGRAD>Y\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>Wsrm2um\n")
                results.write(
                    "<GENERATOR-SUMMARY>Wavelengths range 5.00-28.30 um with a "
                    "resolution of 2400 RP. Molecular radiative-transfer "
                    "enabled; Continuum flux module enabled;\n")
                results.write("<GENERATOR-TRANS>02-01\n")

            elif self.scope == "NIRCam-Grism":
                results.write(
                    "<GENERATOR-INSTRUMENT>JWST_NIRCam-Grism: JWST/NIRCam is "
                    "the Near Infrared Camera onboard the James Webb Space "
                    "Telescope. While mainly designed for imagining, "
                    "this instrument enables slitless spectroscopy by "
                    "employing a grism. Spectroscopy is performed on the long "
                    "wavelength channel with two filters F322W2 (2.5 to 4 um) "
                    "and F444W (3.9 to 5 um), with an effective throughput of "
                    "0.3 and an average resolving power of 1600.\n")
                results.write("<GENERATOR-RANGE1>2.50\n")
                results.write("<GENERATOR-RANGE2>5.00\n")
                results.write("<GENERATOR-RANGEUNIT>um\n")
                results.write("<GENERATOR-RESOLUTION>1600\n")
                results.write("<GENERATOR-RESOLUTIONUNIT>RP\n")
                results.write("<GENERATOR-TELESCOPE>SINGLE\n")
                results.write("<GENERATOR-DIAMTELE>6.5\n")
                results.write("<GENERATOR-BEAM>1.0\n")
                results.write("<GENERATOR-BEAM-UNIT>diffrac\n")
                results.write("<GENERATOR-TELESCOPE1>1\n")
                results.write("<GENERATOR-TELESCOPE2>0.0\n")
                results.write("<GENERATOR-TELESCOPE3>1.0\n")
                results.write("<GENERATOR-NOISE>CCD\n")
                results.write("<GENERATOR-NOISE1>10.96\n")
                results.write("<GENERATOR-NOISE2>0.005\n")
                results.write("<GENERATOR-NOISEOTEMP>50\n")
                results.write("<GENERATOR-NOISEOEFF>0.3\n")
                results.write("<GENERATOR-NOISEOEMIS>0.10\n")
                results.write("<GENERATOR-NOISETIME>{}\n".format(
                    self.exposure_time))
                results.write("<GENERATOR-NOISEFRAMES>{}\n".format(
                    self.exposure_count))
                results.write("<GENERATOR-NOISEPIXELS>8\n")
                results.write("<GENERATOR-TRANS-APPLY>N\n")
                results.write("<GENERATOR-TRANS-SHOW>N\n")
                results.write("<GENERATOR-LOGRAD>Y\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>Wsrm2um\n")
                results.write(
                    "<GENERATOR-SUMMARY>Wavelengths range 2.50-5.00 um"
                    "with a resolution of 1600 RP. Molecular radiative-transfer"
                    " enabled; Continuum flux module enabled;\n")
                results.write("<GENERATOR-TRANS>02-01\n")

            elif self.scope == "NIRISS-SOSS":
                results.write(
                    "<GENERATOR-INSTRUMENT>JWST_NIRISS-SOSS: JWST/NIRISS is "
                    "the Near Infrared Imager and Slitless Spectrograph ("
                    "NIRISS). The single-object slitless-spectroscopy (SOSS) "
                    "configuration enables medium-resolution (R~700)ang "
                    "spectroscopy at 0.6um-2.8um, in 3 cross-dispersed orders "
                    "for a single bright target. Efficiencies / throughputs "
                    "vary substantially (0.1 to 0.6) across the main two "
                    "orders, with an equivalent throughput of 0.4 considered "
                    "for this configuration.\n")
                results.write("<GENERATOR-RANGE1>0.60\n")
                results.write("<GENERATOR-RANGE2>2.80\n")
                results.write("<GENERATOR-RANGEUNIT>um\n")
                results.write("<GENERATOR-RESOLUTION>700\n")
                results.write("<GENERATOR-RESOLUTIONUNIT>RP\n")
                results.write("<GENERATOR-TELESCOPE>SINGLE\n")
                results.write("<GENERATOR-DIAMTELE>6.5\n")
                results.write("<GENERATOR-BEAM>1.0\n")
                results.write("<GENERATOR-BEAM-UNIT>diffrac\n")
                results.write("<GENERATOR-TELESCOPE1>1\n")
                results.write("<GENERATOR-TELESCOPE2>0.0\n")
                results.write("<GENERATOR-TELESCOPE3>1.0\n")
                results.write("<GENERATOR-NOISE>CCD\n")
                results.write("<GENERATOR-NOISE1>11.55\n")
                results.write("<GENERATOR-NOISE2>0.005\n")
                results.write("<GENERATOR-NOISEOTEMP>50\n")
                results.write("<GENERATOR-NOISEOEFF>0.4\n")
                results.write("<GENERATOR-NOISEOEMIS>0.10\n")
                results.write("<GENERATOR-NOISETIME>{}\n".format(
                    self.exposure_time))
                results.write("<GENERATOR-NOISEFRAMES>{}\n".format(
                    self.exposure_count))
                results.write("<GENERATOR-NOISEPIXELS>8\n")
                results.write("<GENERATOR-TRANS-APPLY>N\n")
                results.write("<GENERATOR-TRANS-SHOW>N\n")
                results.write("<GENERATOR-LOGRAD>Y\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>Wsrm2um\n")
                results.write(
                    "<GENERATOR-SUMMARY>Wavelengths range 0.60-2.80 um with a "
                    "resolution of 700 RP. Molecular radiative-transfer "
                    "enabled; Continuum flux module enabled;\n")
                results.write("<GENERATOR-TRANS>02-01\n")

            elif self.scope == "NIRSpec-1000":
                results.write(
                    "<GENERATOR-INSTRUMENT>JWST_NIRSpec-1000: JWST/NIRSpec is "
                    "the infrared spectrometer onboard the James Webb Space "
                    "Telescope. The mid-resolution configuration (RP~1000) "
                    "considers the mid resolving power configuration for each "
                    "of the available instrument gratings/filters "
                    "combinations (G140M+F070LP, G140M+F100LP, G235M+F170LP, "
                    "G395M+F290LP). The three gratings cover the 1.0 to 5.3 "
                    "um spectral region with an average resolving power of "
                    "1000.\n")
                results.write("<GENERATOR-RANGE1>1.00\n")
                results.write("<GENERATOR-RANGE2>5.30\n")
                results.write("<GENERATOR-RANGEUNIT>um\n")
                results.write("<GENERATOR-RESOLUTION>1000\n")
                results.write("<GENERATOR-RESOLUTIONUNIT>RP\n")
                results.write("<GENERATOR-TELESCOPE>SINGLE\n")
                results.write("<GENERATOR-DIAMTELE>6.5\n")
                results.write("<GENERATOR-BEAM>1.0\n")
                results.write("<GENERATOR-BEAM-UNIT>diffrac\n")
                results.write("<GENERATOR-TELESCOPE1>1\n")
                results.write("<GENERATOR-TELESCOPE2>0.0\n")
                results.write("<GENERATOR-TELESCOPE3>1.0\n")
                results.write("<GENERATOR-NOISE>CCD\n")
                results.write("<GENERATOR-NOISE1>16.8\n")
                results.write("<GENERATOR-NOISE2>0.005\n")
                results.write("<GENERATOR-NOISEOTEMP>50\n")
                results.write("<GENERATOR-NOISEOEFF>0.3\n")
                results.write("<GENERATOR-NOISEOEMIS>0.10\n")
                results.write("<GENERATOR-NOISETIME>{}\n".format(
                    self.exposure_time))
                results.write("<GENERATOR-NOISEFRAMES>{}\n".format(
                    self.exposure_count))
                results.write("<GENERATOR-NOISEPIXELS>8\n")
                results.write("<GENERATOR-TRANS-APPLY>N\n")
                results.write("<GENERATOR-TRANS-SHOW>N\n")
                results.write("<GENERATOR-LOGRAD>Y\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>Wsrm2um\n")
                results.write(
                    "<GENERATOR-SUMMARY>Wavelengths range 1.00-5.30 um with a "
                    "resolution of 1000 RP. Molecular radiative-transfer "
                    "enabled; Continuum flux module enabled;\n")
                results.write("<GENERATOR-TRANS>02-01\n")

            elif self.scope == "NIRSpec-2700":
                results.write(
                    "<GENERATOR-INSTRUMENT>JWST_NIRSpec-2700: JWST/NIRSpec is "
                    "the infrared spectrometer onboard the James Webb Space "
                    "Telescope. The high-resolution configuration (RP~2700) "
                    "considers the maximum resolving power for each of the "
                    "available instrument gratings/filters combinations ("
                    "G140H+F070LP, G140H+F100LP, G235H+F170LP, G395H+F290LP). "
                    "The three gratings cover the 1.0 to 5.3 um spectral "
                    "region with an average resolving power of 2700.\n")
                results.write("<GENERATOR-RANGE1>1.00\n")
                results.write("<GENERATOR-RANGE2>5.30\n")
                results.write("<GENERATOR-RANGEUNIT>um\n")
                results.write("<GENERATOR-RESOLUTION>2700\n")
                results.write("<GENERATOR-RESOLUTIONUNIT>RP\n")
                results.write("<GENERATOR-TELESCOPE>SINGLE\n")
                results.write("<GENERATOR-DIAMTELE>6.5\n")
                results.write("<GENERATOR-BEAM>1.0\n")
                results.write("<GENERATOR-BEAM-UNIT>diffrac\n")
                results.write("<GENERATOR-TELESCOPE1>1\n")
                results.write("<GENERATOR-TELESCOPE2>0.0\n")
                results.write("<GENERATOR-TELESCOPE3>1.0\n")
                results.write("<GENERATOR-NOISE>CCD\n")
                results.write("<GENERATOR-NOISE1>16.8\n")
                results.write("<GENERATOR-NOISE2>0.005\n")
                results.write("<GENERATOR-NOISEOTEMP>50\n")
                results.write("<GENERATOR-NOISEOEFF>0.3\n")
                results.write("<GENERATOR-NOISEOEMIS>0.10\n")
                results.write("<GENERATOR-NOISETIME>{}\n".format(
                    self.exposure_time))
                results.write("<GENERATOR-NOISEFRAMES>{}\n".format(
                    self.exposure_count))
                results.write("<GENERATOR-NOISEPIXELS>8\n")
                results.write("<GENERATOR-TRANS-APPLY>N\n")
                results.write("<GENERATOR-TRANS-SHOW>N\n")
                results.write("<GENERATOR-LOGRAD>Y\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>Wsrm2um\n")
                results.write(
                    "<GENERATOR-SUMMARY>Wavelengths range 1.00-5.30 um with a "
                    "resolution of 2700 RP. Molecular radiative-transfer "
                    "enabled; Continuum flux module enabled;\n")
                results.write("<GENERATOR-TRANS>02-01\n")

            elif self.scope == "NIRSpec-Prism":
                results.write(
                    "<GENERATOR-INSTRUMENT>JWST_NIRSpec-Prism: JWST/NIRSpec "
                    "is the infrared spectrometer onboard the James Webb "
                    "Space Telescope. This configuration considers the low "
                    "resolving power configuration achieved using a prism ("
                    "Prism+Clear). This settings covers the 0.7 to 5.0 um "
                    "spectral region with an average resolving power of 100. "
                    "This configuration is only recommended for "
                    "point-sources, and the instrumental throughput drops "
                    "substantially below 1um.\n")
                results.write("<GENERATOR-RANGE1>0.70\n")
                results.write("<GENERATOR-RANGE2>5.00\n")
                results.write("<GENERATOR-RANGEUNIT>um\n")
                results.write("<GENERATOR-RESOLUTION>100\n")
                results.write("<GENERATOR-RESOLUTIONUNIT>RP\n")
                results.write("<GENERATOR-TELESCOPE>SINGLE\n")
                results.write("<GENERATOR-DIAMTELE>6.5\n")
                results.write("<GENERATOR-BEAM>1.0\n")
                results.write("<GENERATOR-BEAM-UNIT>diffrac\n")
                results.write("<GENERATOR-TELESCOPE1>1\n")
                results.write("<GENERATOR-TELESCOPE2>0.0\n")
                results.write("<GENERATOR-TELESCOPE3>1.0\n")
                results.write("<GENERATOR-NOISE>CCD\n")
                results.write("<GENERATOR-NOISE1>16.8\n")
                results.write("<GENERATOR-NOISE2>0.005\n")
                results.write("<GENERATOR-NOISEOTEMP>50\n")
                results.write("<GENERATOR-NOISEOEFF>0.4\n")
                results.write("<GENERATOR-NOISEOEMIS>0.10\n")
                results.write("<GENERATOR-NOISETIME>{}\n".format(
                    self.exposure_time))
                results.write("<GENERATOR-NOISEFRAMES>{}\n".format(
                    self.exposure_count))
                results.write("<GENERATOR-NOISEPIXELS>8\n")
                results.write("<GENERATOR-TRANS-APPLY>N\n")
                results.write("<GENERATOR-TRANS-SHOW>N\n")
                results.write("<GENERATOR-LOGRAD>Y\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>Wsrm2um\n")
                results.write(
                    "<GENERATOR-SUMMARY>Wavelengths range 0.70-5.00 um with a "
                    "resolution of 100 RP. Molecular radiative-transfer "
                    "enabled; Continuum flux module enabled;\n")
                results.write("<GENERATOR-TRANS>02-01\n")

            elif self.scope == "Hubble":
                results.write(
                    "<GENERATOR-INSTRUMENT>HST_WFC3-Grism: HST/WFC3 is the "
                    "wide field camera onboard the Hubble space telescope. "
                    "This versatile instrument enables slitless spectroscopy "
                    "at infrared wavelengths with two grisms (G102: 0.8 to "
                    "1.1 um and RP=210; G141: 1.1 to 1.7 um and RP=130). "
                    "Average throughputs (~0.4) and resolving powers (~160) "
                    "are considered for this configuration across the whole "
                    "infrared coverage of the grisms.\n")
                results.write("<GENERATOR-RANGE1>0.84\n")
                results.write("<GENERATOR-RANGE2>1.65\n")
                results.write("<GENERATOR-RANGEUNIT>um\n")
                results.write("<GENERATOR-RESOLUTION>170\n")
                results.write("<GENERATOR-RESOLUTIONUNIT>RP\n")
                results.write("<GENERATOR-TELESCOPE>SINGLE\n")
                results.write("<GENERATOR-DIAMTELE>2.4\n")
                results.write("<GENERATOR-BEAM>1.0\n")
                results.write("<GENERATOR-BEAM-UNIT>diffrac\n")
                results.write("<GENERATOR-TELESCOPE1>1\n")
                results.write("<GENERATOR-TELESCOPE2>0.0\n")
                results.write("<GENERATOR-TELESCOPE3>1.0\n")
                results.write("<GENERATOR-NOISE>CCD\n")
                results.write("<GENERATOR-NOISE1>20.0\n")
                results.write("<GENERATOR-NOISE2>0.005\n")
                results.write("<GENERATOR-NOISEOTEMP>288\n")
                results.write("<GENERATOR-NOISEOEFF>0.4\n")
                results.write("<GENERATOR-NOISEOEMIS>0.10\n")
                results.write("<GENERATOR-NOISETIME>{}\n".format(
                    self.exposure_time))
                results.write("<GENERATOR-NOISEFRAMES>{}\n".format(
                    self.exposure_count))
                results.write("<GENERATOR-NOISEPIXELS>8\n")
                results.write("<GENERATOR-TRANS-APPLY>N\n")
                results.write("<GENERATOR-TRANS-SHOW>N\n")
                results.write("<GENERATOR-LOGRAD>Y\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>Wsrm2um\n")
                results.write(
                    "<GENERATOR-SUMMARY>Wavelengths range 0.84-1.65 um with a "
                    "resolution of 170 RP. Molecular radiative-transfer "
                    "enabled; Continuum flux module enabled;\n")
                results.write("<GENERATOR-TRANS>02-01\n")

            elif self.scope == "Spitzer-Short-High":
                results.write(
                    "<GENERATOR-INSTRUMENT>Spizter Space Telescope IRS Short "
                    "High\n")
                results.write("<GENERATOR-RANGE1>9.9\n")
                results.write("<GENERATOR-RANGE2>19.6\n")
                results.write("<GENERATOR-RANGEUNIT>um\n")
                results.write("<GENERATOR-RESOLUTION>600\n")
                results.write("<GENERATOR-RESOLUTIONUNIT>RP\n")
                results.write("<GENERATOR-TELESCOPE>SINGLE\n")
                results.write("<GENERATOR-DIAMTELE>0.85\n")
                results.write("<GENERATOR-BEAM>1.0\n")
                results.write("<GENERATOR-BEAM-UNIT>diffrac\n")
                results.write("<GENERATOR-TELESCOPE1>1\n")
                results.write("<GENERATOR-NOISE>CCD\n")
                results.write("<GENERATOR-NOISE1>10.7\n")
                results.write("<GENERATOR-NOISE2>0.03\n")
                results.write("<GENERATOR-NOISEOTEMP>31\n")
                results.write("<GENERATOR-NOISEOEFF>0.4\n")
                results.write("<GENERATOR-NOISEOEMIS>0.10\n")
                results.write("<GENERATOR-NOISETIME>{}\n".format(
                    self.exposure_time))
                results.write("<GENERATOR-NOISEFRAMES>{}\n".format(
                    self.exposure_count))
                results.write("<GENERATOR-NOISEPIXELS>8\n")
                results.write("<GENERATOR-TRANS-APPLY>N\n")
                results.write("<GENERATOR-TRANS-SHOW>N\n")
                results.write("<GENERATOR-LOGRAD>Y\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>Wsrm2um\n")
                results.write(
                    "<GENERATOR-SUMMARY>Wavelengths range 9.9-19.6 um with a "
                    "resolution of 600 RP. Molecular radiative-transfer "
                    "enabled; Continuum flux module enabled;\n")
                results.write("<GENERATOR-TRANS>02-01\n")

            elif self.scope == "ALMA_Band7":
                results.write(
                    "< GENERATOR - INSTRUMENT > ALMA_Band7: ALMA(Atacama "
                    "Large Millimeter / submillimeter Array) is a mm / sub - "
                    "mm interferometer located on the Chajnantor plateau("
                    "Chile) at 5000 meters.High - resolution spectroscopy is "
                    "provided in the 31 - 950 GHz via a set of ten receivers "
                    "/ bands / frontends, combined with a versatile "
                    "correlator backend configuration.Highest spectral "
                    "resolution is 7.6 kHz, and highest spatial resolution is "
                    "10 milliarcseconds.Band-7 samples the 275-373 GHz range "
                    "employing SIS receivers.")
                results.write("< GENERATOR - RANGE1 > 345.7")
                results.write("< GENERATOR - RANGE2 > 345.9")
                results.write("< GENERATOR - RANGEUNIT > GHz")
                results.write("< GENERATOR - RESOLUTION > 500")
                results.write("< GENERATOR - RESOLUTIONUNIT > kHz")
                results.write("< GENERATOR - TELESCOPE > ARRAY")
                results.write("< GENERATOR - DIAMTELE > 12")
                results.write("< GENERATOR - BEAM > 0.2")
                results.write("< GENERATOR - BEAM - UNIT > arcsec")
                results.write("< GENERATOR - TELESCOPE1 > 54")
                results.write("< GENERATOR - TELESCOPE2 > 0.0")
                results.write("< GENERATOR - TELESCOPE3 > 1.0")
                results.write("< GENERATOR - NOISE > TRX")
                results.write("< GENERATOR - NOISE1 > 147")
                results.write("< GENERATOR - NOISE2 > 0")
                results.write("< GENERATOR - NOISEOTEMP > 270")
                results.write("< GENERATOR - NOISEOEFF > 0.85")
                results.write("< GENERATOR - NOISEOEMIS > 0.05")
                results.write("< GENERATOR - NOISETIME > 3600")
                results.write("< GENERATOR - NOISEFRAMES > 1")
                results.write("< GENERATOR - NOISEPIXELS > 8")
                results.write("< GENERATOR - TRANS - APPLY > N")
                results.write("< GENERATOR - TRANS - SHOW > Y")
                results.write("< GENERATOR - TRANS > 02 - 01")
                results.write("< GENERATOR - LOGRAD > N")
                results.write("< GENERATOR - GAS - MODEL > Y")
                results.write("< GENERATOR - CONT - MODEL > Y")
                results.write("< GENERATOR - CONT - STELLAR > N")
                results.write("< GENERATOR - RADUNITS > Jy")
                results.write(
                    "< GENERATOR - SUMMARY > Wavelengths range 345.7 - 345.9 "
                    "GHz with a resolution of 500 kHz.Molecular "
                    "radiative-transfer enabled; Continuum flux module "
                    "enabled; Interferometric observations;")

            else:
                print("Scope input {} not understood. Only JWST instruments, "
                      "Hubble, and Spitzer-Short-High supported here"
                      .format(self.scope))

            # ATMOSPHERE
            results.write("<ATMOSPHERE-DESCRIPTION>{} - CU Boulder GCM (Eric "
                          "Wolf)\n".format(self.planet_data["Name"]))
            results.write("<ATMOSPHERE-STRUCTURE>Equilibrium\n")
            results.write("<ATMOSPHERE-WEIGHT>{}\n".format(
                self.planet_data["MWeightTotal"]))
            results.write("<ATMOSPHERE-PRESSURE>{}\n".format(
                self.planet_data["SurfacePressure"]))
            results.write("<ATMOSPHERE-PUNIT>bar\n")
            # Needs to be modified if there are more gasses,added gasses
            # must be put in manually
            results.write("<ATMOSPHERE-NGAS>3\n")
            results.write("<ATMOSPHERE-GAS>N2,CO2,H2O\n")
            results.write("<ATMOSPHERE-TYPE>HIT[22],HIT[2],HIT[1]\n")
            results.write("<ATMOSPHERE-ABUN>1,1,1\n")
            results.write("<ATMOSPHERE-UNIT>scl,scl,scl\n")
            results.write("<ATMOSPHERE-NAERO>2\n")
            results.write("<ATMOSPHERE-AEROS>WaterIce,Cloud\n")
            results.write(
                "<ATMOSPHERE-ATYPE>CRISM_Wolff[reff=2.0um 0.31-99.57um],"
                "White_GSFC[reff=1.0um 0.10-1000000.00um]\n")
            results.write("<ATMOSPHERE-AABUN>{},{}\n".
                          format(self.planet_data["IceCloudAbundance"],
                                 self.planet_data["LiquidCloudAbundance"]))
            results.write("<ATMOSPHERE-AUNIT>scl,scl\n")
            results.write("<ATMOSPHERE-NMAX>0\n")
            results.write("<ATMOSPHERE-LMAX>0\n")
            results.write("<ATMOSPHERE-ASIZE>{},{}\n".format(
                self.planet_data["MeanIceCloudSize"],
                self.planet_data["MeanLiquidCloudSize"]/2))
            results.write(
                "<ATMOSPHERE-LAYERS-MOLECULES>Altitude,N2,CO2,H2O,WaterIce"
                ",Cloud\n")
            results.write("<ATMOSPHERE-LAYERS>{}\n".format(self.n_layers))
            for i, lvl in enumerate(self.atmosphere):
                results.write(
                    "<ATMOSPHERE-LAYER-{:.0f}>{:.3E},{:.3E},{:.3E},{:.3E},"
                    "{:.3E},{:.3E},{:.3E},{:.3E}\n".format(
                        lvl["Layer"], lvl["Pressure"], lvl["Temperature"],
                        lvl["Height"], lvl["Height"], lvl["CO2"], lvl["H2O"],
                        lvl["IceCloud"], lvl["LiquidCloud"]))

            print("    Successfully added {} layers to the atmosphere".format(
                self.n_uplayers))
            results.write("<SURFACE-TEMPERATURE>{}\n".format(
                self.planet_data["SurfaceTemperature"]))
            results.write("<SURFACE-ALBEDO>{}\n".format(
                self.planet_data["Albedo"]))
            results.write("<SURFACE-EMISSIVITY>{}\n".format(
                1 - self.planet_data["Albedo"]))
            print("    Successfully created PSG Config file from GCM results"
                  " for {}".format(self.planet_data["Name"]))
            print("    The file's name is {}".format(self._psginput_name))

    def send(self, keep_files=("trn", "lyr", "rad", "noi", "log", "atm", "err"),
             run=True):
        """This function send the the file file_name to the NASA GSFC PSG for 
        analysis It will not return anything, but will write files in the 
        current directory.

        Arguments: keep_files: a tuple of 3 letter strings which can include:
        trn, lyr, rad, noi, err, cfg, atm, log, str"""
        # Cumbersome giant file of everything
        alloutputname = self._psginput_name.split("psginput.txt")[
                            0] + "psgoutput_all.txt"
        filestem = self._psginput_name.split("psginput.txt")[0] + "psgoutput_"
        print("Sending to PSG")
        if run:
            command = "curl -d type=all --data-urlencode file@{} " \
                      "https://psg.gsfc.nasa.gov/api.php > {}" \
                .format(self._psginput_name, alloutputname)
            print("    {}".format(command))
            os.system(command)
        print("    Successfully connected to NASA PSG")
        with open(alloutputname, "r+") as allfile:
            sections = 0
            for i, line in enumerate(allfile):
                if "results_" in line:
                    sections += 1
            if sections == 0:  # Happens if the PSG is down or file is incorrect
                print("    No results returned from PSG")
                exit()
            allfile.seek(0)
            filetail = allfile.readline().split("_")[1].split("\n")[0]
            # print(filetail)
            # print(sections)
            self.returned_files = []
            for i in range(sections):
                outname = filestem + filetail
                self.returned_files.append(outname)
                with open(outname, "w") as wri:
                    line = allfile.readline()
                    count = 0
                    while "results_" not in line and count < 10000:
                        wri.write(line)
                        line = allfile.readline()
                        count += 1
                    if "results_" in line:
                        filetail = line.split("_")[1].split("\n")[0]
        print("    {} files created".format(sections))

        for i, fil in enumerate(self.returned_files):
            if fil.split("_")[-1].split(".")[0] not in keep_files:
                os.system("rm -v {}".format(fil))

    def plot_setup(self):
        rad_file = None
        for fil in self.returned_files:
            if "rad" in fil:
                rad_file = fil
        self._file_stem = rad_file.split("_psg")[0]
        name_parts = self._file_stem.split("_")
        self._title_stem = "{} {} {}".format(
            self.planet, name_parts[1], name_parts[2])

        radfil = np.loadtxt(rad_file, unpack=True, skiprows=16)
        self.Wavelengths = radfil[0]
        self.Total = radfil[1]
        self.Noise = radfil[2]
        self.Stellar = radfil[3]
        if len(radfil) == 5:
            self.is_transit = False
            print("    Exoplanet was not transitting")
            self.Thermal = radfil[4]
        else:
            self.is_transit = True
            self.Planet = radfil[4]
            self.Transit = radfil[5]
            self.Thermal = radfil[6]
        self._plot_range = (self.Wavelengths.min(), self.Wavelengths.max())

        trnfil = np.loadtxt(self._file_stem + "_psgoutput_trn.txt",
                            unpack=True, skiprows=18)
        self.tTotalSpec = trnfil[1]
        self.tN2Spec = trnfil[2]
        self.tCO2Spec = trnfil[3]
        self.tH2OSpec = trnfil[4]
        self.tIceSpec = trnfil[5]
        self.tCloudSpec = trnfil[6]
        try:
            self.tCIASpec = trnfil[7]
        except IndexError:
            self.tCIASpec = np.linspace(1, 1, len(self.tWavelengthSpec))

        noi_fil = np.loadtxt(self._file_stem + "_psgoutput_noi.txt",
                             unpack=True, skiprows=9)
        self.nTotal = noi_fil[1]
        self.nSource = noi_fil[2]
        self.nDetector = noi_fil[3]
        self.nTelescope = noi_fil[4]
        self.nBackground = noi_fil[5]

        print("Ready to Make Plots")

    def depth_plot(self):
        fig = plt.figure(figsize=(10, 6))
        ax = plt.gca()
        ax.step(self.Wavelengths, -self.Transit / self.Stellar * 1e+6,
                linewidth=0.25, c="b", label="Transit TransitDepth ",
                where="post")
        ax.set_title("Transit TransitDepth Without Noise\n{}".format(
            self._title_stem))
        ax.set_xlabel(r"Wavelengths ($\mu m$)")
        ax.set_ylabel("Signal (ppm)")
        ax.set_xlim(*self._plot_range)
        ax.xaxis.grid(True)
        fig.savefig("{}_dpth.png".format(self._file_stem))
        plt.cla()

    def depth_height(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        true_depth = -self.Transit / self.Stellar - self.planet_data["DepthEff"]
        den = (2 * self.planet_data["Radius"]
               * self.planet_data["ScaleHeight"])
        ax.step(self.Wavelengths, (true_depth * self.star_data["Radius"]**2
                                   / den),
                linewidth=0.25, c="b", label="Transit TransitDepth ",
                where="post")
        ax.set_title("Atmospheric Height Without Noise\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel(r"Scale Heights (km)")
        ax.set_xlim(*self._plot_range)
        ax.xaxis.grid(True)
        fig.savefig("{}_dpthh.png".format(self._file_stem))
        plt.cla()

    def emission(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.Thermal * 1.0e9, linewidth=0.5,
                c="black", where="post")
        ax.set_title("Planet emission\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel(r"Signal ($\frac{nW}{m^2 \mu m}$)")
        ax.set_xlim(*self._plot_range)
        ax.xaxis.grid(True)
        fig.savefig("{}_emis.png".format(self._file_stem))
        plt.cla()

    def absorption(self):
        print("    Planet blackbody peak wavelength: {:.2f}um".format(
            2897.8 / self.planet_data["SurfaceTemperature"]))
        planck = []
        distfactor = self.star_data["Radius"] ** 2 / (
                self.star_data["Distance"] * 3.086e16) ** 2
        for i, wave in enumerate(self.Wavelengths):
            wave = wave * 1.e-6
            valnum = (2 * 6.6261e-34 * 2.99792e8 ** 2)
            valden = wave ** 5 * (np.exp(
                6.6261e-34 * 2.99792e8
                / (wave * 1.38e-23 * self.star_data["Temperature"])) - 1)
            val = valnum / valden
            val *= distfactor
            planck.append(val)
        planck = np.array(planck)
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.Thermal * 1.0e9, linewidth=0.5,
                c="black", where="post")
        ax.step(self.Wavelengths, planck * 1.0e9, linewidth=0.5, c="blue",
                where="post")
        ax.set_title("Planet Absorption\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel(r"Signal ($\frac{nW}{m^2 \mu m}$)")
        ax.set_xlim(*self._plot_range)
        ax.xaxis.grid(True)
        fig.savefig("{}_abso.png".format(self._file_stem))
        plt.cla()

    def raw(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.Total + self.Noise,
                linewidth=0.5, c="black", where="post")
        ax.set_title("raw Signal\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel(r"Signal ($\frac{W}{m^2}$)")
        ax.set_xlim(*self._plot_range)
        ax.xaxis.grid(True)
        fig.savefig("{}_raw.png".format(self._file_stem))
        plt.cla()

    def star(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.Stellar, linewidth=0.5, c="red",
                where="post")
        ax.set_title("star\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel(r"Signal ($\frac{W}{m^{2}}$)")
        ax.set_xlim(*self._plot_range)
        ax.xaxis.grid(True)
        fig.savefig("{}_star.png".format(self._file_stem))
        plt.cla()

    def signal_and_noise(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths,
                -self.Transit / self.Stellar * 1.0e+6, linewidth=0.5,
                c="g", label="TransitDepth", where="post")
        ax.step(self.Wavelengths, self.Noise / self.Stellar * 1.0e+6,
                linewidth=0.5, c="r", label="Noise", where="post")
        ax.legend(loc=0)
        ax.set_title("Signal and Noise\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Signal (ppm)")
        ax.set_xlim(*self._plot_range)
        ax.xaxis.grid(True)
        fig.savefig("{}_sann.png".format(self._file_stem))
        plt.cla()

    def signal_noise_ratio(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, -self.Transit / self.Noise*100,
                linewidth=0.5, c="g", label="Signal/Noise", where="post")
        ax.axhline(1.0, label="1")
        ax.legend()
        ax.set_title("Signal to Noise\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Ratio (%)")
        ax.set_xlim(*self._plot_range)
        ax.xaxis.grid(True)
        fig.savefig("{}_snr.png".format(self._file_stem))
        plt.cla()

    def trn_total(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.tTotalSpec, linewidth=0.5, c="g",
                where="post")
        ax.set_title("Transmittance\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Transmittance (normalized)")
        ax.set_xlim(*self._plot_range)
        ax.set_ylim(0, 1.1)
        ax.xaxis.grid(True)
        fig.savefig("{}_tttl.png".format(self._file_stem))
        plt.cla()

    def trn_co2(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.tCO2Spec, linewidth=0.5, c="orange",
                where="post")
        ax.set_title("CO2 Transmittance\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Transmittance (normalized)")
        ax.set_xlim(*self._plot_range)
        ax.set_ylim(0, 1.1)
        ax.xaxis.grid(True)
        fig.savefig("{}_tCO2.png".format(self._file_stem))
        plt.cla()

    def trn_n2(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.tN2Spec, linewidth=0.5, c="purple",
                where="post")
        ax.set_title("N2 Transmittance\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Transmittance (normalized)")
        ax.set_xlim(*self._plot_range)
        ax.set_ylim(0, 1.1)
        ax.xaxis.grid(True)
        fig.savefig("{}_tN2.png".format(self._file_stem))
        plt.cla()

    def trn_h2o(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.tH2OSpec, linewidth=0.5, c="b",
                where="post")
        ax.set_title("H2O Transmittance\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Transmittance (normalized)")
        ax.set_xlim(*self._plot_range)
        ax.set_ylim(0, 1.1)
        ax.xaxis.grid(True)
        fig.savefig("{}_tH2O.png".format(self._file_stem))
        plt.cla()

    def trn_ice(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.tIceSpec, linewidth=0.5, c="gray",
                where="post")
        ax.set_title("Ice Transmittance\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Transmittance (normalized)")
        ax.set_xlim(*self._plot_range)
        ax.set_ylim(0, 1.1)
        ax.xaxis.grid(True)
        fig.savefig("{}_tice.png".format(self._file_stem))
        plt.cla()

    def trn_cloud(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.tCloudSpec, linewidth=0.5, c="black",
                where="post")
        ax.set_title("Cloud Transmittance\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Transmittance (normalized)")
        ax.set_xlim(*self._plot_range)
        ax.set_ylim(0, 1.1)
        ax.xaxis.grid(True)
        fig.savefig("{}_tcld.png".format(self._file_stem))
        plt.cla()

    def trn_cia(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.tCIASpec, linewidth=0.5, c="red",
                where="post")
        ax.set_title("tCIASpec Transmittance\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Transmittance (normalized)")
        ax.set_xlim(*self._plot_range)
        ax.set_ylim(0, 1.1)
        ax.xaxis.grid(True)
        fig.savefig("{}_tcia.png".format(self._file_stem))
        plt.cla()

    def trn_species(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.tTotalSpec, linewidth=0.5, c="g",
                where="post")
        ax.step(self.Wavelengths, self.tCO2Spec, linewidth=0.5, c="orange",
                alpha=0.4, where="post")
        ax.step(self.Wavelengths, self.tN2Spec, linewidth=0.5, c="purple",
                where="post")
        ax.step(self.Wavelengths, self.tH2OSpec, linewidth=0.5, c="b",
                alpha=0.2, where="post")
        plt.legend(["Total", "CO2", "N2", "H2O"], loc=4)
        ax.set_title("Species Transmittance\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Transmittance (normalized)")
        ax.set_xlim(*self._plot_range)
        ax.set_ylim(0, 1.1)
        ax.xaxis.grid(True)
        fig.savefig("{}_tspc.png".format(self._file_stem))
        plt.cla()

    def trn_aero(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.tTotalSpec, linewidth=0.5, c="g",
                where="post")
        ax.step(self.Wavelengths, self.tIceSpec, linewidth=0.5, c="gray",
                where="post")
        ax.step(self.Wavelengths, self.tCloudSpec, linewidth=0.5, c="black",
                where="post")
        ax.step(self.Wavelengths, self.tCIASpec, linewidth=0.5, c="red",
                where="post")
        plt.legend(["Total", "Ice", "Cloud", "tCIASpec"], loc=4)
        ax.set_title("tCIASpec Transmittance\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Transmittance (normalized)")
        ax.set_xlim(*self._plot_range)
        ax.set_ylim(0, 1.1)
        ax.xaxis.grid(True)
        fig.savefig("{}_taer.png".format(self._file_stem))
        plt.cla()

    def trn_all(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.tTotalSpec, linewidth=0.5, c="g")
        ax.step(self.Wavelengths, self.tCO2Spec, linewidth=0.5, c="orange",
                alpha=0.4, where="post")
        ax.step(self.Wavelengths, self.tN2Spec, linewidth=0.5, c="purple",
                where="post")
        ax.step(self.Wavelengths, self.tH2OSpec, linewidth=0.5, c="b",
                alpha=0.2, where="post")
        ax.step(self.Wavelengths, self.tIceSpec, linewidth=0.5, c="gray",
                where="post")
        ax.step(self.Wavelengths, self.tCloudSpec, linewidth=0.5, c="black",
                where="post")
        ax.step(self.Wavelengths, self.tCIASpec, linewidth=0.5, c="red",
                where="post")
        ax.legend(["Total", "CO2", "N2", "H2O", "Ice", "Cloud", "CIA"], loc=4)
        ax.set_title("Transmittance\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Transmittance (normalized)")
        ax.set_xlim(*self._plot_range)
        ax.set_ylim(0, 1.1)
        ax.xaxis.grid(True)
        fig.savefig("{}_tall.png".format(self._file_stem))
        plt.cla()

    def noise_components(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.nTotal, linewidth=0.5,
                c="black", label="Total", where="post")
        ax.step(self.Wavelengths, self.nSource, linewidth=0.5, c="b",
                label="nSource", where="post")
        ax.step(self.Wavelengths, self.nDetector, linewidth=0.5, c="g",
                label="nDetector", where="post")
        ax.step(self.Wavelengths, self.nTelescope, linewidth=0.5,
                c="gold", label="nTelescope", where="post")
        ax.step(self.Wavelengths, self.nBackground, linewidth=0.5,
                c="r", label="nBackground", where="post")
        ax.step(self.Wavelengths, self.Transit, linewidth=0.5,
                c="black", label="Signal", where="post")
        ax.legend(markerscale=5, loc=4)
        ax.set_title("Noise\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel(r"Noise ($\log{\frac{W}{m^{2}}}$)")
        ax.set_xlim(*self._plot_range)
        ax.set_yscale("log")
        ax.xaxis.grid(True)
        fig.savefig("{}_nois.png".format(self._file_stem))
        plt.cla()

    def noise_ppm(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, self.nTotal / self.Stellar * 1.0e6,
                linewidth=0.5, c="black", label="Total", where="post")
        ax.legend(markerscale=5, loc=4)
        ax.set_title("Noise\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel(r"Noise ($\log{\frac{W}{m^{2}}}$)")
        ax.set_xlim(*self._plot_range)
        ax.set_yscale("log")
        ax.xaxis.grid(True)
        fig.savefig("{}_noisppm.png".format(self._file_stem))
        plt.cla()


class PSGCompared(object):
    """This function plots all _rad PSG outputs together to compare their
    peaks and features. The most useful part is the depth combined plots,
    which can give you a sense of signal strength."""

    def __init__(self, radfiles):
        super(PSGCompared, self).__init__()
        self.radfiles = glob.glob(radfiles)
        self.colors = ["r", "orange", "green", "yellowgreen", "orange",
                       "yellowgreen", "g", "b", "indigo", "m", "r",
                       "cyan", "blue", "violet"]
        self.alphas = [0.7, 0.9, 0.7, 0.9, 0.9, 0.9, 0.7, 0.5, 0.5, 0.7, 0.6,
                       0.6, 0.4, 0.6]
        self.files = glob.glob(radfiles)
        self.outputfile = self.files[0].split("_")[0]

        self.WavelengthSpec = None

    def N2CombPlot(self):
        fig = plt.figure(figsize=(12, 15))
        plt.tick_params(axis="both", which="both", bottom="on", top="off",
                        labelbottom="on", left="on", right="off",
                        labelleft="on", labelsize=16, length=8, width=2)
        for i, file in enumerate(self.radfiles):
            if "barN2" in file:
                if "barCO2" in self.files[i]:
                    C = float(self.files[i].split("_")[2].split("barCO2")[0])
                else:
                    C = 0
                contents = np.loadtxt(self.files[i], unpack=True, skiprows=16)
                self.WavelengthSpec = contents[0]
                Transit = contents[5]
                Stellar = contents[3]
                sig = -Transit / Stellar * 1.0e6
                ax.step(self.WavelengthSpec, sig, c=self.colors[i],
                         linewidth=0.5, label="{0:6.4f} Bar CO2".format(C),
                         where="post")
                plt.axhline(np.mean(sig[2630:2634]), c=self.colors[i])
                plt.annotate(
                    "{0:6.4f} Bar: {1:4.0f}".format(C, np.mean(sig[2630:2634])),
                    (
                        np.max(self.WavelengthSpec) - 3.8,
                        np.mean(sig[2630:2634]) - 6),
                    fontsize=16)
        ax.set_title(r"\huge Transit Depths for 1 bar N$_{2}$ Planet with Line "
                  r"at CO$_{2}$ 15\si{\micro \meter} Peak")
        plt.xlabel(r"Wavelengths ($\si{\micro \meter}$)", fontsize=18)
        ax.set_xlim(np.min(self.WavelengthSpec) - 2,
                 np.max(self.WavelengthSpec) + 2)
        plt.set_ylabel("{Signal (ppm)", fontsize=20)
        leg = plt.legend(loc=2, fontsize=15)
        for label in leg.get_lines():
            label.set_linewidth(4)
            ax = plt.axes()
            ax.xaxis.grid(True)
        fig.savefig(self.outputfile + "N2CombinedPlots.png", dpi=300)
        plt.cla()

    def NoN2CombPlot(self):
        fig = plt.figure(figsize=(12, 15))
        plt.tick_params(axis="both", which="both", bottom="on", top="off",
                        labelbottom="on", left="on", right="off",
                        labelleft="on", labelsize=16, length=8, width=2)
        for i, file in enumerate(self.radfiles):
            if "barN2" not in file:
                if "barCO2" in file:
                    C = float(file.split("barCO2")[0].split("_")[-1])
                    contents = np.loadtxt(file, unpack=True, skiprows=16)
                    self.WavelengthSpec = contents[0]
                    Transit = contents[5]
                    Stellar = contents[3]
                    sig = -Transit / Stellar * 1.0e6
                    ax.step(self.WavelengthSpec, sig, c=self.colors[i],
                             linewidth=0.5, label="{0:6.4f} Bar CO2".format(C),
                             where="post")
                    plt.axhline(np.mean(sig[2630:2634]), c=self.colors[i])
                    plt.annotate("{0:6.4f} Bar: {1:4.0f}".format(C, np.mean(
                        sig[2630:2634])),
                                 (np.max(self.WavelengthSpec) - 3.8,
                                  np.mean(sig[2630:2634]) - 8), fontsize=16)
        ax.set_title(r"\huge Transit Depths for 0 bar N$_{2}$ Planet with Line "
                  "at CO$_{2}$ 15\si{\micro \meter} Peak")
        plt.xlabel(r"Wavelengths ($\si{\micro \meter}$)", fontsize=18)
        ax.set_xlim(np.min(self.WavelengthSpec) - 2,
                 np.max(self.WavelengthSpec) + 2)
        plt.set_ylabel("{Signal (ppm)", fontsize=20)
        leg = plt.legend(loc=2, fontsize=15)
        for label in leg.get_lines():
            label.set_linewidth(4)
            ax = plt.axes()
            ax.xaxis.grid(True)
        fig.savefig(self.outputfile + "NoN2CombinedPlots.png", dpi=300)
        plt.cla()

    def PeakCompare(self):
        fig = plt.figure(figsize=(8, 6))
        N2peaks = []
        Nopeaks = []
        Ns = []
        Nos = []
        for i in range(len(self.files)):
            contents = np.loadtxt(self.files[i], unpack=True, skiprows=16)
            self.WavelengthSpec = contents[0]
            Transit = contents[5]
            Stellar = contents[3]
            Signal = -Transit / Stellar * 1.0e6
            if "barN2" in self.files[i]:
                N = float(self.files[i].split("_")[1].split("barN2")[0])
                if "barCO2" in self.files[i]:
                    C = float(self.files[i].split("_")[2].split("barCO2")[0])
                else:
                    C = 0
                    Ns.append(C)
                    N2peaks.append(np.max(Signal))
                    if "barN2" not in self.files[i]:
                        C = float(
                            self.files[i].split("barCO2")[0].split("_")[-1])
                        Nos.append(C)
                        Nopeaks.append(np.max(Signal))
        plt.scatter(Ns, N2peaks, c="b", s=100, label="With Nitrogen")
        plt.scatter(Nos, Nopeaks, c="r", s=100, label="Without Nitrogen")
        ax.set_title("Peak Transit TransitDepth")
        plt.xlabel(r"CO$_{2}$ Abundance ($\si{\bar}$)")
        plt.set_ylabel("Peak Signal (ppm)")
        leg = plt.legend(loc=4)
        """for label in leg.get_lines():
            label.set_linewidth(4)"""
        ax = plt.axes()
        ax.xaxis.grid(True)
        fig.savefig(self.outputfile + "PeakCompare.png", dpi=300)
        plt.cla()


def PSGCompare(filename):
    x = PSGCompared(filename)
    x.N2CombPlot()
    x.NoN2CombPlot()
    x.PeakCompare()
