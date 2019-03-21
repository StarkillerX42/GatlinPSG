import os
import sys
import time
from astropy.table import Table
import astropy.units as u
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import warnings
from pathlib import Path
# import starcoder42 as s
warnings.filterwarnings("ignore")
import pandexo.engine.justdoit as jdi
assert sys.version_info.major >= 3, "Only runs on Python 3.6 or higher"
assert sys.version_info.minor >= 6, "Only runs on Python 3.6 or higher"

__author__ = "Dylan_Gatlin"
__version__ = 3.7
psg_scopes = ["MIRI-LRS", "MIRI-MRS", "NIRCam-Grism", "NIRISS-SOSS",
              "NIRSpec-1000", "NIRSpec-2700", "NIRSpec-Prism", "Hubble",
              "Spitzer-Short-High", "ALMA_Band_7"]
pandexo_scopes = ['WFC3 G141', 'MIRI LRS', 'NIRISS SOSS', 'NIRSpec G140M',
                  'NIRSpec G140H', 'NIRSpec G235M', 'NIRSpec G235H',
                  'NIRSpec G395M', 'NIRSpec G395H', 'NIRSpec Prism',
                  'NIRCam F322W2', 'NIRCam F444W']
rad_unit_types = ["rel", "rkm", "Wsrm2um", "raw", "Jy"]
file_types = ["trn", "lyr", "rad", "noi", "log", "atm", "err"]


# noinspection PyShadowingNames
class PSG(object):
    """This class is meant to help interface between NASA's PSG and climate
    models. The goal is to create a versatile PSG object which can be
    interacted with via jupyter notebooks. The proper way to run the PSG is:

    planet = PSG.PSG(planet_name)
    planet.fetch_archive(is_earth)
    planet.from_cdf(cdf_file, phase)
    planet.calculate(atmosphere_ceiling, n_uplayers)
    planet.write(scope, exp_time, n_exposures, rad_units)
    planet.send(keep_files, run)
    planet.plot_setup()
    planet.depth_plot()


    Required Files:

        In the folder you"re running this, one file must be accessible.
        cdf_file, a text file following Eric Wolf's formatting
        cdf_file must be the terminator profile ending in
        _terminator.txt or another similar formatted file. If the observation
        phase is different than 180, the corresponding atmosphere profile
        must be given.

    PSG Init Arguments:

        planet_name: (str) The name of the planet according to NASA"s API.
        For a list of known planets, run this function with a dummy planet name
        and read through exoplanets.csv, which should appear in your directory.

        cdf_file: (str) The name of the file which you would like to
        read in, described above.

        is_earth: (bool)Whether or not this planet is a fake exoplanet with
        the same mass and radius as Earth If is_earth is True, the planet
        will be imagined as if it was around the star of the variable "planet"

        n_uplayers: (int) The number of isothermal layers between the top given
        layer and the top of the atmosphere

        phase: The orbital phase, 180 is a transit, 0 is an occultation

    PSG.calculate Arguments:

        skprow: The argument skiprows from np.loadtxt, telling you how many
        rows to skip from the atmopshere profile file.

        atmosphere_ceiling: (float)The pressure where the atmosphere ends
        The PSG will only produce useful results if there are layers at
        extremely low pressure_profile. If the atmosphere profile isn't high
        enough, imaginary upper layers can be added. Default is 0., but 1e-6
        is recommended

    PSG.write Arguments:

        scope: The scope used in the observation. For a list of possible inputs
        look at PSG.psg_scopes

        exposure time: The length of each exposure

        exposure count: The number of exposures

        rad_units: The type of return from the PSG. For a list of possible
        values, check PSG.rad_unit_types

    PSG.send Arguments

        run: A bool of whether or not to actually send it to PSG. If you ran it
        recently, skipping this step will save time.

        keep_files: A tuple of file types you'd like returned. For a list of
        file types, check PSG.file_types

    PSG.plot_setup Arguments

        None

    PSG.pandexo Arguments

        scope: Different from the scope in PSG.write, this refers to the same
        thing, but it must be formatted according to pandexo. See
        PSG.pandexo_scopes for a list of options.
    """

    def __init__(self, planet_name: str):
        super(PSG, self).__init__()

        self.planet_name = planet_name
        # Planet Variables
        self.is_earth = None
        self.planet_data = {}
        self.star_data = {}
        self.atmosphere = None
        self.n_downlayers = None
        self.n_layers = None
        self.phase = None
        self.cdf_file = None
        self.profile_file = None
        self.cdf_atmosphere = None
        self.netcdf = None
        self.lat_weight = None
        self.lon_weight = None
        self.phase_mask = None
        self.mask = None

        # PSG Variables
        self.scope = None
        self.exposure_time = None
        self.exposure_count = None
        self.rad_units = None
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
        self.nReal = None

        # Pandexo variables
        self.pand_fil = None
        self.pandexo_result = None

        print(f"Starting PSG for {planet_name}")
        self._file_dir = Path("./psg_files")
        if not self._file_dir.is_dir():
            self._file_dir.mkdir()

    def fetch_archive(self, is_earth: bool = False):
        """This function fetches an archive from the NASA exoplanet archive
        which is then used in later steps. It only works in the planet name
        provided when the class is initialized matches an item in the archive.

        Arguments:
            is_earth: (bool) Whether or not this model used is Earth-like, this
            will replace the values found in the exoplanet archive like mass and
             radius with Earth values.
        """
        self.is_earth = is_earth
        exoplanets = self._file_dir.joinpath("exoplanets.csv")
        if not exoplanets.is_file():
            need_file = True
        else:
            st = os.stat(exoplanets)
            age = time.time()-st.st_mtime
            if age > (3600 * 24 * 2.0):
                need_file = True
            else:
                need_file = False
        if need_file:
            print("    Retrieving planet variables from NASA's Exoplanet "
                  "Archive")
            import requests
            r = requests.get("https://exoplanetarchive.ipac.caltech.edu/cgi"
                             "-bin/nstedAPI/nph-nstedAPI?table=exoplanets"
                             "&select=pl_name,pl_masse,pl_rade,pl_orbsmax,"
                             "pl_orbincl,pl_trandep,st_teff,st_rad,st_radv,"
                             "st_dist,st_optmag,pl_insol,st_metfe,st_logg"
                             "&format=csv")
            lines = r.text[2:].splitlines()
            with open(exoplanets.absolute(), "w") as fil:
                for line in lines:
                    fil.write(line + "\n")
        # Extracts necessary details about the exoplanet from NASA"s API
        # Defines
        (planet_name, planet_emass, planet_erad,
         sma, inclination, transit_depth,
         star_temp, star_srad, star_velocity,
         star_distance, star_magnitude,
         insolation, metallicity, star_gravity) = [None] * 14
        exoplanets = np.loadtxt(exoplanets.absolute(), delimiter=",",
                                skiprows=1, dtype=np.str, comments="'")
        if sum(exoplanets[:, 0] == self.planet_name) == 1:
            line = exoplanets[exoplanets[:, 0] == self.planet_name][0]
            (planet_name, planet_emass, planet_erad,
             sma, inclination, transit_depth,
             star_temp, star_srad, star_velocity,
             star_distance, star_magnitude,
             insolation, metallicity, star_gravity) = line
        else:
            print("    Planet not found. Inupts must be given manually.")
        if star_velocity == "":
            star_velocity = 0.
        if star_magnitude == "nan":
            star_magnitude = 10.
        if star_gravity == "":
            star_gravity = 4.
        if self.is_earth:
            planet_name = "ExoEarth like " + self.planet_name
            planet_emass = 1.
            planet_erad = 1.
            star_srad = 1.
            sma = 1.
            inclination = 90.
            insolation = 1361.
            star_velocity = 0.
            star_temp = 5700
            star_distance = 10.
            star_magnitude = 10.
            metallicity = 0.0122
            star_gravity = 4.4
            transit_depth = str(float(planet_erad)
                                / float(star_srad) * 0.009154)

        # Converts NASA"s values to proper units
        self.planet_data["Name"] = planet_name
        self.planet_data["Mass"] = float(planet_emass) * 5.9736e24
        self.planet_data["ERadius"] = float(planet_erad)
        self.planet_data["Radius"] = float(planet_erad) * 6371.0
        self.planet_data["Diameter"] = self.planet_data["Radius"] * 2
        self.planet_data["SemiMajorAxis"] = float(sma)
        self.planet_data["Inclination"] = float(inclination)
        self.planet_data["TransitDepth"] = float(transit_depth) / 100
        self.planet_data["Insolation"] = float(insolation) * 1361.
        self.star_data["SRadius"] = float(star_srad)
        self.star_data["Radius"] = self.star_data["SRadius"] * 1.39e6
        self.star_data["Velocity"] = float(star_velocity)
        self.star_data["Temperature"] = float(star_temp)
        self.star_data["Distance"] = float(star_distance)
        self.star_data["Magnitude"] = float(star_magnitude)
        self.star_data["Metallicity"] = np.log10(float(metallicity))
        self.star_data["Gravity"] = float(star_gravity)

        self.planet_data["DepthEff"] = (self.planet_data["Radius"]
                                        / self.star_data["Radius"]) ** 2
        self.planet_data["Gravity"] = (6.67e-11
                                       * self.planet_data["Mass"]
                                       / (self.planet_data["Radius"]
                                          * 1000) ** 2)
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
        print("    Exoplanet Archive fetched. planet_data and star_data filled")

    def from_cdf(self, cdf_file, phase=180, is_transit=True):
        """This function imports an atmosphere profile from a netcdf file. If
        a netcdf file isn't available, simply specify the PSG attribute
        self.cdf_file before running calculate.

        Args:
            cdf_file: (path-like) The location of a netcdf file to be used
            phase: (number) The orbital phase where 180 is a transit and 0 is
                an occultation. In degrees.
            is_transit: (bool) Whether or not you want the PSG pipeline to
            consider your input as a transit. This may be useful to promote
            continuity for thermal phase curves when phase is close to 180"""

        self.phase = phase
        self.cdf_file = Path(cdf_file)
        self.is_transit = is_transit
        # print(self.cdf_file.absolute())
        self.netcdf = Dataset(cdf_file)
        # print(self.netcdf)
        print("    Accessing netCDF contents")

        # noinspection PyShadowingNames
        def hybrid2pressure(cdf):
            """This function takes a netcdf4 object as input and returns the
            p_pro and ilev_pressures"""
            hyam = cdf["hyam"][:]
            hybm = cdf["hybm"][:]
            p_s = cdf["PS"][0]
            p0 = cdf["P0"][0]
            ip = np.zeros(
                (len(cdf["lat"]), len(cdf["lon"]), len(cdf["lev"])))
            # noinspection PyShadowingNames
            for i, lat in enumerate(cdf["lat"][:]):
                for j, lon in enumerate(cdf["lon"][:]):
                    ip[i, j, :] = hyam * p0 + hybm * p_s[i, j]
            ip = ip[:]
            for i, lon in enumerate(cdf["lat"][:]):
                for j, lat in enumerate(cdf["lon"][:]):
                    ip[i, j, :] = hyam * p0 + hybm * p_s[i, j]
            return ip, ip

        def hybrid2height(cdf, g, molar_mass):
            """Converts netcdf hybrid coefficients to heights in m"""
            ps = cdf["PS"][0]
            p0 = cdf["P0"][0]
            hyai = cdf["hyai"][:]
            hybi = cdf["hybi"][:]
            t = cdf["T"][0]
            hyam = cdf["hyam"][:]
            hybm = cdf["hybm"][:]
            r = 8.3144 / molar_mass
            z = np.zeros((len(cdf["lat"]), len(cdf["lon"]), len(cdf["lev"])+1))
            iz = np.zeros((len(cdf["lat"]), len(cdf["lon"]), len(cdf["lev"])+1))
            for i, lat in enumerate(cdf["lat"][:]):
                for j, lon in enumerate(cdf["lon"]):
                    for k in range(len(cdf["lev"]), 0, -1):
                        # print(k)
                        p1 = hyai[k] * p0 + hybi[k] * ps[i, j]
                        p2 = hyai[k-1] * p0 + hybi[k-1] * ps[i, j]
                        delta_z = r * t[k-1, i, j] / g * np.log(p1 / p2)
                        iz[i, j, k-1] = iz[i, j, k] + delta_z
                        p_prime = hyam[k-1] * p0 + hybm[k-1] * ps[i, j]
                        z_scale = -r * t[k-1, i, j] / g * np.log(p_prime / p1)
                        z[i, j, k-1] = iz[i, j, k] + z_scale
            return z, iz

        earth_facing_long = (180-phase) % 360  # Longitude is opposite of PSG's
        # phase conventions, earth_facing_long refers to the longitude of the
        # cdf file, phase refers to the orbital phase in PSG. Transit is for
        # phase of 180, but it grabs from the terminator at 90 and 270. For
        # non-transits, the whole earth-facing side is averaged.
        if self.is_transit:
            if (phase < 175) or (phase > 185):
                print("    is_transit is True, but was likely not transitting "
                      "because the phase is not close to 180")
            # Use % to account for sub-zero angs
            left_ang = (earth_facing_long - 90.) % 360.
            right_ang = (earth_facing_long + 90.) % 360.
            left_bool = ((left_ang-5. <= self.netcdf["lon"][:])
                         & (self.netcdf["lon"][:] <= left_ang + 5.))
            right_bool = ((right_ang-5. <= self.netcdf["lon"][:])
                          & (self.netcdf["lon"][:] <= right_ang + 5.))
            phase_mask = left_bool ^ right_bool
            lon_weight = np.ones(len(self.netcdf["lon"]))  # If it's a transit,
            # we need a different longitude weight than a disk average
            lat_weight = np.ones(len(self.netcdf["lat"][:]))
            if phase != 180:
                print("    Phase is not 180, but is_transit is True")
        else:  # Else, the Earth-facing side is selected
            left_ang = (earth_facing_long-90.) % 360.
            right_ang = (earth_facing_long+90.) % 360.
            # print(left_ang, right_ang)
            if left_ang < right_ang:
                phase_mask = ((left_ang <= self.netcdf["lon"][:])
                              & (self.netcdf["lon"][:] <= right_ang))
            else:
                phase_mask = ((left_ang <= self.netcdf["lon"][:])
                              ^ (right_ang >= self.netcdf["lon"][:]))
            lon_weight = np.cos(np.deg2rad(self.netcdf["lon"][:] + phase + 180))
            lat_weight = np.round(
                np.cos(np.deg2rad(self.netcdf["lat"][:])) ** 2, 5)
        weight_grid = np.outer(lat_weight, lon_weight)
        self.lon_weight = lon_weight
        self.lat_weight = lat_weight
        self.phase_mask = phase_mask
        self.mask = (weight_grid > 0) & phase_mask
        t_pro = np.average(
                    np.average(
                        self.netcdf["T"][0, :, :, phase_mask],
                        axis=2, weights=lon_weight[phase_mask]),
                    axis=1, weights=lat_weight)
        q_pro = np.average(
            np.average(
                self.netcdf["Q"][0, :, :, phase_mask],
                axis=2, weights=lon_weight[phase_mask]),
            axis=1, weights=lat_weight)
        ch4_vmr = self.netcdf["ch4vmr"][0]
        co2_vmr = self.netcdf["co2vmr"][0]
        n2_vmr = 1-co2_vmr-ch4_vmr
        # All need to be averaged by cos^2 in lat and cos in lon, but np.average
        # Can't handle an average over 2 of the 3 given dimensions
        cldliq = np.average(
            np.average(
                self.netcdf["CLDLIQ"][0, :, :, phase_mask],
                axis=2, weights=lon_weight[phase_mask]),
            axis=1, weights=lat_weight)
        cldice = np.average(
            np.average(
                self.netcdf["CLDICE"][0, :, :, phase_mask],
                axis=2, weights=lon_weight[phase_mask]),
            axis=1, weights=lat_weight)
        reffliq = np.average(
            np.average(
                self.netcdf["REL"][0, :, :, phase_mask],
                axis=2, weights=lon_weight[phase_mask]),
            axis=1, weights=lat_weight)
        reffice = np.average(
            np.average(
                self.netcdf["REI"][0, :, :, phase_mask],
                axis=2, weights=lon_weight[phase_mask]),
            axis=1, weights=lat_weight)
        p_surf = np.average(
            np.average(
                self.netcdf["PS"][0, :, phase_mask],
                axis=1, weights=lon_weight[phase_mask]),
            axis=0, weights=lat_weight) / 100 / 1000
        t_surf = np.average(
            np.average(
                self.netcdf["TS"][0, :, phase_mask],
                axis=1, weights=lon_weight[phase_mask]),
            axis=0, weights=lat_weight)
        albedos = self.netcdf["FUS"][0, -1] / self.netcdf["FDS"][0, -1]
        albedos[~self.mask] = np.nan
        albedo = np.nansum(albedos * (weight_grid/np.nansum(
            weight_grid[self.mask])))
        # print(albedo)

        amounts = [n2_vmr, co2_vmr, ch4_vmr]
        m_weight_dry = np.average([0.0280, 0.0440, 0.01604],
                                  weights=amounts)
        heights, iheights = hybrid2height(self.netcdf,
                                          self.planet_data["Gravity"],
                                          m_weight_dry)
        # print(heights, iheights)

        pressures, ipressures = hybrid2pressure(self.netcdf)
        p_pro = np.average(np.average(pressures[:, phase_mask], axis=0,
                                      weights=lat_weight), axis=0) / 100 / 1000
        h_pro = np.average(np.average(heights[:, phase_mask, :-1],
                                      axis=0, weights=lat_weight), axis=0)
        q_mmr = q_pro / (1.-q_pro)
        q_vmr_dry = q_mmr * 28.0059 / 18.01528
        q_vmr_wet = q_vmr_dry / (1 + q_vmr_dry)
        co2_vmr_wet = co2_vmr / (1 + q_vmr_wet)
        ch4_vmr_wet = ch4_vmr / (1 + q_vmr_wet)
        n2_vmr_wet = 1.0-co2_vmr_wet-ch4_vmr_wet-q_vmr_wet
        self.cdf_atmosphere = Table([h_pro, p_pro, t_pro, n2_vmr_wet,
                                     co2_vmr_wet, ch4_vmr_wet, q_vmr_wet,
                                     cldliq, cldice, reffliq, reffice],
                                    names=["Heights", "Pressures", "Temps",
                                           "N2", "CO2", "CH4", "H2O",
                                           "LiquidCloud", "IceCloud",
                                           "LiquidCloudSize", "IceCloudSize"])
        stem_str = str(self._file_dir.absolute().joinpath(
            self.cdf_file.name)).split("_aqua")[0]
        if self.is_transit:
            self.profile_file = Path(stem_str + "_transit.txt").absolute()
        else:
            self.profile_file = Path(stem_str + "_{:.0f}_pro.txt".format(phase))
        # print(self.profile_file)
        with open(self.profile_file, "w") as out:
            out.write("# Simulation:  {}\n".format(
                str(self.cdf_file).split(".cam")[0]))
            out.write("# {}, {:.2f}% N2, {:.2f}% CO2, {:.2f}% CH4\n".format(
                self.planet_name, n2_vmr, co2_vmr, ch4_vmr))
            out.write("# Mean Vertical Profile\n")
            out.write("# Mass mixing ratios (kg/kg) relative to "
                      "the moist air mass (i.e dry air mass + water vapor "
                      "mass)\n")
            out.write("# Pressure are total pressure (i.e P_dry + "
                      "P_h2o)\n")
            out.write("# Surf Press (mb)    Surf Temp (K)        Surf Albedo   "
                      "       Dry Molar Weight     Phase\n")
            out.write("#{:<20.7f}#{:<20.7f}#{:<20.7f}#{:<20.7f}#{:<20.7f}\n".
                      format(p_surf, t_surf, albedo, m_weight_dry*1000, phase))
            out.write("{:4.20s}{:>20.20s}{:>20.20s}{:>20.20s}{:>20.20s}"
                      "{:>20.20s}{:>20.20s}{:>20.20s}{:>20.20s}{:>20.20s}"
                      "{:>20.20s}{:>20.20s}\n".format("# Lev", "Height (m)",
                                                      "Pressure (bar)",
                                                      "Temperature (K)",
                                                      "N2 (vmr)", "CO2 (vmr)",
                                                      "CH4 (vmr)", "H2O (vmr)",
                                                      "Liquid Clouds",
                                                      "Ice Clouds",
                                                      "Liquid Size",
                                                      "Ice Size"))
            for i, lvl in enumerate(self.cdf_atmosphere):
                out.write(
                    "{:5.0f}{:20.7f}{:20.7f}{:20.7f}{:20.7f}{:20.7E}{:20.7E}"
                    "{:20.7E}{:20.7E}{:20.7E}{:20.3f}{:20.3f}\n".format(
                        i, lvl["Heights"], lvl["Pressures"],
                        lvl["Temps"], lvl["N2"], lvl["CO2"], lvl["CH4"],
                        lvl["H2O"], lvl["LiquidCloud"], lvl["IceCloud"],
                        lvl["LiquidCloudSize"], lvl["IceCloudSize"]))
        print("    Output file written to "
              f"{self.profile_file.relative_to(Path('.').absolute())}")

    def calculate(self, atmosphere_ceiling=0, n_uplayers: int = 0):

        """
        1. Reads model file
        2. Converts to PSG units
        3. Adds layers to TOA

        Args:
            atmosphere_ceiling: (float) The pressure in bars of the top
                of the atmosphere. Recommend 1e-6
            n_uplayers: (int) The number of layers to add to the atmosphere.
                Recommend 7
        """

        # GCM Inputs
        # TODO Integrate pathlib support beyond here
        with open(self.profile_file, "r") as fp:
            lines = fp.readlines()
            line = lines[6]
            c = line.split("#")
            self.planet_data["SurfacePressure"] = float(c[1])
            self.planet_data["SurfaceTemperature"] = float(c[2])
            self.planet_data["Albedo"] = float(c[3])
            self.planet_data["MWeightDry"] = float(c[4])
            self.planet_data["Phase"] = float(c[5])

        profile = np.loadtxt(self.profile_file, skiprows=8, unpack=True)
        layers_profile = profile[0] + 1
        self.n_downlayers = len(layers_profile)
        heights_profile = np.flipud(profile[1])  # m
        pressure_profile = np.flipud(profile[2])  # mbar
        temperature_profile = np.flipud(profile[3])  # K
        n2_profile = np.flipud(profile[4])  # vmr
        co2_profile = np.flipud(profile[5])  # vmr
        ch4_profile = np.flipud(profile[6])
        h2o_profile = np.flipud(profile[7])  # vmr
        liquid_cloud_profile = np.flipud(profile[8])  # vmr
        ice_cloud_profile = np.flipud(profile[9])  # vmr
        liquid_cloud_size_profile = np.flipud(profile[10])  # um
        ice_cloud_size_profile = np.flipud(profile[11])  # um
        self.atmosphere = Table([layers_profile, heights_profile,
                                 pressure_profile, temperature_profile,
                                 n2_profile, co2_profile, ch4_profile,
                                 h2o_profile, liquid_cloud_profile,
                                 ice_cloud_profile, liquid_cloud_size_profile,
                                 ice_cloud_size_profile],
                                names=["Layer", "Height", "Pressure",
                                       "Temperature", "N2", "CO2", "CH4", "H2O",
                                       "LiquidCloud", "IceCloud",
                                       "LiquidCloudSize", "IceCloudSize"])
        self.atmosphere["Height"] /= 1e3
        print("    Atmosphere Profile File Read")
        self.planet_data["MeanLiquidCloudSize"] = np.average(
            self.atmosphere["LiquidCloudSize"],
            weights=self.atmosphere["LiquidCloud"])
        self.planet_data["MeanIceCloudSize"] = np.average(
            self.atmosphere["IceCloudSize"],
            weights=self.atmosphere["IceCloud"])
        self.planet_data["LiquidCloudAbundance"] = int(
            any(i > 0. for i in self.atmosphere["LiquidCloud"]))
        self.planet_data["IceCloudAbundance"] = int(
            any(i > 0. for i in self.atmosphere["IceCloud"]))
        self.planet_data["EffectiveTemp"] = (self.planet_data["Insolation"]
                                             * (1-self.planet_data["Albedo"])
                                             / (4 * 5.67e-8)) ** 0.25
        self.planet_data["ScaleHeight"] = (
                8.3144598
                * np.mean(self.atmosphere["Temperature"])
                / self.planet_data["MWeightDry"]
                / self.planet_data["Gravity"])

        # Adding layers
        if atmosphere_ceiling != 0:
            # This handles the imaginary upper atmosphere levels from
            # the top layer to the specified pressure atmosphere_ceiling
            # Defines upper atmosphere abundances
            ind = self.n_downlayers-1
            (lyr, upheight, uppres, uptemp, upN2, upCO2, upCH4, upH2O, upLiq,
             upIce, uplsize, upisize) = self.atmosphere[ind]
            # Planet radius at top of atmosphere
            # Hypsometric equation to find top of atmosphere
            uptop = (-np.log(atmosphere_ceiling / uppres)
                     * self.planet_data["ScaleHeight"]
                     + upheight)

            up_heights = np.geomspace(upheight + 2, uptop, n_uplayers)
            up_pressures = uppres * np.exp(-(up_heights-upheight)
                                           / self.planet_data["ScaleHeight"])
            up_indexes = range(n_uplayers)
            for i, h, p in zip(up_indexes, up_heights, up_pressures):
                self.atmosphere.add_row(
                    [self.n_downlayers + i + 1, h, p, uptemp,
                     upN2, upCO2, upCH4, upH2O, upLiq, upIce,
                     uplsize, upisize])
            self.n_layers = self.n_downlayers + n_uplayers
            print("    Added {} layers to the atmosphere".format(n_uplayers))

    def zero_atmosphere(self):
        """This is an optional method that zeros out the atmosphere for
        background calculations. It does not work reliably. This function
        should be skipped for normal operations, although it may serve as a
        template for removing a particular atmospheric species such as CO2"""
        # self.atmosphere["Pressure"] *= 0
        # self.atmosphere["Temperature"] *= 0
        self.atmosphere["N2"] *= 0
        self.atmosphere["CO2"] *= 0
        self.atmosphere["H2O"] *= 0
        self.atmosphere["LiquidCloud"] *= 0
        self.atmosphere["IceCloud"] *= 0
        # self.planet_data["MWeightTotal"] = 0
        # self.planet_data["IceCloudAbundance"] = 0
        # self.planet_data["LiquidCloudAbundance"] = 0
        # self.planet_data["SurfacePressure"] = 0

    def write(self, scope: str = "MIRI-MRS", exposure_time=16,
              exposure_count: int = 110, rad_units: str = "rel"):
        """Writes a PSG input file using parameters defined in previous methods.

        Args:
            scope: The name of the scope, see global varaible psg_scopes for a
                list of supported scopes
            exposure_time: (number) The length of a single exposure, used for
                noise calculations only.
            exposure_count: (int) The number of exposures taken, used for noise
                calculations only.
            rad_units: (str) The desired units of the output, for a list of
                possible inputs, see global varaible rad_unit_types. I recommend
                rel, rkm, and Wsrm2um
        """
        self.scope = scope
        self.exposure_time = exposure_time
        self.exposure_count = exposure_count
        self.rad_units = rad_units
        self._psginput_name = ""
        # print(self.profile_file)
        self.profile_file = str(self.profile_file)
        if "aqua" in self.profile_file:
            self._psginput_name += self.profile_file.split("_aqua")[0]
            end = "_psginput.txt"
        elif ("t" in self.profile_file
              and "s" in self.profile_file
              and "p" in self.profile_file):
            self._psginput_name += self.profile_file.split(".txt")[0]
            end = "_psginput.txt"
        else:
            print("    Assuming the file is a dry file")
            self._psginput_name = self.profile_file.split("_dry")[0]
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
            results.write("<OBJECT-OBS-LONGITUDE>180\n")
            results.write("<OBJECT-OBS-LATITUDE>{:.5f}\n".format(
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
                results.write("<GENERATOR-LOGRAD>N\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>{}\n".format(
                    self.rad_units))
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
                results.write("<GENERATOR-LOGRAD>N\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>{}\n".format(
                    self.rad_units))
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
                results.write("<GENERATOR-LOGRAD>N\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>{}\n".format(
                    self.rad_units))
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
                results.write("<GENERATOR-LOGRAD>N\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>{}\n".format(
                    self.rad_units))
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
                results.write("<GENERATOR-LOGRAD>N\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>{}\n".format(
                    self.rad_units))
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
                results.write("<GENERATOR-LOGRAD>N\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>{}\n".format(
                    self.rad_units))
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
                results.write("<GENERATOR-LOGRAD>N\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>{}\n".format(
                    self.rad_units))
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
                results.write("<GENERATOR-LOGRAD>N\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>{}\n".format(
                    self.rad_units))
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
                results.write("<GENERATOR-LOGRAD>N\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>{}\n".format(
                    self.rad_units))
                results.write(
                    "<GENERATOR-SUMMARY>Wavelengths range 9.9-19.6 um with a "
                    "resolution of 600 RP. Molecular radiative-transfer "
                    "enabled; Continuum flux module enabled;\n")
                results.write("<GENERATOR-TRANS>02-01\n")

            elif self.scope == "OST-TRA":
                results.write("<GENERATOR-INSTRUMENT>OST_MISC-TRA: TRA "
                              "is the Transit Spectrometer of the Mid-"
                              "Infrared Imager / Spectrograph / Coronagraph("
                              "MISC) Instrument for the Origins Space "
                              "Telescope (OST) concept.The instrument offers "
                              "densified pupil spectroscopy with high "
                              "photometric precision (1ppm on timescales of "
                              "hours to days, excluding the fluctuation of "
                              "detector gain).\n")
                results.write("<GENERATOR-RANGE1>5\n")
                results.write("<GENERATOR-RANGE2>20\n")
                results.write("<GENERATOR-RANGEUNIT>um\n")
                results.write("<GENERATOR-RESOLUTION>300\n")
                results.write("<GENERATOR-RESOLUTIONUNIT>RP\n")
                results.write("<GENERATOR-TELESCOPE>SINGLE\n")
                results.write("<GENERATOR-DIAMTELE>9.0\n")
                results.write("<GENERATOR-BEAM>1\n")
                results.write("<GENERATOR-BEAM-UNIT>diffrac\n")
                results.write("<GENERATOR-TELESCOPE1>1\n")
                results.write("<GENERATOR-TELESCOPE2>0.0\n")
                results.write("<GENERATOR-TELESCOPE3>1.0\n")
                results.write("<GENERATOR-NOISE>NEP\n")
                results.write("<GENERATOR-NOISE1>2E-20\n")
                results.write("<GENERATOR-NOISEOTEMP>4.5\n")
                results.write("<GENERATOR-NOISEOEFF>0.1\n")
                results.write("<GENERATOR-NOISEOEMIS>0.1\n")
                results.write("<GENERATOR-NOISETIME>{}\n".format(
                    self.exposure_time))
                results.write("<GENERATOR-NOISEFRAMES>{}\n".format(
                    self.exposure_count))
                results.write("<GENERATOR-NOISEPIXELS>8\n")
                results.write("<GENERATOR-TRANS-APPLY>N\n")
                results.write("<GENERATOR-TRANS-SHOW>N\n")
                results.write("<GENERATOR-TRANS>03-01\n")
                results.write("<GENERATOR-LOGRAD>N\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>{}\n".format(
                    self.rad_units))
            
            elif self.scope == "OST-MIRI":
                results.write("<GENERATOR-INSTRUMENT> OST_MISC-MRS: MRS is "
                              "the Medium-Resolution Spectrometer of the "
                              "Mid-Infrared Imager / Spectrograph / "
                              "Coronagraph(MISC) Instrument for the Origins "
                              "Space Telescope (OST) concept.The instrument "
                              "offers slit and IFU spectroscopy at a variety "
                              "of spectral resolutions and slit-widths (0.33, "
                              "0.55 and 1.0 arcsec).\n")
                results.write("<GENERATOR-RANGE1>5\n")
                results.write("<GENERATOR-RANGE2>36\n")
                results.write("<GENERATOR-RANGEUNIT>um\n")
                results.write("<GENERATOR-RESOLUTION>1200\n")
                results.write("<GENERATOR-RESOLUTIONUNIT>RP\n")
                results.write("<GENERATOR-TELESCOPE>SINGLE\n")
                results.write("<GENERATOR-DIAMTELE>9.0\n")
                results.write("<GENERATOR-BEAM>1\n")
                results.write("<GENERATOR-BEAM-UNIT>diffrac\n")
                results.write("<GENERATOR-TELESCOPE1>1\n")
                results.write("<GENERATOR-TELESCOPE2>0.0\n")
                results.write("<GENERATOR-TELESCOPE3>1.0\n")
                results.write("<GENERATOR-NOISE>NEP\n")
                results.write("<GENERATOR-NOISE1>2E-20\n")
                results.write("<GENERATOR-NOISEOTEMP>4.5\n")
                results.write("<GENERATOR-NOISEOEFF>0.1\n")
                results.write("<GENERATOR-NOISEOEMIS>0.1\n")
                results.write("<GENERATOR-NOISETIME>{}\n".format(
                    self.exposure_time))
                results.write("<GENERATOR-NOISEFRAMES>{}\n".format(
                    self.exposure_count))
                results.write("<GENERATOR-NOISEPIXELS>8\n")
                results.write("<GENERATOR-TRANS-APPLY>N\n")
                results.write("<GENERATOR-TRANS-SHOW>N\n")
                results.write("<GENERATOR-TRANS>03-01\n")
                results.write("<GENERATOR-LOGRAD>N\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>Y\n")
                results.write("<GENERATOR-RADUNITS>{}\n".format(
                    self.rad_units))

            elif self.scope == "ALMA_Band7":
                results.write(
                    "<GENERATOR-INSTRUMENT>ALMA_Band7: ALMA(Atacama "
                    "Large Millimeter / submillimeter Array) is a mm / sub-"
                    "mm interferometer located on the Chajnantor plateau("
                    "Chile) at 5000 meters.High-resolution spectroscopy is "
                    "provided in the 31-950 GHz via a set of ten receivers "
                    "/ bands / frontends, combined with a versatile "
                    "correlator backend configuration.Highest spectral "
                    "resolution is 7.6 kHz, and highest spatial resolution is "
                    "10 milliarcseconds.Band-7 samples the 275-373 GHz range "
                    "employing SIS receivers.\n")
                results.write("<GENERATOR-RANGE1>345.7\n")
                results.write("<GENERATOR-RANGE2>345.9\n")
                results.write("<GENERATOR-RANGEUNIT>GHz\n")
                results.write("<GENERATOR-RESOLUTION>500\n")
                results.write("<GENERATOR-RESOLUTIONUNIT>kHz\n")
                results.write("<GENERATOR-TELESCOPE>ARRAY\n")
                results.write("<GENERATOR-DIAMTELE>12\n")
                results.write("<GENERATOR-BEAM>0.2\n")
                results.write("<GENERATOR-BEAM-UNIT>arcsec\n")
                results.write("<GENERATOR-TELESCOPE1>54\n")
                results.write("<GENERATOR-TELESCOPE2>0.0\n")
                results.write("<GENERATOR-TELESCOPE3>1.0\n")
                results.write("<GENERATOR-NOISE>TRX\n")
                results.write("<GENERATOR-NOISE1>147\n")
                results.write("<GENERATOR-NOISE2>0\n")
                results.write("<GENERATOR-NOISEOTEMP>270\n")
                results.write("<GENERATOR-NOISEOEFF>0.85\n")
                results.write("<GENERATOR-NOISEOEMIS>0.05\n")
                results.write("<GENERATOR-NOISETIME>{}\n".format(
                    self.exposure_time))
                results.write("<GENERATOR-NOISEFRAMES>{}\n".format(
                    self.exposure_count))
                results.write("<GENERATOR-NOISEPIXELS>8\n")
                results.write("<GENERATOR-TRANS-APPLY>N\n")
                results.write("<GENERATOR-TRANS-SHOW>Y\n")
                results.write("<GENERATOR-TRANS>02-01\n")
                results.write("<GENERATOR-LOGRAD>N\n")
                results.write("<GENERATOR-GAS-MODEL>Y\n")
                results.write("<GENERATOR-CONT-MODEL>Y\n")
                results.write("<GENERATOR-CONT-STELLAR>N\n")
                results.write("<GENERATOR-RADUNITS>Jy\n")
                results.write(
                    "<GENERATOR-SUMMARY>Wavelengths range 345.7-345.9 "
                    "GHz with a resolution of 500 kHz.Molecular "
                    "radiative-transfer enabled; Continuum flux module "
                    "enabled; Interferometric observations;\n")

            else:
                print("Scope input {} not understood. Only JWST instruments, "
                      "Hubble, and Spitzer-Short-High supported here"
                      .format(self.scope))

            # ATMOSPHERE
            results.write("<ATMOSPHERE-DESCRIPTION>{}-CU Boulder GCM (Eric "
                          "Wolf)\n".format(self.planet_data["Name"]))
            results.write("<ATMOSPHERE-STRUCTURE>Equilibrium\n")
            results.write("<ATMOSPHERE-WEIGHT>{}\n".format(
                self.planet_data["MWeightDry"]))
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
            results.write("<ATMOSPHERE-AEROS>Cloud,WaterIce\n")
            results.write("<ATMOSPHERE-ATYPE>"
                          "White_GSFC[reff=1.0um 0.10-1000000.00um],"
                          "CRISM_Wolff[reff=2.0um 0.31-99.57um]\n")
            results.write("<ATMOSPHERE-AABUN>{},{}\n".
                          format(self.planet_data["LiquidCloudAbundance"],
                                 self.planet_data["IceCloudAbundance"]))
            results.write("<ATMOSPHERE-AUNIT>scl,scl\n")
            results.write("<ATMOSPHERE-NMAX>0\n")
            results.write("<ATMOSPHERE-LMAX>0\n")
            results.write("<ATMOSPHERE-ASIZE>{},{:.3f}\n".format(
                10,
                self.planet_data["MeanIceCloudSize"]))
            results.write(
                "<ATMOSPHERE-LAYERS-MOLECULES>Altitude,N2,CO2,H2O,Cloud"
                ",WaterIce\n")
            results.write("<ATMOSPHERE-LAYERS>{}\n".format(self.n_layers))
            for i, lvl in enumerate(self.atmosphere):
                results.write(
                    "<ATMOSPHERE-LAYER-{:.0f}>{:.3E},{:.3E},{:.3E},"
                    "{:.3E},{:.3E},{:.3E},{:.3E},{:.3E}\n".format(
                        lvl["Layer"], lvl["Pressure"], lvl["Temperature"],
                        lvl["Height"], lvl["N2"], lvl["CO2"], lvl["H2O"],
                        lvl["LiquidCloud"], lvl["IceCloud"]))

            results.write("<SURFACE-TEMPERATURE>{}\n".format(
                self.planet_data["SurfaceTemperature"]))
            results.write("<SURFACE-ALBEDO>{}\n".format(
                self.planet_data["Albedo"]))
            results.write("<SURFACE-EMISSIVITY>{}\n".format(
                1-self.planet_data["Albedo"]))
            print("    Successfully created PSG Config file from GCM results"
                  " for {}".format(self.planet_data["Name"]))
            print("    The file's name is {}".format(self._psginput_name))

    def send(self, keep_files=("trn", "lyr", "rad", "noi", "log", "atm", "err"),
             run=True, key=None):
        """This function send the the file cdf_file to the NASA GSFC PSG for
        analysis It will not return anything, but will write files in the 
        current directory.

        Arguments: keep_files: a tuple of 3 letter strings which can include:
            trn, lyr, rad, noi, err, cfg, atm, log, str
        run: (bool) of whether or not you actually want to send it to the PSG.
            You can set this to false if you have previously made files with
            correct names
        key: If you want to run more than a limited amount of pings to the PSG,
            you may need a key. The PSG will kick you off if you've accessed it
            too many times. The key must be given by Geronimo Villanueva
        """
        # Cumbersome giant file of everything
        alloutputname = self._psginput_name.split("psginput.txt")[
                            0] + "psgoutput_all.txt"
        filestem = self._psginput_name.split("psginput.txt")[0] + "psgoutput_"
        print("    Sending to PSG")
        if run:
            if key:
                command = (f"curl -d key={key} -d type=all --data-urlencode "
                           f"--speed-time 30 file@{self._psginput_name} "
                           f"https://psg.gsfc.nasa.gov/api.php > "
                           f"{alloutputname}")
            else:
                command = f"curl -d type=all --data-urlencode " \
                          f"file@{self._psginput_name} " \
                          f"https://psg.gsfc.nasa.gov/api.php > {alloutputname}"
            # print("    {}".format(command))
            os.system(command)
        print("    Successfully connected to NASA PSG")
        with open(alloutputname, "r+") as allfile:
            sections = 0
            for i, line in enumerate(allfile):
                if i == 0:
                    if "Your other API call is still running" in line:
                        raise Exception("PSG is hung up on a previous run")
                if "results_" in line:
                    sections += 1
            if sections == 0:  # Happens if the PSG is down or file is incorrect
                raise Exception("No results returned from PSG")
            allfile.seek(0)
            filetail = allfile.readline().split("_")[1].split("\n")[0]
            # print(filetail)
            # print(sections)
            self.returned_files = []
            for _ in range(sections):
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
        if filestem + "log.txt" in self.returned_files:
            print("    PSG returned an error or warning:")
            with open(filestem + "log.txt", "r") as log:
                print(log.readline())

        for i, fil in enumerate(self.returned_files):
            if fil.split("_")[-1].split(".")[0] not in keep_files:
                os.system("rm -v {}".format(fil))

    def plot_setup(self):
        """This takes the results from the PSG and parses them to create arrays
        that can easily be accessed for analysis from the PSG parent object.
        It takes no arguments, but it makes available x.Wavelengths, x.Total,
        etc. Variables beginning with n are noise variables. Variables beginning
        with t are transmission variables (atmospheric transparency between 0
        and 1.

        """
        rad_file = None
        # print(self.returned_files)
        for fil in self.returned_files:
            if "rad" in fil:
                rad_file = fil
        if rad_file is None:
            raise Exception("No radfile returned. The PSG likely lost"
                            " connection")
        self._file_stem = rad_file.split("_psg")[0]
        if ("/" in self._file_stem) or ("\\" in self._file_stem):
            self._file_stem = self._file_stem.replace("\\", "/")
            self._file_stem = self._file_stem.split("/")[-1]
        name_parts = self._file_stem.split("_")
        self._title_stem = " ".join(name_parts)
        radfil = np.loadtxt(rad_file, unpack=True)
        if len(radfil) == 7:
            # print("7 items")
            self.is_transit = True
            self.Wavelengths = radfil[0]
            self.Total = radfil[1]
            self.Noise = radfil[2]
            self.Stellar = radfil[3]
            self.Planet = radfil[4]
            self.Transit = radfil[5]
            self.Thermal = radfil[6]
        elif len(radfil) == 6:
            self.is_transit = True
            self.Wavelengths = radfil[0]
            self.Total = radfil[1]
            self.Stellar = radfil[2]
            self.Thermal = radfil[3]
            self.Transit = radfil[4]
            self.Planet = radfil[5]
        elif len(radfil) == 4:
            self.is_transit = False
            self.Wavelengths = radfil[0]
            self.Total = radfil[1]
            self.Stellar = radfil[2]
            self.Thermal = radfil[3]
        else:
            self.Wavelengths = radfil[0]
            self.Total = radfil[1]
            self.Noise = radfil[2]
            self.Stellar = radfil[3]
            self.Thermal = radfil[4]
            print("    Exoplanet wasn't transitting")
            self.is_transit = False

        self._plot_range = (self.Wavelengths.min(), self.Wavelengths.max())

        if os.path.isfile(self._file_stem + "_psgoutput_trn.txt"):
            trnfil = np.loadtxt(self._file_stem + "_psgoutput_trn.txt",
                                unpack=True)
            self.tWavelengthSpec = trnfil[0]
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

        if os.path.isfile(self._file_stem + "_psgoutput_noi.txt"):
            noi_fil = np.loadtxt(self._file_stem + "_psgoutput_noi.txt",
                                 unpack=True)
            self.nWavelengthSpec = noi_fil[0]
            self.nTotal = noi_fil[1]
            self.nSource = noi_fil[2]
            self.nDetector = noi_fil[3]
            self.nTelescope = noi_fil[4]
            self.nBackground = noi_fil[5]
            self.nReal = np.random.normal(0, self.nTotal, len(self.nTotal))

        print("    Ready to Make Plots")

    def depth_plot(self):
        fig = plt.figure(figsize=(10, 6))
        ax = plt.gca()
        ax.step(self.Wavelengths, -self.Transit / self.Stellar * 1e+6,
                linewidth=0.25, c="b", label="Transit Depth",
                where="post")
        ax.set_title("Transit Depth Without Noise\n{}".format(
            self._title_stem))
        ax.set_xlabel(r"Wavelengths ($\mu m$)")
        ax.set_ylabel("Signal (ppm)")
        ax.set_xlim(*self._plot_range)
        ax.xaxis.grid(True)
        fig.savefig("{}_dpth.png".format(self._file_stem))
        plt.close(fig)

    def depth_height(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        true_depth = -self.Transit / self.Stellar-self.planet_data["DepthEff"]
        den = (2 * self.planet_data["Radius"]
               * self.planet_data["ScaleHeight"])
        ax.step(self.Wavelengths, (true_depth * self.star_data["Radius"] ** 2
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
        plt.close(fig)

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
        plt.close(fig)

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
                / (wave * 1.38e-23 * self.star_data["Temperature"]))-1)
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
        plt.close(fig)

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
        plt.close(fig)

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
        plt.close(fig)

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
        plt.close(fig)

    def signal_noise_ratio(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        ax.step(self.Wavelengths, -self.Transit / self.Noise * 100,
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
        plt.close(fig)

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
        plt.close(fig)

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
        plt.close(fig)

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
        plt.close(fig)

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
        plt.close(fig)

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
        plt.close(fig)

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
        plt.close(fig)

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
        plt.close(fig)

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
        ax.legend(["Total", "CO2", "N2", "H2O"], loc=4)
        ax.set_title("Species Transmittance\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Transmittance (normalized)")
        ax.set_xlim(*self._plot_range)
        ax.set_ylim(0, 1.1)
        ax.xaxis.grid(True)
        fig.savefig("{}_tspc.png".format(self._file_stem))
        plt.close(fig)

    def trn_aerosols(self):
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
        ax.legend(["Total", "Ice", "Cloud", "tCIASpec"], loc=4)
        ax.set_title("tCIASpec Transmittance\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Transmittance (normalized)")
        ax.set_xlim(*self._plot_range)
        ax.set_ylim(0, 1.1)
        ax.xaxis.grid(True)
        fig.savefig("{}_taer.png".format(self._file_stem))
        plt.close(fig)

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
        plt.close(fig)

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
        plt.close(fig)

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
        plt.close(fig)

    def depth_noise(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.gca()
        real_depth = (-self.Transit + np.random.normal(0, self.nTotal)
                      ) / self.Stellar
        ax.step(self.Wavelengths, real_depth * 1e6, linewidth=0.5, c="b",
                where="post")
        ax.set_title("Depth Plotted with Noise\n{}".format(
            self._title_stem))
        ax.set_xlabel("Wavelengths ($\mu m$)")
        ax.set_ylabel("Depth (ppm)")
        ax.set_xlim(*self._plot_range)
        ax.xaxis.grid(True)
        fig.savefig("{}_depth_nois.png".format(self._file_stem))
        plt.close(fig)

    def pandexo_noise(self, scope):
        """This function creates a set instrument parameters using pandas
        This is meant to give us noise simulations to more accurately predict
        transits. Check parent class for more details"""

        # Create a pandexo ready transit file
        self.pand_fil = self._file_stem + "_depth.txt"
        with open(self.pand_fil, "w") as fil:
            for wave, dep in zip(self.Wavelengths, self.Transit):
                fil.write("{:10.7e}{:10.7f}\n".format(wave*1e-6, dep))
        print(self.pand_fil)
        # Object Parameters
        exo_dict = jdi.load_exo_dict()
        exo_dict["observation"]["sat_level"] = 80
        exo_dict["observation"]["sat_unit"] = "%"
        exo_dict["observation"]["noccultations"] = 1
        exo_dict["observation"]["R"] = None
        exo_dict["observation"]["baseline"] = 4. * 60. * 60.
        exo_dict["observation"]["baseline_unit"] = "total"
        exo_dict["observation"]["noise_floor"] = 0

        exo_dict["star"]["type"] = "phoenix"
        exo_dict["star"]["mag"] = self.star_data["Magnitude"]
        exo_dict["star"]["ref_wave"] = 1.25
        exo_dict["star"]["temp"] = self.star_data["Temperature"]
        exo_dict["star"]["metal"] = self.star_data["Metallicity"]
        exo_dict["star"]["logg"] = self.star_data["Gravity"]
        exo_dict["star"]["radius"] = self.star_data["SRadius"]
        exo_dict["star"]["r_unit"] = "R_sun"

        exo_dict["planet"]["type"] = "user"
        exo_dict["planet"]["exopath"] = self.pand_fil
        exo_dict["planet"]["w_unit"] = "um"
        exo_dict["planet"]["f_unit"] = "rp^2/r*^2"
        exo_dict["planet"]["transit_duration"] = (self.exposure_count
                                                  * self.exposure_time)
        exo_dict["planet"]["td_unit"] = "s"
        exo_dict["planet"]["temp"] = self.planet_data["SurfaceTemperature"]
        exo_dict["planet"]["radius"] = self.planet_data["Radius"]
        exo_dict["planet"]["r_unit"] = u.meter
        exo_dict["planet"]["mass"] = self.planet_data["Mass"]
        exo_dict["planet"]["m_unit"] = u.kilogram

        # Instrument Parameters
        inst_dict = jdi.load_mode_dict(scope)

        # Runs PandExo
        output = os.path.abspath(self._file_stem + "_pandexo.txt")
        print(output)
        self.pandexo_result = jdi.run_pandexo(exo_dict, inst_dict,
                                              output_path=output)
