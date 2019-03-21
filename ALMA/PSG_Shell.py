import glob

import PSG

gcmfiles = glob.glob('*lon*.txt')
# print (gcm_files)
filecount = 1
for i, fil in enumerate(gcmfiles):
    print ("")
    lon = int(fil.split('lon')[-1].split('.txt')[0])
    if lon < 180:
        phasein = 180 + lon
    elif lon >= 180:
        phasein = lon - 180
    print ("File {}: {}".format(i + 1, fil))
    x = PSG.PSG('TRAPPIST-1 e', fil, scope='ALMA_Band7', exoearth=False, atmosphere_ceiling=1e-6, uplayers=7, exposure_time=3600,
                exposure_count=1, skprow=8, phase=phasein)
    x.calculate()
    x.send(run=True)
    # x.plot_setup()
    # x.depth_plot()
    # x.depth_height()
    # x.Emission()
    # x.Raw()
    # x.Star()
    # x.signal_and_noise()
    # x.signal_noise_ratio()
    # x.Absorbtion()
    print ("    Radiance Plots Complete")
    # x.trn_total()
    # x.trn_co2()
    # x.trn_n2()
    # x.trn_h2o()
    # x.trn_ice()
    # x.trn_cloud()
    # x.trn_cia()
    # x.trn_species()
    # x.trn_aerosols()
    # x.trn_all()
    # x.noise_components()
    # x.noise_ppm()
print ("PSG: Operation complete")
