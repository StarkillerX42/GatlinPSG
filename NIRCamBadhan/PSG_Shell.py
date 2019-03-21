import glob
import PSG

gcm_files = glob.glob('*terminator.txt')
# print (gcm_files)
filecount = 1
for i, fil in enumerate(gcm_files):
    print("")
    print("File {}: {}".format(i + 1, fil))
    x = PSG.PSG("TRAPPIST-1 e", fil, scope='MIRI-MRS',
                atmosphere_ceiling=1e-6, n_uplayers=7, exposure_time=16,
                exposure_count=114*10)
    x.calculate(skprow=11)
    x.write()
    x.send(run=True)
    x.plot_setup()
    x.depth_plot()
    x.depth_height()
    x.emission()
    x.raw()
    x.star()
    x.signal_and_noise()
    x.signal_noise_ratio()
    x.absorption()
    print("    Radiance Plots Complete")
    x.trn_total()
    x.trn_co2()
    x.trn_n2()
    x.trn_h2o()
    x.trn_ice()
    x.trn_cloud()
    x.trn_cia()
    x.trn_species()
    x.trn_aerosols()
    x.trn_all()
    x.noise_components()
    x.noise_ppm()
    x.depth_noise()
print("PSG: Operation complete")
