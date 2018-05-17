import PSG
import numpy as np
import glob

gcmfiles=glob.glob('*terminator.txt')
print (gcmfiles)
filecount=1
for i,file in enumerate(gcmfiles):
    print ("")
    print ("File {}: {}".format(i+1,file))
    x=PSG.PSG('TRAPPIST-1 e', file, scope='MIRI-MRS', exoearth=False,
              atmosphere_ceiling=1e-6, uplayers=7, exposure_time=30, exposure_count=110, skprow=11)
    x.calculate()
    x.send(atm=True, cfg=True, lay=True, noi=True, rad=True, str=True, trn=True, all=True, err=True)
    x.plot_setup()
    x.depth_plot()
    x.emission()
    x.raw()
    x.star()
    x.signal_and_noise()
    x.signal_noise_ratio()
    x.absorption()
    print ("    Radiance Plots Complete")
    x.trn_total()
    x.trn_co2()
    x.trn_n2()
    x.trn_h2o()
    x.trn_ice()
    x.trn_cloud()
    x.trn_cia()
    x.trn_species()
    x.trn_aero()
    x.trn_all()
    x.noise_components()
print ("PSG: Operation complete")
