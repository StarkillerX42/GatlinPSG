import PSG
import numpy as np
import glob

gcmfiles=glob.glob('trappist1e*aqua_mean*lon0.txt')
#print (gcm_files)
for i,file in enumerate(gcmfiles):
    print ("")
    print ("File {}: {}".format(i+1,file))
    lon=int(file.split('lon')[-1].split('.txt')[0])
    if lon<180:
        phasein=180+lon
    elif lon>=180:
        phasein=lon-180
    x=PSG.PSG('TRAPPIST-1 e', file, scope='MIRI-MRS', exoearth=False,
              atmosphere_ceiling=1e-6, uplayers=7, exposure_time=30, exposure_count=110, skprow=8, phase=phasein)
    x.calculate()
    x.send(atm=False, cfg=False, lay=False, noi=True, rad=True, str=False, trn=True, all=False, err=True)
    #x.PlotSetup()
    #x.depth_plot()
    #x.Emission()
    #x.Raw()
    #x.Star()
    #x.signal_and_noise()
    #x.signal_noise_ratio()
    #x.Absorbtion()
    #print ("    Radiance Plots Complete")
    #x.trn_total()
    #x.trn_co2()
    #x.trn_n2()
    #x.trn_h2o()
    #x.trn_ice()
    #x.trn_cloud()
    #x.trn_cia()
    #x.trn_species()
    #x.trn_aerosols()
    #x.trn_all()
    x.noise_components()
print ("PSG: Operation complete")
