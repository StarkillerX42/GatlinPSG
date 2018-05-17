import PSG
import numpy as np
import glob
import starcoder42 as s
import matplotlib.pyplot as plt

gcmfiles=glob.glob('*terminator.txt')
print (gcmfiles)
filecount=1
for i,file in enumerate(gcmfiles):
    print ("")
    print ("File {}: {}".format(i+1,file))
    tran=PSG.PSG('TRAPPIST-1 e', file, atmosphere_ceiling=1e-6, uplayers=7, exposure_time=30, exposure_count=110, skprow=11)
    occult=PSG.PSG('TRAPPIST-1 e', file, atmosphere_ceiling=1e-6, uplayers=7, exposure_time=30, exposure_count=110, skprow=11, phase=0)
    almost=PSG.PSG('TRAPPIST-1 e', file, atmosphere_ceiling=1e-6, uplayers=7, exposure_time=30, exposure_count=110, skprow=11, phase=178)
    tran.calculate()
    occult.calculate()
    almost.calculate()
    tran.send()
    occult.send()
    almost.send()
    tran.plot_setup()
    occult.plot_setup()
    almost.plot_setup()
Wavelengths=almost.Wavelength
DeltaB=almost.Total+almost.Noise-tran.Total-tran.Noise
Star=occult.Total+occult.Noise
s.describe(tran.Stellar)
plt.step(Wavelengths,(tran.Total-occult.Total)/occult.Total)
