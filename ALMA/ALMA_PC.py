import numpy as np
import glob
import starcoder42 as s
import matplotlib.pyplot as plt

radfiles = glob.glob("trappist1e*psgoutput_rad.txt")
print(radfiles)
phases = []
data = []
for i, radfile in enumerate(radfiles):
	phases.append(radfile.split('CO2_')[1].split('_')[0])
	a = np.loadtxt(radfile, unpack=True, skiprows=16)
	data.append(a)
data = np.array(data)
s.describeArr(data[0, 0, :], printIt=True)
intlist = []
for i, dat in enumerate(data):
	freq = dat[0]
	tota = dat[1]
	stel = dat[2]
	trap = dat[3]
	freqdiff = np.diff(freq)
	traptrap = trap[1:] + trap[:-1]
	trapez = traptrap * freqdiff / 2
	integ = np.sum(trapez)
	intlist.append(integ)


s.Describe(intlist)
plt.scatter(range(len(intlist)), intlist)
