#Phase Curve Analysis
import numpy as np
import matplotlib.pyplot as plt
import glob
import matplotlib.animation as ani
import matplotlib.axes as ax
import starcoder42 as s
%matplotlib inline

radfiles=glob.glob("trappist1e*psgoutput_rad.txt")
print (radfiles)
phases=[]
data=[]
for i, radfile in enumerate(radfiles):
    phases.append(radfile.split('CO2_')[1].split('_')[0])
    a=np.loadtxt(radfile,unpack=True,skiprows=16)
    data.append(a)
#files are in phase*attribute*value dim
#13*5*4162
print (phases)
phases=np.array(phases)
data=np.array(data)

s.Describe(data)
testfile=data[4]
plt.figure(figsize=(100,4))
plt.step(testfile[0],testfile[4],where='post')
plt.axvline(5.045)
plt.axvline(7.67)
plt.axvline(8.415)
plt.axvline(9.4)
plt.axvline(10.05)
plt.axvline(10.405)
plt.axvline(12)
#relative values
plt.axvline(19.48)
plt.axvline(20.46)
plt.axvline(20.16)
plt.axvline(22)
plt.savefig("phase_210_blackbody.png")

#this is a padding value
values=[5.045,7.67,8.415,9.4,10.05,10.405,12.,19.48,20.46,20.16,22.]
indlist=[]
for i,j in enumerate(values):
    ind=1
    count=1

    while testfile[0,ind]<j and count<10000:
        ind+=1
        count+=1
    indlist.append(ind)
#Phase curve at specific wavelengths
print (indlist)
for ind in indlist:
    plt.scatter(phases,data[:,4,ind]*1e9)
    plt.title(ind)
    plt.xlabel('Phase (deg)')
    plt.ylabel('Brightness nW/m^2')
    savename='phaseplot_{:.2f}.png'.format(data[0,0,ind])
    plt.savefig(savename)
    plt.cla()
s.Describe(data[:,4,1000:1500]/data[:,3,1000:1500]*1.0e6)


seplist=[]
for i,wave in enumerate(data[0,0,:]):
    sep=(np.percentile(data[:,4,i],95)-np.percentile(data[:,4,i],5))
    seplist.append(sep)
plt.figure(figsize=(10,6))
plt.step(data[0,0,:],np.array(seplist),where='post')
plt.xlabel(r'Wavelength ($\si{\micro \meter}$)')
plt.ylabel(r"Signal ($\si{\nano\watt\per\meter^{2}\per\micro\meter}$)")
plt.title(r"Phase Curve Difference (max-min)")
plt.savefig("Phase_Curve_Difference.png")



#Video plotting
FFMpegWriter=ani.writers['ffmpeg']
writer=FFMpegWriter(fps=10)
fig=plt.figure(figsize=(8,8))
with writer.saving(fig,"Phase_Curves_Star.mp4",100):
    for i,wave in enumerate(data[0,0,:]):
        plt.scatter(phases,data[:,4,i]/data[:,3,i]*1.0e6)
        plt.title("Thermal Phase Curve wrt Parent star for Index {}".format(i))
        plt.xlabel(r'Phase (degrees)')
        plt.ylabel(r"Signal (ppm)")
        plt.ylim(0,100)
        writer.grab_frame()
        if i % 100 ==0:
            print ("Saved frame {}/{}".format(i,len(data[0,0,:])))
        plt.cla()




i=3540
plt.scatter(phases,data[:,4,i]/data[:,3,i]*1.0e6)
plt.title("Thermal Phase Curve wrt Parent star for Index {}".format(i))
plt.xlabel(r'Phase (degrees)')
plt.ylabel(r"Signal (ppm)")
#plt.ylim(0,100)
goodindexes=[3487,]


phaseb=[]
for i, phase in enumerate(phases):
    phaseb.append(np.mean(data[i,4,3007:]/data[i,3,3007:])*1.0e6)
s.Describe(np.array(phaseb))
plt.figure(figsize=(10,6))
plt.scatter(phases,phaseb)
plt.title("Signal to Noise Ratio of Thermal Phase Curve from 17-28um")
plt.xlabel(r'Phase (degrees)')
plt.ylabel(r"Signal (ppm)")
#plt.ylim(0,150)
plt.xlim(0,360)
ax=plt.axes()
ax.xaxis.grid(True)
plt.savefig("TPC_17_28_avg.png")


phaseb=[]
for i, phase in enumerate(phases):
    phaseb.append(np.mean(data[i,4,3007:]/data[i,2,3007:]))
s.Describe(np.array(phaseb))
plt.figure(figsize=(10,6))
plt.scatter(phases,phaseb)
plt.title("Signal to Noise Ratio of Thermal Phase Curve from 17-28um")
plt.xlabel(r'Phase (degrees)')
plt.ylabel(r"Signal (ppm)")
plt.xlim(0,360)
ax=plt.axes()
ax.xaxis.grid(True)
plt.savefig("TPC_17_28_snr.png")
