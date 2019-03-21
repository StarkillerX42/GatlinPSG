#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  2 09:49:34 2018

@author: tfauchez
"""

from netCDF4 import Dataset
import numpy as np


R_EARTH = 1.0  # the area_weighted_average is indepedent of planet radius
fill_value = -999.0  # set fill-value of variable being used
P = 10  # N2 pressure
Tstar = 4500  # star temperature
Mstar = 0.72
Lstar = 0.189
flux = [1999, 2199, 2249]  # , 1575] #stellar flux at TOA
period = [92.59, 86.20, 84.76]  # , 7.65] #orbital period
lat1 = 46
lon1 = 72
alt = 40

# pick up values from CAM outputs
specific_hum = np.zeros((len(flux), 1, alt, lat1, lon1))
TS = np.zeros((len(flux), 1, lat1, lon1))
FSNT = np.zeros((len(flux), 1, lat1, lon1))
FLNT = np.zeros((len(flux), 1, lat1, lon1))
FSNS = np.zeros((len(flux), 1, lat1, lon1))
FLNS = np.zeros((len(flux), 1, lat1, lon1))
CLDTOT = np.zeros((len(flux), 1, lat1, lon1))
lat = np.zeros((len(flux), lat1))
lon = np.zeros((len(flux), lon1))

for item in range(0, len(flux)):
    i = flux[item]
    j = period[item]
    nc = Dataset("outputs/{}bars/T{}/{}F/T{}K_{}M_{}L_{}F_{"
                 "}days_10bars.cam.h0.avg.nc".format(P, Tstar, i, Tstar,
                                                     Mstar, Lstar, i, j))
    specific_hum[item, :, :, :, :] = nc['Q'][:]  # alt lat lon
    TS[item, :, :, :] = nc['TS'][:]  # surface temperature
    FSNT[item, :, :, :] = nc['FSNT'][
                          :]  # net shortwave flux at the top of the model
    FLNT[item, :, :, :] = nc['FLNT'][
                          :]  # net longwave flux at the top of the model
    FSNS[item, :, :, :] = nc['FSNS'][:]  # net shortwave flux at the surface
    FLNS[item, :, :, :] = nc['FLNS'][:]  # net longwave flux at the surface
    CLDTOT[item, :, :, :] = nc['CLDTOT'][:]
    lat[item, :] = nc['lat'][:]
    lon[item, :] = nc['lon'][:]  # degree East

# calculate water mixing ratio and cloud coverage at the terminator
Q_r = np.zeros((len(flux), 46))
Q_l = np.zeros((len(flux), 46))
MMR_r = np.zeros((len(flux), 46))
MMR_l = np.zeros((len(flux), 46))
VMR_DRY_r = np.zeros((len(flux), 46))
VMR_DRY_l = np.zeros((len(flux), 46))
VMR_WET_r = np.zeros((len(flux), 46))
VMR_WET_l = np.zeros((len(flux), 46))

CLDTOT_r = np.zeros((len(flux), 46))
CLDTOT_l = np.zeros((len(flux), 46))

for i in range(0, len(flux)):
    Q_r[i, :] = np.squeeze(np.asarray(
        specific_hum[i, 0, 39, :, 18]))  # right hemisphere terminator
    Q_l[i, :] = np.squeeze(
        np.asarray(specific_hum[i, 0, 39, :, 54]))  # left hemisphere terminator

    MMR_r[i, :] = Q_r[i, :] / (1.0 - Q_r[i, :])
    MMR_l[i, :] = Q_l[i, :] / (1.0 - Q_l[i, :])

    VMR_DRY_r[i, :] = Q_r[i, :] / (1.0 - Q_r[i, :]) * 28.0059 / 18.01528
    VMR_DRY_l[i, :] = Q_l[i, :] / (1.0 - Q_l[i, :]) * 28.0059 / 18.01528

    VMR_WET_r[i, :] = VMR_DRY_r[i, :] / (1. + VMR_DRY_r[i, :])
    VMR_WET_l[i, :] = VMR_DRY_l[i, :] / (1. + VMR_DRY_l[i, :])

    CLDTOT_r[i, :] = np.squeeze(np.asarray(CLDTOT[i, 0, :, 18]))
    CLDTOT_l[i, :] = np.squeeze(np.asarray(CLDTOT[i, 0, :, 54]))

# average value around the terminator
average_WET_r = np.zeros((len(flux)))
average_WET_l = np.zeros((len(flux)))
average_VMR_wet_tot = np.zeros((len(flux)))

average_CLDTOT_r = np.zeros((len(flux)))
average_CLDTOT_l = np.zeros((len(flux)))
average_CLDTOT_tot = np.zeros((len(flux)))

average_VMR_wet_r = np.mean(VMR_WET_r, axis=1)
average_VMR_wet_l = np.mean(VMR_WET_l, axis=1)
average_VMR_wet_tot = (average_VMR_wet_r + average_VMR_wet_l) / 2

average_CLDTOT_r = np.mean(CLDTOT_r, axis=1)
average_CLDTOT_l = np.mean(CLDTOT_l, axis=1)
average_CLDTOT_tot = (average_CLDTOT_r + average_CLDTOT_l) / 2

# print in a .txt file
for item in range(0, len(flux)):
    i = flux[item]
    TSmax[item] = TS[item, :, :, :].max(axis=(1, 2))
    TSmin[item] = TS[item, :, :, :].min(axis=(1, 2))
    print(i)
    output_files = open("outputs/10bars/T4500/{}F/VMR".format(i), "w")
    output_files.write('N2 pressure: ')
    output_files.write(str(P))
    output_files.write('; Star temperature: ')
    output_files.write(str(Tstar))
    output_files.write('; Stellar flux at TOA: ')
    output_files.write(str(flux[item]))
    output_files.write('; Orbital period: ')
    output_files.write(str(period[item]))
    # print the output you want
    #    output_files.write('\nMean Surface temperature [K]: ')
    #    output_files.write(str(TS_avg[item]))
    #    output_files.write('\nMax Surface temperature [K]: ')
    #    output_files.write(str(TSmax[item]))
    #    output_files.write('\nMin Surface temperature [K]: ')
    #    output_files.write(str(TSmin[item]))
    #    output_files.write('\nH2O volume mixing ration [kg/kg]: ')
    #    output_files.write(str(average_VMR_wet_tot[item]))
    #    output_files.write('\nTOA energy balance: ')
    #    output_files.write(str(average_FNT_tot[item]))
    #    output_files.write(' W.m^2')
    #    output_files.write('\nsurface energy balance: ')
    #    output_files.write(str(average_FNS_tot[item]))
    #    output_files.write(' W.m^2')
    #    output_files.write('\ncloud coverage: ')
    #    output_files.write(str(average_CLDTOT_tot[item]))
    output_files.close()
