# createArrayConfig.py
# Alex Hegedus 12/11/19
# Creates configuration for large array near equator of the Moon.
# Currently gives 10 km diameter Logarithmically spaced circular configuration
# also creates plots of the arrays
# saves and prints the longitudes, latitudes, and altitudes for all antenna (1024 currently)

from pylab import *
import numpy as np

import matplotlib.cm as cm


lunarPath = '/mnt/LinuxData/Lunar_Pipeline/SLDEM2015_128_60S_60N_000_360_FLOAT.IMG'


# m per degree of longitude at lunar equator
mpdeg = 29670.597283903604


numAnts = 6

#pts to make path of rover
pts = [[0., 0.], [15., 0], [15 + sqrt(50), -sqrt(50)], [15+2*sqrt(50), 0], [35., 0], [35 + sqrt(50), -sqrt(50)], [35+2*sqrt(50), 0], [80., 0], [80 + sqrt(50), -sqrt(50)], [80+2*sqrt(50), 0], [100.,0.]]

pts = array(pts)
#antenna locations on
ants = [[15 + sqrt(50)/2, -sqrt(50)/2], [15 + 3*sqrt(50)/2, -sqrt(50)/2], [35 + sqrt(50)/2, -sqrt(50)/2], [80 + 3*sqrt(50)/2, -sqrt(50)/2]]
ants = array(ants)



#new L shape configurations

pts = [[0., 0.], [0, -45], [45, -45]]

pts = array(pts)


ants = [[0., -10], [0,-30], [0, -40], [5, -45], [15, -45], [35, -45]]

ants = array(ants)

antEnds = [([0., -5], [0, -15]), \
        ([0, -25], [0, -34.6]), \
        ([0, -35.4], [0, -45]),
        ([0, -45], [9.6, -45]), \
        ([10.4, -45], [20, -45]),\
        ([30, -45], [40, -45])]

antEnds = array(antEnds)

antEnds[:3, :, 0] -= 0.5
antEnds[3:, :, 1] -= 0.5




Antxs = ants[:,0]
Antys = ants[:,1]

# Antxs = zeros(numAnts)
# Antys = zeros(numAnts)
alts = zeros(numAnts)

# normalize to Lunar Longitude and Latitude, mpdeg = 29670.6 m/degree equatorial Longitude
Antxs /= mpdeg
Antys /= mpdeg

lunarRad = 1737.4*1e3 #km


#min rms is (array([261]), array([52])) with offset 11, 11 into 4x4 degrees around 0, 0 long lat
# 0.002817279089066556 rms around this 5x5 km point in km
midLong = -2. + 63/128.
midLat = 2 - 272./128
#
# #approx
midLong = -1.508
midLat = -0.125


#now for min rms in 10x 10 km
# rmss[289,101] with offset 22, 22
# midLong = -2. + 123/128.
# midLat = 2 - 311./128
#
# #approx
# midLong = -1.04
# midLat = -0.43


#min rms 20x20 km
# (array([215]), array([4])) with offset 44, 44
# midLong = -2. + 259/128.
# midLat = 2 - 48./128
#
# midLong = .023
# midLat = 1.625



Antxs += midLong
Antys += midLat



# take care of wraparound if needed
for i in range(len(Antxs)):
    if Antxs[i] > 360.:
        Antxs[i] -= 360.

    if Antxs[i] < -180.:
        Antxs[i] += 360.

    if Antys[i] > 180.:
        Antys[i] = 180 - (Antys[i] - 180.) #360.

    if Antys[i] < -180.:
        Antys[i] = -180 - (Antys[i] + 180.) #360.







fid = open(lunarPath, 'rb')
data = np.fromfile(fid, dtype=np.dtype(np.float32))
fid.close()


shp = ((15360, 46080))

scale = 0.236901 #km per pix
Simage = data.reshape(shp)

#reshape so 0 Long is in the middle

Simage = append(Simage[:, shp[1]/2:], Simage[:, :shp[1]/2], axis = 1)


#use DSMAP.CAT to change to indices in file to get altitudes (in file as km)

centerLong = 0. #180.
centerLat = 0.

lineProjOffset = 7679.5
sampleProjOffset = 23039.5
resPixperDeg = 128

sampXs = zeros(numAnts)
lineYs = zeros(numAnts)

#(LAT, LON) to LINE and SAMPLE indices
for i in range(numAnts):
    sampX = int(np.round(sampleProjOffset + resPixperDeg*(Antxs[i] - centerLong)))
    lineY = int(np.round(lineProjOffset - resPixperDeg*(Antys[i] - centerLat))) #minus bc lower indices are higher latitudes.  line 1 is 60 lat
    alts[i] = Simage[lineY][sampX%46080]
    sampXs[i] = sampX
    lineYs[i] = lineY


#save to file
savetxt("configLongs.txt", Antxs)
savetxt("configLats.txt", Antys)
savetxt("configAlts.txt", alts)


altitudes = ','.join(map(str, alts))
Long = ','.join(map(str, Antxs))
Lat = ','.join(map(str, Antys))

#print nicely to copy over to eqArrOverTime.c
print('Longs')
print(Long)
print('Lats')
print(Lat)
print('Alts')
print(altitudes)

print('now generating map figures')

#print 3d alt of south pole and antenna locs
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import random

fig = pyplot.figure()
ax = fig.add_subplot(111)


buff = 5#20 #1000#128
buff2 = 5#20 #1000
jump=1



elevMap = Simage[int(min(lineYs))-buff:int(max(lineYs))+buff2, int(min(sampXs))-buff:int(max(sampXs))+buff2]

x = linspace(min(Antxs) - float(buff)/resPixperDeg, max(Antxs) + float(buff)/resPixperDeg, (int(max(sampXs)) - int(min(sampXs)) + 2*buff)/jump)  - midLong
y = linspace(max(Antys) + float(buff)/resPixperDeg, min(Antys) - float(buff)/resPixperDeg, (int(max(lineYs)) - int(min(lineYs)) + 2*buff)/jump) - midLat

x = x*mpdeg#*1e-3
y = y*mpdeg#*1e-3


[X,Y] = meshgrid(x,y)

cs=contourf(X, Y, elevMap*1000, levels = np.linspace(-795., -788, 64))#64)  #, vmin=-794, vmax=-789.5)


# m = cm.ScalarMappable(cmap=cm.rainbow, norm=mpl.colors.Normalize(vmin=-791, vmax=-789.5))
#
# m.set_array(elevMap*1000)
#
# m.set_clim(-791, -789.5)

# cbar=colorbar(m, boundaries=np.linspace(-791, -789.5, 64))

cbar=fig.colorbar(cs)#, boundaries=np.linspace(-791, -789.5, 64), norm=mpl.colors.Normalize(vmin=-794, vmax=-789.5))#, cmap=cm.rainbow, norm=mpl.colors.Normalize(vmin=-794, vmax=-789.5))#cmap=cm.rainbow, norm = cs.norm)
cbar.ax.set_ylabel('Elevation from mean Lunar Radius (m)')


ptxs = pts[:,0]
ptys = pts[:,1]

plot(ptxs, ptys, 'k-')

for i in range(numAnts):
    plot((Antxs[i] - midLong)*mpdeg, (Antys[i] - midLat)*mpdeg, 'r*', markersize=10)

for i in range(len(antEnds)):
    plot(antEnds[i][:,0],antEnds[i][:,1], 'r-')

plot([0.],[0.], 'k.', markersize=22) # lander



# #approx

xlabel(r'm from Longitude $-1.508^{\circ}$')
ylabel(r'm from Latitude $-0.125^{\circ}$')
xlim((-10, 55))
ylim((-55,10))

# xlim((-20, 120))
# ylim((-70,70))
title('Elevation Map (km) & Antenna Locations')

savefig('PRIME_Map2.png')

#
# clf()
#
# print('max(alts)*1e3 - min(alts)*1e3 is ' + str(max(alts)*1e3 - min(alts)*1e3))
#
# fig = pyplot.figure()
# ax = fig.add_subplot(111)
#
# elevMap = Simage[15360/2 - resPixperDeg*2:15360/2 + resPixperDeg*2, 46080/2 - resPixperDeg*2:46080/2 + resPixperDeg*2]  #shp = ((15360, 46080))
#
# x = linspace(-2, 2, resPixperDeg*4)
# y = linspace(2, -2, resPixperDeg*4)
# [X,Y] = meshgrid(x,y)
#
#
# #center on 0 longitude
# cs=contourf(X, Y, elevMap, 64) #append(elevMap[:, shp[1]/4:], elevMap[:, :shp[1]/4], axis = 1)
#
# cbar=fig.colorbar(cs)#cmap=cm.rainbow, norm = cs.norm)
# cbar.ax.set_ylabel('Elevation from mean Lunar Radius (km)')#, fontsize=25)
# #cbar.ax.set_tick_params(labelsize=40)
#
#
# plot([-1.508], [-0.125], 'r*', ms=10., label='Lowest rms 5x5 km patch')
#
# plot([-1.04], [-0.43], 'b*', ms=10., label='Lowest rms 10x10 km patch')
#
# plot([0.023], [1.625], 'k*', ms=10., label='Lowest rms 20x20 km patch')
#
# legend(loc=4)
#
# # ax.xaxis.set_tick_params(labelsize=20)
# # ax.yaxis.set_tick_params(labelsize=20)
#
# xlabel(r'Longitude $\circ$')#, fontsize=25)
# ylabel(r'Latitude $\circ$')#, fontsize=25)
# #xlim((-40, 20))
# #ylim((-20,40))
# title("Elevation Map (km) of Lunar Equator Region & Array Location")#, fontsize=25)
#
# savefig('MEDarrayMap0center.png')


#elevMap is 512x512

# rmss = zeros((490, 490))
#
# for i in range(11, 501):
#     for j in range(11, 501):
#         subArr = elevMap[i - 11:i+11, j-11:j+11].flatten()
#         rmss[i-11, j-11] = std(subArr)
#
# print 'min rms is ' + str(where(rmss == rmss.min()))
#min rms is (array([261]), array([52])) with offset 11, 11 into 4x4 degrees around 00
#midLong = -1.508
# midLat = -0.125
# 0.002817279089066556 rms around this 5.1x5.1 km point in km
#max altitude diff for 1024 array is 38.5 m across 10 km array
#max altitude diff for 1024 array is 18.3 m across 6 km array



# rmss = zeros((468, 468))
#
# for i in range(22, 490):
#     for j in range(22, 490):
#         subArr = elevMap[i - 22:i+22, j-22:j+22].flatten()
#         rmss[i-22, j-22] = std(subArr)
#
# print 'min rms is ' + str(where(rmss == rmss.min()))
#min rmss[289,101] with offset 22, 22 into 4x4 around 00
# midLong = -1.04
# midLat = -0.43
# 0.0055846117078073125 rms around this 10.2x10.2 km point in km
#max altitude diff for 1024 array is 43 m across 10 km array

#max(alts)*1e3 - min(alts)*1e3 = 110m for 20 km array, 28m for 6 km array




# rmss = zeros((424, 424))
#
# for i in range(44, 468):
#     for j in range(44, 468):
#         subArr = elevMap[i - 44:i+44, j-44:j+44].flatten()
#         rmss[i-44, j-44] = std(subArr)
#
# print 'min rms is ' + str(where(rmss == rmss.min()))
#min rms is (array([215]), array([4])) with offset 44, 44 into 4x4 deg around 00 long lat

# 0.013486054847264227 km rms around this 20.4x20.4 km point in km
#max altitude diff for 1024 array is 217.713177204 m across 20 km array



#
#make big plot with all -60 to 60 lat data in SLDEM file
# fig = pyplot.figure(figsize=(30, 10))
# ax = fig.add_subplot(111)
#
#
#
#
# elevMap = Simage[::2, ::2]
#
# x = linspace(-180, 180, shp[1]/2)
# y = linspace(60, -60, shp[0]/2)
# [X,Y] = meshgrid(x,y)
#
#
# #center on 0 longitude
# cs=contourf(X, Y, elevMap, 64) #append(elevMap[:, shp[1]/4:], elevMap[:, :shp[1]/4], axis = 1)
#
# cbar=fig.colorbar(cs)#cmap=cm.rainbow, norm = cs.norm)
# cbar.ax.set_ylabel('Elevation from mean Lunar Radius (km)', fontsize=25)
# #cbar.ax.set_tick_params(labelsize=40)
#
#
# plot([midLong], [midLat], 'k*', ms=30.)
#
# ax.xaxis.set_tick_params(labelsize=20)
# ax.yaxis.set_tick_params(labelsize=20)
#
# xlabel(r'Longitude $\circ$', fontsize=25)
# ylabel(r'Latitude $\circ$', fontsize=25)
# #xlim((-40, 20))
# #ylim((-20,40))
# title("Elevation Map (km) of Lunar Equator Region & Array Location", fontsize=25)
#
# savefig('BIGarrayMap0center.png')
