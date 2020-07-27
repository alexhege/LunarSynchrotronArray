# createArrayConfig.py
# Alex Hegedus 12/11/19
# Creates configuration for large array near equator of the Moon.
# Currently gives 10 km diameter Logarithmically spaced circular configuration
# also creates plots of the arrays
# saves and prints the longitudes, latitudes, and altitudes for all antenna (1024 currently)

from pylab import *
import numpy as np

import matplotlib.cm as cm


#returns CDF of a single petal arc length.  defined between t = 0-pi/4, returns 0-1
#expansion of above since it doesn't play nice defined between t = 0-pi/4, returns 0-1
#https://www.wolframalpha.com/input/?i=expand+y%3D+%28Sqrt%5B19+%2B+4+Cos%5B4+t%5D+-+15+Cos%5B8+t%5D%5D+Sec%5B2+t%5D+%28-%28Sqrt%5B30%5D+ArcTan%5BSin%5B2+t%5D%2FSqrt%5B-17%2F30+%2B+Cos%5B4+t%5D%2F2%5D%5D%29+%2B+15+Sqrt%5B-17+%2B+15+Cos%5B4+t%5D%5D+Sin%5B2+t%5D%29%29%2F%2860+Sqrt%5B2%5D+Sqrt%5B-17+%2B+15+Cos%5B4+t%5D%5D%29+
def cdfPetal(t):

    t = t + 0j


    cdf = 1./2.2664*(-1*(49*sin(t)*sqrt(4*cos(4*t) - 15*cos(8*t) + 19)*cos(t))/(4*sqrt(2.)*(cos(t)**2 - sin(t)**2)*(15*sin(t)**4 + 15*cos(t)**4 - 90*sin(t)**2 *cos(t)**2 - 17)) + (45*sin(t)**5*sqrt(4*cos(4*t) - 15*cos(8*t) + 19)*cos(t))/(4*sqrt(2.)*(cos(t)**2 - sin(t)**2)*(15*sin(t)**4 + 15*cos(t)**4 - 90*sin(t)**2* cos(t)**2 - 17)) + (45*sin(t)*sqrt(4*cos(4*t) - 15*cos(8*t) + 19)*cos(t)**5)/(4*sqrt(2.)*(cos(t)**2 - sin(t)**2)*(15*sin(t)**4 + 15*cos(t)**4 - 90*sin(t)**2 *cos(t)**2 - 17)) - (75*sin(t)**3* sqrt(4. *cos(4*t) - 15*cos(8*t) + 19)*cos(t)**3)/(2.*sqrt(2.)*(cos(t)**2 - sin(t)**2)*(15*sin(t)**4 + 15*cos(t)**4 - 90*sin(t)**2*cos(t)**2 - 17)) - (sqrt(15*cos(4*t) - 17)*sqrt(4*cos(4*t) - 15*cos(8*t) + 19)*np.arctan(sin(2*t)/sqrt(1/2.*cos(4*t) - 17./30)))/(4.*sqrt(15.)*(cos(t)**2 - sin(t)**2)*(15*sin(t)**4 + 15*cos(t)**4 - 90*sin(t)**2* cos(t)**2 - 17.)))

    assert cdf.imag < 0.005

    return abs(cdf)


#UPDATE this to where you have downloaded SLDEM data
lunarPath = '/mnt/LinuxData/Lunar_Pipeline/SLDEM2015_128_60S_60N_000_360_FLOAT.IMG'


# m per degree of longitude at lunar equator
mpdeg = 29670.597283903604

#diameter of array, or maximum separation in meters
maxSep = 10000.

#how to partition array, currently explicitly simulating 32*32 = 1024 antenna
numinArm = 16
numArms = 8

dishRadius = maxSep/2.  # in meters

print 'dish Radius is ' + str(dishRadius)

numAnts = numArms * numinArm

Antxs = zeros(numAnts)
Antys = zeros(numAnts)
alts = zeros(numAnts)



#use parametric equation for flower petal configuration
#config1
# for i in range(numAnts):
#     t = float(i)/numAnts*2*pi
#     Antxs[i] = (1 + cos(4*t))*cos(t) * dishRadius/2.
#     Antys[i] = (1 + cos(4*t))*sin(t) * dishRadius/2.


#parametric euation starts at end of petal, where separations are largest, near 300 m, I upped to 400 to fit
dr = ones(numinArm)*400

fh = logspace(log10(200), log10(400), numinArm/2)

dr[:numinArm/2] = fh

#bigger jumps first, since parametric equation starts at end of petal
dr = flipud(dr)

#makes a symmetic 400 distance between 2 at end of petal on opposite sides
dr[0] = 150. #200.

#length of curve, should do integral to figure out length, guess for now, update done, arc length is 2.2664 from below link
# https://www.wolframalpha.com/input/?i=integral+from+0+to+pi%2F4+of+sqrt%28%28sin%28t%29+%28-%281+%2B+cos%284+t%29%29%29+-+4+sin%284+t%29+cos%28t%29%29%5E2+%2B+%28cos%28t%29+%281+%2B+cos%284+t%29%29+-+4+sin%28t%29+sin%284+t%29%29%5E2%29
totDist = dishRadius/2.*2.2664

#sum of distance so far
dt = 0.

#angular progress so far
progressAngle = 0.
progressAngleY = 0.

angStep = .001

for i in range(numinArm):
    # t = float(i)/numinArm*2*pi/8. #only first arm
    dt += dr[i]
    progressFraction = dt/totDist
    cdfFrac = cdfPetal(progressAngle)
    #asymmetric counter
    cdfFracY = cdfPetal(progressAngleY)
    #march through CDF of arc length so far until it matches what we want
    while(cdfFrac < progressFraction):
        progressAngle += angStep
        cdfFrac = cdfPetal(progressAngle)

        #asymmetric
        if progressAngle < pi/8:
            progressAngleY += angStep#*.95#1.05 #.0015 #half time for first half

        else:
            progressAngleY += angStep#*1.05#.95 #.0005 #1.5 time for second half

    # t = dt/totDist*2*pi/8. # assumes phase linearly increases with distance travelled, not true.  CDF indexing?

    t = progressAngle

    Antxs[i] = (1 + cos(4*t))*cos(t) * dishRadius/2.
    Antys[i] = (1 + cos(4*t))*sin(t) * dishRadius/2.

    Antys[i] = (1 + cos(4*progressAngleY))*sin(progressAngleY) * dishRadius/2. * 0.8



#finish other half of first petal
Antxs[numinArm:numinArm*2] = Antxs[0:numinArm]
Antys[numinArm:numinArm*2] = Antys[0:numinArm]*-1.

#asymmetric other half of first petal
Antys[numinArm:numinArm*2] = flipud(Antys[0:numinArm]*-1.)

#mirror this petal on other side
Antxs[numinArm*2:numinArm*4] = Antxs[0:numinArm*2]*-1
Antys[numinArm*2:numinArm*4] = Antys[0:numinArm*2]*-1

#mirror these petals on other axis
Antxs[numinArm*4:numinArm*8] = Antys[0:numinArm*4]
Antys[numinArm*4:numinArm*8] = Antxs[0:numinArm*4]


normsMeters = np.sqrt(Antxs**2 + Antys**2)

Antxsm = Antxs.copy()
Antysm = Antys.copy()


savetxt("configLongsMeters.txt", Antxs)
savetxt("configLatsMeters.txt", Antys)



# normalize to Lunar Longitude and Latitude, mpdeg = 29670.6 m/degree equatorial Longitude
Antxs /= mpdeg
Antys /= mpdeg

lunarRad = 1737.4*1e3 # meters


#min rms is (array([261]), array([52])) with offset 11, 11 into 4x4 degrees around 0, 0 long lat
# 0.0028 rms around this 5x5 km point in km
# midLong = -2. + 63/128.
# midLat = 2 - 272./128
#
# #approx
# midLong = -1.508
# midLat = -0.125


# min rms in 10x 10 km
# rmss[289,101] with offset 22, 22
midLong = -2. + 123/128.
midLat = 2 - 311./128
#
# #approx
midLong = -1.04
midLat = -0.43


#min rms 20x20 km patch location
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





#Import lunar SLDEM data

fid = open(lunarPath, 'rb')
data = np.fromfile(fid, dtype=np.dtype(np.float32))
fid.close()


shp = ((15360, 46080))

scale = 0.236901 #km per pix
Simage = data.reshape(shp)

#reshape so 0 Long is in the middle

Simage = append(Simage[:, shp[1]/2:], Simage[:, :shp[1]/2], axis = 1)


#use DSMAP.CAT to change to indices in file to get altitudes (in file as km)

centerLong = 0.
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
savetxt("configAltsMeters.txt", alts)


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

#buffer of space beyond array to include in frame
buff = 5
buff2 = 5
jump=1



elevMap = Simage[int(min(lineYs))-buff:int(max(lineYs))+buff2, int(min(sampXs))-buff:int(max(sampXs))+buff2]

x = linspace(min(Antxs) - float(buff)/resPixperDeg, max(Antxs) + float(buff)/resPixperDeg, (int(max(sampXs)) - int(min(sampXs)) + 2*buff)/jump)  - midLong
y = linspace(max(Antys) + float(buff)/resPixperDeg, min(Antys) - float(buff)/resPixperDeg, (int(max(lineYs)) - int(min(lineYs)) + 2*buff)/jump) - midLat

x = x*mpdeg*1e-3
y = y*mpdeg*1e-3


[X,Y] = meshgrid(x,y)

cs=contourf(X, Y, elevMap, 64)


cbar=fig.colorbar(cs)
cbar.ax.set_ylabel('Elevation from mean Lunar Radius (km)')


for i in range(numAnts):
    plot((Antxs[i] - midLong)*mpdeg*1e-3, (Antys[i] - midLat)*mpdeg*1e-3, 'k*')

#change for center of array, currently at minimum rms elevation in 10x10 km square
xlabel(r'km from Longitude $-1.04^{\circ}$')
ylabel(r'km from Latitude $-0.43^{\circ}$')
title('Elevation Map (km) & Antenna Locations')

savefig('arrayMap1024-10.png')


clf()

print('max(alts)*1e3 - min(alts)*1e3 is ' + str(max(alts)*1e3 - min(alts)*1e3))

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
# cs=contourf(X, Y, elevMap, 64)
#
#
# #Some options for enlarging the axes labels are commented out
# cbar=fig.colorbar(cs)
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
#
# title("Elevation Map (km) of Lunar Equator Region & Array Location")#, fontsize=25)
#
# savefig('MEDarrayMap0center.png')


# #make big plot with all -60 to 60 lat data in SLDEM file
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
# cs=contourf(X, Y, elevMap, 64)
#
# cbar=fig.colorbar(cs)
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



#code used to calculate rms values of elevation variance on the surface

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
# 0.002817 rms around this 5.1x5.1 km point in km
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
# 0.005584 rms around this 10.2x10.2 km point in km
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

# 0.0134 km rms around this 20.4x20.4 km point in km
#max altitude diff for 1024 array is 217.713 m across 20 km array
