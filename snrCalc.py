# -*- coding: cp1252 -*-
# snrCalc.py
# Alex Hegedus 12/11/19
#snr calculations for lunar near side array for earth's radiation belts
# creates noise budget figures from Hegedus et al.


from pylab import *
from scipy.interpolate import UnivariateSpline

#MHz
freqs = array([0.01, 0.0136, 0.0185, 0.0251, 0.0341, 0.0464, 0.0631, 0.0858, 0.1166, 0.1585, 0.2154, 0.2929, 0.3981, 0.5412, 0.7356, 1.0])

#Nov 1st 2016 Simulated total integrated brightness
jys = array([2.87221036, 3.08760165, 3.28986448, 3.47002421, 3.61532536, 3.71087727, 3.74478608, 3.70769761, 3.59555882, 3.41312112, 3.17395549, 2.89424464, 2.59118944, 2.2806108 , 1.97589993, 1.68965141])



spl = UnivariateSpline(freqs*1e6, jys, k=3, s=0.)

#total brightness from freq[index] to 1.0 MHz.  We assume below some cutoff we won't see usable signal
cumBrightnessRight = zeros(len(freqs))

#from 0.01 to freq[ind]
cumBrightnessLeft = zeros(len(freqs))

for i in range(len(freqs)):
    cumBrightnessRight[i] = spl.integral(1e6*freqs[i], 1e6*freqs[-1])
    cumBrightnessLeft[i] = spl.integral( 1e6*freqs[0], 1e6*freqs[i])

fig, ax = subplots()

plot(freqs, jys, 'o-')
xlabel('Frequency (MHz)')
ylabel('Total Flux Density (Jy)')

# ylim((6, 16))
xscale('log')

title('Integrated Flux Density of Earth Radiation Belts from Lunar Distances')

savefig('freqbright.png')

plot(linspace(0.01, 1.0, 100), spl(linspace(0.01, 1.0, 100)*1e6))
xlabel('Frequency (MHz)')
ylabel('Spline Total Brightness/Hz (Jy)')

title('Spline Integrated Brightness of Earth Radiation Belts from Lunar Distances')

savefig('freqbrightspline.png')

fig, ax = subplots()
plot(freqs, cumBrightnessRight)
xlabel('Frequency (MHz)')
ylabel('Total Brightness from Freq - 1.0 MHz (Jy)')

title('Integrated Brightness from Freq - 1.0 MHz\n of Earth Radiation Belts from Lunar Distances')

# savefig('freqbrightright.png')


fig, ax = subplots()
plot(freqs, cumBrightnessLeft)
xlabel('Frequency (MHz)')
ylabel('Total Brightness from 0.01 MHz - Freq (Jy)')

title('Integrated Brightness from 0.01 MHz - Freq\n of Earth Radiation Belts from Lunar Distances')

# savefig('freqbrightleft.png')



#galactic noise model novacco and brown 1978
B0 = 1.38 * 10**-19 # W/m2/Hz/sr

taus = 3.28 * linspace(0.01, 1.0, 100)**-.64


Bmodelgalactic = B0* linspace(0.01, 1.0, 100)**-.76 * exp(-1*taus)



fig, ax = subplots()
plot(linspace(0.01, 1.0, 100)*1e3, Bmodelgalactic, label='Galactic Noise')
xlabel('Frequency (kHz)')
ylabel('Galctic Brightness (W/m2/Hz/sr)')

xscale('log')
yscale('log')

ylim((1e-22, 1e-19))

xlim((1e2, 1e3))

title('Galactic Brightness')

# savefig('galacticBright.png')

#quasi thermal noise from plasma
#Mhz freq, w/m2/Hz/sr brightness y f2/f1 / x2/x1
y1 = 1e-21 #10e-17
y0 = 3e-20 #2e-16

x1 = .500
x0 = .100


m = log(y1/y0)/log(x1/x0)

a = y0/(x0**m)

#F = F0*(x/x0)**m

qtnPlas = a*(linspace(0.01, 1.0, 100))**m

#physics

#blackbody calculations
h = 6.626e-34  #plank constant
kb = 1.38e-23 #boltzmann constant
c = 3e8 # speed of light

freqsBlack = linspace(0.01, 1.0, 100)*1e6
temp = 288. #288 kelvin earth, may be 288, more like 255 actually


blackBodyEarth = 2*h*freqsBlack**3/c**2*1/(exp(h*freqsBlack/kb/temp) - 1) # returns W/m2/sr/Hz

#scale with https://www.acs.org/content/acs/en/climatescience/energybalance/energyfromsun.html as guide

#Rs = 695510 #km
#astrounit = 1.496e+8 #km



totSigBlackMoon = blackBodyEarth*(6371./384400.)**2 #equivalent brightness at moon distnaces[[]], still W/m2/sr/Hz


totSigBlackMoon *= 4*pi # you are a receiver over the new r2 scaled ground, taking up half the sky, the other half is underground
#now in W/m2/Hz, like jansky e26


#diameter of earth in pic is 50 pixels, pix size is .038 degrees
perpixval = totSigBlackMoon/pi/(25)**2


# totSigBlackMoon = blackBodyEarth* (.038 * 6731000*2/1.9.)**2  *(6371./384400.)**2 *2*pi#deg squared times meters per degree #scaled to lunar distances



#for solar wind
kPereV = 1.160452e4

neqtn = 5. #/cc
teqtn = 1e5 #Kelvin  1.160 452 21 x 104 K per electonVolt 10 eV typical dayside

teqtn = 12.*kPereV
neqtn = 8.

#from lunar literature
# neqtn = 250. #dayside lunar surface enhancement less than 60 deg from zenith Imamura et al 2012 SELENE radio occultation
# neqtn = 1000.  # from Luna 19 22 radio refraction
# neqtn = 8.
#On the lunar night side ne shows a range of 2–0.002 cm^3
#and Te has a range of 15–50 eV 58eV max. Chandran et al 2013 Plasma electron temper variability
#On the day side, ne >= 8cm^3 and Te >= 12 eV

#SEFD is brightness average*4pi in W/m^2/Hz, 1 Jy = 1e-26 W/m^2/Hz

#SEFD con 1.46e-18  1000/cc density Luna measurements
#SEFD mod 4.6e-19   250/cc density  SZA < 60 deg SELENE noise
#SEFD opt 1.38e-19  8/cc density (amplifier dom) night side



fig, ax = subplots()


for neqtn in [8, 250, 1000]:


    gamma = .5 #gamma**2 = 0.5 = -3 dB in SunRISE.  In Zav, seems to be gamma**1= 0.5
    Ldipole = 5.#meters on 1 side of dipole
    qtnPlasFreqs =  linspace(0.01, 1.0, 100)*1e6
    lambdas = (3e8/qtnPlasFreqs)
    qtnPlasEqnv2 = 5e-5*neqtn*teqtn/Ldipole/(qtnPlasFreqs)**3  #*gamma #?
    z0 = 120*pi
    Rrs = 2*pi/3.*z0*(5./lambdas)**2
    qtnPlasEqn = qtnPlasEqnv2/2./Rrs/(lambdas)**2

    plot(linspace(0.01, 1.0, 100)*1e3, qtnPlasEqn, label='Quasithermal Plasma Noise ' + str(neqtn)+ '/cc')

ylabel('(Equivalent) Brightness (W/m2/Hz/sr)')
xlabel('Frequency (kHz)')


xscale('log')
yscale('log')



xlim((1e2, 1e3))

ylim((1e-22, 1e-16))

ampNoisev2gamma = ones(100)*1e-16

ampNoisebright = ampNoisev2gamma/2./Rrs/(lambdas)**2/gamma**2

ampNoise = ones(100)*5e-20  #old 4e-21. but sunRISE more like 5e-20

ampNoise = ampNoisebright

plot(linspace(0.01, 1.0, 100)*1e3, ampNoise, label='Amplifier Noise')

plot(linspace(0.01, 1.0, 100)*1e3, Bmodelgalactic, label='Galactic Noise')


# tots = qtnPlasEqn + ampNoise # + Bmodelgalactic #assumed able to take out
#
# plot(linspace(0.01, 1.0, 100)*1e3, tots, label = 'Total Noise')

# plot(linspace(0.01, 1.0, 100)*1e3, totssqrt, label = 'Total Add Noise')

legend()

title('Constant Noise Terms on Lunar Surface', fontsize='16')




# savefig('galacticPlasBright.png')
# savefig('galacticPlasBrightEqnOpt.png')
savefig('galacticPlasBrightEqnAll.png')

# noisespl = UnivariateSpline(linspace(0.01, 1.0, 100)*1e6, tots, k=3, s=0.)
#
# cumNoiseRight = zeros(len(freqs))
#
# #from 0.01 to freq[ind]
# cumNoiseLeft = zeros(len(freqs))
#
# for i in range(len(freqs)):
#     cumNoiseRight[i] = noisespl.integral(1e6*freqs[i], 1e6*freqs[-1])
#     cumNoiseLeft[i] = noisespl.integral( 1e6*freqs[0], 1e6*freqs[i])
#
#
# avBrightness = 10 * 1e-26 # jy
# # avNoise = noisespl.integral(1e6*0.1, 1e6*1)/(0.9*1e6)
#
# avNoise = noisespl.integral(1e6*0.5, 1e6*1)/(0.5*1e6)
# print('avNoise con is ' + str(avNoise))
#
# # #conservative
# # avNoise = 1e-20
#
#
# avSignal = spl.integral(1e6*.5, 1e6*1.0)/500000. * 1e-26 # jy
#
#
# N = 100
# tau = 360
#
# targetsnr = 5.0
#
# dv = 5e5
#
# snr = avBrightness/(avNoise)*sqrt(2*N*(N-1)*tau*dv)
#
# targetTau = (targetsnr*avNoise/avSignal)**(2)/2/N/(N-1)/dv
#
# targetTaus = zeros(100)
#
# targetTaus[0] = (targetsnr*avNoise/avSignal)**(2)/2/dv
#
# for n in range(2, 101):
#
#     targetTaus[n-1] = (targetsnr*avNoise/avSignal)**(2)/2/n/(n-1)/dv
#
#
# targetsnr = 3.0
#
# targetTaus3 = zeros(100)
#
# targetTaus3[0] = (targetsnr*avNoise/avSignal)**(2)/2/dv
#
# for n in range(2, 101):
#
#     targetTaus3[n-1] = (targetsnr*avNoise/avSignal)**(2)/2/n/(n-1)/dv
#
# clf()
# plot(range(20, 101), targetTaus[19:], label = 'SNR=5')
# plot(range(20, 101), targetTaus3[19:], label = 'SNR=3')
# legend()
#
#
#
# xlabel('Number of Antenna (#N)')
# ylabel('Seconds of Integration Time to Target SNR')
#
# title('Integration Time for 3 & 5 SNR Detection of Radiation Belts')
#
# savefig('snrTime.png')
