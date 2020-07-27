#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
################################################################################

'''
LunarEarthPicFreqIntegration.py
Alex Hegedus 4/9/18
2/6/19
12/11/19
alexhege@umich.edu
This simulation code has been put together to simulate imaging Earth Radiation belts
Get positions with SPICE code in eqArrOverTimeEarth.c

Reads in the .dat files of the synchrotron emission and converts them into a CASA ground
truth file that is given to the virtual array.  Currently, only the ~0.75 MHz image is used as an average for the 0.5 – 1.0 MHz range.

It uses CASA routines to create a synthetic array using the J2000 data calculated
from the SPICE kernel in the previous script.  It then feeds the truth image into the
simulated array to create the visibility data and puts it into a Measurement Set (MS) file.
This MS file is the main output of this script.


instrution on how to run the program
within casa run this command:

%run LunarEarthPicFreqIntegration.py  -outDir . -correlate True -numSC 128

or in terminal
casa --nologger --nologfile --nogui --agg -c LunarEarthPicFreqIntegration.py  -outDir . -correlate True -numSC 1024| tee out.out
nohup casa --nologger --nologfile --nogui --agg -c LunarEarthPicFreqIntegration.py  -outDir . -correlate True -numSC 128 | tee mags.out &

'''
#frame of reference is J2000/eme2000
#Unless you are doing milliarcsecond astronomy, you can ignore that bias and rotation.
#For most applications, J2000/FK5=EME2000=ICRF=ICRF2.

#import the required python modules
import time
import os,sys
import shutil
from collections import defaultdict
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import random
import math
import scipy.interpolate as interpolate

### import casa modules
from rmtables import rmtables
from importvla import importvla
from importuvfits import importuvfits
from listobs import listobs
from applycal import applycal
from setjy import setjy
from gaincal import gaincal
from fluxscale import fluxscale
from accum import accum
#from plotxy import plotxy
from plotcal import plotcal
from clean import clean
from split import split
from clearcal import clearcal
from imstat import imstat
from imregrid import imregrid
from exportfits import exportfits
import casac
from taskinit import *
import numpy as np
import random
import math
from casa import imhead
from casa import *
from casac import *
from concat import concat
import signal
import time

import matplotlib.image as plimg
import scipy.ndimage.interpolation as spndint

from pylab import *

import ephem
import datetime

################################################################################

##define input options for running the code

def _getParser():
    argParser = argparse.ArgumentParser(description= __doc__,
                                        formatter_class=argparse. \
                                        RawDescriptionHelpFormatter )

    #
    # argParser.add_argument( '-componentFile', required = True, type=str,\
    #                           help="text file which contains the info to build the truth image.")
    #

    argParser.add_argument( '-correlate', required = False,
                            action='store_true', \
                            help = """ #True: simulate manual correlation to form visibilities (slow)
                            # False: use CASA sm.predict to form visibilities (much faster)""" )

    argParser.add_argument( '-thermalNoise', required = False,
                            action='store_true', \
                            help = """ turns on galactic thermal noise """ )

    argParser.add_argument( '-phaseNoise', required = False,
                            action='store_true', \
                            help = """ turns on positional and clock bias noise """ )

    argParser.add_argument( '-numSC', required = False,
                            default = 6, type = int,\
                            help = """ randomly removes sc if less than 6 """ )


    argParser.add_argument( '-startPos', required = False,
                            default = 0, type = int,\
                            help = """ initial index into gps pos files """ )

    argParser.add_argument( '-angPos', required = False,
                            default = -1, type = int,\
                            help = """ index into phase of propogation angle """ )

    argParser.add_argument( '--nologger', required = False,
                            action='store_true', \
                            help = """ For casa no log window""" )

    argParser.add_argument( '--nologfile', required = False,
                            action='store_true', \
                            help = """ For casa no log window""" )

    argParser.add_argument( '--nogui', required = False,
                            action='store_true', \
                            help = """ For casa no log window""" )

    argParser.add_argument( '-c', required = False,
                            action='store_true', \
                            help = """ For casa run as script""" )

    argParser.add_argument( '--agg', required = False,
                            action='store_true', \
                            help = """ For casa run as script""" )


    args = argParser.parse_known_args()

    return args[0]



###############################################################################

################################################################################

def arraySim(lunarMS,freqMHz,truthImage, img, Ant_x, Ant_y, Ant_z, args, baselines, dishSize=10, timestep = 0):
    '''
    This function will simulate the observations and creates an MS file
    '''
    print 'Now in arraySim '
    correlate = args.correlate
    phaseNoise = args.phaseNoise
    thermalNoise = args.thermalNoise

    numSC = len(Ant_x)
    numbl = numSC*(numSC-1)/2

    #one time things to compute
    dishDiameter = np.ones(numSC)*dishSize

    mounts = []
    antnames = []

    for i in range(numSC):
        mounts.append('X-Y')
        antnames.append('S'+str(i))

    sm.open(lunarMS)
    refloc=me.observatory('VLA')
    #set the antenna configuration
    sm.setconfig(telescopename='VLA',
                 x=Ant_x,y=Ant_y,z=Ant_z,
                 dishdiameter=dishDiameter,
                 mount=mounts,
                 antname=antnames,
                 coordsystem='local',referencelocation=refloc)



    sm.setspwindow(spwname='HF',
                   freq=str(freqMHz)+'MHz',
                   deltafreq='6100Hz', freqresolution='6100Hz', nchannels=1)
    # Set the feeds to be X Y.
    print 'Setting the feed polarization ...'
    sm.setfeed('perfect X Y')


    print 'Setting the source ...'

    srcdir = CMEDir
    sm.setfield(sourcename='CME1', sourcedirection=srcdir)


    print 'Setting integration times ...'
    rday1 = 60017.83
    rday1 += 0.1*timestep/60./60./24. #add seconds to day
    newday=me.epoch(rf='UTC',v0=str(rday1)+'d')
    #integration time doesn't factor meaningfully since we calculate
    #different noise levels in a later script noiseCopies.py
    sm.settimes(integrationtime='.066s',usehourangle=False,referencetime=newday)

    sm.setauto(0.0);
    print 'Observing ...'
    sm.observe(sourcename='CME1',
               spwname='HF',
               observer='SunRISE',
               starttime='0s', stoptime='.066s')
    # Fourier transform model image into data set

    freq = freqMHz*1e6

    au = 1.496e11 #AU in m


    #insert new baselines into MS file here

    baselines2 = np.reshape(baselines[:,timestep,:], (numbl,3), order='F').T

    # print 'Ant x: ' + str(Ant_x)
    # print 'Ant y: ' + str(Ant_y)
    # print 'Ant z: ' + str(Ant_z)
    # print 'baselines timestep inp to F order ' + str(baselines[:,timestep,:])
    # print 'output baselines ' + str(baselines2)

    ia.open(truthImage)
    cs=ia.coordsys()
    print(truthImage)

    pix = ia.getchunk()
    imSize = np.shape(pix)[0]
    pix = pix.reshape((imSize, imSize))
    img = pix

    # print 'img ' + str(img)

    # print 'max pos im ' + str(np.where(img == np.amax(img)))

    tb.open(lunarMS, nomodify=False)
    actual_uvw = tb.getcol("UVW")
    tb.putcol("UVW", baselines2)


    #insert man made correlation visibilities here, from time series code

    print args.correlate
    # correlate = False


    if correlate:


        print 'now correlating'
        res1=400
        res=400
        width = res1*.038 * np.pi/180.
        dres = width/res
        cell_rad=dres
        #qa.convert(qa.quantity("1arcmin"),"rad")['value'] # for 1 arcmin
        imSize = res
        print 'res ra dec ', res, ra, dec
        #real world sky width
        #assuming eash point is armin squared
        c = 3e8
        kwavnum = 2*np.pi*freq
        wavelen = c/freq
        kb = 1.380648e-23
        #calculate width of each pixel
        dres = width/res

        #make cs here
        ia.fromshape('temp.truth',shape=[res,res,1,1],overwrite=True)
        cs=ia.coordsys()
        cs.setunits(['rad','rad','','Hz'])
        arcminDiv = 2.28 * 400/float(400) #from resolution of truth images
        cell_rad=qa.convert(qa.quantity(str(arcminDiv)+"arcmin"),"rad")['value']
        cs.setincrement([-cell_rad,cell_rad],'direction')
        cs.setreferencevalue([CMEDir['m0']['value'], CMEDir['m1']['value']],
                              type="direction")
        cs.setrestfrequency(freqMHz)
        cs.setreferencevalue(str(freqMHz)+'MHz', 'spectral')

        ####

        #https://casaguides.nrao.edu/index.php/N891_simdata_(CASA_3.4)


        #which s/c is farthest down along pz?
        minZ = np.where(phasescpos[:,timestep,2] == min(phasescpos[:,timestep,2]))[0][0]

        #get delays in seconds, from diff in z from minZ
        delays = np.zeros(numSC)
        for i in range(numSC):
            delays[i] = phasescpos[i][timestep][2] - phasescpos[minZ][timestep][2]
            delays[i] = delays[i]/c

        #change later for more frequencies
        fourierCoeff = np.zeros(numSC, dtype=np.complex128)

        visibilities = np.zeros(numbl, dtype=np.complex128)

        print np.shape(img)
        d1, d2 = np.shape(img)

        #loop over other pixels, adding in apporpriate delays for each sc fourierCoeff
        #like calculation of phase center, but t in front of variables for temporary, don't keep
        for i in range(d1):
            for j in range(d2):
                if j == 0:
                    # print('now on corr pixel ' + str(d1*i + j))
                    pass
                if img[i][j] == 0.0:
                    continue


                #try casa way

                s = cs.toworld([i,j,0,0],'s')['string']

                tra = qa.convert(qa.toangle(s[0]), 'rad')['value']
                tdec = qa.convert(qa.toangle(s[1]), 'rad')['value']

                # print 'pixel at '
                # print tra, tdec

                #print tra, tdec, i, j

                #unit direction pointing in direction of source/phase center
                tx = np.cos(tdec) * np.cos(tra)
                ty = np.cos(tdec) * np.sin(tra)
                tz = np.sin(tdec)


                ts = np.array([tx, ty, tz])


                # #get x and y axes perp to s

                tx = np.cross(ts, [0, 0, 1]) #, ts)
                ty = np.cross(ts, tx)

                #normalize since cross product isn't already 1 for some reason
                # this is
                tx = tx/np.linalg.norm(tx)
                ty = ty/np.linalg.norm(ty)
                tz = ts/np.linalg.norm(ts)

                tphasescpos = np.zeros((numSC,3))

                for k in range(numSC):
                    tphasescpos[k][0] = np.dot(allscpos[k][timestep], tx)
                    tphasescpos[k][1] = np.dot(allscpos[k][timestep], ty)
                    tphasescpos[k][2] = np.dot(allscpos[k][timestep], tz)

                #which s/c is farthest down along pz?
                tminZ = np.where(tphasescpos[:,2] == min(tphasescpos[:,2]))[0][0]

                tdelays = np.zeros(numSC)
                for k in range(numSC):
                    tdelays[k] = tphasescpos[k][2] - tphasescpos[tminZ][2]
                    tdelays[k] = tdelays[k]/c


                #phase of furthest down z in this projection, delay from this one
                phase = np.exp(1j*2*np.pi*np.random.random())
                C = img[i][j]*phase#*prefact
                for k in range(numSC):
                    fourierCoeff[k] = C*np.exp(-1j*kwavnum*(tdelays[k] - delays[k]))#*(1+(delta/360.*2*np.pi)))

                k=0
                for m in range(numSC):
                  for n in range(m+1, numSC):
                    toAdd = fourierCoeff[n]*fourierCoeff[m].conj()
                    toAdd = toAdd/np.sqrt(abs(toAdd))
                    visibilities[k] += toAdd
                    k+=1



        ############## do baselines & cross correlation, multiply fourierCoeff


        ia.close()



    print("iFence['model'",freqMHz)

    if not correlate:
        print 'Predicting the visibilities ...'

        sm.predict(imagename=truthImage)
    #sm.predict(truthImage)

    qs = zeros(numbl, dtype=np.complex128)
    dqs = zeros(numbl, dtype=np.complex128)

    pure = zeros(numbl, dtype=np.complex128)

    if correlate:

        print 'Predicting the visibilities ...'

        sm.predict(imagename=truthImage)

        #prepare plots

        data=tb.getcol("DATA")
        cdata=tb.getcol("CORRECTED_DATA")


        print("predict vis vals")

        print data[0][0][:]

        absq = zeros(numbl)
        angq = zeros(numbl, dtype=np.complex128)



        ##

        print "Now replacing CASA data with self made correlated data, comparing it"
        for i in range(numbl):
            #"check casa uvw and calculated"
            #print actual_uvw[:,i], baselines[:,i]

            q = data[0][0][i]/visibilities[i]

            pure[i] = data[0][0][i]

            qs[i]=  q
            absq[i] = abs(q)
            angq[i] = angle(q) #arcsin(q.imag)

            #print data[0][0][i], visibilities[i], q, abs(q), q/abs(q)
            #print data[1][0][i], visibilities[i]

            data[0][0][i] = visibilities[i]
            data[1][0][i] = visibilities[i]

            cdata[0][0][i] = visibilities[i]
            cdata[1][0][i] = visibilities[i]


        tb.putcol("DATA", data)
        tb.putcol("CORRECTED_DATA", cdata)



    sm.close()


    return


################################################################################

def CMEConcatModel(inputMS, outlunarMS):
    '''
    this function concats nearby frequency channels
    inputMS: a list of CME files which will be concatenated
    outCME: the output file
    '''
    concat(vis=inputMS,
           concatvis=outlunarMS,
           freqtol='',
           dirtol='',
           respectname=False,
           timesort=False,
           copypointing=True,
           visweightscale=[])



    return



#adds noise to CORRECTED_DATA column of MS, with SEFD given or calculated, dt per visibility (seconds), dv(taken from MS)
#calculates average galactic brightness, plasma quasithermal noise, and amplifier noise, add for SEFD
def addThermalNoise(MS, SEFD=-1., dv = 6100., dt = 60.):
    '''
    addThermalNoise(MS, dt=60., SEFD=-1.), estimates proper low frequency noise components and corrupts the
    CORRECTED_DATA column of MS. Frequency bandwidth and integration time are taken directly from the MS file.
    It is currently assumed that dt and dv is identical for each
    measurement.  When left to -1., System Equivalent Flux Density SEFD is calculated using Novacco & Brown 1973
    Galactic brightness models, a long wavelength approximation (Meyer-Vernet 2000) for quasithermal electron noise for typical
    solar wind density and temperature values, and a rough estimate of amplifier noise.  SEFD is in units of Jy,
    or 10^-26W/m2/Hz, and may be compared directly to the 'signal' data.  This is used with the integration time
    and the frequency bandwidth (taken automatically from the MS) to calculate the noise added to the visibilities.

    NOTE: SunRISE currently has an integration fraction of 6.6 ms / second, and this is hard coded in.
    '''

    ##calculate SEFD
    ms.open(MS)
    spws = ms.getspectralwindowinfo()
    sumData = ms.getscansummary()
    ms.close()

    tb.open(MS, nomodify=False)
    specInd = tb.getcol("DATA_DESC_ID")
    tb.close()

    hz = spws[str(specInd[0])]['Chan1Freq']

    # dv = spws[str(specInd[0])]['ChanWidth']
    #
    # dt = sumData['1']['0']['IntegrationTime']


    freqMHz = hz/(1.e6)

    #calculate frequency dependent brightness from novacco model
    B0 = 1.38 * 10**-19 # W/m2/Hz/sr
    taus = 3.28 * freqMHz**-.64
    Bmodelgalactic = B0* freqMHz**-.76 * exp(-1*taus)
    Bmodelgalactic = Bmodelgalactic*1.e26 #total jansky/sr, 1.42e7 Jy at 10MHz

    print('Bmodelgalactic ' + str(Bmodelgalactic) + ' Jy/sr')

    #solar wind temp and density
    kPereV = 1.160452e4
    teqtn = 12.*kPereV
    neqtn = 8.

    #brightness from quasithermal noise
    gamma = .5 #gamma**2 = 0.5 = -3 dB in SunRISE.  In Zav, seems to be gamma**1= 0.5
    Ldipole = 50. #2.5#meters on 1 side of dipole
    qtnPlasFreqs =  freqMHz*1e6
    lambdas = (3e8/qtnPlasFreqs)
    qtnPlasEqnv2 = 5e-5*neqtn*teqtn/Ldipole/(qtnPlasFreqs)**3  #*gamma #?
    z0 = 120*pi
    Rrs = 2*pi/3.*z0*(5./lambdas)**2
    qtnPlasEqn = qtnPlasEqnv2/2./Rrs/(lambdas)**2

    print('qtnPlasEqn ' + str(qtnPlasEqn)+ ' Jy/sr')

    ampNoise = Bmodelgalactic*0.1 # from farside, above .2 MHz, sky noise dominates by 10x or more , before for SunRISE 5e-21*1e26

    print('ampNoise ' + str(ampNoise)+ ' Jy/sr')

    SEFDCalc = 4*pi*(ampNoise + qtnPlasEqn + Bmodelgalactic*0.5) #total noise in Jy

    #if SEFD unset, use calculated SEFD, else, use provided SEFD
    if SEFD == -1.:
        SEFD = SEFDCalc

    print('SEFD is ' + str(SEFD) + ' Jy')

    print('dv is ' + str(dv) + ' Hz')
    print('dt is ' + str(dt) + ' seconds')

    tb.open(MS, nomodify=False)

    data=tb.getcol("DATA")
    cdata=tb.getcol("CORRECTED_DATA")

    numbl = shape(data)[2]

    #0.9 correlator efficiency, 2 antenna
    noise = SEFD/0.9/sqrt(2*(dt)*dv) # 6.6 ms / second integration for SunRISE

    print('Noise rms per polarization added is ' + str(noise) + ' Jy')

    for k in range(numbl):
        toAdd = np.random.normal(0, noise)*1.j + np.random.normal(0, noise)

        #leave noise free data in data column, clean checks/uses corrected data first
        data[0][0][k] += toAdd
        data[1][0][k] += toAdd


        cdata[0][0][k] += toAdd
        cdata[1][0][k] += toAdd

    tb.putcol("DATA", data)
    tb.putcol("CORRECTED_DATA", cdata)

    tb.close()

    return


#create truth image from component file describing gaussian position off phase center
def CMEform(comp, truthName, CMEDir, fMHz=-1.):
    '''
    Form the truth image, We assume that the CME moves outward in frequency.

    comp is a list of doubles with structure = (MHz, delta, offangle, major, minor, PA, flux, resolution(arcmin), imSize)

    if fMHz is defined, it takes the place of comp[0]

    truthName (string) is the basename of the .truth created.  CMEDir is the center of the image.

    Puts a single gaussian source on a truth model to be fed through the simulation pipeline.

    '''
    ##
    if fMHz == -1.:
        freq0 = comp[0]
    else:
        freq0 = fMHz
    delta = comp[1]
    delta1 = qa.toangle(str(delta)+'deg')
    offangle = comp[2]
    # offangle1 = initoffangle #qa.toangle(str(offangle)+'deg') #for random
    offangle1 = qa.toangle(str(offangle)+'deg')
    NewDir=me.shift(CMEDir, offset=delta1, pa=offangle1)
    major = comp[3]
    SizeMajor=qa.toangle(str(major)+'deg')
    minor = comp[4]
    SizeMinor=qa.toangle(str(minor)+ 'deg')
    PA= comp[5] #(90 + initoffangle['value'])%360. #comp[5]
    PA=qa.toangle(str(PA)+'deg')
    Flux= comp[6]
    resolution = comp[7] #in arcmin
    imSize = int(comp[8])

    #Construct an empty casa image from a shape
    ia.fromshape(truthName,shape=[imSize,imSize,1,1],overwrite=True)
    #adding components to  the empty image file
    cl.addcomponent(dir=NewDir, \
               flux=Flux, fluxunit='Jy', freq=str(freq0)+'MHz', \
                   shape="Gaussian", majoraxis=SizeMajor, minoraxis=SizeMinor, positionangle=PA)

    # cl.addcomponent(dir=CMEDir, flux=1e10, fluxunit='Jy', freq=comp[0],
    #                 shape='point')
    #
    # cl.addcomponent(dir=NewDir, flux=1e9, fluxunit='Jy', freq=comp[0],
    #                 shape='point')
    #
    # cl.addcomponent(dir=NewDir2, flux=7e9, fluxunit='Jy', freq=comp[0],
    #                  shape='point')

    cs=ia.coordsys()
    cs.setunits(['rad','rad','','Hz'])
    cell_rad=qa.convert(qa.quantity(str(resolution)+"arcmin"),"rad")['value']
    cs.setincrement([-cell_rad,cell_rad],'direction')
    cs.setreferencevalue([CMEDir['m0']['value'], CMEDir['m1']['value']],
                          type="direction")
    cs.setrestfrequency(comp[0])
    # add important header keywords
    imhead(imagename=truthName,mode="put",hdkey="object",hdvalue="Model CME")
    imhead(imagename=truthName,mode="put",hdkey="imtype",hdvalue='Intensity')
    imhead(imagename=truthName,
          mode="put",hdkey="observer",hdvalue="simulation")
    imhead(imagename=truthName,
          mode="put",hdkey="date-obs",hdvalue="2023/03/15/00:00:00")
    imhead(imagename=truthName,mode="put",hdkey="reffreqtype",hdvalue='TOPO')
    imhead(imagename=truthName,
           mode="put",hdkey="restfreq",hdvalue=str(comp[0])+'MHz')
    imhead(imagename=truthName,mode='list')
    cs.setreferencevalue(str(comp[0])+'MHz', 'spectral')
    cs.settelescope('VLA')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(),subtract=False)

    # pix = ia.getchunk()
    # pix = pix.reshape((imSize, imSize))

    # p = np.flipud(pix.T)

    ia.close()
    cl.done()
    qa.done()
    me.done()

    return #pix


def imageForm(filename, truthImage, freq, dataArr):

    brightness = sum(dataArr)

    print 'tot brightness is ' + str(brightness)
    imSize = 400
    zoomimg = spndint.zoom(dataArr,float(imSize)/400)
    zdims = np.shape(zoomimg)

    for i in range(zdims[0]):
        for j in range(zdims[1]):
            if zoomimg[i][j] < 0.0:
                zoomimg[i][j] = 0.0

    newbrightness = sum(zoomimg)
    zoomimg *= brightness/newbrightness

    print 'tot brightness is ' + str(sum(zoomimg))


    z = zoomimg.copy().T
    z= np.fliplr(z)  #these operations flip to CASA style of storing data
    #which starts at lower left corner, and goes up by columns left to right


    casaArr = z.reshape((imSize, imSize, 1, 1))

    toutfile=os.path.join('.', truthImage)
    #if the truth image already exists remove it

    if os.path.exists(toutfile):
        shutil.rmtree(toutfile)

    ia.fromarray(truthImage, pixels=casaArr, overwrite=True)


    cs=ia.coordsys()
    cs.setunits(['rad','rad','','Hz'])
    res1=400
    res = 400
    width = res1*.038 * np.pi/180.
    dres = width/res
    cell_rad=dres #qa.convert(qa.quantity("1arcmin"),"rad")['value']
    imSize = res
    cs.setincrement([-cell_rad,cell_rad],'direction')
    cs.setreferencevalue([CMEDir['m0']['value'], CMEDir['m1']['value']],
                          type="direction")
    cs.setrestfrequency(freq)


    # add important header keywords
    imhead(imagename=truthImage,mode="put",hdkey="object",hdvalue="DRAGN")
    imhead(imagename=truthImage,mode="put",hdkey="imtype",hdvalue='Intensity')
    imhead(imagename=truthImage,
          mode="put",hdkey="observer",hdvalue="simulation")
    imhead(imagename=truthImage,
          mode="put",hdkey="date-obs",hdvalue="2023/03/15/00:00:00")
    imhead(imagename=truthImage,mode="put",hdkey="reffreqtype",hdvalue='TOPO')
    imhead(imagename=truthImage,
           mode="put",hdkey="restfreq",hdvalue=str(freq))
    imhead(imagename=truthImage,mode='list')
    cs.setreferencevalue(str(freq)+'Hz', 'spectral')
    Telescope='VLA' #or else it breaks and whines
    cs.settelescope(Telescope)
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")

    pix = ia.getchunk()
    pix = pix.reshape((imSize, imSize))

    zoomimg = pix


    ia.close()

    return zoomimg


def imageData(MS, imageName, imWidthDeg=2., nIter=10, Threshold = '3000000.0Jy'):
    '''
    imageData(MS, imageName, imWidthDeg=10., nIter=10, Threshold = '30000000.0Jy')
    Takes in ms file MS (string), creates image with basename of imageName (string),
    imWidthDeg is width of image desired, in degrees.  nIter is number of model update & deconvolution cleaning cycles
    to halt after, Threshold is the alternate way of deciding when to stop cleaning, when no unmodeled sources are above Threshold
    '''
    print('cleaning dirty image')

    # hz=float(MS[11:16])*10**6
    # hz=float(MS[11:14])*10**6 #for 3 subs

    #grab frequency from MS
    ms.open(MS)
    spws = ms.getspectralwindowinfo()
    ms.close()

    tb.open(MS, nomodify=False)
    specInd = tb.getcol("DATA_DESC_ID")
    uvws = tb.getcol("UVW")

    hz = spws[str(specInd[0])]['Chan1Freq']

    uvs = uvws.T[:,:2]
    uvnorms = np.sum(np.abs(uvs)**2,axis=-1)**(1./2)
    largestBL = amax(uvnorms)


    wavelen = 3e8/hz
    cellsizerad = wavelen/largestBL/16. #radians, oversample by factor of 4, 16 may be better, but takes longer.
    csas = cellsizerad*180./np.pi*3600 #to arcseconds


    imsize = int(imWidthDeg*3600./csas) #int(1024/60.*3600./csas)
    print('calculated imsize is ' + str(imsize))
    if imsize < 64:
        imsize = 64

    #get next power of 2
    imsize = int(pow(2, math.ceil(math.log(imsize, 2))))

    print('imsize used is ' + str(imsize))

    print('cell size for ' + MS + ' is ' + str(csas) +' arcsec')

    #change flags, being stupid
    f = tb.getcol("FLAG")
    f2 = tb.getcol("FLAG_ROW")
    flaglen = np.shape(f)[2]
    for i in range(flaglen):
        f[0][0][i] = False
        f[1][0][i] = False
        f2[i] = False


    tb.putcol("FLAG", f)
    tb.putcol("FLAG_ROW", f2)

    tb.close()

    try:
        tclean(vis=MS,imagename=imageName,
                    outlierfile='',
                    field='',spw='',
                    selectdata=False,
                    nterms=1,
                    gridder='widefield', wprojplanes = 1, facets = 1,
                    niter=0,gain=0.8,threshold=Threshold,
                    deconvolver='hogbom',
                    interactive=False,
                    mask=[],
                    imsize=[imsize, imsize],
                    cell=[str(csas)+'arcsec', str(csas)+'arcsec'] ,#cell=['1arcsec', '1arcsec'],
                    phasecenter='',
                    stokes='I',
                    startmodel='', #truthim,
                    weighting='briggs',robust=0.
                    )

        imName = imageName + '.image'
        psfName = imageName + '.psf'
        modelName = imageName + '.model'

        # imview(raster={'file': 'EarthSynch_CONCAT_0.760869565217MHz_1_timestep.ms', 'colorwedge': True}, zoom={'blc': [175,175], 'trc': [225,225]}, out='psf128.png')
        imview(raster={'file': psfName, 'colorwedge': True}, axes={'y':'Declination'}, out=psfName+'BARpsf.png')
        imview(raster={'file': imName, 'colorwedge': True},  axes={'y':'Declination'}, out=imName+'BARimDirty.png')

        tclean(vis=MS,imagename=imageName,
                    outlierfile='',
                    field='',spw='',
                    selectdata=False,
                    nterms=1,
                    gridder='widefield', wprojplanes = 1, facets = 1,
                    niter=nIter,gain=0.8,threshold=Threshold,
                    deconvolver='hogbom',
                    interactive=False,
                    mask=[],
                    imsize=[imsize, imsize],
                    cell=[str(csas)+'arcsec', str(csas)+'arcsec'] ,#cell=['1arcsec', '1arcsec'],
                    phasecenter='',
                    stokes='I',
                    startmodel='', #truthim,
                    weighting='briggs',robust=0.
                    )

        imview(raster={'file': imName, 'colorwedge': True},  axes={'y':'Declination'}, out=imName+'BARimClean.png')
        imview(raster={'file': modelName, 'colorwedge': True},  axes={'y':'Declination'}, out=modelName+'BARmodel.png')

    except:
        print('Clean wanted to crash on ' + MS)
        pass

    return


###########################################################################

if __name__ == "__main__":

    """
    main section
    """

    import argparse
    ##define the input argument, only need numSC currently
    args = _getParser()


    #in index=1 we have the sun mostly overhead the array
    xyzs = np.loadtxt('eqXYZ_EarthCentered.txt')
    xyzs = xyzs.reshape((args.numSC, 48, 3))
    xyzs *= 1e3

    ras = np.loadtxt('RAs.txt')
    decs = np.loadtxt('Decs.txt')

    sunRAdeg = ras[0]
    sunDecdeg = decs[0]


    #This points to Earth from Moon
    CMEDir = me.direction('J2000', str(sunRAdeg)+'deg', str(sunDecdeg) + 'deg')

    print 'CMEDir'
    print CMEDir

    #compute baselines to put into ms

    ra = CMEDir['m0']['value'] #returns radians
    dec = CMEDir['m1']['value'] #returns radians

    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)

    s = np.array([x, y, z])

    #get x and y axes perp to s
    #get correct ones though!!!!

    x = np.cross(s, [0, 0, 1])
    y = np.cross(s, x)

    #normalize since cross product isn't already 1 for some reason
    # this is
    px = x/np.linalg.norm(x)
    py = y/np.linalg.norm(y)
    pz = s/np.linalg.norm(s)

    srcdir = [px, py, pz]
    length = 1

    #use argument startPos, defaults to 0 in not in command line
    phasescpos = xyzs[:,args.startPos,:]

    phasescpos = np.expand_dims(phasescpos, axis=1)
    allscpos = phasescpos.copy()

    numSC = args.numSC

    #process data
    for l in range(length):

        #reproject pos and err into new frames
        for i in range(numSC):
            phasescpos[i][l][0] = np.dot(allscpos[i][l], px)
            phasescpos[i][l][1] = np.dot(allscpos[i][l], py)
            phasescpos[i][l][2] = np.dot(allscpos[i][l], pz)



    phasescerr = np.zeros(np.shape(phasescpos))

    #now make baselines from phase centered  projected positions
    numSC = np.shape(phasescpos)[0] #args.numSC
    numbl = numSC*(numSC - 1)/2

    baselines = np.zeros((numbl, length, 3), dtype=np.float64)
    baselineserr = np.zeros((numbl, length, 3), dtype=np.float64)



    for l in range(length):
        #calc baselines for this timestep
        k=0
        for i in range(numSC):
          for j in range(i+1, numSC):
            baselines[k][l] =  phasescpos[j][l] - phasescpos[i][l]
            baselineserr[k][l] =  phasescerr[j][l] - phasescerr[i][l]
            k+=1


    #baselines2 = np.reshape(baselines[:,l,:], (numbl,3), order='F').T

    bigBL = 0.
    for i in range(numbl):
        for j in range(length):
            cont = np.sqrt(baselines[i][j][0]**2 + baselines[i][j][1]**2)
            if cont > bigBL:
                bigBL = cont

    zBLErrs = baselineserr[:,:,2] #z component in new frame, in meters
    zErrs = phasescerr[:,0,2] #get just first timestep, is the same throughout anyway, prints shorter
	#print out errors
    print 'total z errors in meters'
    print str(zErrs)
    print 'errors relative to first SC'
    print str(zErrs - zErrs[0])
    print 'total z errors in seconds'
    print str(zErrs/(3e8))
    print 'errors relative to first SC'
    print str((zErrs - zErrs[0])/3e8)
    print 'total relative errors in phase = seconds * 2pi'
    print str(2*np.pi*(zErrs  - zErrs[0])/3e8)

    print args

    #close previosuly exising components
    cl.done()
    qa.done()
    me.done()
    countForm = 0

    length=1


    freqs=[.702, .5, 3.001, 2.001, 4.002, .098]
    truthImages = ['AKR_Blob_0.702MHz.truth', 'Hiss_Blob_0.500MHz.truth', 'AR3_Blob_3.001MHz.truth', 'MF2_Blob_2.001MHz.truth',  'MF4_Blob_4.002MHz.truth', 'TCE_Blob_0.098MHz.truth']

    freqs=[4.002, .098]
    truthImages = ['MF4_Blob_4.002MHz.truth', 'TCE_Blob_0.098MHz.truth']

    freqs=[.702, .5, 3.001, 4.002]
    truthImages = ['AKR_Blob_0.702MHz.truth', 'Hiss_Blob_0.500MHz.truth', 'AR3_Blob_3.001MHz.truth',  'MF4_Blob_4.002MHz.truth']


    for f in range(len(freqs)): #currently only using 1 frequency band, may expand later

        freq = freqs[f]

        #Always referece same antenna locations in sm.setconfig to allow MS consistency after we swap stuff later on
        Ant_x = phasescpos[:,0,0]
        Ant_y = phasescpos[:,0,1]
        Ant_z = phasescpos[:,0,2]

        truthImage = truthImages[f]

        # ia.open(truthImage)
        # ia.fft(amp=truthImage + 'Amp.im', phase=truthImage + 'Phase.im')
        # viewer(infile=truthImage + 'Amp.im',displaytype='raster',
        #       outfile=truthImage + 'Amp.im'+'.jpg',outformat='jpg',
        #       gui=False)
        # viewer(infile=truthImage + 'Phase.im',displaytype='raster',
        #     outfile=truthImage + 'Phase.im'+'.jpg',outformat='jpg',
        #     gui=False)
        # ia.close()
        ## define the list of nearby frequency channels to image


        #iterate over time
        for l in range(length):


            counter = 0
            #not changing baselines within our 1 day averaging, so only one set of baselines/visibilities necessary
            for i in range(1):


                print 'tic freq'

                sm.close()
                me.done()
                vp.reset()
                #start up the empty measurement set
                lunarMS='EarthSynch128-sim1-%06.3fMHz_timestep_%d.ms'%(freq,l)
                print("components[i,0] (MHz), timestep", (freq,l))



                #call the arraySim function to simulate the observations

                arraySim(lunarMS, freq, truthImage, 0., Ant_x, Ant_y, Ant_z, args, baselines, dishSize=10, timestep= l)

                sm.close()



                print('Now fitting Gaussian model to recovered data')


        #image
        imageData(lunarMS, truthImage+'NOISELESS', imWidthDeg=10., nIter=0, Threshold = '3000000.0Jy')
        #add noise
        addThermalNoise(lunarMS, SEFD=-1., dv = 10000., dt = 60.)

        #image
        imageData(lunarMS, truthImage, imWidthDeg=10., nIter=0, Threshold = '3000000.0Jy')



        #fit data

        truthData = imstat(truthImage)

        brightPoint = truthData['maxposf'].split(',')
        bFlux = truthData['max'][0]
        rastring = brightPoint[0]
        decstring = brightPoint[1]
        brightDirTruth = me.direction('J2000', rastring, decstring)
        ravaltruth = brightDirTruth['m0']['value'] #returns radians
        decvaltruth = brightDirTruth['m1']['value'] #returns radians

        print('True Input Params (flux (Jy), ra (rad), dec (rad), major (arcsec), minor (arcsec), PA (degrees))')
        print('%2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f'%(bFlux, ravaltruth, decvaltruth, 0., 0., 0.))#components[0,3]*3600., components[0,4]*3600., components[0,5]))


        try:
            fitCleanData = imfit(truthImage+'.image')

            majguess = fitCleanData['results']['component0']['shape']['majoraxis']['value']
            minguess = fitCleanData['results']['component0']['shape']['minoraxis']['value']
            raguess = fitCleanData['results']['component0']['shape']['direction']['m0']['value']
            decguess = fitCleanData['results']['component0']['shape']['direction']['m1']['value']
            fluxguess = fitCleanData['results']['component0']['flux']['value'][0]
            paguess = fitCleanData['results']['component0']['shape']['positionangle']['value']

            print('Dirty Image Determined Params (flux (Jy), ra (rad), dec (rad), major (arcsec), minor (arcsec), PA (degrees))')
            print('%2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f'%(fluxguess, raguess, decguess, majguess, minguess, paguess))

            raoff = (raguess - ra)/pi*180.*3600.  #distance in arcsec from center of image, expected format of uvmodelfit
            decoff = (decguess - dec)/pi*180.*3600.

        except:
            print('Failed imfit on ' + truthImage+'.image')


        try:

            uvclname = lunarMS+'.uv.cl'

            if os.path.exists(uvclname):
                shutil.rmtree(uvclname)

            uvmodelfit(lunarMS, niter=5, comptype='G', sourcepar=[fluxguess,raoff,decoff,majguess, minguess/majguess, paguess], outfile=uvclname) # Output component list file

            cl.open(uvclname)
            fit = cl.getcomponent(0)            # stores component information

            fluxuv = fit['flux']['value'][0]     #   to store the I,Q,U,V, flux
            rauv = fit['shape']['direction']['m0']['value']
            decuv =fit['shape']['direction']['m1']['value']
            bmajuv = fit['shape']['majoraxis']['value']*60. #    to get major axis arcsec
            bminuv = fit['shape']['minoraxis']['value']*60.  #   to get minor axis arcsec
            pauv = fit['shape']['positionangle']['value']


            print('UVModelFit Determined Params (flux (Jy), ra (rad), dec (rad), major (arcsec), minor (arcsec), PA (degrees))')
            print('%2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f'%(fluxuv, rauv, decuv, bmajuv, bminuv, pauv))

        except:
            print('Failed uvmodelfit on ' + lunarMS)
