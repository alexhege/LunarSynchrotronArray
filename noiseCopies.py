# noiseCopies.py
# Alex Hegedus 12/11/19
# This file is to be run after the ms file has been generated.
# This script runs tclean with different settings
# First it creates a noiseless image
# then adds noise for different integration times or noise levels and images
# multiple robustness weighting parameters used in each case (-1, -0.5, 0, 0.5, 1.0)

#import the required python modules
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

import subprocess
################################################################################

#np.sum(np.abs(x)**2,axis=-1)**(1./2) #for uv norm

MS='EarthSynch_CONCAT_0.760869565217MHz_1_timestep.ms'

# MS='EarthSynchNoNoise2048_CONCAT_0.760869565217MHz_1_timestep.ms'

print 'saving UV stats'
plotms(MS, xaxis='u', yaxis='v', spw='', timerange='', plotfile='UVConcat.jpg', expformat = 'jpg', showgui=False, overwrite=True)
plotms(MS, xaxis='uwave', yaxis='vwave', spw='', timerange='', plotfile='UVwave.jpg', expformat = 'jpg', showgui=False, overwrite=True)
tb.open(MS)
uvws = tb.getcol('UVW')
tb.close()
uvws = uvws.T
uvs = uvws[:,:2]
uvNorms = np.zeros(np.shape(uvs)[0])
for n in range(len(uvNorms)):
    uvNorms[n] = np.linalg.norm(uvs[n,:])

np.savetxt('UVWDataEarth.out', uvws)
np.savetxt('UVNormEarth.out', uvNorms)


figure()
hist(uvNorms/450.)
savefig('uvNormWaveHist.png')
clf()


#copy MS for noise and no noise copies, copy data into both columns, noise will be added to DATA
returned_value = subprocess.call('cp -rf '+MS+' EarthSynchNoise_CONCAT_0.760869565217MHz_1_timestep.ms', shell=True)  # returns the exit code in unix
print('returned value:', returned_value)

returned_value = subprocess.call('cp -rf '+MS+' EarthSynchNoNoise_CONCAT_0.760869565217MHz_1_timestep.ms', shell=True)  # returns the exit code in unix
print('returned value:', returned_value)

tb.open('EarthSynchNoise_CONCAT_0.760869565217MHz_1_timestep.ms', nomodify=False)

cdata = tb.getcol("CORRECTED_DATA")

tb.putcol("DATA", cdata)

tb.close()

tb.open('EarthSynchNoNoise_CONCAT_0.760869565217MHz_1_timestep.ms', nomodify=False)

data = tb.getcol("DATA")


tb.putcol("CORRECTED_DATA", data)

tb.putcol("DATA", data)

tb.close()


#make 2048 version

# concat(vis=['EarthSynchNoNoise_CONCAT_0.760869565217MHz_1_timestep.ms', 'EarthSynchNoNoise_CONCAT_0.760869565217MHz_1_timestep.ms', 'EarthSynchNoNoise_CONCAT_0.760869565217MHz_1_timestep.ms', 'EarthSynchNoNoise_CONCAT_0.760869565217MHz_1_timestep.ms'],
#        concatvis='EarthSynchNoNoise1024_CONCAT_0.760869565217MHz_1_timestep.ms',
#        freqtol='',
#        dirtol='',
#        respectname=False,
#        timesort=False,
#        copypointing=True,
#        visweightscale=[])
#
# vis2048 = []
# for i in range(4):
#     vis2048.append('EarthSynchNoNoise_CONCAT_0.760869565217MHz_1_timestep.ms')
# #
# #
# concat(vis=vis2048,
#        concatvis='EarthSynchNoNoise2048_CONCAT_0.760869565217MHz_1_timestep.ms',
#        freqtol='',
#        dirtol='',
#        respectname=False,
#        timesort=False,
#        copypointing=True,
#        visweightscale=[])


##now clean them!!
MS='EarthSynchNoNoise_CONCAT_0.760869565217MHz_1_timestep.ms'
ff=16384
imsize = 64 # 128 # for 10 km,  64# for 6 km, 20 km is 256#for big, small is 256 #2048
csas = 829. #498. # for 10 km, 829.# for 6 km, 20km is 250. # 165 #14.7

tclean(vis=MS,imagename=MS+'.dirty'+str(ff),
            outlierfile='',
            field='',spw='',
            selectdata=False,
            nterms=1,
            gridder='widefield', wprojplanes = 1, facets = 1,
            niter=0,gain=0.1,threshold='0.0mJy',
            deconvolver='hogbom',
            interactive=False,
            mask=[],
            imsize=[imsize, imsize],
            cell=[str(csas)+'arcsec', str(csas)+'arcsec'] ,#cell=['1arcsec', '1arcsec'],
            phasecenter='',
            stokes='I',
            startmodel='', #truthim,
            weighting='briggs',robust=-0.5
            )

imName = MS+'.dirty'+str(ff) + '.image'
psfName = MS+'.dirty'+str(ff) + '.psf'
#     outfile=os.path.join('.', imName+'.jpg')
# #if the truth image already exists remove it
#
#     if os.path.exists(outfile):
#         os.remove(outfile)
# viewer(infile=imName,displaytype='raster',
#       outfile=imName+'.jpg',outformat='jpg',
#       gui=False)
#
# viewer(infile=psfName,displaytype='raster',
#       outfile=psfName+'.jpg',outformat='jpg',
#       gui=False)

# imview(raster={'file': 'EarthSynch_CONCAT_0.760869565217MHz_1_timestep.ms', 'colorwedge': True}, zoom={'blc': [175,175], 'trc': [225,225]}, out='psf128.png')
# imview(raster={'file': 'EarthSynch-00.736MHz.truthPhase.im', 'colorwedge': True},  axes={'y':'Declination'}, out='EarthSynch-00.736MHz.truthPhase.im.png')
imview(raster={'file': psfName, 'colorwedge': True}, axes={'y':'Declination'}, out=psfName+'BARpsf.png')
imview(raster={'file': imName, 'colorwedge': True},  axes={'y':'Declination'}, out=imName+'BARim.png')

print 'saving UV stats'

plotms(MS, xaxis='uvwave', yaxis='amp', spw='', timerange='', plotfile='UVwaveampNoNoise'+str(ff)+'.jpg', expformat = 'jpg', showgui=False, overwrite=True)



###now add diff noise levels into UVwaveampNoNoise


avNoise = 1.38e7 #Jy,

nmode = 'Opt' #put in file names to indicate noise level

print('no integration raw avNoise is ' + str(avNoise) + ' Jy, 1 Jy = 1e-26 W/m^2/Hz')


#SEFD is brightness average*4pi in W/m^2/Hz, 1 Jy = 1e-26 W/m^2/Hz

#SEFD con 1.46e-18  1000/cc density Luna measurements
#SEFD mod 4.6e-19   250/cc density  SZA < 60 deg SELENE noise
#SEFD opt 1.38e-19  8/cc density (amplifier dom) night side


#make copies using unix command line

returned_value = subprocess.call('cp -rf EarthSynchNoNoise_CONCAT_0.760869565217MHz_1_timestep.ms EarthSynchNoise1hr_CONCAT_0.760869565217MHz_1_timestep.ms', shell=True)  # returns the exit code in unix
print('returned value:', returned_value)

returned_value = subprocess.call('cp -rf EarthSynchNoNoise_CONCAT_0.760869565217MHz_1_timestep.ms EarthSynchNoise2hr_CONCAT_0.760869565217MHz_1_timestep.ms', shell=True)  # returns the exit code in unix
print('returned value:', returned_value)


returned_value = subprocess.call('cp -rf EarthSynchNoNoise_CONCAT_0.760869565217MHz_1_timestep.ms EarthSynchNoise4hr_CONCAT_0.760869565217MHz_1_timestep.ms', shell=True)  # returns the exit code in unix
print('returned value:', returned_value)


returned_value = subprocess.call('cp -rf EarthSynchNoNoise_CONCAT_0.760869565217MHz_1_timestep.ms EarthSynchNoise6hr_CONCAT_0.760869565217MHz_1_timestep.ms', shell=True)  # returns the exit code in unix
print('returned value:', returned_value)


returned_value = subprocess.call('cp -rf EarthSynchNoNoise_CONCAT_0.760869565217MHz_1_timestep.ms EarthSynchNoise12hr_CONCAT_0.760869565217MHz_1_timestep.ms', shell=True)  # returns the exit code in unix
print('returned value:', returned_value)


returned_value = subprocess.call('cp -rf EarthSynchNoNoise_CONCAT_0.760869565217MHz_1_timestep.ms EarthSynchNoise24hr_CONCAT_0.760869565217MHz_1_timestep.ms', shell=True)  # returns the exit code in unix
print('returned value:', returned_value)


MSs = ['EarthSynchNoise1hr_CONCAT_0.760869565217MHz_1_timestep.ms', 'EarthSynchNoise2hr_CONCAT_0.760869565217MHz_1_timestep.ms', 'EarthSynchNoise4hr_CONCAT_0.760869565217MHz_1_timestep.ms', 'EarthSynchNoise6hr_CONCAT_0.760869565217MHz_1_timestep.ms', 'EarthSynchNoise12hr_CONCAT_0.760869565217MHz_1_timestep.ms', 'EarthSynchNoise24hr_CONCAT_0.760869565217MHz_1_timestep.ms']


hours = [1, 2, 4, 6, 12, 24]

#all in muinutes actually
#loop through, add noise for different integration times, and for each time
# make multiple images with different robustness weighting values  -0.5 works best for log circular array
for i in range(len(hours)):
    hourIntegration = hours[i]
    noise = avNoise/0.8/sqrt(2*hourIntegration*60*60*500000) #sensitivity of single polarization, 2 antenna, 0.8 efficiency 500kHz bandwidth
    noise /= 16. #change 1024 to 16k antenna
    print avNoise
    print noise

    MS = MSs[i]

    tb.open(MS, nomodify=False)



# noise = avNoise/sqrt(2*2*24*60*60*500000/3)

    data=tb.getcol("DATA")
    cdata=tb.getcol("CORRECTED_DATA")

    # print("pre noise vals")

    # print data[0][0][:]

    print 'single visibility noise rms for ' + str(hourIntegration) + ' hours integration'
    print noise*sqrt(2)



    oldData = data[0][0][:].copy()

    numbl = shape(data)[2]

    for k in range(numbl):
        toAdd = np.random.normal(0, noise)*1.j + np.random.normal(0, noise)

        #leave noise free data in data column, clean checks/uses corrected data first
        # data[0][0][k] += toAdd
        # data[1][0][k] += toAdd


        cdata[0][0][k] += toAdd
        cdata[1][0][k] += toAdd

    tb.putcol("DATA", data)
    tb.putcol("CORRECTED_DATA", cdata)

    # print("post noise vals")
    #
    # print data[0][0][:]
    #
    # print 'noise diff'
    # print str(data[0][0][:] - oldData)

    #
    # tclean(vis=MS,imagename=MS+'.dirty'+nmode+str(ff),
    #             outlierfile='',
    #             field='',spw='',
    #             selectdata=False,
    #             nterms=1,
    #             gridder='widefield', wprojplanes = 1, facets = 1,
    #             niter=200,gain=0.1,threshold='0.0mJy',
    #             deconvolver='hogbom',
    #             interactive=False,
    #             mask=[],
    #             imsize=[imsize, imsize],
    #             cell=[str(csas)+'arcsec', str(csas)+'arcsec'] ,#cell=['1arcsec', '1arcsec'],
    #             phasecenter='',
    #             stokes='I',
    #             startmodel='', #truthim,
    #             weighting='briggs',robust=-0.5
    #             )
    #
    # imName = MS+'.dirty'+nmode+str(ff) + '.image'
    # # psfName = MS+'.dirtyMod'+str(ff) + '.psf'
    #
    # # imview(raster={'file': psfName, 'colorwedge': True}, axes={'y':'Declination'}, out=psfName+'BARpsf'+str(hours[i])+'.png')
    # imview(raster={'file': imName, 'colorwedge': True},  axes={'y':'Declination'}, out=imName+'BARim'+str(hours[i])+'.png')


    tclean(vis=MS,imagename=MS+'.n5dirty'+nmode+str(ff),
                outlierfile='',
                field='',spw='',
                selectdata=False,
                nterms=1,
                gridder='widefield', wprojplanes = 1, facets = 1,
                niter=0,gain=0.1,threshold='0.0mJy',
                deconvolver='hogbom',
                interactive=False,
                mask=[],
                imsize=[imsize, imsize],
                cell=[str(csas)+'arcsec', str(csas)+'arcsec'] ,#cell=['1arcsec', '1arcsec'],
                phasecenter='',
                stokes='I',
                startmodel='', #truthim,
                weighting='briggs',robust=-0.5
                )

    imName = MS+'.n5dirty'+nmode+str(ff) + '.image'
    # psfName = MS+'.dirtyMod'+str(ff) + '.psf'

    # imview(raster={'file': psfName, 'colorwedge': True}, axes={'y':'Declination'}, out=psfName+'BARpsf'+str(hours[i])+'.png')
    imview(raster={'file': imName, 'colorwedge': True},  axes={'y':'Declination'}, out=imName+'BARim'+str(hours[i])+'.png')



    # print 'saving UV stats'
    plotms(MS, xaxis='uvwave', yaxis='amp', spw='', timerange='', plotfile=MS+'UVwaveampNoise'+str(ff)+'-'+str(hours[i])+'.jpg', expformat = 'jpg', showgui=False, overwrite=True)

    tclean(vis=MS,imagename=MS+'.5dirty'+nmode+str(ff),
                outlierfile='',
                field='',spw='',
                selectdata=False,
                nterms=1,
                gridder='widefield', wprojplanes = 1, facets = 1,
                niter=0,gain=0.1,threshold='0.0mJy',
                deconvolver='hogbom',
                interactive=False,
                mask=[],
                imsize=[imsize, imsize],
                cell=[str(csas)+'arcsec', str(csas)+'arcsec'] ,#cell=['1arcsec', '1arcsec'],
                phasecenter='',
                stokes='I',
                startmodel='', #truthim,
                weighting='briggs',robust=0.5
                )

    imName = MS+'.5dirty'+nmode+str(ff) + '.image'
    # psfName = MS+'.dirtyMod'+str(ff) + '.psf'

    # imview(raster={'file': psfName, 'colorwedge': True}, axes={'y':'Declination'}, out=psfName+'BARpsf'+str(hours[i])+'.png')
    imview(raster={'file': imName, 'colorwedge': True},  axes={'y':'Declination'}, out=imName+'BARim'+str(hours[i])+'.png')



    # print 'saving UV stats'
    plotms(MS, xaxis='uvwave', yaxis='amp', spw='', timerange='', plotfile=MS+'UVwaveampNoise'+str(ff)+'-'+str(hours[i])+'.jpg', expformat = 'jpg', showgui=False, overwrite=True)

    tclean(vis=MS,imagename=MS+'.0dirty'+nmode+str(ff),
                outlierfile='',
                field='',spw='',
                selectdata=False,
                nterms=1,
                gridder='widefield', wprojplanes = 1, facets = 1,
                niter=0,gain=0.1,threshold='0.0mJy',
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

    imName = MS+'.0dirty'+nmode+str(ff) + '.image'
    # psfName = MS+'.dirtyMod'+str(ff) + '.psf'

    # imview(raster={'file': psfName, 'colorwedge': True}, axes={'y':'Declination'}, out=psfName+'BARpsf'+str(hours[i])+'.png')
    imview(raster={'file': imName, 'colorwedge': True},  axes={'y':'Declination'}, out=imName+'BARim'+str(hours[i])+'.png')



    # print 'saving UV stats'
    plotms(MS, xaxis='uvwave', yaxis='amp', spw='', timerange='', plotfile=MS+'UVwaveampNoise'+str(ff)+'-'+str(hours[i])+'.jpg', expformat = 'jpg', showgui=False, overwrite=True)

    tclean(vis=MS,imagename=MS+'.n1dirty'+nmode+str(ff),
                outlierfile='',
                field='',spw='',
                selectdata=False,
                nterms=1,
                gridder='widefield', wprojplanes = 1, facets = 1,
                niter=0,gain=0.1,threshold='0.0mJy',
                deconvolver='hogbom',
                interactive=False,
                mask=[],
                imsize=[imsize, imsize],
                cell=[str(csas)+'arcsec', str(csas)+'arcsec'] ,#cell=['1arcsec', '1arcsec'],
                phasecenter='',
                stokes='I',
                startmodel='', #truthim,
                weighting='briggs',robust=-1.
                )

    imName = MS+'.n1dirty'+nmode+str(ff) + '.image'
    # psfName = MS+'.dirtyMod'+str(ff) + '.psf'

    # imview(raster={'file': psfName, 'colorwedge': True}, axes={'y':'Declination'}, out=psfName+'BARpsf'+str(hours[i])+'.png')
    imview(raster={'file': imName, 'colorwedge': True},  axes={'y':'Declination'}, out=imName+'BARim'+str(hours[i])+'.png')



    # print 'saving UV stats'
    plotms(MS, xaxis='uvwave', yaxis='amp', spw='', timerange='', plotfile=MS+'UVwaveampNoise'+str(ff)+'-'+str(hours[i])+'.jpg', expformat = 'jpg', showgui=False, overwrite=True)

    tclean(vis=MS,imagename=MS+'.1dirty'+nmode+str(ff),
                outlierfile='',
                field='',spw='',
                selectdata=False,
                nterms=1,
                gridder='widefield', wprojplanes = 1, facets = 1,
                niter=0,gain=0.1,threshold='0.0mJy',
                deconvolver='hogbom',
                interactive=False,
                mask=[],
                imsize=[imsize, imsize],
                cell=[str(csas)+'arcsec', str(csas)+'arcsec'] ,#cell=['1arcsec', '1arcsec'],
                phasecenter='',
                stokes='I',
                startmodel='', #truthim,
                weighting='briggs',robust=1.
                )

    imName = MS+'.1dirty'+nmode+str(ff) + '.image'
    # psfName = MS+'.dirtyMod'+str(ff) + '.psf'

    # imview(raster={'file': psfName, 'colorwedge': True}, axes={'y':'Declination'}, out=psfName+'BARpsf'+str(hours[i])+'.png')
    imview(raster={'file': imName, 'colorwedge': True},  axes={'y':'Declination'}, out=imName+'BARim'+str(hours[i])+'.png')



    # print 'saving UV stats'
    plotms(MS, xaxis='uvwave', yaxis='amp', spw='', timerange='', plotfile=MS+'UVwaveampNoise'+str(ff)+'-'+str(hours[i])+'.jpg', expformat = 'jpg', showgui=False, overwrite=True)
