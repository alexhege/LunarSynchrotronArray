This is the README for the LunarSynchrotronArray package, maintained by Dr. Alex Hegedus alexhege@umich.edu

This set of codes should guide you through making the figures in the paper, as well as hopefully being accessible enough for changing the code for your own array.  I would encourage you to please reach out to collaborate if that is the case!

Requirements:  CASA 4.7.1 (or greater?) built on python 2.7
Example link for Red Hat 7 
https://casa.nrao.edu/download/distro/casa/release/el7/casa-release-4.7.1-el7.tar.gz

gcc 4.8.5 or above (or below?)

SPICE (I use cspice here)
https://naif.jpl.nasa.gov/naif/toolkit_C.html
As seen in lunar_furnsh.txt which loads the SPICE kernels, you also must download 
KERNELS_TO_LOAD = ( '/home/alexhege/SPICE/LunarEph/moon_pa_de421_1900-2050.bpc'
                   '/home/alexhege/SPICE/LunarEph/moon_080317.tf'
                   '/home/alexhege/SPICE/LunarEph/moon_assoc_me.tf'
                   '/home/alexhege/SPICE/LunarEph/pck00010.tpc'
                   '/home/alexhege/SPICE/LunarEph/naif0008.tls'
                   '/home/alexhege/SPICE/LunarEph/de430.bsp' )

All of which can be found at 
https://naif.jpl.nasa.gov/pub/naif/generic_kernels/

SLDEM2015_128_60S_60N_000_360_FLOAT.IMG for the lunar surface data by LRO LOLA
Found at 
http://imbrium.mit.edu/DATA/SLDEM2015/GLOBAL/FLOAT_IMG/




Codes within this package in order of use in pipeline:
0. snrCalc.py 
This regular python script creates the noise budget plots with frequency, not necessary for pipeline.  Information about the frequency dependence of the November 1st 2016 stormy synchrotron brightness as a function of frequency is also contained here.

1. createArrayConfig.py
This regular python script sets your array configurations (currently log spaced circular for 1024 elements).  Reads the SLDEM data and saves the longitude, latitude, and altitude of each antenna.
Also you must update the variable lunarPath to the current location of the SLDEM file

2. eqArrOverTimeEarth.c 
This is runs the SPICE kernel, taking in the Long Lat and Alt data from the first script, as well as the period of time you want to observe.  This will align the frame of the antenna on the Lunar Surface to the area of the sky containing the targeted Earth overhead.  This saves a number of variables as output files, most important being the XYZ position of each antenna in J2000 coordinates, and the J2000 coordinates of the Earth from the Moon’s perspective.  
You also must update lunar_furnsh.txt with the new path names for the required frame and ephemeris files. Update below in lunar_furnsh.txt with the directory you installed the files to.
KERNELS_TO_LOAD = ( '/home/alexhege/SPICE/LunarEph/moon_pa_de421_1900-2050.bpc'
                   '/home/alexhege/SPICE/LunarEph/moon_080317.tf'
                   '/home/alexhege/SPICE/LunarEph/moon_assoc_me.tf'
                   '/home/alexhege/SPICE/LunarEph/pck00010.tpc'
                   '/home/alexhege/SPICE/LunarEph/naif0008.tls'
                   '/home/alexhege/SPICE/LunarEph/de430.bsp' )

Currently, you also must copy paste any new Long Lat Alt data into the C script manually, as I have not implemented a dynamic routine for C that can read in variable length arrays from a file.  Currently has 10 km logarithmic circular array loaded in.
Then compile the script with
gcc eqArrOverTimeEarth.c -o eqArrOverTimeEarth -I/home/alexhege/SPICE/cspice/include /home/alexhege/SPICE/cspice/lib/cspice.a -lm -std=c99
update the included directories to wherever you have installed the SPICE libraries
and run with ./eqArrOverTimeEarth 
which will create the output files.  


3. Run LunarEarthPicFreqIntegration.py
This is the CASA script that does the work of simulating the radio array now that everything is defined.  It reads in the array position data and the position of the target in the sky.
To run within casa run this command:

%run LunarEarthPicFreqIntegration.py  -outDir . -correlate True -numSC 1024

or in terminal
casa --nologger --nologfile --nogui --agg -c LunarEarthPicFreqIntegration.py  -outDir . -correlate True -numSC 1024| tee out.out
nohup casa --nologger --nologfile --nogui --agg -c LunarEarthPicFreqIntegration.py  -outDir . -correlate True -numSC 1024 | tee earth.out &


It reads in the .dat files of the synchrotron emission maps in the EarthSynchrotronMapVsFrequency subfolder and converts them into a CASA ground truth file that is given to the virtual array.  Currently, only the ~0.75 MHz image is used as an average for the 0.5 – 1.0 MHz range.

It uses CASA routines to create a synthetic array using the J2000 data calculated from the SPICE kernel in the previous script.  It then feeds the truth image into the simulated array to create the visibility data and puts it into a Measurement Set (MS) file.  This MS file is the main output of this script.

4. Run noiseCopies.py
This file takes the newly created MS file and creates a noiseless image of the recovered emission.  It then adds noise according to the set integration time and SEFD noise level.  Currently makes images for a range of integration times up to 24 hours and over several robust weighting scheme values.  



