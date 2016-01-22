#!/usr/bin/env python

import numpy as np
import pylab
import math

# see ./calibrate_distortion.py
refrac_wave=np.array([3000,3500,4000,5000,5400,6000,7000,8000])
refrac_arcsec=np.array([44.166347,43.365612,42.8640697818,42.292551282,42.1507465805,41.990386,41.811009,41.695723])

# see ~/software/platedesign/trunk/pro/plate/apo_rdistort.pro
# at r=145.18 
distort_wave=np.array([4000,5500,6000,8000,10000])
distort_arcsec=np.array([-0.005,0.,0.001,0.004,0.005])/217.*3600 # arcsec

if True :
    pylab.figure()
    pylab.plot(refrac_wave,refrac_arcsec,"o-")
    pylab.xlabel("Wavelength (A)")
    pylab.ylabel(r"$\alpha$ (arcsec)")
    pylab.figure()
    pylab.plot(distort_wave,distort_arcsec,"o-")
    pylab.xlabel("Wavelength (A)")
    pylab.ylabel(r"$\delta$ (arcsec)")

design_guide_wave=5400.
obs_guide_wave=5000.
obs_wave_b=4900.
obs_wave_r=7500.

distortion_scale=1.5
sig=0.86 # aperture loss coef, for a seeing of 1.5 arcsec and fiber diam of 2 arcsec
#sig=0.78 # aperture loss coef, for a seeing of 1.5 arcsec and fiber diam of 1.8 arcsec
alt=70. # deg

print "Design guide wave    = %d A"%design_guide_wave
print "Obs. guide wave      = %d A"%obs_guide_wave
print "Obs. target wave (b) = %d A"%obs_wave_b
print "Obs. target wave (r) = %d A"%obs_wave_r
print "Aperture loss sigma  = %f arcsec"%sig
print "Distortion scaling   = %f"%distortion_scale
print "Altitude             = %f deg"%alt

for design_wave in [4000.,5400.] :
    
# the total adr=adr_prefact*sqrt(airmass**2-1)


    adr_prefact_b=-np.interp(obs_wave_b,refrac_wave,refrac_arcsec)+np.interp(design_wave,refrac_wave,refrac_arcsec)+np.interp(obs_guide_wave,refrac_wave,refrac_arcsec)-np.interp(design_guide_wave,refrac_wave,refrac_arcsec)  # arcsec
    adr_prefact_r=-np.interp(obs_wave_r,refrac_wave,refrac_arcsec)+np.interp(design_wave,refrac_wave,refrac_arcsec)+np.interp(obs_guide_wave,refrac_wave,refrac_arcsec)-np.interp(design_guide_wave,refrac_wave,refrac_arcsec)  # arcsec

    disto_b=distortion_scale*np.interp(obs_wave_b,distort_wave,distort_arcsec)-np.interp(design_wave,distort_wave,distort_arcsec)-distortion_scale*np.interp(obs_guide_wave,distort_wave,distort_arcsec)+np.interp(design_guide_wave,distort_wave,distort_arcsec) # arcsec
    disto_r=distortion_scale*np.interp(obs_wave_r,distort_wave,distort_arcsec)-np.interp(design_wave,distort_wave,distort_arcsec)-distortion_scale*np.interp(obs_guide_wave,distort_wave,distort_arcsec)+np.interp(design_guide_wave,distort_wave,distort_arcsec) # arcsec

    disto_b*=distortion_scale
    disto_r*=distortion_scale

    print ""
    print "Design target wave   = %d A"%design_wave
    print "------------------------------"
    print "ADR prefactor (b)    = %f arcsec"%adr_prefact_b
    print "Distortion (b)       = %f arcsec"%disto_b
    print "ADR prefactor (r)    = %f arcsec"%adr_prefact_b
    print "Distortion (r)       = %f arcsec"%disto_b

    # ratio = exp(-(adr*cpa+distort)**2/2*sig**2)/exp(-(+adr*cpa-distort)**2/2*sig**2)
    #       = exp( -2*adr*cpa*distort/sig**2 )
    # max delta ratio = exp( 4*adr*cpa*distort/sig**2 )

    d2r=math.pi/180.
    delta_r1r2=np.exp(4*adr_prefact_r/np.tan(alt*d2r)*disto_r/sig**2)
    delta_b1b2=np.exp(4*adr_prefact_b/np.tan(alt*d2r)*disto_b/sig**2)

    print "max delta b1/b2      = %f"%delta_b1b2
    print "max delta r1/r2      = %f"%delta_r1r2




if False :
    adr_prefact = adr_prefact_r
    disto       = disto_r
    pylab.figure()
    a=pylab.subplot(1,2,1)
    airmass=np.linspace(1.,1.4,50)
    pa=0.
    r1r2=np.exp(-np.cos(pa*d2r)*2*adr_prefact*np.sqrt(airmass**2-1)*disto/sig**2)
    a.plot(airmass,r1r2)
    a.set_xlabel("airmass")
    a.set_ylabel("r1/r2 for SPA=0")

    a=pylab.subplot(1,2,2)
    pa=np.linspace(-180,180,100)
    for airmass in [1.2] :
        r1r2=np.exp(-np.cos(pa*d2r)*2*adr_prefact*np.sqrt(airmass**2-1)*disto/sig**2)
        a.plot(pa,r1r2,c="b")
    a.set_xlabel("PA (deg)")
    a.set_ylabel("r1/r2 for airmass=1.2")

    pylab.show()


