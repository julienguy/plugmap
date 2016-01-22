#!/usr/bin/env python

import numpy as np
import pylab
import math



# the total adr=adr_prefact*sqrt(airmass**2-1)
adr_prefact=1.1 # arcsec (42.8640697818-41.7591303, between 4000A and 7450A)

disto=0.13 # arcsec (at r=200m, from image-heights.txt, between 4000A and 7450A)

# relative aperture loss is ~ exp(-x**2/2/sig**2) for x = offset in arcsec
sig=0.87 # aperture loss coef, for a seeing of 1.5 arcsec and fiber diam of 2 arcsec

d2r=math.pi/180.

pylab.figure()
a=pylab.subplot(1,2,1)
airmass=np.linspace(1.,1.4,50)
spa=0.
r1r2=np.exp(-np.cos(spa*d2r)*4*adr_prefact*np.sqrt(airmass**2-1)*disto/2./sig**2)
a.plot(airmass,r1r2)
a.set_xlabel("airmass")
a.set_ylabel("r1/r2 for SPA=0")

a=pylab.subplot(1,2,2)
spa=np.linspace(-180,180,100)
for airmass in [1.2] :
    r1r2=np.exp(-np.cos(spa*d2r)*4*adr_prefact*np.sqrt(airmass**2-1)*disto/2./sig**2)
    a.plot(spa,r1r2,c="b")
a.set_xlabel("SPA (deg)")
a.set_ylabel("r1/r2 for airmass=1.2")

pylab.show()

