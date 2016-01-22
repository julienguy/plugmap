#!/usr/bin/env python

import numpy as np
import numpy.random
import pylab
import math
import sys
import string

class OpticalDistortion() :
    def __init__(self,platescale) :

        self.platescale=platescale # has units
        
        # see ~/software/platedesign/trunk/pro/plate/ad2xyfocal.pro
        coef=np.array([-0.000137627, -0.00125238, 1.5447e-09, 
                        8.23673e-08, -2.74584e-13, -1.53239e-12, 
                        6.04194e-18, 1.38033e-17, -2.97064e-23, 
                        -3.58767e-23])
        self.achromatic_distortion_pol=np.poly1d(coef[::-1])

        # see ~/software/platedesign/trunk/pro/plate/apo_rdistort.pro
        mm_per_rad =platescale*180/math.pi
        self.chromatic_distort_radii=np.arcsin(np.linspace(0,90,10)*math.pi/(60*180))*mm_per_rad
        print "RADII=",self.chromatic_distort_radii
        
        self.chromatic_distort_wave=np.array([5300,4000,5500,6000,8000,10000,15350,15950,16550])
        nw=self.chromatic_distort_wave.size
        nr=self.chromatic_distort_radii.size
        
        self.chromatic_distort=np.array([
                [0.,36.26,72.53,108.84,145.18,181.53,217.90,254.29,290.77,327.44],
                [0.,-0.002,-0.003,-0.004,-0.005,-0.005,-0.005,-0.004,-0.002,0.003],
                [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
                [0.,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,-0.001],
                [0.,0.001,0.003,0.003,0.004,0.004,0.004,0.003,0.002,-0.003],
                [0.,0.002,0.004,0.005,0.005,0.005,0.005,0.005,0.003,-0.004],
                [0.,0.003,0.006,0.007,0.008,0.008,0.008,0.008,0.004,-0.006],
                [0.,0.003,0.006,0.008,0.008,0.009,0.009,0.008,0.004,-0.006],
                [0.,0.004,0.006,0.008,0.009,0.009,0.009,0.008,0.004,-0.007]])
        
        # apply scaling
        scale=np.zeros((nr))
        scale[1:]=self.chromatic_distort_radii[1:]/self.chromatic_distort[0,1:]
        self.chromatic_distort[1:] *= scale
        self.chromatic_distort[0]=0.
        # sort wave
        ii=np.argsort(self.chromatic_distort_wave)
        self.chromatic_distort_wave=self.chromatic_distort_wave[ii]
        for j in range(nr) :
            self.chromatic_distort[:,j]=self.chromatic_distort[ii,j]
        # in ad2xyfocal, a reference wavelength of 5000A instead of 5500A is used !!
        ref_distort = np.zeros((nr))
        for j in range(nr) :
            ref_distort[j]=np.interp(5000,self.chromatic_distort_wave,self.chromatic_distort[:,j])
        self.chromatic_distort -= ref_distort
        
        """
        pylab.plot(self.chromatic_distort_wave,self.chromatic_distort[:,-1],"o-")
        
        ww=np.linspace(4000,8000,200)*u.angstrom
        r=self.chromatic_distort_radii[-1]
        dd=np.zeros((ww.size))
        for i in range(ww.size) :
            dd[i]=self.chromatic_distortion(r,ww[i]).to(u.mm).value
        pylab.plot(ww,dd,c="r")
        pylab.show()
        """


    def chromatic_distortion(self,radius,wavelength) : # with radius and wave with units , returns delta r to be added
        i=np.where(self.chromatic_distort_wave>=wavelength)[0]
        if i.size == 0 :
            i=1
        else :
            i=min(max(1,i[0]),self.chromatic_distort_radii.size-1)
        dist1=np.interp(radius,self.chromatic_distort_radii,self.chromatic_distort[i-1])
        dist2=np.interp(radius,self.chromatic_distort_radii,self.chromatic_distort[i])
        dist=np.interp(wavelength,[self.chromatic_distort_wave[i-1],self.chromatic_distort_wave[i]],[dist1,dist2])
        return dist
    
    def distortion(self,radius,wavelength) : 
        return self.achromatic_distortion_pol(radius) + self.chromatic_distortion(radius,wavelength)


what=0

if len(sys.argv)>1 :
    what=string.atoi(sys.argv[1])

# see ./calibrate_distortion.py
refrac_wave=np.array([3000,3500,4000,5000,5400,6000,7000,8000])
refrac_arcsec=np.array([44.166347,43.365612,42.8640697818,42.292551282,42.1507465805,41.990386,41.811009,41.695723])

# see ~/software/platedesign/trunk/pro/plate/apo_rdistort.pro
# at r=145.18 
#distort_wave=np.array([4000,5500,6000,8000,10000])
#distort_arcsec=np.array([-0.005,0.,0.001,0.004,0.005])/217.*3600 # arcsec

numpy.random.seed(12)
rmax=300.
n=60
x=rmax*(2*numpy.random.uniform(size=n)-1.)
y=rmax*(2*numpy.random.uniform(size=n)-1.)
r=np.sqrt(x**2+y**2)
x=x[r<rmax]
y=y[r<rmax]
r=r[r<rmax]
xx=np.linspace(-rmax,rmax,100)
n=x.size

refwave=5400.
d2r=math.pi/180.
alt=85.
pa=-20.
platescale=217.7358 
scalefact=3000. # for display
scale=scalefact/3600.*platescale
pylab.figure(figsize=(8,7.5))
pylab.plot(x,y,"o",markersize=16,color="white")
pylab.plot(xx,np.sqrt(rmax**2-xx**2),"-",color="gray")
pylab.plot(xx,-np.sqrt(rmax**2-xx**2),"-",color="gray")
pylab.plot(xx,0*xx,"--",color="gray")
pylab.xlabel("XFOCAL ~ RA")
pylab.ylabel("YFOCAL ~ Dec")

wave=np.linspace(3600,9000,10)
adr=np.interp(wave,refrac_wave,refrac_arcsec)-np.interp(refwave,refrac_wave,refrac_arcsec)
adrx=scale/np.tan(alt*d2r)*adr*np.sin(pa*d2r)
adry=scale/np.tan(alt*d2r)*adr*np.cos(pa*d2r) 
distortion            = OpticalDistortion(platescale)

pylab.text(200,-250,"SPECTRO #1",fontsize=16)
pylab.text(200,+250,"SPECTRO #2",fontsize=16)

pylab.arrow(0, rmax, 0, rmax*.1, head_width=rmax*0.05, head_length=rmax*0.1, fc='k', ec='k')
pylab.arrow(rmax*np.sin(pa*d2r), rmax*np.cos(pa*d2r), rmax*np.sin(pa*d2r)*0.1, rmax*np.cos(pa*d2r)*0.1, head_width=rmax*0.05, head_length=rmax*0.1, fc='k', ec='k')
pylab.text(0,rmax*1.25,"North",horizontalalignment='center', fontsize=16)
pylab.text(rmax*1.25*np.sin(pa*d2r),rmax*1.25*np.cos(pa*d2r),"Zenith",horizontalalignment='center', fontsize=16)

colors=[]
for w in range(wave.size) :
    colors.append(pylab.cm.RdYlBu(float(wave.size-w)/wave.size))

if what==0 :
    for w in range(wave.size) :
        pylab.plot(x+adrx[w],y+adry[w],"o",markersize=8,color=colors[w],markeredgewidth=0)
if what==1 :
    for w in range(wave.size) :
        rscale = np.zeros((n))
        for i in range(n) :
            rscale[i]=1.+scalefact*(distortion.distortion(r[i],wave[w])-distortion.distortion(r[i],refwave))/r[i]
        
        pylab.plot(x*rscale,y*rscale,"o",markersize=8,color=colors[w],markeredgewidth=0)
if what==2 :
    for w in range(wave.size) :
        rscale = np.zeros((n))
        for i in range(n) :
            rscale[i]=1.+scalefact*(distortion.distortion(r[i],wave[w])-distortion.distortion(r[i],refwave))/r[i]
        
        pylab.plot(x*rscale+adrx[w],y*rscale+adry[w],"o",markersize=8,color=colors[w],markeredgewidth=0)

rscale=1.5
pylab.xlim([-rscale*rmax,rscale*rmax])
pylab.ylim([-rscale*rmax,rscale*rmax])

    
pylab.show()
