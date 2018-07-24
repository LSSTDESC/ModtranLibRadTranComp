from lsst.sims.photUtils import SignalToNoise
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.photUtils import Bandpass,Sed
from Throughputs import Throughputs

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.constants import *

class Telescope(Throughputs):
    def __init__(self,airmass=1,atmos=True,aerosol=True,**kwargs):
        Throughputs.__init__(self,**kwargs)
        
        #self.filters=filterlist 

        params=['mag_sky','m5','FWHMeff','Tb','Sigmab','zp','counts_zp','Skyb','flux_sky']

        self.data={}
        for par in params:
            self.data[par]={}
        
        self.data['FWHMeff']=dict(zip('ugrizy',[0.92,0.87,0.83,0.80,0.78,0.76]))
         
        self.atmos=atmos

        self.Load_Atmosphere(airmass)

        self.Inputs()
        self.Sky()
        self.ZP()
        
        
    @property
    def FWHMeff(self):
        return self.data['FWHMeff']

    @property
    def mag_sky(self):
        return self.data['mag_sky']

    @property
    def m5(self):
        return self.data['m5']

    @property
    def Tb(self):
        return self.data['Tb']

    @property
    def Sigmab(self):
        return self.data['Sigmab']
    @property
    def zp(self):
        return self.data['zp']

    @property
    def ADU_zp(self):
        return self.data['counts_zp']

    @property
    def flux_sky(self):
        return self.data['flux_sky']

    def Inputs(self):
        
        for filtre in self.filterlist:    
            myup=self.Calc_Integ_Sed(self.darksky,self.system[filtre])
            self.data['Tb'][filtre]=self.Calc_Integ(self.atmosphere[filtre])
            self.data['Sigmab'][filtre]=self.Calc_Integ(self.system[filtre])
            self.data['mag_sky'][filtre]=-2.5*np.log10(myup/(3631.*self.Sigmab[filtre]))

    def Sky(self):

        for filtre in self.filterlist:
            self.Calc_m5(filtre)

    def Calc_m5(self,filtre):

        filtre_trans=self.system[filtre]
        wavelen_min, wavelen_max, wavelen_step=filtre_trans.getWavelenLimits(None,None,None)
                    
        bandpass=Bandpass(wavelen=filtre_trans.wavelen, sb=filtre_trans.sb)

        flatSedb = Sed()
        flatSedb.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
        flux0b=np.power(10.,-0.4*self.mag_sky[filtre])
        flatSedb.multiplyFluxNorm(flux0b)
        photParams = PhotometricParameters(bandpass=filtre)
        norm=photParams.platescale**2/2.*photParams.exptime/photParams.gain
        if self.atmos:
            self.data['m5'][filtre]=SignalToNoise.calcM5(flatSedb,self.atmosphere[filtre],self.system[filtre],photParams=photParams,FWHMeff=self.FWHMeff[filtre])
            adu_int= flatSedb.calcADU(bandpass=self.atmosphere[filtre], photParams=photParams)
            self.data['flux_sky'][filtre]=adu_int*norm
        else:
            self.data['m5'][filtre]=SignalToNoise.calcM5(flatSedb,self.system[filtre],self.system[filtre],photParams=photParams,FWHMeff=self.FWHMeff[filtre])
            adu_int= flatSedb.calcADU(bandpass=self.system[filtre], photParams=photParams)
            self.data['flux_sky'][filtre]=adu_int*norm


    def ZP(self):
        
        for filtre in self.filterlist:
            self.ZP_filtre(filtre)
        #print 'zeropoints',self.data['zp'],self.data['counts_zp']
        #self.data['zp']=dict(zip(['u','g','r','i','z','y'],[27.03,28.53,28.27,27.91,27.49,26.78]))

    def ZP_filtre(self,filtre):

        photParams=PhotometricParameters(bandpass=filtre)
        Diameter=2.*np.sqrt(photParams.effarea*1.e-4/np.pi) # diameter in meter
        Cte=3631.*np.pi*Diameter**2*2.*photParams.exptime/4/h/1.e36
        #print('hello Cte',Cte,Diameter,h,photParams.exptime)

        self.data['Skyb'][filtre]=Cte*np.power(Diameter/6.5,2.)*np.power(2.*photParams.exptime/30.,2.)*np.power(photParams.platescale,2.)*np.power(10.,0.4*(25.-self.mag_sky[filtre]))*self.Sigmab[filtre]
            
        Zb=181.8*np.power(Diameter/6.5,2.)*self.Tb[filtre]
        mbZ=25.+2.5*np.log10(Zb) 
        filtre_trans=self.system[filtre]
        wavelen_min, wavelen_max, wavelen_step=filtre_trans.getWavelenLimits(None,None,None)
        bandpass=Bandpass(wavelen=filtre_trans.wavelen, sb=filtre_trans.sb)
        flatSed = Sed()
        flatSed.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
        flux0=np.power(10.,-0.4*mbZ)
        flatSed.multiplyFluxNorm(flux0)
        photParams=PhotometricParameters(bandpass=filtre)
        counts = flatSed.calcADU(bandpass, photParams=photParams) #number of counts for exptime
        self.data['zp'][filtre]=mbZ
        #print 'hello',counts/self.photParams.exptime
        self.data['counts_zp'][filtre]=counts/2.*photParams.exptime

    def Calc_Integ(self,bandpass):
        resu=0.
        dlam=0
        for i,wave in enumerate(bandpass.wavelen):
            if i < len(bandpass.wavelen)-1:
                dlam=bandpass.wavelen[i+1]-wave
                resu+=dlam*bandpass.sb[i]/wave
            #resu+=dlam*bandpass.sb[i]

        return resu

    def Calc_Integ_Sed(self,sed,bandpass,wavelen=None, fnu=None):
      
        use_self = sed._checkUseSelf(wavelen, fnu)
        # Use self values if desired, otherwise use values passed to function.
        if use_self:
            # Calculate fnu if required.
            if sed.fnu is None:
                # If fnu not present, calculate. (does not regrid).
                sed.flambdaTofnu()
            wavelen = sed.wavelen
            fnu = sed.fnu
        # Make sure wavelen/fnu are on the same wavelength grid as bandpass.
        wavelen, fnu = sed.resampleSED(wavelen, fnu, wavelen_match=bandpass.wavelen)
    
        # Calculate the number of photons.
        nphoton = (fnu / wavelen * bandpass.sb).sum()
        dlambda = wavelen[1] - wavelen[0]
        return nphoton * dlambda

    def flux_to_mag(self, flux, band, zp=None):
        if zp is None:
            zp = self.zero_points(band)
        #print 'zp',zp,band
        m = -2.5 * np.log10(flux) + zp
        return m

    def mag_to_flux(self, mag, band, zp=None):
        if zp is None:
            zp = self.zero_points(band)
        return np.power(10., -0.4 * (mag-zp))


    def zero_points(self, band):
        return np.asarray([self.zp[b] for b in band])

    def mag_to_flux_e_sec(self,mag,band,trans,sed):

        #this should be debugged at some point
        photrams=PhotometricParameters(bandpass=band)
        E_per_sec = sed.calcADU(bandpass=trans, photParams=photParams)
        e_per_sec/=exptime/photParams.gain
        return e_per_sec
        
