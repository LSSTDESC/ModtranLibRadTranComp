from lsst.sims.photUtils import SignalToNoise
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.photUtils import Bandpass,Sed
from MyThroughputs import Throughputs

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.constants import *

#=======================================================================================
class Telescope(Throughputs):
    def __init__(self,airmass=1,atmos=True,aerosol=False,libradtran=False,**kwargs):
        Throughputs.__init__(self,**kwargs)
        print("**** Telescope.__init__******")
        #self.filters=filterlist 

        params=['mag_sky','m5','FWHMeff','Tb','Sigmab','zp','counts_zp','Skyb','flux_sky']

        self.data={}
        for par in params:
            self.data[par]={}
        
        self.data['FWHMeff']=dict(zip('ugrizy',[0.92,0.87,0.83,0.80,0.78,0.76]))
         
        self.atmos=atmos
        self.libradtran=libradtran

        self.Load_Atmosphere(airmass)
        self.sed=None
        self.sedAB0=None
        
        self.Inputs()
        self.Sky()
        self.ZP()
        self.Set_SED_AB0()  # set a AB reference source
        
    #-------------------------------------------------------------------------    
    @property
    def FWHMeff(self):
        return self.data['FWHMeff']
    #-------------------------------------------------------------------------
    @property
    def mag_sky(self):
        return self.data['mag_sky']
    #-------------------------------------------------------------------------
    @property
    def m5(self):
        return self.data['m5']
    #-------------------------------------------------------------------------
    @property
    def Tb(self):
        return self.data['Tb']
    #-------------------------------------------------------------------------
    @property
    def Sigmab(self):
        return self.data['Sigmab']
    #-------------------------------------------------------------------------
    @property
    def zp(self):
        return self.data['zp']
    #-------------------------------------------------------------------------
    @property
    def ADU_zp(self):
        return self.data['counts_zp']
    #-------------------------------------------------------------------------
    @property
    def flux_sky(self):
        return self.data['flux_sky']
    #-------------------------------------------------------------------------
    def Inputs(self):
        
        for filtre in self.filterlist:    
            myup=self.Calc_Integ_Sed(self.darksky,self.system[filtre])
            self.data['Tb'][filtre]=self.Calc_Integ(self.atmosphere[filtre])
            self.data['Sigmab'][filtre]=self.Calc_Integ(self.system[filtre])
            self.data['mag_sky'][filtre]=-2.5*np.log10(myup/(3631.*self.Sigmab[filtre]))
    #-------------------------------------------------------------------------
    def Sky(self):

        for filtre in self.filterlist:
            self.Calc_m5(filtre)
    #-------------------------------------------------------------------------
    def Calc_m5(self,filtre):
        """
        Calc_m5(filtre):
        
        Compute m5 or SNR at five sigma
        Tool function implemented by Phillie Gris (IN2P3)
        
        """
        # get telescope passband (no atmosphere)
        filtre_trans=self.system[filtre]
        wavelen_min, wavelen_max, wavelen_step=filtre_trans.getWavelenLimits(None,None,None)
                    
        bandpass=Bandpass(wavelen=filtre_trans.wavelen, sb=filtre_trans.sb)
        # create a Flat sed S_nu from the sky brightness magnitude
        flatSedb = Sed()
        flatSedb.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
        flux0b=np.power(10.,-0.4*self.mag_sky[filtre])
        flatSedb.multiplyFluxNorm(flux0b)
        
        # Get LSST photometric parameters
        photParams = PhotometricParameters(bandpass=filtre)
        norm=photParams.platescale**2/2.*photParams.exptime/photParams.gain
        
        # Use LSST sims (SignalToNoise) to calculate M5 with atmosphere or without atmosphere
        if self.atmos:
            self.data['m5'][filtre]=SignalToNoise.calcM5(flatSedb,self.atmosphere[filtre],self.system[filtre],photParams=photParams,FWHMeff=self.FWHMeff[filtre])
            adu_int= flatSedb.calcADU(bandpass=self.atmosphere[filtre], photParams=photParams)
            self.data['flux_sky'][filtre]=adu_int*norm
        else:
            self.data['m5'][filtre]=SignalToNoise.calcM5(flatSedb,self.system[filtre],self.system[filtre],photParams=photParams,FWHMeff=self.FWHMeff[filtre])
            adu_int= flatSedb.calcADU(bandpass=self.system[filtre], photParams=photParams)
            self.data['flux_sky'][filtre]=adu_int*norm

    #-------------------------------------------------------------------------
    def ZP(self):
        
        for filtre in self.filterlist:
            self.ZP_filtre(filtre)
        #print 'zeropoints',self.data['zp'],self.data['counts_zp']
        #self.data['zp']=dict(zip(['u','g','r','i','z','y'],[27.03,28.53,28.27,27.91,27.49,26.78]))
    #-------------------------------------------------------------------------
    def ZP_filtre(self,filtre):
        """
        ZP_filtre() : Compute zero point in filter band
        - platescale is 0.2 arcsec per pixel
        
        Tool function implemented by Phillie Gris (IN2P3)
        
        """
        # extract parameters from lsst_sims
        photParams=PhotometricParameters(bandpass=filtre)
        # compute Diameter in meters : D=2*R = 2*sqrt(S/pi)
        Diameter=2.*np.sqrt(photParams.effarea*1.e-4/np.pi) # diameter in meter
        # AB flux is 3.6307805477e-20 erg/cm2/s/Hz or 3.6307805477e-16 erg/m2/s/Hz
        #   or 3.6307e-23 J/m2/s/Hz, h=6.626e-34 J.s
        # What is the meaning of this Cte ????
        # especcialy what is 1e36 ?
        Cte=3631.*np.pi*Diameter**2*2.*photParams.exptime/4/h/1.e36
        #print('Telescope::ZP_filtre: hello Cte=',Cte, ' Diam=',Diameter, 'h=',h,' exptime=',photParams.exptime)
        
        # What is the meaning of Skyb ?????
        self.data['Skyb'][filtre]=Cte*np.power(Diameter/6.5,2.)*np.power(2.*photParams.exptime/30.,2.)*np.power(photParams.platescale,2.)*np.power(10.,0.4*(25.-self.mag_sky[filtre]))*self.Sigmab[filtre]
        
        #What is the meaning of Sigmab, Tb, Zb   such Zb is used to calculate zero point 
        Zb=181.8*np.power(Diameter/6.5,2.)*self.Tb[filtre]
        mbZ=25.+2.5*np.log10(Zb) 
        #filter without atmosphere
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
        #print('Telescope::ZP_filtre hello',counts/self.photParams.exptime)
        self.data['counts_zp'][filtre]=counts/2.*photParams.exptime
    #-------------------------------------------------------------------------
    def Calc_Integ(self,bandpass):
        """
        Calc_Integ():
        Compute sum  F(lambda).dlambda/lambda in band (no unit)
             - F(lamba) : pass band
             - lambda   : wavelength
             
             
        Tool function implemented by Phillie Gris (IN2P3)
          
        """
        resu=0.
        dlam=0
        for i,wave in enumerate(bandpass.wavelen):
            if i < len(bandpass.wavelen)-1:
                dlam=bandpass.wavelen[i+1]-wave
                resu+=dlam*bandpass.sb[i]/wave
            #resu+=dlam*bandpass.sb[i]

        return resu
   
    #-------------------------------------------------------------------------
    def Calc_Integ_Sed(self,sed,bandpass,wavelen=None, fnu=None):
        """
        Calc_Integ_Sed(self,sed,bandpass,wavelen=None, fnu=None)
        
        Compute sum of  S_nu*F(lambda).dlambda/lambda in band units in erg/cm2/s/Hz
             - S_nu     : SED in erg/cm2/s/Hz
             - F(lamba) : pass band
             - lambda   : wavelength
             - force to use the SED S_nu (erg/cm2/s/Hz) instead of S_lamba (erg/cm2/s/nm)
             
        Tool function implemented by Phillie Gris (IN2P3)
        """
      
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
        fnu=np.nan_to_num(fnu)  # SDC(29/06/18) reset to 0 out of band where there are nan
        
        # Calculate the number of photons.
        nphoton = (fnu / wavelen * bandpass.sb).sum()
        dlambda = wavelen[1] - wavelen[0]
        
        return nphoton * dlambda
    
    #---------------------------------------------------------------
    def CalcMyMagnitudes(self):
        """
        CalcMyMagnitudes(sed)
        
        - author : Sylvie Dagoret-Campagne
        - affiliation : LAL/IN2P3/CNRS/FRANCE
        - date   : July 4th 2018
        
        Check how LSST Sim compute the magnitudes.
        Compute here with self.Calc_Integ_Sed() and Sed.calcMag()
        Used just for debug purpose
        
        """
        all_mag1=[]
        all_mag2=[]
        #sed.flambdaTofnu()
       
        for i,band in enumerate(self.filterlist):
            filter=self.lsst_atmos[band]
            #phinorm=filter.sbTophi()
            # resample the wavelength each time for the filter
            wl,fnu=self.sed.getSED_fnu()
            wavelen, fnu = self.sed.resampleSED(wl, fnu, wavelen_match=filter.wavelen)
            fnu=np.nan_to_num(fnu)  # SDC(29/06/18) reset to 0 out of band where there are nan
            
            self.sed=Sed(wavelen=wavelen,fnu=fnu,name=self.sed.name)
            mag1=self.sed.calcMag(bandpass=filter,wavelen=wavelen,fnu=fnu) 
            mag2=-2.5*np.log10(self.Calc_Integ_Sed(self.sed,filter))
            all_mag1.append(mag1)
            all_mag2.append(mag2)
            print('CalcMyMagnitudes :: band = {}, mag1= {} , mag2= {}'.format(i,mag1,mag2))
        return np.array(all_mag1), np.array(all_mag2)

    #---------------------------------------------------------------
    def CalcMyZP(self):
        """
        CalcMyZP()
        
        - author : Sylvie Dagoret-Campagne
        - affiliation : LAL/IN2P3/CNRS/FRANCE
        - date   : July 5th 2018
        
        Calculate the Zero Points for all bands.
        
        This calculation should assume 1 second exposure and unit electronic gain
        """   
       
        for i,band in enumerate(self.filterlist):
        #for filtre in self.filterlist:
            filtre=self.lsst_atmos[band]
            
            # parameters
            photParams = PhotometricParameters(bandpass=band)
            Diameter=2.*np.sqrt(photParams.effarea*1.e-4/np.pi) # diameter in meter
            exptime=2*photParams.exptime
            gain=photParams.gain
            
            # lsst sim calculation
            # by definition Zero point is defined for unit gain and unit exposure
            zp1=filtre.calcZP_t(photParams)+2.5*np.log10(gain/exptime)
            
            # my calculation : Zero point should be calculated for unit gain and per second of exposure
    
            # in Jansky divided by J (photon energy E=hc/lambda)
            Snu_Tl_dldivl_AB0=self.Calc_Integ_Sed(self.sedAB0,filtre)
            
            # in photoelectron per meter squared per meters per second
            # h is the Planck constant h=6.626x 10^-34 J.s
            dN_PhEl_AB0=Snu_Tl_dldivl_AB0/h*1e-26*np.pi*Diameter**2/4.
            
            zp2=+2.5*np.log10(dN_PhEl_AB0)
            
            print("CalcMyZP :: band = {}, zp1(lsst_sim) = {}, zp2(me)= {}, deltaZP= {}".format(i,zp1,zp2,zp1-zp2))
           
    #---------------------------------------------------------------
    def CalcMyPhElMagnitudes(self):
        """
        CalcMyElectronagnitudes(sed)
        
        - author : Sylvie Dagoret-Campagne
        - affiliation : LAL/IN2P3/CNRS/FRANCE
        - date   : July 5th 2018
        
        Calculate the instrumental magnitude (Photoelectrn unit) for all bands.
        
        """
       
        all_magPhEl=[]
       
       
        for i,band in enumerate(self.filterlist):
            filter=self.lsst_atmos[band]
            
            #typical parameters of the band
            photParams = PhotometricParameters(bandpass=band)
            Diameter=2.*np.sqrt(photParams.effarea*1.e-4/np.pi) # diameter in meter
            exptime=2*photParams.exptime
            
            
            # resample the wavelength each time for the filter
            wl,fnu=self.sed.getSED_fnu()
            wavelen, fnu = self.sed.resampleSED(wl, fnu, wavelen_match=filter.wavelen)
            fnu=np.nan_to_num(fnu)  # SDC(29/06/18) reset to 0 out of band where there are nan
            
            #this SED_nu is now in Jansky, units of 10-23 erg/cm2/s/Hz
            # 1 erg=10-7 J
            # 1 cm^-2 = 10^4 m^-2
            # we have to multiply the SED_nu by 10-26 to be in J/m2/s/Hz
            self.sed=Sed(wavelen=wavelen,fnu=fnu,name=self.sed.name)
            
            # in Jansky divided by J (photon energy E=hc/lambda)
            Snu_Tl_dldivl=self.Calc_Integ_Sed(self.sed,filter)
            
            # in photoelectron per meter squared per meters per second
            # h is the Planck constant h=6.626x 10^-34 J.s
            dN_el=Snu_Tl_dldivl/h*1e-26*np.pi*Diameter**2/4.*exptime

            mag_el=-2.5*np.log10(dN_el)
            print('CalcMyPhElMagnitudes :: band = {}, mag= {}'.format(i,mag_el))
            
            all_magPhEl.append(mag_el)
            
        return np.array(all_magPhEl)
    #---------------------------------------------------------------
    def CalcMyADUMagnitudes(self):
        """
        CalcMyADUMagnitudes(sed)
        
        - author : Sylvie Dagoret-Campagne
        - affiliation : LAL/IN2P3/CNRS/FRANCE
        - date   : July 5th 2018
        
        Calculate the instrumental magnitude (ADU unit) for all bands.
        
        
        """
       
        all_magADU=[]
       
       
        for i,band in enumerate(self.filterlist):
            filter=self.lsst_atmos[band]
            
            #typical parameters of the band
            photParams = PhotometricParameters(bandpass=band)
            Diameter=2.*np.sqrt(photParams.effarea*1.e-4/np.pi) # diameter in meter
            exptime=2*photParams.exptime
            gain=photParams.gain
            
            
            # resample the wavelength each time for the filter
            wl,fnu=self.sed.getSED_fnu()
            wavelen, fnu = self.sed.resampleSED(wl, fnu, wavelen_match=filter.wavelen)
            fnu=np.nan_to_num(fnu)  # SDC(29/06/18) reset to 0 out of band where there are nan
            
            #this SED_nu is now in Jansky, units of 10-23 erg/cm2/s/Hz
            # 1 erg=10-7 J
            # 1 cm^-2 = 10^4 m^-2
            # we have to multiply the SED_nu by 10-26 to be in J/m2/s/Hz
            self.sed=Sed(wavelen=wavelen,fnu=fnu,name=self.sed.name)
            
            
            # in Jansky divided by J (photon energy E=hc/lambda)
            Snu_Tl_dldivl=self.Calc_Integ_Sed(self.sed,filter)
            
            # in photoelectron per meter squared per meters per second
            # h is the Planck constant h=6.626x 10^-34 J.s
            dN_ADU=Snu_Tl_dldivl/h*1e-26*np.pi*Diameter**2/4.*exptime/gain

            mag_ADU=-2.5*np.log10(dN_ADU)
            mag_ADU2=-2.5*np.log10(self.sed.calcADU(bandpass=filter,photParams=photParams,wavelen=wavelen,fnu=fnu))
            
            print('CalcMyADUMagnitudes :: band = {}, mag1= {}, mag2={}, deltaM={}'.format(i,mag_ADU,mag_ADU2,mag_ADU-mag_ADU2))
            
            all_magADU.append(mag_ADU)
            
        return np.array(all_magADU)
    #---------------------------------------------------------------
    def CalcMyADUMagnitude_filter(self,band):
        """
        CalcMyADUMagnitude_filter(band)
        
        - author : Sylvie Dagoret-Campagne
        - affiliation : LAL/IN2P3/CNRS/FRANCE
        - date   : July 5th 2018
        
        Calculate the instrumental magnitude (ADU unit) for one band.
        
        """
       
        filter=self.lsst_atmos[band]
            
        #typical parameters of the band
        photParams = PhotometricParameters(bandpass=band)
        Diameter=2.*np.sqrt(photParams.effarea*1.e-4/np.pi) # diameter in meter
        exptime=2*photParams.exptime
        gain=photParams.gain
            
            
        # resample the wavelength each time for the filter
        wl,fnu=self.sed.getSED_fnu()
        wavelen, fnu = self.sed.resampleSED(wl, fnu, wavelen_match=filter.wavelen)
        fnu=np.nan_to_num(fnu)  # SDC(29/06/18) reset to 0 out of band where there are nan
            
        #this SED_nu is now in Jansky, units of 10-23 erg/cm2/s/Hz
        # 1 erg=10-7 J
        # 1 cm^-2 = 10^4 m^-2
        # we have to multiply the SED_nu by 10-26 to be in J/m2/s/Hz
        self.sed=Sed(wavelen=wavelen,fnu=fnu,name=self.sed.name)
            
            
        # in Jansky divided by J (photon energy E=hc/lambda)
        Snu_Tl_dldivl=self.Calc_Integ_Sed(self.sed,filter)
            
        # in photoelectron per meter squared per meters per second
        # h is the Planck constant h=6.626x 10^-34 J.s
        dN_ADU=Snu_Tl_dldivl/h*1e-26*np.pi*Diameter**2/4.*exptime/gain

        mag_ADU=-2.5*np.log10(dN_ADU)
        mag_ADU2=-2.5*np.log10(self.sed.calcADU(bandpass=filter,photParams=photParams,wavelen=wavelen,fnu=fnu))
            
        #print('CalcMyADUMagnitude_filter :: band = {}, mag1(me)= {}, mag2(lsst_sim)={}, deltaM={}'.format(band,mag_ADU,mag_ADU2,mag_ADU-mag_ADU2))
                    
        return mag_ADU
    
    
    #---------------------------------------------------------------------
    def CalcMyABMagnitudes(self):
        """
        CalcMyABMagnitudes()
        
        - author : Sylvie Dagoret-Campagne
        - affiliation : LAL/IN2P3/CNRS/FRANCE
        - date   : July 4th 2018
        
        Calculate the magnitude in AB system unit for all bands.
        """
       
        all_magAB=[]
       
       
        for i,band in enumerate(self.filterlist):
            filter=self.lsst_atmos[band]
            
            # resample the wavelength each time for the filter
            wl,fnu=self.sed.getSED_fnu()
            wavelen, fnu = self.sed.resampleSED(wl, fnu, wavelen_match=filter.wavelen)
            fnu=np.nan_to_num(fnu)  # SDC(29/06/18) reset to 0 out of band where there are nan
            
            self.sed=Sed(wavelen=wavelen,fnu=fnu,name=self.sed.name)
            
            mag1=-2.5*np.log10(self.Calc_Integ_Sed(self.sed,filter))
            mag2=-2.5*np.log10(self.Calc_Integ_Sed(self.sedAB0,filter))
            all_magAB.append(mag1-mag2)

            
            mag3=self.sed.calcMag(bandpass=filter,wavelen=wavelen,fnu=fnu)
            
            print('CalcMyABMagnitudes :: band = {}, mag1={} , mag2={} , deltaM(me)={}, mag3(lsst)={}'.format(i,mag1,mag2,mag1-mag2,mag3))
        return np.array(all_magAB)
    #---------------------------------------------------------------------
    def CalcMyABMagnitude_filter(self,band):
        """
        CalcMyABMagnitudes_filter()
        
        - author : Sylvie Dagoret-Campagne
        - affiliation : LAL/IN2P3/CNRS/FRANCE
        - date   : July 4th 2018
        
        Calculate the magnitude in AB system unit for all bands.
        """
       
        
        filter=self.lsst_atmos[band]
            
        # resample the wavelength each time for the filter
        wl,fnu=self.sed.getSED_fnu()
        wavelen, fnu = self.sed.resampleSED(wl, fnu, wavelen_match=filter.wavelen)
        fnu=np.nan_to_num(fnu)  # SDC(29/06/18) reset to 0 out of band where there are nan
            
        self.sed=Sed(wavelen=wavelen,fnu=fnu,name=self.sed.name)
            
        mag1=-2.5*np.log10(self.Calc_Integ_Sed(self.sed,filter))
        mag2=-2.5*np.log10(self.Calc_Integ_Sed(self.sedAB0,filter))
        mag3=self.sed.calcMag(bandpass=filter,wavelen=wavelen,fnu=fnu)
            
        #print('CalcMyABMagnitude_filter :: band = {}, mag1={} , mag2={} , deltaM(me)={}, mag3(lsst)={}'.format(band,mag1,mag2,mag1-mag2,mag3))
        return mag3   
    
    #---------------------------------------------------------------
    def CalcMyABMagnitudesErrors(self):
        """
        CalcMyABMagnitudesErrors(self)
        
        - author : Sylvie Dagoret-Campagne
        - affiliation : LAL/IN2P3/CNRS/FRANCE
        - date   : July 4th 2018
        
        Calculate magnitude errors for all bands
        
        """
        all_magABErr=[]
        
        for i,band in enumerate(self.filterlist):
           
            
            filtre_atm=self.lsst_atmos[band]
            filtre_syst=self.lsst_system[band]
            
            wavelen_min, wavelen_max, wavelen_step=filtre_syst.getWavelenLimits(None,None,None)
            
            photParams = PhotometricParameters(bandpass=band)
            FWHMeff=self.data['FWHMeff'][band]
           
            
            # create a Flat sed S_nu from the sky brightness magnitude
            skysed = Sed()
            skysed.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
            flux0b=np.power(10.,-0.4*self.mag_sky[band])
            skysed.multiplyFluxNorm(flux0b)
                  
            #calcMagError filled according doc    
            magerr=SignalToNoise.calcMagError_sed( \
                self.sed,filtre_atm,skysed,filtre_syst,photParams,FWHMeff,verbose=False)
                                    
            all_magABErr.append(magerr)
        return np.array(all_magABErr)        
    #---------------------------------------------------------------
    def CalcMyABMagnitudesError_filter(self,band,SkyBrightnessMag,FWHMGeom):
        """
        CalcMyABMagnitudesError_filter(self,band,SkyBrightnessMag,FWHMGeom)
        
        - author : Sylvie Dagoret-Campagne
        - affiliation : LAL/IN2P3/CNRS/FRANCE
        - date   : July 5th 2018
        
        Calculate magnitude errors for one band.
        
        Input args:
        - band : filter band
        - SkyBrighnessMag : Sky Brighness Magnitude in the band
        - FWHMGeom : Geometrical PSF in the band
        
        """      
            
        filtre_atm=self.lsst_atmos[band]
        filtre_syst=self.lsst_system[band]
            
        wavelen_min, wavelen_max, wavelen_step=filtre_syst.getWavelenLimits(None,None,None)
       
        #calculation of effective PSF
        FWHMeff=SignalToNoise.FWHMgeom2FWHMeff(FWHMGeom)
        
        photParams = PhotometricParameters(bandpass=band)
       
        # create a Flat sed S_nu from the sky brightness magnitude
        skysed = Sed()
        skysed.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
        flux0b=np.power(10.,-0.4*SkyBrightnessMag)
        skysed.multiplyFluxNorm(flux0b)
                  
        #calcMagError filled according doc    
        mag_err=SignalToNoise.calcMagError_sed(self.sed,filtre_atm,skysed,filtre_syst,photParams,FWHMeff,verbose=False)
                                    
        return mag_err         
            
    #---------------------------------------------------------------
    def Plot_Filter(self):
        plt.figure(figsize=(5,4))
        for i,band in enumerate(self.filterlist):
            filter=self.lsst_atmos[band]
            #phinorm=filter.sbTophi()
            #print('phinorm',phinorm)
            plt.plot(filter.wavelen, filter.sb, 'k:')
            #plt.plot(filter.wavelen, phinorm, 'r.')
        plt.show()
    
    #-------------------------------------------------------------------------
    def flux_to_mag(self, flux, band, zp=None):
        if zp is None:
            zp = self.zero_points(band)
        print('Telescope::flux_to_mag: zp',zp,band)
        m = -2.5 * np.log10(flux) + zp
        return m
    #-------------------------------------------------------------------------
    def mag_to_flux(self, mag, band, zp=None):
        if zp is None:
            zp = self.zero_points(band)
        return np.power(10., -0.4 * (mag-zp))

    #-------------------------------------------------------------------------
    def zero_points(self, band):
        return np.asarray([self.zp[b] for b in band])
    #-------------------------------------------------------------------------
    def mag_to_flux_e_sec(self,mag,band,trans,sed):

        #this should be debugged at some point
        photrams=PhotometricParameters(bandpass=band)
        E_per_sec = sed.calcADU(bandpass=trans, photParams=photParams)
        e_per_sec/=exptime/photParams.gain
        return e_per_sec
    #-------------------------------------------------------------------------
    def Set_SEDAB(self):
        """
        Set AB source :
        Enter the SED in erg/cm2/s/nm,
        """
        
        M0=48.6           # magnitude of a AB source
        S_nu0=10**(-M0/2.5) # flux in erg/cm2/s/Hz : 3.630780547701003e-20
        c=2.99792458e10    # speed of light in cm/s in CGS
        nm_to_cm=1e-7      # conversion nm to cm 
        wavelength=np.arange(300.,1151.,1)
        S_lambda0=S_nu0*c/(nm_to_cm)/wavelength**2   # in erg/cm2/s/nm
        
        
        self.Set_SED(wavel=wavelength,newsed=S_lambda0,name='AB-source')

    #-------------------------------------------------------------------------
    def Set_SED_AB0(self):
        """
        Set AB source :
        Enter the SED in erg/cm2/s/nm,
        """
        
        M0=48.6           # magnitude of a AB source
        S_nu0=10**(-M0/2.5) # flux in erg/cm2/s/Hz : 3.630780547701003e-20
        c=2.99792458e10    # speed of light in cm/s in CGS
        nm_to_cm=1e-7      # conversion nm to cm 
        wavelength=np.arange(300.,1151.,1)
        S_lambda0=S_nu0*c/(nm_to_cm)/wavelength**2   # in erg/cm2/s/nm
        
        
        self.sedAB0=Sed(wavelen=wavelength,flambda=S_lambda0,fnu=None,name='AB0-source')
        self.sedAB0.flambdaTofnu()
        
    #-------------------------------------------------------------------------
    def Set_SED(self,wavel,newsed,name='Pickles'):
        """
        Set_SED(self,wavel,sed):
        Enter the SED in erg/cm2/s/nm,
        """

        self.sed=Sed(wavelen=wavel,flambda=newsed,fnu=None,name=name)
        self.sed.flambdaTofnu()
    #-------------------------------------------------------------------------
    def Plot_SED(self):
        wl,fnu=self.sed.getSED_fnu()
        plt.plot(wl,fnu,'b-')
        plt.xlabel("$\lambda$ (nm)")
        plt.ylabel("Flux in Jansky ($10^{-23}erg/cm^2/s/Hz$")
        title='SED '
        plt.title(title)
        plt.grid(True)
    #------------------------------------------------------------------------
