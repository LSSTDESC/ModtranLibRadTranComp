import os
import matplotlib.pyplot as plt
from lsst.sims.photUtils import Bandpass
from lsst.sims.photUtils import Sed
import numpy as np
import sys
#from astropy.table import Table
from scipy.interpolate import interp1d

# to simulate atmosphere with libradtran
PATH_LIBRADTRAN='./libradtran'
sys.path.append(PATH_LIBRADTRAN)
from libsimulateTranspLSSTScattAbsAer3 import *

#------------------------------------------------------------------------
# Definition of data format for the atmospheric grid
#-----------------------------------------------------------------------------
WLMIN=300.  # Minimum wavelength : PySynPhot works with Angstrom
WLMAX=1150. # Minimum wavelength : PySynPhot works with Angstrom
WL=np.arange(WLMIN,WLMAX,1) # Array of wavelength in Angstrom
NBWL=len(WL)

atm_pressure = 782.5
atm_pwv=4.0
atm_aer=0
atm_ozone=300.
atm_wl0=500.



#==========================================================================================
class Throughputs(object):
    def __init__(self,**kwargs):
        print("**** Throughputs.__init__****")
        params={}
        params['through_dir']='LSST_THROUGHPUTS_BASELINE'
        params['atmos_dir']='THROUGHPUTS_DIR'
        params['atmos']=True
        params['aerosol']=False
        params['telescope_files']=['detector.dat', 'lens1.dat', 'lens2.dat', 'lens3.dat', 'm1.dat', 'm2.dat', 'm3.dat']
        params['filterlist']='ugrizy'
        params['wave_min']=300.
        params['wave_max']=1150.
        params['libradtran']=False
        for par in ['through_dir','atmos_dir','atmos','aerosol=','telescope_files','filterlist','libradtran']:
            if par in kwargs.keys():
                    params[par]=str(kwargs[par])


        self.throughputsDir = os.getenv(params['through_dir'])
        if os.path.exists(os.path.join(os.getenv(params['atmos_dir']), 'atmos')):
            self.atmosDir = os.path.join(os.getenv(params['atmos_dir']), 'atmos')
        else:
            self.atmosDir = os.getenv(params['atmos_dir'])

        self.telescope_files=params['telescope_files']
        self.filter_files=['filter_'+f+'.dat' for f in params['filterlist']]
        self.wave_min=params['wave_min']
        self.wave_max=params['wave_max']
        
        self.filterlist = params['filterlist']
        self.filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}
        self.filternum = {'u':1, 'g':2, 'r':3, 'i':4, 'z':5, 'y':6}

        self.lsst_std = {}
        self.lsst_system = {}
        self.mean_wavelength={}
        self.lsst_detector = {}
        self.lsst_atmos={}
        self.lsst_atmos_aerosol={}
        self.airmass=-1.
        self.aerosol=params['aerosol']
        self.libradtran=params['libradtran']
        self.Load_System()
        self.Load_DarkSky()
        
        if params['atmos']:
            self.Load_Atmosphere()
        else:
            for f in self.filterlist:  
                self.lsst_atmos[f]= self.lsst_system[f]
                self.lsst_atmos_aerosol[f]=self.lsst_system[f]

        self.lsst_telescope={}
        
        self.Load_Telescope()
        
        self.Mean_Wave()
        
    #-------------------------------------------------------------------------
    @property
    def system(self):
        return self.lsst_system
    #-------------------------------------------------------------------------
    @property
    def telescope(self):
        return self.lsst_telescope
    #-------------------------------------------------------------------------
    @property
    def atmosphere(self):
        return self.lsst_atmos
    #-------------------------------------------------------------------------
    def aerosol(self):
        return self.lsst_atmos_aerosol
    #-------------------------------------------------------------------------
    def Load_System(self):
        
        for f in self.filterlist:
            self.lsst_std[f] = Bandpass()
            #self.lsst_std[f].readThroughput(os.path.join(self.throughputsDir, 'total_'+f+'.dat'))
            self.lsst_system[f] = Bandpass()
            index = [i for i,x in enumerate(self.filter_files) if f+'.dat' in x]
            self.lsst_system[f].readThroughputList(self.telescope_files+[self.filter_files[index[0]]], 
                                                   rootDir=self.throughputsDir,wavelen_min=self.wave_min,wavelen_max=self.wave_max)

            

        print('System loaded from',self.throughputsDir)
    #-------------------------------------------------------------------------
    def Load_Telescope(self):

        for system in self.telescope_files+self.filter_files:
            self.lsst_telescope[system]=Bandpass()
            self.lsst_telescope[system].readThroughputList([system], 
                                                   rootDir=self.throughputsDir,wavelen_min=self.wave_min,wavelen_max=self.wave_max)
    #-------------------------------------------------------------------------
    def Load_DarkSky(self):

        self.darksky = Sed()
        self.darksky.readSED_flambda(os.path.join(self.throughputsDir, 'darksky.dat'))
    #-------------------------------------------------------------------------    
    def Load_Atmosphere(self, airmass=1.2):

        self.airmass=airmass
        if self.airmass > 0.:
            atmosphere = Bandpass()
            
            #official atmospheric simulation
            if not self.libradtran:
                path_atmos=os.path.join(self.atmosDir, 'atmos_%d.dat' %(self.airmass*10))
                if os.path.exists(path_atmos):
                    atmosphere.readThroughput(os.path.join(self.atmosDir, 'atmos_%d.dat' %(self.airmass*10)))
                else:
                    atmosphere.readThroughput(os.path.join(self.atmosDir, 'atmos.dat'))
                self.atmos= Bandpass(wavelen=atmosphere.wavelen, sb=atmosphere.sb)
            else: # with libradtran
                path,thefile=ProcessSimulationaer(airmass,atm_pwv,atm_ozone,atm_wl0,atm_aer)
                #path,thefile=ProcessSimulation(am,pwv,ozone)
                atmdata = np.loadtxt(os.path.join(path,thefile))
                wl = atmdata[:,0]
                atm = atmdata[:,1]  
                func=interp1d(wl,atm,kind='linear') 
                transm=func(WL)
                atmosphere.wavelen=WL
                atmosphere.sb=transm
                self.atmos= Bandpass(wavelen=WL, sb=transm)

            for f in self.filterlist:
                wavelen, sb = self.lsst_system[f].multiplyThroughputs(atmosphere.wavelen, atmosphere.sb)
                self.lsst_atmos[f]= Bandpass(wavelen=wavelen, sb=sb)
                

            if self.aerosol:
                atmosphere_aero = Bandpass()
                atmosphere_aero.readThroughput(os.path.join(self.atmosDir, 'atmos_%d_aerosol.dat' %(self.airmass*10)))
                self.atmos_aerosol= Bandpass(wavelen=atmosphere_aero.wavelen, sb=atmosphere_aero.sb)

                for f in self.filterlist:
                    wavelen, sb = self.lsst_system[f].multiplyThroughputs(atmosphere_aero.wavelen, atmosphere_aero.sb)
                    self.lsst_atmos_aerosol[f]= Bandpass(wavelen=wavelen, sb=sb)
        else:
          for f in self.filterlist:  
              self.lsst_atmos[f]= self.lsst_system[f]
              self.lsst_atmos_aerosol[f]=self.lsst_system[f]
    #---------------------------------------------------------------------------- 
    
    #-------------------------------------------------------------------------    
    def Set_Atmosphere(self, airmass,wavelenatm,transmatm):

        self.airmass=airmass
        if self.airmass > 0.:
            atmosphere = Bandpass()
           
            # interpolation
            func=interp1d(wavelenatm,transmatm,kind='linear') 
            transm=func(WL)
            
            atmosphere.wavelen=WL
            atmosphere.sb=transm
            self.atmos= Bandpass(wavelen=WL, sb=transm)
            
            # loop on filters to update effective throughput
            for f in self.filterlist:
                wavelen, sb = self.lsst_system[f].multiplyThroughputs(atmosphere.wavelen, atmosphere.sb)
                self.lsst_atmos[f]= Bandpass(wavelen=wavelen, sb=sb)
                
    #---------------------------------------------------------------------------- 
    
    
    
    
    
    #-----------------------------------------------------------------------------
    def Plot_Throughputs(self):

        #colors=['b','g','r','m','c',[0.8,0,0]]
        #self.filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}
        style=[',',',',',',',']
        for i,band in enumerate(self.filterlist):
            #plt.plot(self.lsst_std[band].wavelen,self.lsst_std[band].sb,linestyle='-',color=self.filtercolors[band], label='%s - std' %(band))
            plt.plot(self.lsst_system[band].wavelen,self.lsst_system[band].sb,linestyle='--',color=self.filtercolors[band], label='%s - syst' %(band))
            plt.plot(self.lsst_atmos[band].wavelen,self.lsst_atmos[band].sb,linestyle='-.',color=self.filtercolors[band], label='%s - syst+atm' %(band))
            if len(self.lsst_atmos_aerosol) > 0:
                plt.plot(self.lsst_atmos_aerosol[band].wavelen,self.lsst_atmos_aerosol[band].sb,linestyle='-',color=self.filtercolors[band], label='%s - syst+atm+aero' %(band))

        plt.plot(self.atmos.wavelen, self.atmos.sb, 'k:', label='X =%.1f atmos' %(self.airmass),linestyle='-')
        if len(self.lsst_atmos_aerosol) > 0:
            plt.plot(self.atmos_aerosol.wavelen, self.atmos_aerosol.sb, 'k:', label='X =%.1f atm+aero' %(self.airmass),linestyle='--')
        plt.legend(loc=(0.85, 0.1), fontsize='smaller', fancybox=True, numpoints=1)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Sb (0-1)')
        plt.title('System throughput')
        plt.grid(True)
        plt.ylim(0.,1.)
        plt.show()
    #-------------------------------------------------------------------------
    def Plot_DarkSky(self):

        plt.plot(self.darksky.wavelen,self.darksky.flambda,'k:',linestyle='-')
        plt.xlabel('Wavelength (nm)')                                                                                   
        plt.ylabel('flambda (ergs/cm$^2$/s/nm)')                                                                                          
        plt.title('Dark Sky SED')    
        plt.grid(True)
        plt.show()
    #-------------------------------------------------------------------------
    def Plot_Throughputs_Spectrum(self,wavelength,fluxes,z):
        
        fig, ax1 = plt.subplots()
        style=[',',',',',',',']

        for i,band in enumerate(self.filterlist):
            #print 'Mean wave',band,np.mean(self.lsst_std[band].wavelen)
            #plt.plot(self.lsst_std[band].wavelen,self.lsst_std[band].sb,linestyle='-',color=self.filtercolors[band], label='%s - std' %(band))
            #plt.plot(self.lsst_system[band].wavelen,self.lsst_system[band].sb,linestyle='--',color=self.filtercolors[band], label='%s - system' %(band))
            plt.plot(self.lsst_system[band].wavelen,self.lsst_system[band].sb,linestyle='-',color=self.filtercolors[band], label='%s - system' %(band))
            #plt.plot(self.lsst_atmos[band].wavelen,self.lsst_atmos[band].sb,linestyle='-.',color=self.filtercolors[band], label='%s - system+atmos' %(band))
            #ax1.plot(self.lsst_atmos_aerosol[band].wavelen,self.lsst_atmos_aerosol[band].sb,linestyle='-',color=self.filtercolors[band], label='%s - system+atmos+aerosol' %(band))
        plt.grid(True)
    #-------------------------------------------------------------------------
    def Mean_Wave(self):

        """
        thedir='../Cosmo_Maf/Instruments/Landolt/'
        names=dict(zip(['B','I','R','U','V'],['sb_-41A.dat','si_-25A.dat','sr_-21A.dat','sux_modified.dat','sv_-27A.dat']))
        for band in 'BRIUV':
            t = Table.read(thedir+names[band],format='ascii')
            print band,np.min(t['col1']),np.max(t['col1']),np.mean(np.max(t['col1'])-np.min(t['col1']))
        """
        
        if 'g' in self.filterlist:
            for band in self.filterlist:
                self.mean_wavelength[band]=np.sum(self.lsst_atmos[band].wavelen*self.lsst_atmos[band].sb)/np.sum(self.lsst_atmos[band].sb)
        else:
            for band in self.filterlist:
                self.mean_wavelength[band]=0. 
    #-------------------------------------------------------------------------