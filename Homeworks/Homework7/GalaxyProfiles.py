'''
Author: Thomas Joyce

Class: ASTR 400B - Galaxies and Cosmology

Description: Script Containing the MassProfile class for a given galaxy at a specific snap number
'''

### Importations ###
# Numpy and Astropy
import numpy as np
import astropy.units as u
import astropy.constants as const
# External (must be in same directory as Homework5)
from ReadFile import Read 
from CenterOfMass import CenterOfMass


### Classes ### 
class MassProfile:
    ''' Class to compute the mass profile of a galaxy at a specific timestep '''

    def __init__(self, galaxy, snap):
        ''' 
        Description:
            
        Inputs:
            - galaxy: string name for galaxy (MW, M31, M33, etc)
            - snap: interger snapshot number representing a point in time
        '''
        
        # Reconstructing filename #
        ilbl = '000' + str(snap) # adding snap number to the value '000'
        ilbl = ilbl[-3:] # only the last 3 digits being used
        
        self.filename = '%s_'%(galaxy) + ilbl + '.txt' # reconstructing filename
        
        # Defining global variable for galaxy name #
        self.gname = galaxy
        
        # Obtaining Galaxy data #
        self.time, self.total, self.data = Read(self.filename) 
    ### END __init__
    
    
    def MassEnclosed(self, Ptype, R_m):
        ''' 
        Description:
            
        Inputs:
            - Ptype: Particle ype as an interger. 1 = dark matter halo, 2 = disk stars, 3 = bulge stars
            - R_m: array of radai in kpc (NON astropy quantity)
            
        Returns:
            - M_enc: array of enclosed mass points, each index corresponding to the indeces of R_m in mass units of 1e10 M_sun (Astropy Quantity)
        '''
        
        # Defining Output Variable #
        M_enc = np.zeros(len(R_m)) # initialized zero array to keep track of enclosed mass
        
        ### Obtaining Ptype Particle Data #
        self.index = np.where(self.data['type'] == Ptype) # using indeces for a specific particle type
        
        # Reading data to obtain corresponding quantities (with units)
        self.m = (self.data['m'][self.index]) # unitless for now (but in units of 1e10 M_sun)
        
        # Positions (in kpc)
        self.x = (self.data['x'][self.index] * u.kpc)
        self.y = (self.data['y'][self.index] * u.kpc)
        self.z = (self.data['z'][self.index] * u.kpc)
        
        # velocities (in km/s)
        self.vx = (self.data['vx'][self.index] * (u.km / u.second))
        self.vy = (self.data['vy'][self.index] * (u.km / u.second))
        self.vz = (self.data['vz'][self.index] * (u.km / u.second))
        
        
        # Galaxy Center of mass relative to Disk Particles #
        Gal_COM = CenterOfMass(self.filename, 2) # using only disk particles for COM determination
        
        Gal_COM_Pos = Gal_COM.COM_P(0.1) # COM position within (+/- 0.1 Kpc)
        
        # Determining new position vectors for particle elements (with respect to the cetner of mass frame) #
        x_comframe = self.x - Gal_COM_Pos[0]
        y_comframe = self.y - Gal_COM_Pos[1]
        z_comframe = self.z - Gal_COM_Pos[2]
        
        r_comframe = np.sqrt( (x_comframe ** 2) + (y_comframe ** 2) + (z_comframe ** 2) ) # radius vector for all particle elements
        
        # For loop to iterate through radai list (R_m) #
        for rad_i in range(0, len(R_m)):
            LessthanIndex = np.where(r_comframe <= (R_m[rad_i] * u.kpc))
            
            Masses = self.m[LessthanIndex] # all masses less than radius
            SummedMasses = 0 # Initializing summed variable
            for mass in Masses:
                SummedMasses += mass
            
            # Storing Mass Profile at a given radius #
            M_enc[rad_i] = (SummedMasses * 1e10)
        ###
        
        return (M_enc * u.Msun)
    ### END MassEnclosed
            
    
    def MassEnclosedTotal(self, R_m):
        ''' 
        Description:
            
        Inputs:
            - R_m: array of radai in kpc (NON Astropy quantity)
            
        Returns:
            - MassTotal: array of total mass evaluated at each point of R_m in units of 1e10 M_sun (Astropy Quantity)
        '''
        
        # Non-Compatible Galaxies #
        # Galaxies that lack particle type 1, 2, or 3
        # Store as Galaxy name (as used in gname), and a list of the types absent from text file --> ex: {"M33": [3]}
        BadGalaxies = {"M33": [3]}
        
        
        ### Computing Total Mass ###
        if self.gname not in BadGalaxies:
            
            # DarkMatter Masses
            Mass1 = MassProfile.MassEnclosed(self, 1, R_m)
            # Disk Star Masses
            Mass2 = MassProfile.MassEnclosed(self, 2, R_m)
            # Bulge Star Masses 
            Mass3 = MassProfile.MassEnclosed(self, 3, R_m)
            
            MassTotal = Mass1 + Mass2 + Mass3
            
            
        else: # if Galaxy in Bad Galaxies Dictionary
            ParticleTypes = [1,2,3]
            MassTotal = (np.zeros(len(R_m)) * u.Msun) # initilizing output variable
            for t in ParticleTypes:
                if t not in BadGalaxies[self.gname]:
                    MassTotal += MassProfile.MassEnclosed(self, t, R_m)
        ###
        
        return MassTotal
    ### END MassEnclosedTotal
        
        
    def HernquistMass(self,r, a, M_halo):
        ''' 
        Description: Function that defines the Hernquist 1990 Mass
            
        Inputs:
            - r: Galactocentric distance in kpc (Astropy Quantity)
            - a: Scale radius of the Hernquist Profile in kpc (Astropy Quantity)
            - M_halo: Total Halo mass in units of 1e13 M_sun (NON Astropy Quantity)
            
        Returns:
            - Hmass: Hernquist mass derived at the specific radius point r in units of 1e12 Msun (Astropy Quantity)
        '''
    
        Const = M_halo * (1e12) * u.Msun # Correcting units of constant
       
        Var = (r ** 2) / ((a + r) ** 2)
       
        Hmass = Const * Var 
       
        return Hmass
       ### END HernquistMass
    
    
    def CircularVelocity(self, Ptype, R_v):
        ''' 
        Description: Computes the circular speed using a spherical symmetry approximation. 
        V_circ = sqrt(G * M / r)
            
        Inputs:
            - Ptype: Particle ype as an interger. 1 = dark matter halo, 2 = disk stars, 3 = bulge stars
            - R_v: array of radai in kpc (NON Astropy quantity)
            
        Returns:
            - V_circ: The circular speed of evaluated a specific radius given the mass enclosed from MassEnclosed in km/s (Astropy Quantity)
        '''
        
        M = MassProfile.MassEnclosed(self, Ptype, R_v) # determining mass using Mass Enclosed
        
        G = const.G.to(u.kpc * ((u.km ** 2) / (u.s ** 2)) / u.Msun) # Converting G constant to appropriate units
        
        V_circ = np.sqrt((G * M) / R_v) # circular velocity in km/s
        
        return V_circ
    ### END CircularVelocity
    
    
    def TotalCircularVelocity(self, R_v):
        ''' 
        Description: Computes the circular speed using a spherical symmetry approximation FOR ALL PARTICLE TYPES. 
        V_circ = sqrt(G * M / r)
            
        Inputs:
            - R_v: array of radai in kpc (NON Astropy quantity)
            
        Returns:
            - V_circ: The circular speed of evaluated a specific radius given the mass enclosed from MassEnclosed in km/s (Astropy Quantity)
        '''
        
        M = MassProfile.MassEnclosedTotal(self, R_v) # determining mass using Mass Enclosed
        
        G = const.G.to(u.kpc * ((u.km ** 2) / (u.s ** 2)) / u.Msun) # Converting G constant to appropriate units
        
        V_circ = np.sqrt((G * M) / R_v) # circular velocity in km/s
        
        return V_circ
    ### END TotalCircularVelocity
    
    
    def HernquistVCirc(self, r, a, M_halo):
        ''' 
        Description: Computes the circular speed using the hernquist approximation 
        V_circ = sqrt(G * M / r)
            
        Inputs:
            - r: Galactocentric distance in kpc (Astropy Quantity)
            - a: Scale radius of the Hernquist Profile in kpc (Astropy Quantity)
            - M_halo: Total Halo mass in units of 1e13 M_sun (NON Astropy Quantity)
            
        Returns:
            - V_circ: The circular speed of evaluated a specific radius given the mass enclosed from MassEnclosed in km/s (Astropy Quantity)
        '''
        
        M = MassProfile.HernquistMass(self, r, a, M_halo)
        
        G = const.G.to(u.kpc * ((u.km ** 2) / (u.s ** 2)) / u.Msun) # Converting G constant to appropriate units
        
        V_circ = np.sqrt((G * M) / r) # circular velocity in km/s
        
        return V_circ
    ### END HernquistVCirc

### END MassProfile