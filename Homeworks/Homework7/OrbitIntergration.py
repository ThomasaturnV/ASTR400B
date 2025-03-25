'''
Author: Thomas Joyce

Class: ASTR 400B - Galaxies and Cosmology

Description: Script Containing the M33 Analytic Orbit class, to determine the galaxies utilizing
the leapfrog method
'''

### Importations ###
# Numpy and Astropy
import numpy as np
import astropy.units as u
import astropy.constants as const
# Plotting
import matplotlib.pyplot as plt
# Latex Equation Display
from IPython.display import Latex
# External (must be in same directory as Homework7)
from ReadFile import Read 
from CenterOfMass import CenterOfMass
from GalaxyProfiles import MassProfile
from GalaxyMass import ComponentMass



### Classes ### 
class M33AnalyticOrbit:
    '''
    Description: Calculates the analytical orbit of M33 around M31
    '''
    
    def __init__(self, OutputFileName): 
        ''' 
        Description: Initializes and instance of the class by defining many self parameters
        that will be used in calculations later on based upon the initial conditions at
        timestep (snpashot) 0. 
            
        Inputs:
            - OutputFileName: string, filename for the output file, describing the orbital positions
            of M33 relative to M31. Suggested that you use a '.txt' file. 
        '''

        # Gravitational Constant Conversion (the value is 4.498502151575286e-06 kpc^3 * Msun / Gyr^2) #
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        self.OutputFileName = OutputFileName
        
        
        ### Initial Position and Velocity Vectors for Disk (Snapshot 0) ###
        # M33 #
        DiskCOM_M33 = CenterOfMass('M33_000.txt', 2)

        DiskCOM_M33_pos = DiskCOM_M33.COM_P(0.1, 2) # in kpc
        DiskCOM_M33_vec = DiskCOM_M33.COM_V(DiskCOM_M33_pos[0], DiskCOM_M33_pos[1], DiskCOM_M33_pos[2]).value # unitless
        
        # M31 #
        DiskCOM_M31 = CenterOfMass('M31_000.txt', 2)
        
        DiskCOM_M31_pos = DiskCOM_M31.COM_P(0.1, 2) # in kpc

        DiskCOM_M31_vec = DiskCOM_M31.COM_V(DiskCOM_M31_pos[0], DiskCOM_M31_pos[1], DiskCOM_M31_pos[2]).value # unitless
        
        
        # Initial Relative Position of M33 with respect to M31 #
        SeperationX = DiskCOM_M33_pos[0] - DiskCOM_M31_pos[0]
        SeperationY = DiskCOM_M33_pos[1] - DiskCOM_M31_pos[1]
        SeperationZ = DiskCOM_M33_pos[2] - DiskCOM_M31_pos[2]
        
        self.r0 = np.array([SeperationX.value, SeperationY.value, SeperationZ.value]) # unitless
        
        # Initial Relative Velocity of M33 with respect to M31 #
        RelVecX = DiskCOM_M33_vec[0] - DiskCOM_M31_vec[0]
        RelVecY = DiskCOM_M33_vec[1] - DiskCOM_M31_vec[1]
        RelVecZ = DiskCOM_M33_vec[2] - DiskCOM_M31_vec[2]
        
        self.v0 = np.array([RelVecX, RelVecY, RelVecZ]) # unitless
        
        
        ### Mass Components for each galaxy (using snapshot 0 as a reference) ###
        self.rdisk = 5 # disk scale length in kpc
        self.Mdisk = (ComponentMass('M31_000.txt', 2) * 1e12).value
        
        
        self.rbulge = 1 # bulge scale length in kpc
        self.Mbulge = (ComponentMass('M31_000.txt', 3) * 1e12).value
        

        self.rhalo = 45 # halo scale length in kpc
        self.Mhalo = (ComponentMass('M31_000.txt', 1) * 1e12).value
    ### END __init__
    
    
    def HernquistAccel(self, M, r_a, r):
        '''
        Description: Determines the Hernquist acceleration for the halo and bulge if we 
        assume a hernquist profile
            
        Inputs:
            - M: float, Mass of halo/bulge component
            - r_a: float, scale length of halo/bulge component
            - r: list of floats, 3 vector of position of M33 relative to M31 (ex: [x_pos, y_pos, z_pos])
        
        Returns:
            - Hern_a: list of floats, 3 vector hernquist acceleration derived at a specific r (ex: [x_acc, y_acc, z_acc])
        '''
        
        # Position Magnitude #
        rmag = np.sqrt((r[0] ** 2) + (r[1] ** 2) + (r[2] ** 2))
        
        # Hernquist Acceleration #
        Hern_a =  -self.G*M/(rmag * (r_a + rmag) ** 2) * r 
        
        return Hern_a
    ### END HernquistAccel
    
    
    
    def MiyamotoNagaiAccel(self, M, r_d, r):
        '''
        Description: Determines the acceleration for the disk using a disk profile derived
        from the Miyamoto-Nagai 1975 article. 
            
        Inputs:
            - M: float, Mass of disk component
            - r_d: float, scale length of disk component
            - r: list of floats, 3 vector of position of M33 relative to M31 (ex: [x_pos, y_pos, z_pos])
        
        Returns:
            - MN_a: list of floats, 3 vector Miyamoto-Nagai acceleration derived at a specific r (ex: [x_acc, y_acc, z_acc])
        '''
        
        R = np.sqrt((r[0] ** 2) + (r[1] ** 2)) # R in cylindrical coordinates, radial component
        B = r_d + np.sqrt((r[2] ** 2) + ((self.rdisk / 5) ** 2)) # self consistent constant for Miyamoto-Nagai acceleration
        
        # Z component has an additional scaling factor in z that we must account for (doing so as an array)
        ZComponentScaling = np.array([1, 1, (B / np.sqrt((r[2] ** 2) + ((self.rdisk / 5) ** 2)))])
        
        # Determining Acceleration
        MN_a = ((-self.G * M) / (((R ** 2) + (B ** 2)) ** 1.5)) * r * ZComponentScaling
        
        return MN_a    
    ### END MiyamotoNagaiAccel
     
    
    def M31Accel(self, Position): # input should include the position vector, r
        ''''
        Description: Determines the acceleration contribution from each component of the galaxy
        by factoring in the hernquist acceleration for the halo and bulge as well as the 
        Miyamoto-Nagai acceleration from the disk. The total acceleration is the resulting vector 
        addition of each component.
            
        Inputs:
            - Position: list of floats, 3 vector of position of M33 relative to M31 (ex: [x_pos, y_pos, z_pos])
        
        Returns:
            - Acceleration: list of floats, 3 vector total acceleration derived at a specific position (ex: [x_acc, y_acc, z_acc])
        '''

        BulgeAcc = M33AnalyticOrbit.HernquistAccel(self, self.Mbulge, self.rbulge, Position)
        
        HaloAcc = M33AnalyticOrbit.HernquistAccel(self, self.Mhalo, self.rhalo, Position)
        
        DiskAcc = M33AnalyticOrbit.MiyamotoNagaiAccel(self, self.Mdisk, self.rdisk, Position)
        
        Acceleration = BulgeAcc + HaloAcc + DiskAcc
        
        return Acceleration
    ### END M31Accel
    
    
    def LeapFrog(self, Position, Velocity, dt): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        ''''
        Description: 
            
        Inputs:
            - Position: list of floats, 3 vector of position of M33 relative to M31 (ex: [x_pos, y_pos, z_pos])
            - Velocity: list of floats, 3 vector of velocity of M33 relative to M31 (ex: [x_vec, y_vec, z_vec])
            - dt: float, timestep increment in Gyr
        
        Returns:
            - Acceleration: list of floats, 3 vector total acceleration derived at a specific position (ex: [x_acc, y_acc, z_acc])
        '''
        
        # Predicting position at next half timestep #
        rhalf = Position + (Velocity * (dt / 2))
        
        # predicting the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = Velocity + (self.M31Accel(rhalf) * dt)
        
        # Implementing Leadpfrog intregration method to find the position at the full timestep
        # by implementing half timestep knowledge to discretize movements
        rnew = Position + (0.5 * (vnew + Velocity) * dt)
        
        return rnew, vnew
    ### END LeapFrog
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        ''''
        Description: Performs repeated leapfrog integration iterations to 
        calculate the orbit of M33 relative to M31. Keep in mind that this method will
        consider M33 as a point source, not an extended body with rotation!
        
        Saves the resulting orbit determination in an output txt file with the 
        name of self.OutputFileName, initialized in __init__. 
            
        Inputs:
            - t0: float, initial time in intergration series in Gyr
            - dt: float, timestep increment in Gyr
            - tmax: float, final time in intergration series in Gyr
        '''

        # Initializing iterable time as the first point in the series #
        t = t0
        
        # Initializing an empty array of size :  rows: int(tmax/dt)+1  , columns: 7
        orbit = np.zeros( ((int(tmax/dt) + 1), 7) )
        
        # First row of orbit output
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        # Initializing iterable position and velocity as the first point in the series #
        r_i, v_i = self.r0, self.v0
        
        # Initializing a counter for the orbit index point.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # Integration #
        while t < tmax:
            t += dt # advancing timestep
        
            r_i, v_i = M33AnalyticOrbit.LeapFrog(self, r_i, v_i, dt)
         
            orbit[i] = t, *tuple(r_i), *tuple(v_i)
            
            i += 1
        ###
        
        
        # Saving Data as an output file #
        np.savetxt(self.OutputFileName, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
    ### END OrbitIntegration
    
### END M33AnalyticOrbit

