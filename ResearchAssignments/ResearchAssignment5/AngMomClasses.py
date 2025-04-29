''' 
Author: Thomas Joyce

Class: ASTR 400B - Galaxies and Cosmology

Research Topic: Modelling Dynamical Friction-based Angular Momemtum Evolution

Version: 0.25 - Working Basic Version of Code

Description: Python outline script to determine the angular momentum evolution of the M31 and M33 galactic merger. 
This script will aim to determine the angular momentum transfer between the M33 and M31 galaxies as a result of 
M33 traveling through M31's dark matter halo. This will essentially modelling the angular momentum effects brought 
upon through the dynamical friction process being simulated in this study. The goal of this study is to create a gif
of M31's angular momentum as a function of radius iterating in time. Ideally you would be able to see a "wake" of angular momentum
brough on upon by the increase of DM particle density behind M33 as it travels in time. The second plot will be how the orbital angular momentum
of M33 changes in time from the dynamical friction slowing the orbit. This study will also take into account the direction of angular momentum as well, but for the purpose
of this outline and assignment I will focus soley on the magnitude. Implementation of direction based analysis will come in a later version.  
'''

### Importations ###
# Numpy and Astropy
import numpy as np
import astropy.units as u
import astropy.constants as const
# FileManagement 
import os
import glob
# Plotting 
import matplotlib.pyplot as plt
# Image Management - Used for Gifs
from PIL import Image
# External (must be in same directory as ResearchAssignment3)
from ReadFile import Read
from CenterOfMassLimited import CenterOfMass


### Classes ###
class HostGalaxy_AngMomEvolution:
    ''' 
    Description: Class to determine the angular momentum evolution of a host galaxy. This class will create plot 1: A gif
    of angular moentum as a function of radius iterating in time. The angular Momentum is calculated relative to radial
    bins that is user definined in the initialize function. 

    Note: Only the dark matter particles are being considered in this class, as it is the dominant effect, overshadowing the other contributions. 
    '''

    def __init__(self, GalaxyName, a, NShells, RadialLimit, SatelliteGalaxy):
        '''
        Description: Initializes an instance of the class. 

        Inputs:
            - GalaxyName: string, name of the galaxy (ex: 'MW', 'M31', 'M33')
            - a: float, scale factor of host gaaxy (in kpc)
            - NShells: int, number of radial shells to consider for angular momentum bins (ex: 10)
            - RadialLimit: float, maximum radial limit of particles to be considered (in kpc)
            - SatelliteGalaxy: object, Satellite Galaxy object (an instance of the SatelliteGalaxy_AngMomEvolution class)
        '''

        # Establishing self parameters from class inputs #
        self.GalaxyName = GalaxyName
        self.a = a 
        self.NShells = NShells
        self.RadialLimit = RadialLimit
        self.SatGal = SatelliteGalaxy # satellite galaxy object

        # Defining relevant COM parameters #
        self.Delta = 0.1 # kpc (tolerance for COM)
        self.VolDec = 1.5 # volume divison factor for iterative COM calculation (set to 1.5 becuase of the Dark Matter Correction, since halo is very large)
    ### END __init__


    def AngularMomentum(self, ParticlePositions, ParticleMass, ParticleVelocity):
        '''
        Description: Computes the angular momentum of one particle given its mass, position vector, and velocity vector. 
        Each vector is relatve to the center of mass frame. Can operate on a list of particles as well. For example:

        Particle Positions = [[x1, y1, z1], [x2, y2, z2], ...]
        Particle Velocity = [[vx1, vy1, vz1], [vx2, vy2, vz2], ...]
        Particle Mass = [m1, m2, ...]

        Inputs:
            - ParticlePositions: list of floats, of particle positions in kpc (3D vector)
            - ParticleMass: list of floats, of particle masses in (1e10 Msun)
            - ParticleVelocity: list of floats, particle velocities in km/s (3D vector)

        Returns:
            - L: list of floats, angular momentum of the particles in 1e10 kpc * Msun * km/s (3D vector)
        '''
        L = ParticleMass[:, np.newaxis] * np.cross(ParticlePositions, ParticleVelocity)

        return L
    ### END AngularMomentum


    def AngularMomentumEvolution(self, start, end, n):
        ''' 
        Description: Determines the radial angular momentum contribution for each shell of the galaxy at each specific time. Saves
        the data to a text file for later referencing. 
            
        Inputs:
            - start: interger, the number corresponding to the first snapshot to be read in
            - end: interger, the number correponding to the last snapshot to be read in
            - n: interger, a number corresponding to the amount of intervals over which the COM will be returned
        '''
        
        
        # Defining Output Filename #
        OutputName =  self.GalaxyName + '_DMAngularMomentum.txt'
        
        snap_ids = np.array(np.arange(start, end, n)) # array of snap designations to be read in
        
        # Check to end function if equal divisions of n are not achieved #
        if abs(end - start) % n != 0:
            print('Interval is incompatible!!!')
            print('Function Ending ...')
            return
        
        AngMom = np.zeros([len(snap_ids), (self.NShells + 1)])

        SatelliteGalaxyPositions = np.zeros([len(snap_ids)]) 

        BinAvergaeAmount = np.zeros([len(snap_ids), self.NShells])

        # Looping through SnapShots #
        for i, snap_id in enumerate(snap_ids):
            # Determining File Name #
            ilbl = '000' + str(snap_id) # adding snap number to the value '000'
            ilbl = ilbl[-3:] # only the last 3 digits being used
            
            HostFileName = '%s_'%(self.GalaxyName) + ilbl + '.txt' # reconstructing filename
            SatFileName = '%s_'%(self.SatGal.SatGalaxyName) + ilbl + '.txt' # reconstructing filename
            
            # Creating Center of Mass Objects #
            gal_COM = CenterOfMass(HostFileName, 1, self.RadialLimit) 
            Satgal_COM = CenterOfMass(SatFileName, 1, self.SatGal.SatRadialLimit) 

            # Determining Satellite Galaxy Position in Host Galaxy Reference Frame #
            SatPos, SatVel = self.SatGal.SatGalaxyVectors(Satgal_COM, gal_COM)
            SatelliteGalaxyPositions[i] = (np.sqrt((SatPos[0]**2 + SatPos[1]**2 + SatPos[2]**2))) # Saving Position magnitude
            
            # Computing center of mass for host galaxy #
            gal_P = gal_COM.COM_P(self.Delta, self.VolDec)
            gal_V = gal_COM.COM_V(gal_P[0], gal_P[1], gal_P[2])
            
            Time = gal_COM.time / 1000 # Time in Gyr of Snapshot

            AngMom[i][0] = Time.value # time is the first index (column)


            ### --- Modifications --- ####

            # COM Reference Frame Vectors (position and velocity) #
            # Position:
            gal_COMX = gal_COM.x - gal_P[0].value # x position within the center of mass frame
            gal_COMY = gal_COM.y - gal_P[1].value # y position within the center of mass frame 
            gal_COMZ = gal_COM.z - gal_P[2].value # z position within the center of mass frame
            # Velocity:
            gal_COMVx = gal_COM.vx - gal_V[0].value # x velocity within the center of mass frame
            gal_COMVy = gal_COM.vy - gal_V[1].value # y velocity within the center of mass frame
            gal_COMVz = gal_COM.vz - gal_V[2].value # z velocity within the center of mass frame

            gal_masses = gal_COM.m # masses of the particles in the galaxy

            # Reorganizing the position and velocity vectors into a 3D vector format for angular momentum calculation #
            COMPositions = np.zeros([len(gal_COMX), 3]) # initializing the position array
            COMVelocities = np.zeros([len(gal_COMX), 3]) # initializing the velocity array
            
            for index in range(len(gal_COMX)):
                COMPositions[index] = [gal_COMX[index], gal_COMY[index], gal_COMZ[index]]
                COMVelocities[index] = [gal_COMVx[index], gal_COMVy[index], gal_COMVz[index]] 
            ###

            COMPositions = np.array(COMPositions) # defining numpy array for easier use
            COMVelocities = np.array(COMVelocities) # defining numpy array for easier use

            # Finding Position Magnitudes #
            COMPositionMagnitudes = np.sqrt((COMPositions[:,0]**2) + (COMPositions[:,1]**2) + (COMPositions[:,2]**2))

            # Looping through each radial shell #
            for shell in range(1, self.NShells + 1):
                thickness = self.RadialLimit / self.NShells # Shell thickness
                # Shell limits
                shell_min = (shell - 1) * thickness
                shell_max = shell * thickness

                # Selecting particles within the current shell #
                shell_indices = np.where((COMPositionMagnitudes >= shell_min) & (COMPositionMagnitudes < shell_max))[0]
                shell_COMPositions = COMPositions[shell_indices]
                shell_COMVelocities = COMVelocities[shell_indices]
                shell_gal_masses = gal_masses[shell_indices]

                # Calculate the angular momentum for the particles in the current shell
                gal_L_shell = self.AngularMomentum(np.array(shell_COMPositions), np.array(shell_gal_masses), np.array(shell_COMVelocities))
                L_tot_shell = np.sum(gal_L_shell, axis=0)  # Total Angular Momentum for the shell

                # Store the total angular momentum for the shell in the AngMom array
                AngMom[i][shell] = np.sqrt((L_tot_shell[0] ** 2) + (L_tot_shell[1] ** 2) + (L_tot_shell[2] ** 2))  # Store the magnitude of the total angular momentum for the shell
                BinAvergaeAmount[i][shell - 1] = len(shell_indices)
            ###

            print('Iteration: ' + str(i+1) + ' / ' + str(len(snap_ids)))
        ###

        # Printing Median number of particles in each shell (for optimizing trends) #
        print(f'Median number of particles in each shell: {np.median(BinAvergaeAmount, axis=0)}')

        # Saving to txt file #
        header = "{:>20s}".format("t") + "".join(["{:>20s}".format(f"L_tot_shell_{i}") for i in range(1, self.NShells + 1)])
        fmt = "%20.3f" * (self.NShells + 1)
        np.savetxt(OutputName, AngMom, fmt=fmt, comments='#', header=header)

        self.SatelliteGalaxyPositions = np.array(SatelliteGalaxyPositions) # Saving Satellite Positions into a self variable
    ### END AngularMomentumEvolution


    def ReadAngMomFile(self, FileName):
        ''' 
        Description: Reads the angular momentum file and creates a data array as a self parameter for later use. 
        A shell Radaii array is also created for plotting purposes.

        Inputs:
            - FileName: string, name of the file to be read in (ex: self.GalaxyName + '_DMAngularMomentum.txt')
        '''
        self.data = np.genfromtxt(FileName, dtype = None, names = True, skip_header = 0) # data array creation

        # Reading the Time and L_tot values from the file #
        self.Time = self.data['t']
        self.L_tot_Distribution = np.array(list(zip(*[self.data['L_tot_shell_' + str(i)] for i in range(1, self.NShells + 1)]))) # reading the entire shell distribution of angular momentum

        # Creating Radial Array (from shell ranges) #
        shell_thickness = self.RadialLimit / self.NShells
        self.ShellRadaii = np.arange((shell_thickness / 2), (self.RadialLimit + (shell_thickness / 2)), shell_thickness) # 1D array of shell radii (kpc)
    ### END ReadAngMomFile


    def PlotAngMom(self, iteration):
        '''
        Description: Plots the angular momentum distribution of the galaxy at a specific time. The plot is saved as a png file.

        Inputs:
            - iteration: int, the iteration number to be plotted. Represents the (time+1) point. (1-indexed)

        Returns:
            - FigName: string, the name of the figure file that was saved, to be combined into a gif.
        '''


        # Figure #
        plt.figure((f'{self.GalaxyName} Evolution-iteration{iteration}'), figsize = (12, 8))

        # Plotting #
        plt.scatter(self.ShellRadaii, self.L_tot_Distribution[(iteration - 1)], color = 'black', s = 150, label = f'{self.GalaxyName} Angular Momenumtum Distribution')
        plt.vlines(self.SatelliteGalaxyPositions[(iteration - 1)], ymin = -10, ymax = max(self.L_tot_Distribution[(iteration - 1)]) + 20, color = 'gray', linestyle = '--', linewidth = 5, label = f'{self.SatGal.SatGalaxyName} Position') # M33 position in the host galaxy reference frame

        # Plot Formatting #
        plt.xlabel('Radial Shell (kpc)')
        plt.ylabel(r'Angular Momentum $(10^{10} kpc \, \dot{M}_{\odot} \, \mathrm{km/s})$')
        plt.title(f'Angular Momentum Distribution of {self.GalaxyName} at {self.Time[(iteration - 1)]} Gyr') # f'Angular Momentum Distribution of {self.GalaxyName} at {time} Gyr'

        plt.grid(axis = 'both', linestyle = '--', linewidth = 0.5, color = 'gray')

        plt.legend()

        # Figure Saving #
        FigName = f'{self.GalaxyName}_AngularMomentumEvolution-T{self.Time[(iteration - 1)]}.png'
        plt.savefig(FigName, dpi = 250)

        return FigName
    ### END PlotAngMom


    def GifGeneration(self):
        '''
        Description: Generates a gif of the angular momentum evolution of the galaxy over time. The gif is created by iterating through the
        saved png files and combining them into a gif. 
        '''

        # Initializing Frames List (To save each frame of the gif) #
        Frames = np.zeros([len(self.Time)]) # empty list to save each frame of the gif

        # Looping through each frame #
        for t in range(len(self.Time)):
            FigName = self.PlotAngMom(t + 1)
            Frames[t] = Image.open(FigName) # opening the image and saving it to the frames list
        ###

        # Saving Gif #
        Frames[0].save((f'{self.GalaxyName}_AngularMomentumEvolution.gif'), format='GIF', append_images=Frames[1:], save_all=True, duration=400, loop=0) # saves the whole thing as a gif
    ### END GifGeneration
### END HostGalaxy_AngMomEvolution



class SatelliteGalaxy_AngMomEvolution:
    ''' 
    Description: Class to determine the orbital angular momentum evolution of a satellite galaxy relative to the host's center of mass frame. 
    This class will create plot 2: A graph of the orbital angular momentum in time. The satellite galaxy will be approximated as a point mass
    made up of the mass from the dark matter particles. 
        
    Note: Only the dark matter particles are being considered in this class, as it is the dominant effect, overshadowing the other contributions. 
    '''

    def __init__(self, SatGalaxyName, Sat_a, SatRadialLimit, HostGalaxyName, Host_a, HostRadialLimit):
        '''
        Description: Initializes and instance of the class. 
        
        Inputs:
            - SatGalaxyName: string, name of the satellite galaxy (ex: 'MW', 'M31', 'M33')
            - Sat_a: float, satellite galaxy scale factor (in kpc)
            - SatRadialLimit: float, maximum radial limit of particles to be considered (in kpc)
            - HostGalaxyName: string, name of the host galaxy (ex: 'MW', 'M31', 'M33')
            - HostradialLimit: float, maximum radial limit of particles to be considered (in kpc)
        '''

        # Establishing self parameters from class inputs #
        self.SatGalaxyName = SatGalaxyName
        self.Sat_a = Sat_a
        self.SatRadialLimit = SatRadialLimit

        self.HostGalaxyName = HostGalaxyName    
        self.Host_a = Host_a
        self.HostRadialLimit = HostRadialLimit

        # Defining relevant COM parameters #
        self.Delta = 0.1 # kpc (tolerance for COM)
        self.VolDec = 1.5 # volume divison factor for iterative COM calculation (set to 1.5 becuase of the Dark Matter Correction, since halo is very large)
    ### END __init__


    def AngularMomentum(self, ParticlePositions, ParticleMass, ParticleVelocity):
        '''
        Description: Computes the angular momentum of one particle given its mass, position vector, and velocity vector. 
        Each vector is relatve to the center of mass frame. THIS FUNCTION ONLY OPERATES ON ONE PARTICLE (THE
        SATELLITE GALAXY IS TREATED AS A POINT PARTICLE)

        Inputs:
            - ParticlePositions: list of floats, of particle position in kpc (3D vector)
            - ParticleMass: float, prticle mass in (1e10 Msun)
            - ParticleVelocity: list of floats, particle velocity in km/s (3D vector)

        Returns:
            - L: list of floats, angular momentum of the particle in 1e10 kpc * Msun * km/s (3D vector)
        '''
        L = ParticleMass * np.cross(ParticlePositions, ParticleVelocity)

        return L
    ### END AngularMomentum


    def SatGalaxyVectors(self, Satgal_COM, Hostgal_COM):
        '''
        Description: Determines the Satellite's position and velocity vector relative to the host's COM reference
        frame. 
        
        Inputs:
            - Satgal_COM: CenterOfMass object for the satellite galaxy
            - Hostgal_COM: CenterOfMass object for the host galaxy
        
        Returns:
            - satPositions: list of floats, position vector of the satellite galaxy in the host galaxy's center of mass frame (3D vector)
            - satVelocities: list of floats, velocity vector of the satellite galaxy in the host galaxy's center of mass frame (3D vector)
        '''

        # Computing center of mass for object #
        Satgal_P = Satgal_COM.COM_P(self.Delta, self.VolDec)
        Satgal_V = Satgal_COM.COM_V(Satgal_P[0], Satgal_P[1], Satgal_P[2])

        Hostgal_P = Hostgal_COM.COM_P(self.Delta, self.VolDec)
        Hostgal_V = Hostgal_COM.COM_V(Hostgal_P[0], Hostgal_P[1], Hostgal_P[2])
            

        # COM Reference Frame Vectors (position and velocity) #
        # Position:
        gal_COMX = Satgal_P[0].value - Hostgal_P[0].value  # x position within the host's center of mass frame
        gal_COMY = Satgal_P[1].value - Hostgal_P[1].value # y position within the host's center of mass frame 
        gal_COMZ = Satgal_P[2].value - Hostgal_P[2].value # z position within the host's center of mass frame
        # Velocity:
        gal_COMVx = Satgal_V[0].value - Hostgal_V[0].value # x velocity within the host's center of mass frame
        gal_COMVy = Satgal_V[1].value - Hostgal_V[1].value # y velocity within the host's center of mass frame
        gal_COMVz = Satgal_V[2].value - Hostgal_V[2].value # z velocity within the host's center of mass frame

        # Reorganizing the position and velocity vectors into a 3D vector format for angular momentum calculation #
        satPositions = np.array([gal_COMX, gal_COMY, gal_COMZ])
        satVelocities = np.array([gal_COMVx, gal_COMVy, gal_COMVz])

        return satPositions, satVelocities
    ### END SatGalaxyVectors


    def OrbitalAngularMomentumEvolution(self, start, end, n):
        ''' 
        Description: Determines the orbital angular momentum of the satellite galaxy relative to the host's COM frame. Saves
        the data to a text file for later referencing. 
            
        Inputs:
            - start: the interger number corresponding to the first snapshot to be read in
            - end: the interget number correponding to the last snapshot to be read in
            - n: an interger number corresponding to the amout of intervals over which the COM will be returned
        '''
        
        
        # Defining Output Filename #
        OutputName =  f'{self.SatGalaxyName}_OrbitalAngularMomentum.txt'
        
        snap_ids = np.array(np.arange(start, end, n)) # array of snap designations to be read in
        
        # Check to end function if equal divisions of n are not achieved #
        if abs(end - start) % n != 0:
            print('Interval is incompatible!!!')
            print('Function Ending ...')
            return
        
        AngMom = np.zeros([len(snap_ids), 5]) # defining output list to store values
        


        # Looping through SnapShots #
        for i, snap_id in enumerate(snap_ids):
            # Determining File Name #
            ilbl = '000' + str(snap_id) # adding snap number to the value '000'
            ilbl = ilbl[-3:] # only the last 3 digits being used
            
            SatFileName = '%s_'%(self.SatGalaxyName) + ilbl + '.txt' # reconstructing filename
            HostFileName = '%s_'%(self.HostGalaxyName) + ilbl + '.txt' # reconstructing filename
            
            # Creating Center of Mass Object #
            Satgal_COM = CenterOfMass(SatFileName, 1, self.SatRadialLimit) # Center of mass object using only DM particles (majority of mass) - Satellite Galaxy
            Hostgal_COM = CenterOfMass(HostFileName, 1, self.HostRadialLimit) # Center of mass object using only DM particles (majority of mass) - Host Galaxy

            Time = Satgal_COM.time / 1000 # Time in Gyr of Snapshot
            SatMass = np.sum(Satgal_COM.m) # mass of the entire galaxy
            
            COMPositions, COMVelocities = self.SatGalaxyVectors(Satgal_COM, Hostgal_COM)

            # Orbital Angular Momentum Calculation #
            Orbital_L = self.AngularMomentum(np.array(COMPositions), np.array(SatMass), np.array(COMVelocities))
            L_mag = np.sqrt((Orbital_L[0] ** 2) + (Orbital_L[1] ** 2) + (Orbital_L[2] ** 2))

            # Saving Output values to data array #
            AngMom[i][0] = Time.value

            AngMom[i][1] = Orbital_L[0]
            AngMom[i][2] = Orbital_L[1]
            AngMom[i][3] = Orbital_L[2]

            AngMom[i][4] = L_mag

            print('Iteration: ' + str(i+1) + ' / ' + str(len(snap_ids)))
        ###

        # Saving to txt file #
        np.savetxt(OutputName, AngMom, fmt = "%18.3f"*5, comments = '#', header = "{:18s}{:18s}{:18s}{:18s}{:18s}".format('t', 'L_x', 'L_y', 'L_z', 'L_mag'))
    ### END AngularMomentumEvolution


    def ReadAngMomFile(self, FileName):
        ''' 
        Description: Reads the orbital angular momentum file and creates data arrays for later referencing. 

        Inputs:
            - FileName: string, name of the file to be read in (ex: self.GalaxyName + '_DMAngularMomentum.txt')
        '''
        self.data = np.genfromtxt(FileName, dtype = None, names = True, skip_header = 0) # data array creation

        # Reading the Time and L_tot values from the file #
        self.Time = self.data['t']
        self.L_x = self.data['L_x']
        self.L_y = self.data['L_y']
        self.L_x = self.data['L_x']
        self.L_mag = self.data['L_mag']
    ### END ReadAngMomFile


    def PlotAngMom(self):
        '''
        Description: Plots the Orbital Angular Momentum as a Function of time and saves the 
        plot as a png image. Time is on the x-axis and Orbital Angular Momentum is on the y-axis. 
        '''

        # Figure #
        plt.figure((f'{self.SatGalaxyName} Orbital Angular Momentum Evolution'), figsize = (12, 8))

        # Plotting the Angular Momentum Magnitudes #
        plt.plot(self.Time, self.L_mag, color = 'black', linewidth = 5)

        # Plot Formatting #
        plt.xlabel('Time (Gyr)')
        plt.ylabel(r'Orbital Angular Momentum $(10^{10} kpc \, \dot{M}_{\odot} \, \mathrm{km/s})$')
        plt.title(f'Orbital Angular Momentum Evolution of {self.SatGalaxyName} Orbiting {self.HostGalaxyName}')

        plt.grid(axis = 'both', linestyle = '--', linewidth = 0.5, color = 'gray')

        plt.legend()

        # Saving Figure #
        FigName = f'{self.SatGalaxyName}OrbitalAngMomEvolution.png'
        plt.savefig(FigName, dpi = 250)
    ### END PlotAngMom
### END SatelliteGalaxy_AngMomEvolution




### MAIN Function ###
def MAIN():
    ''' Executes the Class Functions'''

    
    M33 = SatelliteGalaxy_AngMomEvolution('M33', 20, (100 * 20), 'M31', 45, (10 * 45)) # Creating M33 object

    M31 = HostGalaxy_AngMomEvolution('M31', 45, 30, (15 * 45), M33) # Creating M31 object



    M31.AngularMomentumEvolution(0, 800, 80)
    M31.ReadAngMomFile('M31_DMAngularMomentum.txt')
    M31.GifGeneration()

    #M33.OrbitalAngularMomentumEvolution(0, 800, 40)
    #M33.ReadAngMomFile('M33_OrbitalAngularMomentum.txt')
    #M33.PlotAngMom()
### END MAIN

### Execution ###
MAIN()
# if __name__ == __MAIN__
#   MAIN()


### Optional Add-ons ###
# Make a folder to store all the lil images
# switch into that folder for all the saving
# then switch back out of it to save the gif in the main folder
# then delete the folder and all the images in it (optional) --> make an option to...



# decrease binning for the gif to preoperly capture trend
###