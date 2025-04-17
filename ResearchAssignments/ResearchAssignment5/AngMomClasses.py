''' 
Author: Thomas Joyce

Class: ASTR 400B - Galaxies and Cosmology

Research Topic: Modelling Dynamical Friction-based Angular Momemtum Evolution

Version: 0.01 - Planning

Description: Python outline script to determine the angular momentum evolution of the MW and M31 galactic merger. 
This script will aim to determine the angular moentum transfer for each particle type (1 = dark matter, 
2 = disk stars, and 3 = bulge stars), as a function of time throughout the merging process. This will essentially 
model the angular momentum effects brought upon through the dynamical friction process being simulated in this study. 
The goal of this code is to produce a 6 panel plot (or a 2 panel plot) displaying the angular momentum magnitude over time
for each particle type of each galaxy. For the 6 panel version it would be 2 x 3 (2 rows and 3 columns) where each row is a 
galaxy (MW or M31) and each column is a particle type (1,2 or 3). In the 2 panel version it has each compnent (1,2 and 3) over plotted for 
the specific galaxy in question (MW or M31). This study will also take into account the direction of angular momentum as well, but for the purpose
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
# External (must be in same directory as ResearchAssignment3)
from ReadFile import Read 
from CenterOfMass import CenterOfMass


### Classes ###
class HostGalaxy_AngMomEvolution:
    ''' 
    
    '''

    def __init__(self, GalaxyName, a, NShells):
        '''
        Description: Initializes an instance of the class. 

        yada yada yada... only DM particles being considered here...

        a = scale factor
        
        '''

        self.GalaxyName = GalaxyName
        self.a = a # Scale factor
        self.NShells = NShells

        # Defining relevant COM parameters #
        self.Delta = 0.1 # kpc (tolerance for COM)
        self.VolDec = 1.5 # volume divison factor for iterative COM calculation (set to 1.5 becuase of the Dark Matter Correction, since halo is very large)

    
    ### END __init__

    def AngularMomentumEvolution(self, start, end, n):
        ''' 
        Description: Designs a txt file in a 7 column format, 
        tracking the time, position, and velocity of a galaxies center of mass
            
        Inputs:
            - galaxy: string name of galxy (ex: "MW")
            - start: the interger number corresponding to the first snapshot to be read in
            - end: the interget number correponding to the last snapshot to be read in
            - n: an interger number corresponding to the amout of intervals over which the COM will be returned
            
        '''
        
        
        # Defining Output Filename #
        OutputName =  self.GalaxyName + '_DMAngularMomentum.txt'
        
        snap_ids = np.array(np.arange(start, end, n)) # array of snap designations to be read in
        
        # Check to end function if equal divisions of n are not achieved #
        if abs(end - start) % n != 0:
            print('Interval is incompatible!!!')
            print('Function Ending ...')
            return
        
        AngMom = np.zeros([len(snap_ids), 4]) # defining output list to store values
        


        # For loop to iterate through snap shots #
        for i, snap_id in enumerate(snap_ids):
            # Determining File Name #
            ilbl = '000' + str(snap_id) # adding snap number to the value '000'
            ilbl = ilbl[-3:] # only the last 3 digits being used
            
            FileName = '%s_'%(self.GalaxyName) + ilbl + '.txt' # reconstructing filename
            
            # Creating Center of Mass Object #
            gal_COM = CenterOfMass(FileName, 1, (self.a * 3)) # Center of mass object using user defined particle type
            
            # Computing center of mass for object #
            gal_P = gal_COM.COM_P(self.Delta, self.VolDec)
            gal_V = gal_COM.COM_V(gal_P[0], gal_P[1], gal_P[2])
            
            Time = gal_COM.time / 1000 # Time in Gyr of Snapshot

            AngMom[i][0] = Time.value # time is the first index (column)


            ### --- Modifications --- ####

            # COM Reference Fram Vectors (position and velocity) #
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
            COMPositions = []
            COMVelocities = []
            for index in range(len(gal_COMX)):
                COMPositions.append([gal_COMX[index], gal_COMY[index], gal_COMZ[index]])
                COMVelocities.append([gal_COMVx[index], gal_COMVy[index], gal_COMVz[index]])

            COMPositionMagnitudes = np.sqrt((COMPositions[:,0]**2) + (COMPositions[:,1]**2) + (COMPositions[:,2]**2))


            for shell in range(1, self.NShells + 1):
                # Thickness of each shell = 3 * a / NShells 
                thickness = (3 * self.a) / self.NShells
                # Shell limits
                shell_min = (shell - 1) * thickness
                shell_max = shell * thickness

                # Select particles within the current shell
                shell_indices = np.where((COMPositionMagnitudes >= shell_min) & (COMPositionMagnitudes < shell_max))[0]
                shell_COMPositions = COMPositions[shell_indices]
                shell_COMVelocities = COMVelocities[shell_indices]
                shell_gal_masses = gal_masses[shell_indices]

                # Calculate the angular momentum for the particles in the current shell
                gal_L_shell = AngularMomentum(np.array(shell_COMPositions), np.array(shell_gal_masses), np.array(shell_COMVelocities))

                L_tot_shell = np.sum(gal_L_shell, axis=0)  # Total Angular Momentum for the shell

                # Store the total angular momentum for the shell in the AngMom array
                AngMom[i][shell] = np.sqrt((L_tot_shell[0] ** 2) + (L_tot_shell[1] ** 2) + (L_tot_shell[2] ** 2))  # Store the magnitude of the total angular momentum for the shell


            print('Iteration: ' + str(i+1) + ' / ' + str(len(snap_ids)))

            # Saving to txt file #
            header = "{:>15s}".format("t") + "".join(["{:>15s}".format(f"L_tot_shell_{i}") for i in range(1, Nshells + 1)])
            fmt = "%15.3f" * (Nshells + 1)

            np.savetxt(OutputName, AngMom, fmt=fmt, comments='#', header=header)

            ''' 
            NOW WE NEED TO IMPLEMENT THE SHELL METHOD!!! 
            Thickness of each shell = 3 * a / NShells 


            Heres probably what is going to hapen... we are going to have to use the galcompositions and determine L_tot for each shell
            then we need to save that to a unique variable somehow... 


            Basically:

            we will hve to make a shell loop and then iterate the Gal_L on the values that fit within eahc shell doing an iteration in i (starting at 1)
            and ((i * thickness) - thickness) to (i * thickness)
                we then grab those values and pass them to gal_L and then to L_tot and save it to the AngMom array, where AngMom have to have NShells + 1

                then we can slap that whole thing into the savetext method

                

                header = "{:>15s}".format("t") + "".join(["{:>15s}".format(f"L_tot_shell_{i}") for i in range(1, Nshells + 1)])
                fmt = "%15.3f" * (Nshells + 1)

                np.savetxt(OutputName, AngMom, fmt=fmt, comments='#', header=header)


            '''



            # Galaxy Angular Momentum Calculation #
            gal_L = AngularMomentum(np.array(COMPositions), np.array(gal_masses), np.array(COMVelocities))

            L_tot = np.sum(gal_L, axis = 0) # Total Angular Momentum



            
            # Storing values to AngMom array #
            AngMom[i][0] = Time.value # time is the first index (column)
            
            AngMom[i][1] = L_tot[0] # x position is the second index (column)
            AngMom[i][2] = L_tot[1] # y position is the third index (column)
            AngMom[i][3] = L_tot[2] # z position is the fourth index (column)
            
            # Print statement for code running visibility #
            print(L_tot)
            print('Iteration: ' + str(i+1) + ' / ' + str(len(snap_ids)))
            
        # Saving to txt file #
            np.savetxt(OutputName, AngMom, fmt = "%15.3f"*4, comments = '#', 
                    header = "{:15s}{:15s}{:15s}{:15s}".format('t', 'L_tot_x', 'L_tot_y', 'L_tot_z'))
    ### END AngularMomentumEvolution

### END HostGalaxy_AngMomEvolution



class SatelliteGalaxy_AngMomEvolution:
    ''' 
    
    '''

    def __init__(self, SatGalaxyName, HostGalaxyName):
        '''
        here i was to treat the satellite as a point particle intilialiing its position 

        remember to treat this as a point mass
        '''

    ### END __init__

### END SatelliteGalaxy_AngMomEvolution
































### Functions ###
def AngularMomentum(ParticlePositions, ParticleMass, ParticleVelocity):
    '''
    Description: Computes the angular momentum of one particle given its mass, position vector, and velocity vector. 
    Each vector is relatve to the center of mass frame. Can operate on a list of particles as well. For example:

    Particle Positions = [[x1, y1, z1], [x2, y2, z2], ...]
    Particle Velocity = [[vx1, vy1, vz1], [vx2, vy2, vz2], ...]
    Particle Mass = [m1, m2, ...]

    Inputs:
        - ParticlePositions: list of floats, of particle positions in kpc (3D vector)
        - ParticleMass: list of floats, of particle masses in Msun
        - ParticleVelocity: list of floats, particle velocities in km/s (3D vector)
    '''
    L = ParticleMass[:, np.newaxis] * np.cross(ParticlePositions, ParticleVelocity)

    return L
### END AngularMomentum
     


# HW 6 - Modified #
def AngularMomentumEvolution(galaxy, ParticleType, start, end, n):
        ''' 
        Description: Designs a txt file in a 7 column format, 
        tracking the time, position, and velocity of a galaxies center of mass
            
        Inputs:
            - galaxy: string name of galxy (ex: "MW")
            - ParticleType: interger, specifying which particle type is used when finding COM (1, 2 or 3)
            - start: the interger number corresponding to the first snapshot to be read in
            - end: the interget number correponding to the last snapshot to be read in
            - n: an interger number corresponding to the amout of intervals over which the COM will be returned
            
        '''
        # Defining relevant COM parameters #
        Delta = 0.1
        VolDec = 2
        if galaxy == 'M33':
            VolDec = 4
        if ParticleType == 1: # Dark Matter Correction, since halo is very large
            VolDec = 1.5
        
        
        # Defining Output Filename #
        OutputName =  galaxy + '_AngularMomentum_P' + str(ParticleType) + '.txt'
        
        snap_ids = np.array(np.arange(start, end, n)) # array of snap designations to be read in
        
        # Check to end function if equal divisions of n are not achieved #
        if abs(end - start) % n != 0:
            print('Interval is incompatible!!!')
            print('Function Ending ...')
            return
        
        AngMom = np.zeros([len(snap_ids), 4]) # defining output list to store values
        
        # For loop to iterate through snap shots #


        for i, snap_id in enumerate(snap_ids):
            # Determining File Name #
            ilbl = '000' + str(snap_id) # adding snap number to the value '000'
            ilbl = ilbl[-3:] # only the last 3 digits being used
            
            FileName = '%s_'%(galaxy) + ilbl + '.txt' # reconstructing filename
            
            # Creating Center of Mass Object #
            gal_COM = CenterOfMass(FileName, ParticleType) # Center of mass object using user defined particle type
            
            # Computing center of mass for object #
            gal_P = gal_COM.COM_P(Delta, VolDec)
            gal_V = gal_COM.COM_V(gal_P[0], gal_P[1], gal_P[2])
            
            Time = gal_COM.time / 1000 # Time in Gyr of Snapshot


            ### --- Modifications --- ####

            # COM Reference Fram Vectors (position and velocity) #
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
            COMPositions = []
            COMVelocities = []
            for index in range(len(gal_COMX)):
                COMPositions.append([gal_COMX[index], gal_COMY[index], gal_COMZ[index]])
                COMVelocities.append([gal_COMVx[index], gal_COMVy[index], gal_COMVz[index]])


            # Galaxy Angular Momentum Calculation #
            gal_L = AngularMomentum(np.array(COMPositions), np.array(gal_masses), np.array(COMVelocities))

            L_tot = np.sum(gal_L, axis = 0) # Total Angular Momentum



            
            # Storing values to AngMom array #
            AngMom[i][0] = Time.value # time is the first index (column)
            
            AngMom[i][1] = L_tot[0] # x position is the second index (column)
            AngMom[i][2] = L_tot[1] # y position is the third index (column)
            AngMom[i][3] = L_tot[2] # z position is the fourth index (column)
            
            # Print statement for code running visibility #
            print(L_tot)
            print('Iteration: ' + str(i+1) + ' / ' + str(len(snap_ids)))
            
        # Saving to txt file #
            np.savetxt(OutputName, AngMom, fmt = "%15.3f"*4, comments = '#', 
                    header = "{:15s}{:15s}{:15s}{:15s}".format('t', 'L_tot_x', 'L_tot_y', 'L_tot_z'))
### END OrbitCOM
    
def ReadAngMomFile(FileName):
    ''' 
    Description: reads orbit txt file generated from OrbitCOM and saves the data
    into variables Time, Gal_pos, and Gal_vec 
    
    Inputs:
        - FileName: File name to be read in (string)
        
    Returns:
        - Time: a list of time values in Gyr (NON Astropy Quantity)
        - L_tot_Vector: an array of total angular momentum vectors in time  ([L_tot_x1, L_tot_y1, L_tot_z1], [L_tot_x2, L_tot_y2, L_tot_z2], ...)
        - L_tot_Magnitude: an array of total angular momentum magnitudes in time ([L_tot_magnitude1], [L_tot_magnitude2], ...) 
    '''
    data = np.genfromtxt(FileName, dtype = None, names = True, skip_header = 0) # data array creation
    
    Time = data['t']
    L_tot_Vector = [ data['L_tot_x'], data['L_tot_y'], data['L_tot_z'] ]
    
    L_tot_Magnitude = []
    for index in range(len(Time)):
        L_tot_Magnitude.append(np.sqrt((L_tot_Vector[0][index] ** 2) + (L_tot_Vector[1][index] ** 2) + (L_tot_Vector[2][index] ** 2)))
    
    return np.array(Time), np.array(L_tot_Vector), np.array(L_tot_Magnitude)
### END ReadOrbitFile
            

# Plotting Function #
def PlotAngMom(Time, MWL_tot_Magnitude, M31L_tot_Magnitude, ParticleType):
    '''
    Description:
    
    Inputs:
        - Time: a list of time values in Gyr
        - L_tot_Magnitude: an array of total angular momentum magnitudes in time ([L_tot_magnitude1], [L_tot_magnitude2], ...) 
        - ParticleType: interger, particle type (1 = dark matter, 2 = disk stars, 3 = bulge stars)
    '''

    # Figure #
    FigureNames = ['DarkMatterAngularMomentum', 'DiskStarAngularMomentum', 'BulgeStarAngularMomentum', 'TotalAngularMomentum']
    plt.figure(FigureNames[ParticleType - 1], figsize = (12, 8))

    # Plotting the Angular Momentum Magnitudes #
    Labels = ['Dark Matter', 'Disk Stars', 'Bulge Stars', 'All Particles']
    plt.plot(Time, MWL_tot_Magnitude, label = ('MW ' + Labels[ParticleType - 1]), color = 'blue')
    plt.plot(Time, M31L_tot_Magnitude, label = ('M31 ' + Labels[ParticleType - 1]), color = 'red')

    plt.xlabel('Time (Gyr)')
    plt.ylabel('Angular Momentum (kpc * M_sun * km/s)')
    plt.title('Angular Momentum of ' + Labels[ParticleType - 1])

    plt.grid(axis = 'both', linestyle = '--', linewidth = 0.5, color = 'gray')

    plt.legend()

    plt.savefig(FigureNames[ParticleType - 1] + '.png')
### END PlotAngMom




### MAIN Function ###
def MAIN():
    ''' '''

    #################### USER INPUTS ####################
    GalaxyName = 'MW' # string, name of galaxy ('MW' or 'M31')
    ParticleType = 1 # interger, particle type (1 = dark matter, 2 = disk stars, 3 = bulge stars)
    #################### ----------- ####################

    
    # Creating Total Angular Momentum Files #
    #AngularMomentumEvolution(GalaxyName, ParticleType, 0, 800, 5)


    # Reading in the Angular Momentum Files #
    MW1_T, L_MW_1v, L_MW_1m = ReadAngMomFile('MW_AngularMomentum_P1.txt')
    #MW2_T, L_MW_2v, L_MW_2m = ReadAngMomFile('MW_AngularMomentum_P2.txt')
    #MW3_T, L_MW_3v, L_MW_3m = ReadAngMomFile('MW_AngularMomentum_P3.txt')

    M31_1_T, L_M31_1v, L_M31_1m = ReadAngMomFile('M31_AngularMomentum_P1.txt')
    M31_2_T, L_M31_2v, L_M31_2m = ReadAngMomFile('M31_AngularMomentum_P2.txt')
    M31_3_T, L_M31_3v, L_M31_3m = ReadAngMomFile('M31_AngularMomentum_P3.txt')


    # Plotting the Angular Momentum #
    PlotAngMom(MW1_T, L_MW_1m, L_M31_1m, 1) # Dark Matter 
    PlotAngMom(MW2_T, L_MW_2m, L_M31_2m, 2) # Disk Stars
    PlotAngMom(MW3_T, L_MW_3m, L_M31_3m, 3) # Bulge Stars

    PlotAngMom(MW1_T, (L_MW_1m + L_MW_2m + L_MW_3m), (L_M31_1m + L_M31_2m + L_M31_3m), 4) # Total
### END MAIN

### Execution ###
MAIN()
# if __name__ == __MAIN__
#   MAIN()


### Questions ###
'''
1) Each particle type likely has a different center of mass point, so comparing the effects of DM versus disk stars 
will not be relative to the same point in space

2) I use a voldec of 1.5 for dark matter particles to avoid a division by zero error, since the DM halo is so big. 
Is this a valid approach or should I jsut increase my tolerence (delta)?

3) Overall just really confused about the data analysis...

--> limit it as a function of scale length, maybe like 3 or so (we did this in the hw)


'''

'''
What I learned:

- Basically I need to find the angular momentum evoltuion of the wake particles (the particles being moved around as a result of M31 traveling through the DM halo). 
--> I can do this through a shell method, maybe 3 shells each of size scale radius a. Then I find the ang mom relative to each shell from the COM point, and I can see how each shells ang mom evolves
in theory the shell containing M31, should be expierencing som slowing or reversal, decreasing overall ang mom. 

- Then I need to consider this DM wake's effect on the satellite galaxy as a whole. 
--> first we put a limiting shell on all the particles of the satellite (say 3a cutoff). 
--> then I find the ang mom contribution of ll those particles relative to the host COM frame.
--> this lets me see how the satellite is slowing and decaying over time. 

Note: we could consider M31 as a point mass, and everything should still be valid to simply everything, but we need to check and compare this approximation. 

'''

'''
We need to augment center of mass to only take into account the limited particles to find COM (cut off after 3 scale lengths)
Setting up the shells is only a few extra bits of code within the loop
'''