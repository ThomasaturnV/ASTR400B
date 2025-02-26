### Importations ###
# Numpy and Astropy
import numpy as np
# External (must be in same directory as OrbitCOM)
from CenterOfMass import CenterOfMass

def OrbitCOM(galaxy, start, end, n):
        ''' 
        Description: Designs a txt file in a 7 column format, 
        tracking the time, position, and velocity of a galaxies center of mass
            
        Inputs:
            - galaxy: string name of galxy (ex: "MW")
            - start: the interger number corresponding to the first snapshot to be read in
            - end: the interget number correponding to the last snapshot to be read in
            - n: an interger number corresponding to the amout of intervals over which the COM will be returned
            
        '''
        # Defining relevant COM parameters #
        Delta = 0.1
        VolDec = 2
        if galaxy == 'M33':
            VolDec = 4
        
        
        # Defining Output Filename #
        OutputName = 'Orbit_' + galaxy + '.txt'
        
        snap_ids = np.array(np.arange(start, end, n)) # array of snap designations to be read in
        
        # Check to end function if equal divisions of n are not achieved #
        if abs(end - start) % n != 0:
            print('Interval is incompatible!!!')
            print('Function Ending ...')
            return
        
        orbit = np.zeros([len(snap_ids), 7]) # defining output list to store values
        
        # For loop to iterate through snap shots #
        for i, snap_id in enumerate(snap_ids):
            # Determining File Name #
            ilbl = '000' + str(snap_id) # adding snap number to the value '000'
            ilbl = ilbl[-3:] # only the last 3 digits being used
            
            FileName = '%s_'%(galaxy) + ilbl + '.txt' # reconstructing filename
            
            # Creating Center of Mass Object #
            gal_COM = CenterOfMass(FileName, 2) # Center of mass object using disk particles
            
            # Computing center of mass for object #
            gal_P = gal_COM.COM_P(Delta, VolDec)
            gal_V = gal_COM.COM_V(gal_P[0], gal_P[1], gal_P[2])
            
            Time = gal_COM.time / 1000 # Time in Gyr of Snapshot
            
            # Storing values to orbit array #
            orbit[i][0] = Time.value # time is the first index (column)
            
            orbit[i][1] = gal_P[0].value # x position is the second index (column)
            orbit[i][2] = gal_P[1].value # y position is the third index (column)
            orbit[i][3] = gal_P[2].value # z position is the fourth index (column)
            
            orbit[i][4] = gal_V[0].value # x velocity is the fifth index (column)
            orbit[i][5] = gal_V[1].value # y velocity is the sixth index (column)
            orbit[i][6] = gal_V[2].value # z velocity is the seventh index (column)
            
            # Print statmeent for code running visibility #
            print('Iteration: ' + str(i) + ' / ' + str(len(snap_ids) - 1))
            
        # Saving to txt file #
        np.savetxt(OutputName, orbit, fmt = "%11.3f"*7, comments = '#', 
                   header = "{:10s}{:11s}{:11s}{:11s}{:11s}{:11s}{:11s}".format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
### END OrbitCOM
