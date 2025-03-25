###
'''
Author: Thomas Joyce

Class: ASTR 400B - Galaxies and Cosmology

Description: Reads the MW_000.txt file and returns the snapshot (time), amount of particles,
and the data encoded in the txt file. 
'''
###


### Importations ###
import numpy as np
import astropy.units as u

### Functions ###
def Read(FileName):
    '''
    Description: Reads the txt file and outputs the returns the data encoded.
        
    Inputs: 
        - FileName: the name of the file (MW_000.txt) where the data is located 
    
    Returns: 
        - Time: snapshot time of model
        - NParticles: total number of particles in the data file
        - data: data array containing the information from the file. Format is data['type][particle number]
    '''
    
    # Opening and retrieving Data #
    file = open(FileName, 'r')
    
    # First two lines of data (time and particle count)
    Time = float(file.readline().split()[1]) * u.Myr
    NParticles = int(file.readline().split()[1])
    
    # Closing File
    file.close()
    
    data = np.genfromtxt(FileName, dtype = None, names = True, skip_header = 3) # data array creation
    
    return Time, NParticles, data
### END Read