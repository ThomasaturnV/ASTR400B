###
'''
Author: Thomas Joyce

Class: ASTR 400B - Galaxies and Cosmology

Description: Determines the total mass of all galaxy components of a given
type by utilizing the ComponentMass function. 
'''
###


### Importations ###
# Numpy and Astropy
import numpy as np
import astropy.units as u
# External (must be in same directory as GalaxyMass)
from ReadFile import Read 

### Functions ###
def ComponentMass(FileName, ParticleType):
    '''
    Description: Iterates through all particles of a given type to find the total
    mass of all particles. Returns the total mass of a gien type in units
    of 1e12 solar masses. 
    
    Inputs:
        - FileName: string, Filename to be read in by the Read function
        - ParticleType: int, particle type to be calculated (1 = DM Halo, 2 = Disk, 3 = Bulge)
        
    Returns:
        - TotalMass: float, total mass of a given type in units of 1e12 solar masses (astropy quantity)
    '''
    
    # Defining Output variable #
    TotalMass = 0 * u.M_sun
    
    # Obtaining Data from Input File #
    Time, NParticles, data = Read(FileName) # read function to obtain data
       
    # Organizing by type #
    TypeIndex = np.where(data['type'] == ParticleType)[0][0] # finding beginning index of a particular type
    
    index = TypeIndex # iterative variable initilized 
    
    # Iterating to find the total mass #
    while (data['type'][index] == ParticleType): 
        M = (data['m'][index] * 1e10 * u.M_sun) # mass of particle in solar masses
        TotalMass += M # Adding to total mass
        
        index += 1
        
        if index == NParticles: # special condition for type 3 particles reaching the end of the text file
            break # breaks out of while loop
    
    # Formatting Total Mass #
    TotalMass /= 1e12 # Converting to a base unit of 1e12 Solar masses
    TotalMass = np.round(TotalMass, 3) # rounding answer to 3 decimal places
    
    return TotalMass
### END ComponentMass