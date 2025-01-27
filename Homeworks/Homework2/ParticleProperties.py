###
'''
Author: Thomas Joyce

Class: ASTR 400B - Galaxies and Cosmology

Description: Obtains the Distance, Velocity magnitudes, and the mass of a specified 
particle of a specified type using the MW_000.txt input file, read by the Read function
of the ReadFile python script.
'''
###


### Importations ###
# Numpy and Astropy
import numpy as np
import astropy.units as u
# External (must be in same directory as ParticleProperties)
from ReadFile import Read 

### Functions ###
def ParticleInfo(FileName, ParticleType, ParticleNumber):
    ''' 
    Description: This function obtains the distance, veloicty magnitude and mass of a desired particle.
    
    Inputs: 
        - FileName: the name of the file (MW_000.txt) to be read by ReadFile.py
        - ParticleType: particle type (1, 2, or 3) specifing wether it is darkmatter, disk star, 
        or a bulge star
        - ParticleNumber: particle number (index) of the desired particle in human indexing (100th)
    
    Returns:
        - Distance: the magnitude of the distance (r = sqrt(x^2 + y^2 + z^2)) from the
        center of mass of the milky way in kpc units. 
        - VelocityMagnitude: the magnitude of the velocity (v_r = sqrt(vx^2 + vy^2 + vz^2)) of the 
        particle in km/s units 
        - M: mass of the particle un the units of solar masses 
    '''

    # Obtaining Data from Input File #
    Time, NParticles, data = Read(FileName)
       
    # Organizing by type 
    TypeIndex = np.where(data['type'] == ParticleType)[0][0]
    
    index = TypeIndex + (ParticleNumber - 1) # Index begins from the type desired + (python indexing correction)
    
    # Data extraction (obtaining variables)
    x = data['x'][index]
    y = data['y'][index]
    z = data['z'][index]
    
    Vx = data['vx'][index]
    Vy = data['vy'][index]
    Vz = data['vz'][index]
    
    # Output Variable Construction #
    M = (data['m'][index] * u.M_sun) # mass of sun in solar masses
    
    Distance = np.around((np.sqrt((x ** 2) + (y ** 2) + (z ** 2)) * u.kpc), 3) # 3D magnitude of distance in kpc
    
    VelocityMagnitude = np.around((np.sqrt((Vx ** 2) + (Vy ** 2) + (Vz ** 2)) * (u.km / u.second)), 3) # 3D magnitude of velocity in km/s
    
    return Distance, VelocityMagnitude, M
### END ParticleInfo

