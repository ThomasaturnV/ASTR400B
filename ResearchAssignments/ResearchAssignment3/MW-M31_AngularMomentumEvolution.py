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
# External (must be in same directory as ResearchAssignment3)
from ReadFile import Read 






# Steps:
'''
1) find the center of mass reference frame (velocity and position) of each particle type of M31 and Mw at snap 0
2) design a way to loop over each snapshot and record this information to a file (I believe we did a version of this but maybe just for one particle type in the Hw)
3) With the COM r and v for each particle type in each snapshot for M31 and Mw, then we basically have:
    --> the orbital position of both bodies as a function of time
4) with this known we can compute the angular momentum of each particle (as a point source relative to the COM) and add all these contributions together (do this for each type at each snapshot)
    --> L = r x p = m(r x v)

At the end: I will basically have the ang momentum contribution of each particle type at each moment in time and I can plot all these guys together and see how everything evolves!
'''