"""
@author: MoBueh
"""

#%% Libraries
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np
import math
plt.rcParams["figure.dpi"] = 400
plt.rc('axes', axisbelow=True)


#%%
def radprocessor(i_dir, i_dif, sun_h, sun_az, pos_h, pos_az):
# i_dir   -->|            |
# i_dif   -->|            | 
# sun_h   -->| Radiation  |--> i_sol
# sun_az  -->| Processor  |
# pos_h   -->|            | 
# pos_az  -->|            |
    # Required packages: math
    # Based on the height of the sun and the positioning of a component, indicates the resulting radiation on the component 
    # sun_h -> Sun hight in °
    # sun_az -> Sun azimuth in °
    # i_dir -> Direct irradiation on horizontal plane 
    # i_dif -> Diffuse irradiation on horizontal plane
    # pos_az -> azimuthal orientation of the component, 0° -> north, 90° -> east, 180° -> south, 270° -> west
    # pos_h -> alignment in axis of sun elevation, 0° -> horizontal up, 90°-> vertical, 180° -> horizontal down    
    
    # Direct radiation
    # Angle of Incidence cos(zeta)
    aoi = math.sin(math.radians(sun_h))*math.cos(math.radians(pos_h)) + math.cos(math.radians(sun_h))*math.sin(math.radians(pos_h)) * math.cos(math.radians(abs(pos_az-sun_az)))
    if aoi < 0:
        aoi = 0
    if sun_h<=5:
        i_f = 0
    else:
        i_f = i_dir/math.sin(math.radians(sun_h))*aoi
    # Strahlung durch Bodenreflexion
    i_umg = (i_dir+i_dif)*0.5*0.2*(1-math.cos(math.radians(pos_h)))
    # Diffuse Strahlung
    i_d = i_dif * 0.5 *(1+math.cos(math.radians(pos_h)))
    # Summe
    i_sol = i_f + i_umg + i_d
    return i_sol