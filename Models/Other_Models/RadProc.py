# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 15:20:28 2022

@author: MoBueh
"""

#%%
import math
#%%
def RadProc(i_dir, i_dif, sun_h, sun_az, pos_h, pos_az):
    'Computates the resulting radiation in W/m² on the component based on the height of the sun and the positioning of a component \n \n i_dir  -> direct radiation from TRY in W/m² \n i_dif  -> diffuse radiation from TRY in W/m² \n sun_h  -> hight of sun in ° \n sun_az -> sunazimuth in ° \n pos_h  -> alignment of the component in the axis of the sun elevation \n pos_az -> azimuthal orientation of the component \n for azimuth: 0° -> north, 90° -> east, 180° -> south, 270° -> west \n for hight: 0° -> horizontal top, 90°-> vertical, 180° -> horizontal down'
    # Direkte Strahlung
    # Angle of Incidence cos(zeta)
    aoi = math.sin(math.radians(sun_h))*math.cos(math.radians(pos_h)) + math.cos(math.radians(sun_h))*math.sin(math.radians(pos_h)) * math.cos(math.radians(abs(pos_az-sun_az)))
    if aoi < 0:
        aoi = 0
    if sun_h<=5:
        i_f = 0
    else:
        i_f = i_dir/math.sin(math.radians(sun_h))*aoi
    if i_f > 1200:
        i_f = 1200
    # Strahlung durch Bodenreflexion
    i_umg = (i_dir+i_dif)*0.5*0.2*(1-math.cos(math.radians(pos_h)))
    # Diffuse Strahlung
    i_d = i_dif * 0.5 *(1+math.cos(math.radians(pos_h)))
    # Summe
    i_sol = i_f + i_umg + i_d
    return i_sol
    