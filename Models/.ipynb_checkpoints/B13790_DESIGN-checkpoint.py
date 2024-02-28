# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 12:17:56 2023

@author: 49157
"""

#%% Directory
import os
os.chdir('C:/Users/49157/Desktop/Modellvergleich')

#%% Library
import pandas as pd
import matplotlib.pyplot as plt
from M13790.Modell.B13790_op_air import building_op
from Modelle.RadProc import RadProc

#%%
def building_de(theta_set_heat, theta_set_cool, phi_int, theta_e, V_pkt, theta_sup, theta_m_prev, i_sol,shading, param, operation):
    # Schritt 1
    mass, surface, air_0, operative = operation(phi_int = phi_int, theta_e = theta_e, theta_sup = theta_sup, phi_HC = 0, theta_m_prev = theta_m_prev ,V_pkt = V_pkt, i_sol = i_sol, shading = shading, param = param)
    if theta_set_heat <= air_0 and theta_set_cool >= air_0:
        phi_HC = 0
        theta_air = air_0
        theta_m = mass
        theta_s = surface
        theta_operativ = operative
    # Schritt 2
    else:
        if air_0 > theta_set_cool:
            air_set = theta_set_cool
        else:
            air_set = theta_set_heat
        mass, surface, air_10, operative = operation(phi_int = phi_int, theta_e = theta_e, theta_sup = theta_sup, phi_HC = 10*param['building']['A_f'], theta_m_prev = theta_m_prev ,V_pkt = V_pkt, shading = shading, i_sol = i_sol, param = param)
        phi_HC = 10*param['building']['A_f']*(air_set - air_0) / (air_10 - air_0)
        mass, surface, air_final, operative = operation(phi_int = phi_int, theta_e = theta_e, theta_sup = theta_sup, phi_HC = phi_HC, theta_m_prev = theta_m_prev ,V_pkt = V_pkt, shading = shading, i_sol = i_sol, param = param)
        theta_air = air_final
        theta_m = mass
        theta_s = surface
        theta_operativ = operative
    return phi_HC, theta_m, theta_s, theta_air, theta_operativ


#%% Testing
# param = {}
# param['building'] = {
#         'A_f': 20.96,                                                               # m²    conditioned area
#         'Coe_Am': 2.5,                                                                # -     coefficient for the determination of mass-related area according to table 12
#         'Coe_Cm': 150000,                                                          # J/K*m² coefficiant for the determination of Internal heat storage capacity according to table 12
#         'V': 34.56,                                                                 # m²    Volume of building
#     'transparent components': {
#         'A'          : [0,      0,      0,      0,      0,      0],                     # m²    Area of the component
#         'U'          : 0.8,                    # W/m²/K Heat transfer coefficient determined according to ISO 6946
#         'g_tot'      :    0.4
#         },
#     'opaque components': {
#         'A'          : [6.4,    15.76,  6.4,    15.76,  14.96,  14.96],         # m²    Area of the component
#         'U'          : 0.4,           # W/m²/K Heat transfer coefficient determined according to ISO 6946
#         'abs_coe'    : 0.75,
#         'R_se'       : 0.19
#         }    
#     }

# test = building_de(theta_set_heat = 20, theta_set_cool = 26, phi_int = 0, 
#                        theta_e = 8, 
#                        V_pkt = 75, 
#                        theta_sup = 8, 
#                        theta_m_prev = 20, 
#                        i_sol = [20]*5,
#                        shading = [1]*5, 
#                        param = param, 
#                        operation = building_op)




