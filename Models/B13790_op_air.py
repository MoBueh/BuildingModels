# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 15:25:59 2022

@author: 49157
"""
#%% Directory
import os
os.chdir('C:/Users/49157/Desktop/PuBetA/92_Daten')

#%% Imports und Funktionen
import math
import matplotlib.pyplot as plt


#%%
param = {}
param['building'] = {
        'A_f': 20.96,                                                               # m²    conditioned area
        'Coe_Am': 2.5,                                                                # -     coefficient for the determination of mass-related area according to table 12
        'Coe_Cm': 150000,                                                          # J/K*m² coefficiant for the determination of Internal heat storage capacity according to table 12
        'V': 34.56,                                                                 # m²    Volume of building
    'transparent components': {
        'A'          : [0,      0,      0,      0,      0, 0],                     # m²    Area of the component
        'U'          : 0.8,                    # W/m²/K Heat transfer coefficient determined according to ISO 6946
        'g_tot'      :    0.4
        },
    'opaque components': {
        'A'          : [6.4,    15.76,  6.4,    15.76,  14.96,  14.96],         # m²    Area of the component
        'U'          : 0.4,           # W/m²/K Heat transfer coefficient determined according to ISO 6946
        'abs_coe'    : 0.75,
        'R_se'       : 0.19
        }    
    }

#%% Modell
def building_op(phi_int, theta_e, V_pkt, theta_sup, phi_HC, theta_m_prev, i_sol, shading, param):
    'Building Model based in DIN EN ISO 13790. \ncomputes the temperature at the air, surface and mass nodes based on the following inputs: \n \n phi_int: internal heat gain in W \n theta_e: ambient Temperature in °C \n V_pkt: Infiltration volume flow in m³/h \n theta_sup: Supply temperature of ventilation (no ventilation: theta_sup = theta_e) \n phi_HC: Heating and/or Cooling of zone in W \n theta_m_prev: previous temperature in mass node in °C \n i_sol: solar radiation on walls and roof in W/m². Input has to be a list. \n shading: shading of windows. Input has to be a list. The order must correspond to i_sol'
    
# Parameter
    A_m = param['building']['A_f']*param['building']['Coe_Am']
    A_tot = 4.5 * param['building']['A_f']
    c_m = param['building']['Coe_Cm']*param['building']['A_f']
    H_tr_w = sum([param['building']['transparent components']['U']*A for A in param['building']['transparent components']['A']])
    H_ve = 1200*V_pkt/3600         # Infiltration
    # H_ve = (1-param['ventilation']['etha_wrg'])*1200*V_pkt/3600         # Für Wärmerückgewinnung
    if V_pkt == 0:
        H_ve = 0.00000001
    else:
        H_ve = 1*1200*V_pkt/3600
    H_tr_ms = 9.1*A_m
    H_op = sum([param['building']['opaque components']['U']*A for A in param['building']['opaque components']['A']])
    H_tr_em = 1/(1/H_op+1/H_tr_ms)
    H_tr_is = 3.45*A_tot
    
# solar heat gain inside zone due to transparent components
    phi_sol_trans = []
    for shading, A, i in zip(shading, param['building']['transparent components']['A'], i_sol):
        phi_sol_trans.append(shading* param['building']['transparent components']['g_tot']*0.9*A*i)
# solar heat gain inside zone due to opaque components
    phi_sol_op = []
    for A, i in zip(param['building']['opaque components']['A'], i_sol):
        phi_sol_op.append(param['building']['opaque components']['abs_coe']*param['building']['opaque components']['R_se']*param['building']['opaque components']['U']*A*i)   
# total solar heat gain
    phi_sol = sum(phi_sol_trans) + sum(phi_sol_op)    
    print(phi_sol)
# Equation for Model regarding Appendix C
    H_tr_1 = 1 / (1/H_ve + 1/H_tr_is)                                           # C.6
    H_tr_2 = H_tr_1 + H_tr_w                                                    # C.7
    H_tr_3 = 1 / (1/H_tr_2 + 1/ H_tr_ms)                                        # C.8
    phi_ia = 0.5*phi_int                                                        # C.1
    phi_m = A_m/A_tot * (0.5 * phi_int + phi_sol)                               # C.2
    phi_st = (1-A_m/A_tot-H_tr_w/9.1/A_tot)*(0.5 * phi_int + phi_sol)           # C.3
    phi_mtot = phi_m + H_tr_em * theta_e + H_tr_3 *(phi_st +H_tr_w * theta_e + H_tr_1 * ((phi_ia + phi_HC)/H_ve + theta_sup))/ H_tr_2   # C.5
    theta_m_t = (theta_m_prev * (c_m/3600-0.5 * (H_tr_3 +H_tr_em))+phi_mtot)/(c_m/3600+0.5*(H_tr_3+H_tr_em))  # C.4
    # Output
    theta_m = (theta_m_t + theta_m_prev)/2                                      # C.9
    theta_s = (H_tr_ms * theta_m + phi_st + H_tr_w * theta_e + H_tr_1 *(theta_sup + (phi_ia + phi_HC)/H_ve))/(H_tr_ms + H_tr_w + H_tr_1)     # C.10
    theta_air = (H_tr_is * theta_s + H_ve * theta_sup + phi_ia + phi_HC)/(H_tr_is + H_ve)          # C.11
    theta_operativ = 0.3 * theta_air + 0.7 * theta_s                            # C.12
    return theta_m, theta_s, theta_air, theta_operativ

#%% Testing

    GB = building_op(phi_int = 0, 
                     theta_e = 10, 
                     V_pkt = 75, 
                     theta_sup = 20, 
                     phi_HC = 0,
                     theta_m_prev = 10,
                     i_sol = [0]*5,
                     shading = [1]*5,
                     param = param)



