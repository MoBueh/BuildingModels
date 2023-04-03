# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 16:03:45 2023

@author: 49157
"""

#%% Directory
import os
os.chdir('C:/Users/49157/Desktop/Modellvergleich')

#%% Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Vergleich.B1R1C import B1R1C
from Vergleich.B5R1C_simplified import B5R1C
from Modelle.RadProc import RadProc
plt.rc('axes', axisbelow=True)
plt.rcParams['figure.dpi'] = 800

#%% input
inp = pd.read_csv('Input.csv', index_col = 0)

#%% Parameter-Set
param = {}

# 1R1C
param['B1R1C'] = {
    # 'UA' : 179.2,
    'UA' : 134.4,
    # 'C_m' : 11000000.0,
    # 'C_m' : 16500000.0,
    'C_m' : 26000000.0,
    }

# 5R1C
param['B5R1C'] = {
        'A_f': 100,                                                          # m²    conditioned area
        'Coe_Am': 2.5,                                                         # -     coefficient for the determination of mass-related area according to table 12
        # 'Coe_Cm': 110000,                                                      # J/K*m² coefficiant for the determination of Internal heat storage capacity according to table 12
        # 'Coe_Cm': 165000,                                                      # J/K*m² coefficiant for the determination of Internal heat storage capacity according to table 12
        'Coe_Cm': 260000,                                                     # J/K*m² coefficiant for the determination of Internal heat storage capacity according to table 12
        'V': 750,                                                            # m²    Volume of building
    'transparent': {
        'A'          : [12, 12, 12, 12, 0],                 # m²    Area of the component
        # 'U'          : 0.8,                                                    # W/m²/K Heat transfer coefficient determined according to ISO 6946
        'U'          : 0.6,                                                    # W/m²/K Heat transfer coefficient determined according to ISO 6946
        'g_tot'      : 0.4
        },
    'opaque': {
        'A'          : [63, 63, 63, 63, 100],         # m²    Area of the component
        # 'U'          : 0.4,                                                     # W/m²/K Heat transfer coefficient determined according to ISO 6946
        'U'          : 0.3,                                                     # W/m²/K Heat transfer coefficient determined according to ISO 6946
        'R'          : 0.05,                                                    # Wärmeübergang Strahlung (Absorptionskoeffizient *Oberflächenwärme-Durchlasswiderstand)
        }    
    }
#%% Simulation
days = 60

# 1R1C
GB_theta_i = [0]
# 5R1C
GB_theta_m = [0]
GB_theta_s = [0]
GB_theta_air = [0]
GB_theta_operative = [0]

phi_int = [0]*24 + [400]*(24*days-24)

for time in range(0,24*days,1):
    # 1R1C
    GB_1R1C = B1R1C(phi_int = phi_int[time],
                    theta_e = 0,
                    phi_HC = 0,
                    theta_prev = GB_theta_i[-1],
                    i_sol = 0,
                    param = param['B1R1C'])
    GB_theta_i.append(GB_1R1C)
    # 5R1C
    GB_5R1C = B5R1C(phi_int = phi_int[time],
                    theta_e = 0,
                    V_pkt = 1,
                    theta_sup = inp.loc[time,'t'],
                    phi_HC = 0,
                    theta_m_prev = GB_theta_m[-1],
                    i_sol = [0]*5,
                    shading = [1]*5,
                    param = param['B5R1C'])
    GB_theta_m.append(GB_5R1C[0])
    GB_theta_s.append(GB_5R1C[1])
    GB_theta_air.append(GB_5R1C[2])
    GB_theta_operative.append(GB_5R1C[3])

# df_1R1C = pd.DataFrame({'mittel':GB_theta_i[1:]})
# df_1R1C['leicht'] = GB_theta_i[1:]
df_1R1C['schwer'] = GB_theta_i[1:]

# df_5R1C = pd.DataFrame({'mittel' : GB_theta_air[1:],'surface' : GB_theta_s[1:], 'mass' : GB_theta_m[1:], 'operative' : GB_theta_operative[1:]})
# df_5R1C['leicht'] = GB_theta_air[1:]
df_5R1C['schwer'] = GB_theta_air[1:]
del GB_theta_m, GB_theta_s, GB_theta_air, GB_theta_operative, time, GB_1R1C, GB_5R1C, GB_theta_i

#%%
bis = 24*15

plt.plot(df_5R1C.index[0:bis], [i/100 for i in phi_int][0:bis], c = 'grey', label = 'phi_int/100')
plt.plot([0,1],[-20,-30],linestyle = '-.', c = 'white', label = ' ')

plt.plot([0,1],[-20,-30], c = 'black', label = 'leicht')
plt.plot([0,1],[-20,-30],linestyle = '--', c = 'black', label = 'mittel')
plt.plot([0,1],[-20,-30],linestyle = '-.', c = 'black', label = 'schwer')
plt.plot([0,1],[-20,-30],linestyle = '-.', c = 'white', label = ' ')

plt.plot(df_1R1C.index[0:bis], df_1R1C['mittel'][0:bis], label = '1R1C', c = 'tab:green')
plt.plot(df_1R1C.index[0:bis], df_1R1C['leicht'][0:bis],linestyle = '--', c = 'tab:green', alpha = 0.3)
plt.plot(df_1R1C.index[0:bis], df_1R1C['schwer'][0:bis],linestyle = '-.', c = 'tab:green', alpha = 0.7)


plt.plot(df_5R1C.index[0:bis], df_5R1C['mittel'][0:bis], label = '5R1C', c = 'tab:blue')
plt.plot(df_5R1C.index[0:bis], df_5R1C['leicht'][0:bis],linestyle = '--', c = 'tab:blue', alpha = 0.3)
plt.plot(df_5R1C.index[0:bis], df_5R1C['schwer'][0:bis],linestyle = '-.', c = 'tab:blue', alpha = 0.7)

plt.ylim(-1,4.5)
plt.grid()
plt.xlabel('Zeitschritt in h')
plt.ylabel('Temperatur in °C')
plt.legend(fontsize = 8)






















