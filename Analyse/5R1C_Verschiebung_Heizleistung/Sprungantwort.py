# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 14:29:42 2023

@author: MoBueh
"""

#%% Directory
import os
os.chdir('C:/Users/49157/Documents/Python/BuildingModels')

#%% Library
import pandas as pd
import matplotlib.pyplot as plt
from Models.B13790_op_air import building_op as building_op_air
from Models.B13790_op_mass import building_op as building_op_mass
from Models.B13790_op_surface import building_op as building_op_surface

plt.rc('axes', axisbelow=True)
plt.rcParams['figure.dpi'] = 800

#%% Parameterset
param = {}
param['building'] = {
        'A_f': 100,                                                               # m²    conditioned area
        'Coe_Am': 2.5,                                                                # -     coefficient for the determination of mass-related area according to table 12
        'Coe_Cm': 165000,                                                          # J/K*m² coefficiant for the determination of Internal heat storage capacity according to table 12
        'V': 750,                                                                 # m²    Volume of building
    'transparent components': {
        'A'          : [12, 12, 12, 12, 0],                     # m²    Area of the component
        'U'          : 0.8,                    # W/m²/K Heat transfer coefficient determined according to ISO 6946
        'g_tot'      : 0.4
        },
    'opaque components': {
        'A'          : [63, 63, 63, 63, 63],         # m²    Area of the component
        'U'          : 0.4,           # W/m²/K Heat transfer coefficient determined according to ISO 6946
        'abs_coe'    : 0.3,
        'R_se'       : 0.19
        }    
    }

#%% Air
GB_theta_m = [10]
GB_theta_s = [10]
GB_theta_air = [10]
GB_theta_operative = [10]

days = 20
# phi_HC = [0]*24+[1000]*(24*days-24)
phi_HC = [0]+[1000]*(24*days-1)

for time in range(0,24*days,1):
# Design
    GB = building_op_air(phi_int = 0,
                         theta_e = 10,
                         V_pkt = 20,
                         theta_sup = 10,
                         phi_HC = phi_HC[time],
                         theta_m_prev = GB_theta_m[-1],
                         i_sol = [0]*5,
                         shading = [1]*5,
                         param = param)
    

    
    
    # GB_phi_HC.append(GB[0])
    GB_theta_m.append(GB[0])
    GB_theta_s.append(GB[1])
    GB_theta_air.append(GB[2])
    GB_theta_operative.append(GB[3])
    
air = pd.DataFrame({'air' : GB_theta_air,'surface' : GB_theta_s, 'mass' : GB_theta_m, 'operative' : GB_theta_operative})
del GB_theta_m, GB_theta_s, GB_theta_air, GB_theta_operative, GB, days, time

#%% Surface
GB_theta_m = [10]
GB_theta_s = [10]
GB_theta_air = [10]
GB_theta_operative = [10]

days = 20
# phi_HC = [0]*24+[1000]*(24*days-24)

for time in range(0,24*days,1):
# Design
    GB = building_op_surface(phi_int = 0,
                             theta_e = 10,
                             V_pkt = 20,
                             theta_sup = 10,
                             phi_HC = phi_HC[time],
                             theta_m_prev = GB_theta_m[-1],
                             i_sol = [0]*5,
                             shading = [1]*5,
                             param = param)
    

    
    
    # GB_phi_HC.append(GB[0])
    GB_theta_m.append(GB[0])
    GB_theta_s.append(GB[1])
    GB_theta_air.append(GB[2])
    GB_theta_operative.append(GB[3])
    
surface = pd.DataFrame({'air' : GB_theta_air,'surface' : GB_theta_s, 'mass' : GB_theta_m, 'operative' : GB_theta_operative})
del GB_theta_m, GB_theta_s, GB_theta_air, GB_theta_operative, GB, days, time

#%% Mass
GB_theta_m = [10]
GB_theta_s = [10]
GB_theta_air = [10]
GB_theta_operative = [10]

days = 20
# phi_HC = [0]*24+[1000]*(24*days-24)

for time in range(0,24*days,1):
# Design
    GB = building_op_mass(phi_int = 0,
                         theta_e = 10,
                         V_pkt = 20,
                         theta_sup = 10,
                         phi_HC = phi_HC[time],
                         theta_m_prev = GB_theta_m[-1],
                         i_sol = [0]*5,
                         shading = [1]*5,
                         param = param)
    

    
    
    # GB_phi_HC.append(GB[0])
    GB_theta_m.append(GB[0])
    GB_theta_s.append(GB[1])
    GB_theta_air.append(GB[2])
    GB_theta_operative.append(GB[3])
    
mass = pd.DataFrame({'air' : GB_theta_air,'surface' : GB_theta_s, 'mass' : GB_theta_m, 'operative' : GB_theta_operative})
del GB_theta_m, GB_theta_s, GB_theta_air, GB_theta_operative, GB, days, time


#%% Plots
#%%% Each
# Masse
linewidth = 1
plt.plot(mass.index, mass['air'], label = 'air', linewidth = linewidth, c = 'tab:green')
plt.plot(mass.index, mass['surface'], label = 'surface', linewidth = linewidth, c = 'tab:orange')
plt.plot(mass.index, mass['mass'], label = 'mass', linewidth = linewidth, c = 'tab:blue')
plt.plot(range(0,len(phi_HC)), [phi/50 for phi in phi_HC], label = 'phi_HC/50', c = 'tab:red')
plt.ylabel('Temperatur in °C')
plt.title('mass')
plt.ylim(-2.5,22.5)
# plt.xlim(0, 200)
plt.legend(loc = 'lower right')
plt.grid()
plt.show()

# # Surface
plt.plot(surface.index, surface['air'], label = 'air', linewidth = linewidth, c = 'tab:green')
plt.plot(surface.index, surface['surface'], label = 'surface', linewidth = linewidth, c = 'tab:orange')
plt.plot(surface.index, surface['mass'], label = 'mass', linewidth = linewidth, c = 'tab:blue')
plt.plot(range(0,len(phi_HC)), [phi/50 for phi in phi_HC], label = 'phi_HC/50', c = 'tab:red')
plt.ylabel('Temperatur in °C')
plt.title('surface')
plt.ylim(-2.5,22.5)
plt.legend(loc = 'lower right')
plt.grid()
plt.show()

# Air
plt.plot(air.index, air['air'], label = 'air', linewidth = linewidth, c = 'tab:green')
plt.plot(air.index, air['surface'], label = 'surface', linewidth = linewidth, c = 'tab:orange')
plt.plot(air.index, air['mass'], label = 'mass', linewidth = linewidth, c = 'tab:blue')
plt.plot(range(0,len(phi_HC)), [phi/50 for phi in phi_HC], label = 'phi_HC/50', c = 'tab:red')
plt.ylabel('Temperatur in °C')
plt.title('air')
plt.ylim(-2.5,22.5)
plt.legend(loc = 'lower right')
plt.grid()
plt.show()

#%%% Combined
plt.scatter(2,2, c = 'white', label ='Phi_HC auf:')
plt.plot(air.index, air['operative'], label = 'Luftknoten', linewidth = linewidth, c = 'tab:green')
plt.plot(air.index, surface['operative'], label = 'Oberflächenknoten', linewidth = linewidth, c = 'tab:orange')
plt.plot(air.index, mass['operative'], label = 'Masseknoten', linewidth = linewidth, c = 'tab:blue')


plt.plot([0, 1, 1] + list(range(2,len(phi_HC))), [phi/50 for phi in [0] + phi_HC], label = 'phi_HC/50', c = 'tab:red')

plt.ylabel('Temperatur des Luftknotens in °C')
plt.title('air')
plt.ylim(-2.5,22.5)
plt.xlim(0,10)
plt.legend(loc = 'lower right')
plt.grid()
plt.show()




