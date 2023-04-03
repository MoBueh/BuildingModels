# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 14:29:42 2023

@author: 49157
"""

#%% Directory
import os
os.chdir('C:/Users/49157/Desktop/Modellvergleich')

#%% Library
import pandas as pd
import matplotlib.pyplot as plt
from M13790.Modell.B13790_op_air import building_op as building_op_air
from M13790.Modell.B13790_op_mass import building_op as building_op_mass
from M13790.Modell.B13790_op_surface import building_op as building_op_surface
from M13790.Modell.B13790_DESIGN import building_de
from Modelle.RadProc import RadProc

plt.rc('axes', axisbelow=True)
plt.rcParams['figure.dpi'] = 800

#%% Input
inp = pd.read_csv('Input.csv', index_col = 0)


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
        'A'          : [63, 63, 63, 63, 100],         # m²    Area of the component
        'U'          : 0.4,           # W/m²/K Heat transfer coefficient determined according to ISO 6946
        'abs_coe'    : 0.3,
        'R_se'       : 0.19
        }    
    }

#%% Air
GB_phi_HC = [0]
GB_theta_m = [20]
GB_theta_s = [20]
GB_theta_air = [20]
GB_theta_operative = [20]

days = 1

for time in range(0,24*days,1):
    i_sol = []
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 0))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 90))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 180))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 270))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 0, pos_az = 270))

# Design
    GB = building_de(theta_set_heat = 20, theta_set_cool = 26, phi_int = 0, 
                            theta_e = inp.loc[time, 't'], 
                            V_pkt = 20, 
                            theta_sup = inp.loc[time, 't'], 
                            theta_m_prev = GB_theta_m[-1], 
                            i_sol = i_sol,
                            shading = [1]*5, 
                            param = param, 
                            operation = building_op_air)
    GB_phi_HC.append(GB[0])
    GB_theta_m.append(GB[1])
    GB_theta_s.append(GB[2])
    GB_theta_air.append(GB[3])
    GB_theta_operative.append(GB[4])
    
air = pd.DataFrame({'air' : GB_theta_air,'surface' : GB_theta_s, 'mass' : GB_theta_m, 'operative' : GB_theta_operative, 'phi_HC' : GB_phi_HC})
del GB_phi_HC, GB_theta_m, GB_theta_s, GB_theta_air, GB_theta_operative, GB, days, i_sol, time
    

#%% Surface
GB_phi_HC = [0]
GB_theta_m = [20]
GB_theta_s = [20]
GB_theta_air = [20]
GB_theta_operative = [20]

days = 365

for time in range(0,24*days,1):
    i_sol = []
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 0))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 90))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 180))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 270))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 0, pos_az = 270))
    
# Design
    GB = building_de(theta_set_heat = 20, theta_set_cool = 26, phi_int = 0, 
                            theta_e = inp.loc[time, 't'], 
                            V_pkt = 20, 
                            theta_sup = inp.loc[time, 't'], 
                            theta_m_prev = GB_theta_m[-1], 
                            i_sol = i_sol,
                            shading = [1]*5, 
                            param = param, 
                            operation = building_op_surface)
    GB_phi_HC.append(GB[0])
    GB_theta_m.append(GB[1])
    GB_theta_s.append(GB[2])
    GB_theta_air.append(GB[3])
    GB_theta_operative.append(GB[4])
    
surface = pd.DataFrame({'air' : GB_theta_air,'surface' : GB_theta_s, 'mass' : GB_theta_m, 'operative' : GB_theta_operative, 'phi_HC' : GB_phi_HC})
del GB_phi_HC, GB_theta_m, GB_theta_s, GB_theta_air, GB_theta_operative, GB, days, i_sol, time

#%% Mass
GB_phi_HC = [0]
GB_theta_m = [20]
GB_theta_s = [20]
GB_theta_air = [20]
GB_theta_operative = [20]

days = 365

for time in range(0,24*days,1):
    i_sol = []
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 0))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 90))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 180))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 90, pos_az = 270))
    i_sol.append(RadProc(i_dir = inp.loc[time, 'B'], i_dif = inp.loc[time, 'D'], 
                         sun_h = inp.loc[time, 'sun_h'], sun_az = inp.loc[time, 'sun_az'], 
                         pos_h = 0, pos_az = 270))
    
# Design
    GB = building_de(theta_set_heat = 20, theta_set_cool = 26, phi_int = 0, 
                            theta_e = inp.loc[time, 't'], 
                            V_pkt = 20, 
                            theta_sup = inp.loc[time, 't'], 
                            theta_m_prev = GB_theta_m[-1], 
                            i_sol = i_sol,
                            shading = [1]*5, 
                            param = param, 
                            operation = building_op_mass)
    GB_phi_HC.append(GB[0])
    GB_theta_m.append(GB[1])
    GB_theta_s.append(GB[2])
    GB_theta_air.append(GB[3])
    GB_theta_operative.append(GB[4])
    
mass = pd.DataFrame({'air' : GB_theta_air,'surface' : GB_theta_s, 'mass' : GB_theta_m, 'operative' : GB_theta_operative, 'phi_HC' : GB_phi_HC})
del GB_phi_HC, GB_theta_m, GB_theta_s, GB_theta_air, GB_theta_operative, GB, days, i_sol, time


#%% Kühlleistung
print('Masse Kühlung in W')
print(round(mass[mass['phi_HC']<0]['phi_HC'].sum(),2))
print('')

print('Surface Kühlung in W')
print(round(surface[surface['phi_HC']<0]['phi_HC'].sum(),2))
print('')

print('Luft Kühlung in W')
print(round(air[air['phi_HC']<0]['phi_HC'].sum(),2))

#%% Heizleistung
print('Masse Heizung in W')
print(round(mass[mass['phi_HC']>0]['phi_HC'].sum(),2))
print('')

print('Surface Heizung in W')
print(round(surface[surface['phi_HC']>0]['phi_HC'].sum(),2))
print('')

print('Luft Heizung in W')
print(round(air[air['phi_HC']>0]['phi_HC'].sum(),2))

#%% Gesamt
# Masse
linewidth = 1
plt.plot(air.index, air['mass'], label = 'mass', linewidth = linewidth, c = 'tab:blue')
plt.plot(air.index, air['surface'], label = 'surface', linewidth = linewidth, c = 'tab:orange')
plt.plot(air.index, air['air'], label = 'air', linewidth = linewidth, c = 'tab:green')
plt.ylabel('Temperatur in °C')
plt.title('air')
plt.ylim(15,30)
plt.legend()
plt.grid()
plt.show()

# Surface
linewidth = 1
plt.plot(air.index, surface['mass'], label = 'mass', linewidth = linewidth, c = 'tab:blue')
plt.plot(air.index, surface['surface'], label = 'surface', linewidth = linewidth, c = 'tab:orange')
plt.plot(air.index, surface['air'], label = 'air', linewidth = linewidth, c = 'tab:green')
plt.ylabel('Temperatur in °C')
plt.title('surface')
plt.ylim(15,30)
plt.legend()
plt.grid()
plt.show()

# Air
plt.plot(mass.index, mass['air'], label = 'air', linewidth = linewidth, alpha = 0.5, c = 'tab:green')
plt.plot(mass.index, mass['surface'], label = 'surface', linewidth = linewidth, c = 'tab:orange')
plt.plot(mass.index, mass['mass'], label = 'mass', linewidth = linewidth, c = 'tab:blue')
plt.ylabel('Temperatur in °C')
plt.title('mass')
plt.ylim(15,30)
plt.legend()
plt.grid()
plt.show()

#%% OneDay
# Mass
linewidth = 2
plt.plot(range(0,48,1), air['mass'][24:72], label = 'mass', linewidth = linewidth, c = 'tab:blue')
plt.plot(range(55,103,1), air['mass'][24+4000:72+4000], linewidth = linewidth, c = 'tab:blue')
plt.plot(range(0,48,1), air['surface'][24:72], label = 'surface', linewidth = linewidth, c = 'tab:orange')
plt.plot(range(55,103,1), air['surface'][24+4000:72+4000], linewidth = linewidth, c = 'tab:orange')
plt.plot(range(0,48,1), air['air'][24:72], label = 'air', linewidth = linewidth, c = 'tab:green')
plt.plot(range(55,103,1), air['air'][24+4000:72+4000], linewidth = linewidth, c = 'tab:green')
plt.ylabel('Temperatur in °C')
plt.title('air')
plt.ylim(15,30)
plt.text(15,22,'Winter')
plt.text(72.5,22,'Sommer')
plt.legend()
plt.grid()
plt.show()

# Surface
linewidth = 2
plt.plot(range(0,48,1), surface['mass'][24:72], label = 'mass', linewidth = linewidth, c = 'tab:blue')
plt.plot(range(55,103,1), surface['mass'][24+4000:72+4000], linewidth = linewidth, c = 'tab:blue')
plt.plot(range(0,48,1), surface['surface'][24:72], label = 'surface', linewidth = linewidth, c = 'tab:orange')
plt.plot(range(55,103,1), surface['surface'][24+4000:72+4000], linewidth = linewidth, c = 'tab:orange')
plt.plot(range(0,48,1), surface['air'][24:72], label = 'air', linewidth = linewidth, c = 'tab:green')
plt.plot(range(55,103,1), surface['air'][24+4000:72+4000], linewidth = linewidth, c = 'tab:green')
plt.ylabel('Temperatur in °C')
plt.title('Surface')
plt.ylim(15,30)
plt.text(15,22,'Winter')
plt.text(72.5,22,'Sommer')
plt.legend()
plt.grid()
plt.show()

# Air
linewidth = 2
plt.plot(range(0,48,1), mass['mass'][24:72], label = 'mass', linewidth = linewidth, c = 'tab:blue')
plt.plot(range(55,103,1), mass['mass'][24+4000:72+4000], linewidth = linewidth, c = 'tab:blue')
plt.plot(range(0,48,1), mass['surface'][24:72], label = 'surface', linewidth = linewidth, c = 'tab:orange')
plt.plot(range(55,103,1), mass['surface'][24+4000:72+4000], linewidth = linewidth, c = 'tab:orange')
plt.plot(range(0,48,1), mass['air'][24:72], label = 'air', linewidth = linewidth, c = 'tab:green')
plt.plot(range(55,103,1), mass['air'][24+4000:72+4000], linewidth = linewidth, c = 'tab:green')
plt.ylabel('Temperatur in °C')
plt.title('mass')
plt.ylim(15,30)
plt.text(15,22,'Winter')
plt.text(72.5,22,'Sommer')
plt.legend(loc = 'upper left')
plt.grid()
plt.show()

#%% Komfort
# Air
plt.scatter(inp['t'], air['operative'][1:], s = 1, c = 'tab:orange')
plt.scatter(inp['t'], air['air'][1:], s = 1, c = 'tab:green', alpha = 0.2)
plt.scatter(200,200, s = 15, c = 'tab:orange', label = 'operative')
plt.scatter(200,200, s = 15, c = 'tab:green', label = 'air')
plt.title('air')
plt.ylim(15,30)
plt.xlim(-10,38)
plt.legend(loc = 'upper left')
plt.grid()
plt.ylabel('Temperatur in °C')
plt.xlabel('Außentemperatur in °C')
plt.show()

# Air
plt.scatter(inp['t'], surface['operative'][1:], s = 1, c = 'tab:orange')
plt.scatter(inp['t'], surface['air'][1:], s = 1, c = 'tab:green', alpha = 0.2)
plt.scatter(200,200, s = 15, c = 'tab:orange', label = 'operative')
plt.scatter(200,200, s = 15, c = 'tab:green', label = 'air')
plt.title('surface')
plt.ylim(15,30)
plt.xlim(-10,38)
plt.legend(loc = 'upper left')
plt.grid()
plt.ylabel('Temperatur in °C')
plt.xlabel('Außentemperatur in °C')
plt.show()

# Mass
plt.scatter(inp['t'], mass['operative'][1:], s = 1, c = 'tab:orange')
plt.scatter(inp['t'], mass['air'][1:], s = 1, c = 'tab:green', alpha = 0.2)
plt.scatter(200,200, s = 15, c = 'tab:orange', label = 'operative')
plt.scatter(200,200, s = 15, c = 'tab:green', label = 'air')
plt.title('mass')
plt.ylim(15,30)
plt.xlim(-10,38)
plt.ylabel('Temperatur in °C')
plt.xlabel('Außentemperatur in °C')
plt.legend(loc = 'upper left')
plt.grid()
plt.show()



