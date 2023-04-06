# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 16:03:45 2023

@author: MoBueh
"""

#%% Directory
import os
os.chdir('C:/Users/49157/Documents/Python/BuildingModels')

#%% Libraries
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt

from Models.min_B5R1C_simplified import B5R1C
from Models.min_B1R1C import B1R1C
# from Vergleich.min_B1R1C_JP import B1R1C_JP

# from Vergleich.B1R1C import B1R1C
# from Vergleich.B5R1C_simplified import B5R1C

plt.rc('axes', axisbelow=True)
plt.rcParams['figure.dpi'] = 800

#%% Parameter-Set
param = {}

# 1R1C
param['B1R1C'] = {
    'UA' : 179.2,
    # 'C_m' : 11000000.0,
    'C_m' : 16500000.0,
    # 'C_m' : 26000000.0,
    }

# 5R1C
param['B5R1C'] = {
        'A_f': 100,                                                          # m²    conditioned area
        'Coe_Am': 2.5,                                                         # -     coefficient for the determination of mass-related area according to table 12
        # 'Coe_Cm': 110000,                                                      # J/K*m² coefficiant for the determination of Internal heat storage capacity according to table 12
        'Coe_Cm': 165000,                                                      # J/K*m² coefficiant for the determination of Internal heat storage capacity according to table 12
        # 'Coe_Cm': 260000,                                                      # J/K*m² coefficiant for the determination of Internal heat storage capacity according to table 12
        'V': 750,                                                            # m²    Volume of building
    'transparent': {
        'A'          : [12, 12, 12, 12, 0],                 # m²    Area of the component
        'U'          : 0.8,                                                    # W/m²/K Heat transfer coefficient determined according to ISO 6946
        'g_tot'      : 0.4
        },
    'opaque': {
        'A'          : [63, 63, 63, 63, 100],         # m²    Area of the component
        'U'          : 0.4,                                                     # W/m²/K Heat transfer coefficient determined according to ISO 6946
        'R'          : 0.05,                                                    # Wärmeübergang Strahlung (Absorptionskoeffizient *Oberflächenwärme-Durchlasswiderstand)
        }    
    }

#%% Test
T = 60*24    #min

theta_e = [0]
for i in range(1,(T*10)+1):
    theta_e.append(100*np.sin(i*2*np.pi/T))
plt.plot(range(0, len(theta_e)), theta_e)

# 1R1C
GB_theta_i = [0]
GB_theta_i_JP = [0]

GB_theta_m = [0]
GB_theta_s = [0]
GB_theta_air = [0]
GB_theta_operative = [0]

for time in range(0,T*10,1):
    # 1R1C
    GB_1R1C = B1R1C(phi_int = 0,
                    # theta_e = inp.loc[time,'t'],
                    theta_e = theta_e[time],
                    phi_HC = 0,
                    theta_prev = GB_theta_i[-1],
                    i_sol = 0,
                    param = param['B1R1C'])
    GB_theta_i.append(GB_1R1C)
    # GB_1R1C_JP = B1R1C_JP(phi_int = 0,
    #                    # theta_e = inp.loc[time,'t'],
    #                    theta_e = theta_e[time],
    #                    phi_HC = 0,
    #                    theta_prev = GB_theta_i[-1],
    #                    i_sol = 0,
    #                    param = param['B1R1C'])
    # GB_theta_i_JP.append(GB_1R1C)
    GB_5R1C = B5R1C(phi_int = 0,
                    # theta_e = inp.loc[time,'t'],
                    theta_e = theta_e[time],
                    V_pkt = 1,
                    theta_sup = theta_e[time],
                    phi_HC = 0,
                    theta_m_prev = GB_theta_m[-1],
                    i_sol = [0]*5,
                    shading = [1]*5,
                    param = param['B5R1C'])
    GB_theta_m.append(GB_5R1C[0])
    GB_theta_s.append(GB_5R1C[1])
    GB_theta_air.append(GB_5R1C[2])
    GB_theta_operative.append(GB_5R1C[3])
    
# print(abs(min(GB_theta_air))/100)
# damp_5R1C.append(abs(min(GB_theta_air)))

plt.plot(range(0, len(GB_theta_i)), theta_e)
plt.plot(range(0, len(GB_theta_i)), GB_theta_i)
# plt.plot(range(0, len(GB_theta_i)), GB_theta_i_JP)
plt.plot(range(0, len(GB_theta_i)), GB_theta_air)
plt.xlabel('t in h')
plt.grid()
# print(abs(min(GB_theta_i)))

# print(sum([abs(a-b) for a,b in zip(GB_theta_i, GB_theta_i_JP)]))


#%% Phaseshift
steps = []
a = 23
while a<=24*3:
    a = a+1
    steps.append(a*60)
while a<=24*6:
    a = a+10
    steps.append(a*60)
while a<=24*10:
    a = a+100
    steps.append(a*60)
while a<=8560:
    a = a+200
    steps.append(a*60)
del a
steps.append(8760*60)

phaseshift = pd.DataFrame()
phaseshift['T / h'] = [i/60 for i in steps]

for cm, kat in zip([110000, 165000, 260000], ['leicht', 'mittel', 'schwer']):
    param['B5R1C']['Coe_Cm'] = cm
    param['B1R1C']['C_m'] = cm*param['B5R1C']['A_f']
    ps_1R1C = []
    ps_5R1C = []

    for h in steps:
        print(kat + ': ' + str(h/60))
        T = h
        theta_e = [0]
        for i in range(1,(T*5)+1):
            theta_e.append(100*np.sin(i*2*np.pi/T))
        # plt.plot(range(0, len(theta_e)), theta_e)

        # 1R1C
        GB_theta_i = [0]

        for time in range(0,T*5,1):
            # 1R1C
            GB_1R1C = B1R1C(phi_int = 0,
                            # theta_e = inp.loc[time,'t'],
                            theta_e = theta_e[time],
                            phi_HC = 0,
                            theta_prev = GB_theta_i[-1],
                            i_sol = 0,
                            param = param['B1R1C'])
            GB_theta_i.append(GB_1R1C)
        df_1R1C = pd.DataFrame({'input':theta_e, 'output':GB_theta_i})
        max_in = df_1R1C[(0+4)*T: (1+4)*T]['input'].idxmax()
        max_out = df_1R1C[(0+4)*T: (1+4)*T]['output'].idxmax()
        # print(abs(360/T*max_in-360/T*max_out))
        ps_1R1C.append(abs(360/T*max_in-360/T*max_out))

            
        
        # 5R1C
        GB_theta_m = [0]
        GB_theta_s = [0]
        GB_theta_air = [0]
        GB_theta_operative = [0]

        for time in range(0,T*5,1):
            # 5R1C
            GB_5R1C = B5R1C(phi_int = 0,
                            # theta_e = inp.loc[time,'t'],
                            theta_e = theta_e[time],
                            V_pkt = 1,
                            theta_sup = theta_e[time],
                            phi_HC = 0,
                            theta_m_prev = GB_theta_m[-1],
                            i_sol = [0]*5,
                            shading = [1]*5,
                            param = param['B5R1C'])
            GB_theta_m.append(GB_5R1C[0])
            GB_theta_s.append(GB_5R1C[1])
            GB_theta_air.append(GB_5R1C[2])
            GB_theta_operative.append(GB_5R1C[3])
        df_5R1C = pd.DataFrame({'input':theta_e, 'output':GB_theta_air})
        max_in = df_5R1C[(0+4)*T: (1+4)*T]['input'].idxmax()
        max_out = df_5R1C[(0+4)*T: (1+4)*T]['output'].idxmax()
        # print(abs(360/T*max_in-360/T*max_out))
        ps_5R1C.append(abs(360/T*max_in-360/T*max_out))
    
    phaseshift['1R1C_' + kat] = ps_1R1C
    phaseshift['5R1C_' + kat] = ps_5R1C
    
del df_1R1C, df_5R1C, GB_1R1C, GB_5R1C, GB_theta_air, GB_theta_i, GB_theta_m, GB_theta_operative, GB_theta_s, i, h, kat, max_in, max_out, ps_1R1C, ps_5R1C, steps, T, theta_e, time, cm

#%%% Storing
phaseshift.to_csv('Analyse/Vergleich_5R1C_1R1C/phaseshift.csv')  
    
#%%% Plotting
phaseshift = pd.read_csv('Analyse/Vergleich_5R1C_1R1C/phaseshift.csv', index_col = 0)

plt.plot([0,1],[-20,-30], c = 'black', label = 'leicht')
plt.plot([0,1],[-20,-30],linestyle = '--', c = 'black', label = 'mittel')
plt.plot([0,1],[-20,-30],linestyle = '-.', c = 'black', label = 'schwer')
plt.plot([0,1],[-20,-30],linestyle = '-.', c = 'white', label = ' ')

plt.plot(phaseshift['T / h'], phaseshift['1R1C_mittel'], label = '1R1C', c = 'tab:green')
plt.plot(phaseshift['T / h'], phaseshift['1R1C_leicht'],linestyle = '--', c = 'tab:green', alpha = 0.3)
plt.plot(phaseshift['T / h'], phaseshift['1R1C_schwer'],linestyle = '-.', c = 'tab:green', alpha = 0.3)

plt.plot(phaseshift['T / h'], phaseshift['5R1C_mittel'], label = '5R1C', c = 'tab:blue')
plt.plot(phaseshift['T / h'], phaseshift['5R1C_leicht'], linestyle = '--', c = 'tab:blue', alpha = 0.3)
plt.plot(phaseshift['T / h'], phaseshift['5R1C_schwer'], linestyle = '-.', c = 'tab:blue', alpha = 0.3)
plt.xlim(1,111180/60)
# plt.xlim(0,100)
plt.ylabel('Phasenverschiebung in °')
plt.xlabel('Periodendauer T in h')
plt.ylim(-2,102)
plt.grid()
plt.legend()
plt.show()

#%% Damping
steps = []
a = 23
while a<=24*3:
    a = a+1
    steps.append(a*60)
while a<=24*6:
    a = a+10
    steps.append(a*60)
while a<=24*10:
    a = a+100
    steps.append(a*60)
while a<=8560:
    a = a+200
    steps.append(a*60)
del a
steps.append(8760*60)

damping = pd.DataFrame()
damping['T / h'] = [i/60 for i in steps]


for cm, kat in zip([110000, 165000, 260000], ['leicht', 'mittel', 'schwer']):
    param['B5R1C']['Coe_Cm'] = cm
    param['B1R1C']['C_m'] = cm*param['B5R1C']['A_f']
    
    damp_1R1C = []
    damp_5R1C = []
    for h in steps:
        print(kat + ': ' + str(h/60))
        T = h
        theta_e = [0]
        for i in range(1,(T*20)+1):
            theta_e.append(100*np.sin(i*2*np.pi/T))
            
        # 1R1C
        GB_theta_i = [0]
        for time in range(0,T*20,1):
            # 1R1C
            GB_1R1C = B1R1C(phi_int = 0,
                            # theta_e = inp.loc[time,'t'],
                            theta_e = theta_e[time],
                            phi_HC = 0,
                            theta_prev = GB_theta_i[-1],
                            i_sol = 0,
                            param = param['B1R1C'])
            GB_theta_i.append(GB_1R1C)
        damp_1R1C.append(abs(min(GB_theta_i))/100)
        
        # 5R1C
        GB_theta_m = [0]
        GB_theta_s = [0]
        GB_theta_air = [0]
        GB_theta_operative = [0]
        for time in range(0,T*20,1):
            # 5R1C
            GB_5R1C = B5R1C(phi_int = 0,
                            # theta_e = inp.loc[time,'t'],
                            theta_e = theta_e[time],
                            V_pkt = 1,
                            theta_sup = theta_e[time],
                            phi_HC = 0,
                            theta_m_prev = GB_theta_m[-1],
                            i_sol = [0]*5,
                            shading = [1]*5,
                            param = param['B5R1C'])
            GB_theta_m.append(GB_5R1C[0])
            GB_theta_s.append(GB_5R1C[1])
            GB_theta_air.append(GB_5R1C[2])
            GB_theta_operative.append(GB_5R1C[3])
            
            
            
        damp_5R1C.append(abs(min(GB_theta_air))/100)
        
    damping['1R1C_' + kat + '_G'] = damp_1R1C
    damping['5R1C_' + kat + '_G'] = damp_5R1C

del damp_1R1C, damp_5R1C, GB_1R1C, GB_5R1C, GB_theta_air, GB_theta_i, GB_theta_m, GB_theta_operative, GB_theta_s, i, h, kat, steps, T, theta_e, time, cm

#%%% Storing
damping.to_csv('Analyse/Vergleich_5R1C_1R1C/damping.csv')   
    
#%%% Plotting
damping = pd.read_csv('Analyse/Vergleich_5R1C_1R1C/damping.csv', index_col = 0)

plt.plot([0,1],[-20,-30], c = 'black', label = 'leicht')
plt.plot([0,1],[-20,-30],linestyle = '--', c = 'black', label = 'mittel')
plt.plot([0,1],[-20,-30],linestyle = '-.', c = 'black', label = 'schwer')
plt.plot([0,1],[-20,-30],linestyle = '-.', c = 'white', label = ' ')

plt.plot(damping['T / h'], damping['1R1C_mittel_G'], label = '1R1C', c = 'tab:green')
plt.plot(damping['T / h'], damping['1R1C_leicht_G'],linestyle = '--', c = 'tab:green', alpha = 0.3)
plt.plot(damping['T / h'], damping['1R1C_schwer_G'],linestyle = '-.', c = 'tab:green', alpha = 0.3)

plt.plot(damping['T / h'], damping['5R1C_mittel_G'], label = '5R1C', c = 'tab:blue')
plt.plot(damping['T / h'], damping['5R1C_leicht_G'], linestyle = '--', c = 'tab:blue', alpha = 0.3)
plt.plot(damping['T / h'], damping['5R1C_schwer_G'], linestyle = '-.', c = 'tab:blue', alpha = 0.3)
plt.xlim(0,8760)
plt.ylabel('Dämpfung G')
plt.xlabel('Periodendauer T in h')
plt.ylim(-0.02,1.02)
plt.grid()
plt.legend()
plt.show()

#%% Bode-Diagramm
phaseshift = pd.read_csv('Analyse/Vergleich_5R1C_1R1C/phaseshift.csv', index_col = 0)
damping = pd.read_csv('Analyse/Vergleich_5R1C_1R1C/damping.csv', index_col = 0)

for i in phaseshift.index:
   phaseshift.loc[i, 'f'] = 2*math.pi/(phaseshift.loc[i, 'T / h'])
   
for i in damping.index:
   damping['f'] = 2*math.pi/phaseshift['T / h'] 
   
for typ in damping.columns[1:-1]:
    for i in damping.index:
        damping.loc[i, typ[:-2]+'_dB'] = 20*math.log10(damping.loc[i, typ])
del i, typ

#%%% PLot Damping

plt.plot([0,1],[20,30],linestyle = '--', c = 'black', label = 'leicht')
plt.plot([0,1],[20,30], c = 'black', label = 'mittel')
plt.plot([0,1],[20,30],linestyle = '-.', c = 'black', label = 'schwer')
plt.plot([0,1],[20,30],linestyle = '-.', c = 'white', label = ' ')

plt.plot(damping['f'], damping['1R1C_mittel_dB'], label = '1R1C', c = 'tab:green')
plt.plot(damping['f'], damping['1R1C_leicht_dB'],linestyle = '--', c = 'tab:green', alpha = 0.3)
plt.plot(phaseshift['f'], damping['1R1C_schwer_dB'],linestyle = '-.', c = 'tab:green', alpha = 0.3)
# plt.plot([10**-3, 10**-1], [-3, -3])

plt.plot(damping['f'], damping['5R1C_mittel_dB'], label = '5R1C', c = 'tab:blue')
plt.plot(damping['f'], damping['5R1C_leicht_dB'], linestyle = '--', c = 'tab:blue', alpha = 0.3)
plt.plot(damping['f'], damping['5R1C_schwer_dB'], linestyle = '-.', c = 'tab:blue', alpha = 0.3)
plt.text(7*10**-2, -3.8, 'Schnittpunkt bei \n Zeitkonstanten', fontsize = 6)


for cm in [110000, 165000, 260000]:
    T = 1/(cm*param['B5R1C']['A_f']/param['B1R1C']['UA'])*3600
    plt.scatter(T, -3, c = 'tab:green', s = 5)    
del T


plt.xscale('log')
plt.ylabel('Dämpfung in dB')
plt.xlabel('Frequenz in rad/h')
plt.xticks([0.0008, 0.0009]+[i*0.001 for i in range(1,10)]+[i*0.01 for i in range(1,10)]+[i*0.1 for i in range(1,10)]+[1])
plt.yticks([i/2 for i in list(range(-60,1,5))], [-30, '', -25, '', -20, '', -15, '', -10, '', -5, '', -0])
plt.xlim(math.pi*2/8760, math.pi*2/24)
plt.ylim(-32,2)
plt.grid(alpha = 0.3)
# plt.legend(loc = 'lower left', fontsize = 8)
plt.show()

#%%% Plot Shift
plt.plot([0,1],[20,30],linestyle = '--', c = 'black', label = 'leicht')
plt.plot([0,1],[20,30], c = 'black', label = 'mittel')
plt.plot([0,1],[20,30],linestyle = '-.', c = 'black', label = 'schwer')
plt.plot([0,1],[20,30],linestyle = '-.', c = 'white', label = ' ')

plt.plot(phaseshift['f'], [i*(-1) for i in phaseshift['1R1C_mittel']], label = '1R1C', c = 'tab:green')
plt.plot(phaseshift['f'], [i*(-1) for i in phaseshift['1R1C_leicht']],linestyle = '--', c = 'tab:green', alpha = 0.3)
plt.plot(phaseshift['f'], [i*(-1) for i in phaseshift['1R1C_schwer']],linestyle = '-.', c = 'tab:green', alpha = 0.3)

plt.plot(phaseshift['f'], [i*(-1) for i in phaseshift['5R1C_mittel']], label = '5R1C', c = 'tab:blue')
plt.plot(phaseshift['f'], [i*(-1) for i in phaseshift['5R1C_leicht']], linestyle = '--', c = 'tab:blue', alpha = 0.3)
plt.plot(phaseshift['f'], [i*(-1) for i in phaseshift['5R1C_schwer']], linestyle = '-.', c = 'tab:blue', alpha = 0.3)
plt.xscale('log')
plt.ylabel('Phasenverschiebung in °')

plt.xticks([0.0008, 0.0009]+[i*0.001 for i in range(1,10)]+[i*0.01 for i in range(1,10)]+[i*0.1 for i in range(1,10)]+[1])
plt.xlim(math.pi*2/8760, math.pi*2/24)
plt.xlabel('Frequenz in rad/h')


plt.ylim(-92,2)
plt.yticks([0, -15, -30, -45, -60, -75, -90])
plt.grid(alpha = 0.3)
plt.legend(loc = 'lower left', fontsize = 8)
plt.show()



















