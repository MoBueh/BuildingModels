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

#%% Parameter
prop = {}
prop['wasser'] = {
    'rho'       : 980,                                                          # Density in kg/m³
    'cp'        : 4180                                                          # Specific heat capacity in J/kg/K
    }
prop['Antifrogen'] = {
    'rho'       : 1026,                                                         # Density in kg/m³
    'cp'        : 3500                                                          # Specific heat capacity in J/kg/K     
    }

param = {}
param['TABS'] = {
    'A'         : 15,                                                           # Area of floor
    'da'        : 0.03,                                                        # outside diameter of pipe
    'dp'        : 0.001,                                                        # thickness of pipe
    'dx'        : 0.125,                                                        # Laying distance of pipes
    # 'dx'        : 0.12,                                                        # Laying distance of pipes
    'lambda_p'  : 340,                                                          # Heat transfer pipe (copper)
    'lambda_s'  : 2,                                                          # Heat transfer floor
    'd'         : 0.01,                                                         # Laying depth
    }

#%% Model
def TABS(TABS_theta_in, TABS_Vpkt, B_theta_air):
    #                        _______________
    # theta_in      [°C]    |               | theta_out [°C]
    # --------------------->|               |--------------------->
    # Vpkt          [m³/h]  |               | 
    # --------------------->|     TABS      |
    # B_theta_air   [°C]    |               | Qpkt      [W]
    # --------------------->|               |--------------------->
    #                       |_______________|
    if TABS_Vpkt > 0:
        R_w = (param['TABS']['dx']**0.13*(param['TABS']['da']-2*param['TABS']['dp'])/(TABS_Vpkt*prop['wasser']['rho']/3600*param['TABS']['A']/param['TABS']['dx'])**0.87)/(8*math.pi)
        R_p = (param['TABS']['dx']*math.log(param['TABS']['da']/(param['TABS']['da']-2*param['TABS']['dp'])))/(2*math.pi*param['TABS']['lambda_p'])
        R_x = 1**param['TABS']['dx'] * math.log(param['TABS']['dx']/(math.pi*param['TABS']['da']))/(2*math.pi*param['TABS']['lambda_s'])
        R_z = 1/(2*(TABS_Vpkt*prop['wasser']['rho']/3600)*prop['wasser']['cp'])
        R_d = 1/(1/(1/8 + param['TABS']['d']/param['TABS']['lambda_s']))
        k = 1/(R_w + R_p + R_x + R_z + R_d)
        TABS_theta_out= ((TABS_theta_in - B_theta_air)/math.exp(k*param['TABS']['A']/prop['wasser']['rho']/(TABS_Vpkt/3600)/prop['wasser']['cp']))+B_theta_air
        TABS_Qpkt=prop['wasser']['rho']*TABS_Vpkt/3600*prop['wasser']['cp']* (TABS_theta_in-TABS_theta_out)
    else:
        TABS_Qpkt=np.nan
        TABS_theta_out=np.nan
    TABS_theta_in = TABS_theta_in
    TABS_Vpkt = TABS_Vpkt
    return TABS_theta_in, TABS_theta_out, TABS_Vpkt, TABS_Qpkt

#%% Single Step
# tabs = TABS(TABS_theta_in= 10, TABS_Vpkt = 0.3, B_theta_air = 20)
# print(tabs)
# del tabs

#%% Winter
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_winter.csv', index_col = 0)
df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')
df['HM4_Qpkt'] = [vpkt/1000/3600 * prop['wasser']['rho'] *prop['wasser']['cp'] * (SUP-RET) for vpkt, SUP, RET in zip(df['HM4_Vpkt'], df['HM4_T_SUP'],df['HM4_T_RET'])]
df.loc[df['HM4_Vpkt'] < 40, 'HM4_Vpkt'] = 0
df.loc[df['HM4_Vpkt'] == 0, ['HM4_T_RET', 'HM4_T_SUP']] = np.nan




TABS_theta_in = []
TABS_theta_out = []
TABS_Vpkt = []
TABS_Qpkt = []

for i in df.index:
    tabs = TABS(TABS_theta_in = df.loc[i, 'HM4_T_SUP'], 
            TABS_Vpkt = df.loc[i, 'HM4_Vpkt']/1000,
            B_theta_air = df.loc[i, 'Zone_T'])
    TABS_theta_in.append(tabs[0])
    TABS_theta_out.append(tabs[1])
    TABS_Vpkt.append(tabs[2])
    TABS_Qpkt.append(tabs[3])

df['TABS_theta_in'] = TABS_theta_in
df['TABS_theta_out'] = TABS_theta_out
df['TABS_Vpkt'] = TABS_Vpkt
df['TABS_Qpkt'] = TABS_Qpkt
del tabs, TABS_theta_in, TABS_theta_out, TABS_Vpkt, TABS_Qpkt 

#%%% PLotting Timeseries
for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:], [df.index[0] + timedelta(days = i+1) for i in range(7)][:]):
    # ----------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(6.4, 4.8/1.3))
    # fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]})
    ax2 = plt.twinx()
    # ------------------------------------------------Y1------------------------------------------------
    # Data
    ax1.scatter(df.index, df['HM4_T_SUP'], label = 'Messwert: Vorlauftemperatur', c = 'tab:red', s = 12, alpha = 0.5)
    ax1.plot(df.index, df['HM4_T_SUP'], c = 'tab:red', linestyle = '--', alpha = 0.5)
    ax1.scatter(df.index, df['HM4_T_RET'], label = 'Messwert: Rücklauf', c = 'tab:blue', s = 12, alpha = 0.5)
    ax1.plot(df.index, df['HM4_T_RET'], c = 'tab:blue', linestyle = '--', alpha = 0.5)
    
    # Simulation
    # ax1.scatter(df.index, df['TABS_theta_in'], label = 'Vorlauftemperatur Sim', c = 'tab:red', s = 12, marker = 'x')
    # ax1.plot(df.index, df['TABS_theta_in'], c = 'tab:red', linestyle = '--')
    ax1.scatter(df.index, df['TABS_theta_out'], label = 'Simulation: Rücklauftemperatur ', c = 'tab:blue', s = 12, marker = 'x')
    ax1.plot(df.index, df['TABS_theta_out'], c = 'tab:blue', linestyle = '--')

    # labels
    ax1.set_ylabel('Temperatur in °C')
    ax1.yaxis.set_major_formatter('{:.0f}'.format)
    ax1.set_ylim(15,35)
    ax1.set_yticks(range(15,36,2), minor = True)
    ax1.set_yticks(range(15,36,4))
    # ------------------------------------------------Y2------------------------------------------------
    # Data
    ax2.scatter(df.index, df['TABS_Qpkt'], label = 'Simulation: Wärmeleistung', c = 'tab:green', s = 12, marker = 'x')
    ax2.plot(df.index, df['TABS_Qpkt'], c = 'tab:green', linestyle = '--')
    ax2.scatter(df.index, df['HM4_Qpkt'], label = 'Messwert: Wärmeleistung', c = 'tab:green', s = 12, alpha = 0.5)
    ax2.plot(df.index, df['HM4_Qpkt'], c = 'tab:green', linestyle = '--', alpha = 0.5)
    # labels
    ax2.set_ylabel('Wärmeleistung in W')
    ax2.yaxis.set_major_formatter('{:.0f}'.format)
    ax2.set_ylim(0,7000)
    ax2.set_yticks(range(0,7000,800))
    # ------------------------------------------------X------------------------------------------------
    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))],
                    rotation = 90)
    # ax1.set_xticks([i / 2 for i in range(0,7*2,1)], minor = True)
    ax1.set_xlim(von, bis)
    # ------------------------------------------------Global------------------------------------------------
    ax1.set_title(r'$\bf{Messwerte\ und\ Simulation\ des\ Wärmemengenzählers\ 4}$' + '\nZeitraum: ' + str(von.strftime('%Y-%m-%d')) + ' bis ' + str(bis.strftime('%Y-%m-%d')))
    ax1.grid(True)
    ax1.grid(which='minor', alpha = 0.3)
    ax1.legend(ax1.get_legend_handles_labels()[0]+ax2.get_legend_handles_labels()[0], 
                ax1.get_legend_handles_labels()[1]+ax2.get_legend_handles_labels()[1],
                loc = 'upper left')
    ax1.legend(fontsize = 8, loc = 'upper left')
    ax2.legend(fontsize = 8, loc = 'upper right') 
    
    # text = '\n'.join(f'{key}: {value}' for key, value in param['TABS'].items())
    # ax1.text(von+timedelta(hours = 2), 22, text, ha='center', va='center', fontsize=7)
    plt.show()

#%%% Plotting Energybalance
Q_mess = []
Q_sim = []

for i in df.index.day.unique():
    a = df[df.index.day == i]
    Q_mess.append(a['HM4_Qpkt'].sum()/3600)
    Q_sim.append(a['TABS_Qpkt'].sum()/3600)

# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------SINGLE-PLOT-------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(6.4, 4.8/2))

# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.bar([i-0.2 for i in range(1,8,1)], Q_mess, width = 0.4, label = 'Messwert')
ax1.bar([i+0.2 for i in range(1,8,1)], Q_sim, width = 0.4, label = 'Simulationswert')
# labels
ax1.set_ylabel('Wärmemenge in Wh', fontsize = 8)
# ax1.yaxis.set_major_formatter('{:.1f}'.format)
ax1.set_yticks(range(0,25,2), minor = True)
ax1.set_yticks(range(0,25,4))
ax1.set_ylim(0,25)

# ------------------------------------------------X------------------------------------------------
ax1.set_xticks(range(1,8,1), ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], rotation = 60,fontsize = 8)

# ax1.set_xticks([i / 2 for i in range(0,7*2,1)], minor = True)

# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Tagesbilanzen\ der\ abgegebenen\ Wärmemenge}$' + '\nZeitraum: 27.12.2021 - 02.01.2022')

ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
ax1.legend(fontsize = 8, loc = 'upper left')

ax1.text(6.1, 18.2,  r'$\bf{Gesamte\ Woche:}$' + '\n' + str(round(df['HM4_Qpkt'].sum()/3600,2)) + ' Wh gemessen'
          + '\n' + str(round(df['TABS_Qpkt'].sum()/3600,2)) + ' Wh simuliert' + '\n' +
          str(round((df['TABS_Qpkt'].sum()/3600-df['HM4_Qpkt'].sum()/3600)/(df['HM4_Qpkt'].sum()/3600)*100,2)) + ' % Abweichung',
          fontsize = 7, bbox=dict(facecolor='white', edgecolor='black'))
plt.show()






#%% Sommer
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_summer.csv', index_col = 0)
df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')
df['HM4_Qpkt'] = [vpkt/1000/3600 * prop['wasser']['rho'] *prop['wasser']['cp'] * (SUP-RET) for vpkt, SUP, RET in zip(df['HM4_Vpkt'], df['HM4_T_SUP'],df['HM4_T_RET'])]
df.loc[(df['HM4_Vpkt'] < 40) | (df['HM4_T_RET'] > 20), 'HM4_Vpkt'] = 0
df.loc[df['HM4_Vpkt'] == 0, ['HM4_T_RET', 'HM4_T_SUP']] = np.nan


TABS_theta_in = []
TABS_theta_out = []
TABS_Vpkt = []
TABS_Qpkt = []

for i in df.index:
    tabs = TABS(TABS_theta_in = df.loc[i, 'HM4_T_SUP'], 
            TABS_Vpkt = df.loc[i, 'HM4_Vpkt']/1000,
            B_theta_air = df.loc[i, 'Zone_T'])
    TABS_theta_in.append(tabs[0])
    TABS_theta_out.append(tabs[1])
    TABS_Vpkt.append(tabs[2])
    TABS_Qpkt.append(tabs[3])

df['TABS_theta_in'] = TABS_theta_in
df['TABS_theta_out'] = TABS_theta_out
df['TABS_Vpkt'] = TABS_Vpkt
df['TABS_Qpkt'] = TABS_Qpkt

#%%% PLotting Timeseries
for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:], [df.index[0] + timedelta(days = i+1) for i in range(7)][:]):
    # ----------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(6.4, 4.8/1.2))
    ax2 = plt.twinx()
    # ------------------------------------------------Y1------------------------------------------------
    # Data
    ax1.scatter(df.index, df['HM4_T_SUP'], label = 'Messwert: Vorlauftemperatur', c = 'tab:red', s = 12, alpha = 0.5)
    ax1.plot(df.index, df['HM4_T_SUP'], c = 'tab:red', linestyle = '--', alpha = 0.5)
    ax1.scatter(df.index, df['HM4_T_RET'], label = 'Messwert: Rücklauftemperatur', c = 'tab:blue', s = 12, alpha = 0.5)
    ax1.plot(df.index, df['HM4_T_RET'], c = 'tab:blue', linestyle = '--', alpha = 0.5)
    # Simulation
    # ax1.scatter(df.index, df['TABS_theta_in'], label = 'Simulation: Vorlauftemperatur', c = 'tab:red', s = 12, marker = 'x')
    # ax1.plot(df.index, df['TABS_theta_in'], c = 'tab:red', linestyle = '--')
    ax1.scatter(df.index, df['TABS_theta_out'], label = 'Simulation: Rücklauftemperatur', c = 'tab:blue', s = 12, marker = 'x')
    ax1.plot(df.index, df['TABS_theta_out'], c = 'tab:blue', linestyle = '--')

    # labels
    ax1.set_ylabel('Temperatur in °C')
    ax1.yaxis.set_major_formatter('{:.0f}'.format)
    ax1.set_ylim(0,30)
    ax1.set_yticks(range(0,30,2), minor = True)
    ax1.set_yticks(range(0,30,4))
    # ------------------------------------------------Y2------------------------------------------------
    # Data
    ax2.scatter(df.index, abs(df['TABS_Qpkt']), label = 'Simulation: Kälteleistung', c = 'tab:green', s = 12, marker = 'x')
    ax2.plot(df.index, abs(df['TABS_Qpkt']), c = 'tab:green', linestyle = '--')
    ax2.scatter(df.index, abs(df['HM4_Qpkt']), label = 'Messwert: Kälteleistung', c = 'tab:green', s = 12, alpha = 0.5)
    ax2.plot(df.index, abs(df['HM4_Qpkt']), c = 'tab:green', linestyle = '--', alpha = 0.5)
    # labels
    ax2.set_ylabel('Kälteleistung in -W')
    ax2.yaxis.set_major_formatter('{:.0f}'.format)
    ax2.set_ylim(0,6000)
    ax2.set_yticks(range(0,6001,1600))
    # ------------------------------------------------X------------------------------------------------
    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))],
                    rotation = 90)
    # ax1.set_xticks([i / 2 for i in range(0,7*2,1)], minor = True)
    ax1.set_xlim(von, bis)
    # ------------------------------------------------Global------------------------------------------------
    ax1.set_title(r'$\bf{Messwerte\ und\ Simulation\ des\ Wärmemengenzählers\ 4}$' + '\nZeitraum: ' + str(von.strftime('%Y-%m-%d')) + ' bis ' + str(bis.strftime('%Y-%m-%d')))
    ax1.grid(True)
    ax1.grid(which='minor', alpha = 0.3)
    ax1.legend(ax1.get_legend_handles_labels()[0]+ax2.get_legend_handles_labels()[0], 
                ax1.get_legend_handles_labels()[1]+ax2.get_legend_handles_labels()[1],
                loc = 'upper left')
    ax1.legend(fontsize = 8, loc = 'upper left')
    ax2.legend(fontsize = 8, loc = 'upper right') 
    
    # text = '\n'.join(f'{key}: {value}' for key, value in param['TABS'].items())
    # ax1.text(von+timedelta(hours = 2), 22, text, ha='center', va='center', fontsize=7)
    plt.show()

#%%% Plotting Energybalance
Q_mess = []
Q_sim = []
for i in df.index.day.unique():
    a = df[df.index.day == i]
    Q_mess.append(abs(a['HM4_Qpkt'].sum()/3600))
    Q_sim.append(abs(a['TABS_Qpkt'].sum()/3600))
# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------SINGLE-PLOT-------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(6.4, 4.8/2))
# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.bar([i-0.2 for i in range(1,8,1)], Q_mess, width = 0.4, label = 'Messwert')
ax1.bar([i+0.2 for i in range(1,8,1)], Q_sim, width = 0.4, label = 'Simulationswert')
# labels
ax1.set_ylabel('Kältemenge in -Wh', fontsize = 8)
# ax1.yaxis.set_major_formatter('{:.1f}'.format)
ax1.set_yticks(range(0,25,2), minor = True)
ax1.set_yticks(range(0,25,4))
ax1.set_ylim(0,25)
# ------------------------------------------------X------------------------------------------------
ax1.set_xticks(range(1,8,1), ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], rotation = 60,fontsize = 8)
# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Tagesbilanzen\ der\ abgegebenen\ Wärmemenge}$' + '\nZeitraum: 27.12.2021 - 02.01.2022')
ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
ax1.legend(fontsize = 8, loc = 'upper left')
ax1.text(6.1, 18.2,  r'$\bf{Gesamte\ Woche:}$' + '\n' + str(round(df['HM4_Qpkt'].sum()/3600,2)) + ' Wh gemessen'
          + '\n' + str(round(df['TABS_Qpkt'].sum()/3600,2)) + ' Wh simuliert' + '\n' +
          str(round((df['TABS_Qpkt'].sum()/3600-df['HM4_Qpkt'].sum()/3600)/(df['HM4_Qpkt'].sum()/3600)*100,2)) + ' % Abweichung',
          fontsize = 7, bbox=dict(facecolor='white', edgecolor='black'))
plt.show()




    


