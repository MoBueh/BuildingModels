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

#%% Analyse Parameter
#%%% Delta T
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Modelle/02_WP/ExWirkungsgrad.csv', index_col = 0)
df = df[df['48:TempWP1COP'] >0]
df.index = pd.to_datetime(df.index)
for i in df.index:
    df.loc[i, 'delta_t'] = df.loc[i,'44:TempWP1Vorlauf'] - df.loc[i,'45:TempWP1Ruecklauf']

fig, ax1 = plt.subplots(figsize=(6.4, 4.8/3))
# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.boxplot(df['delta_t'], vert = False)
# labels
ax1.set_yticks([])
ax1.set_ylim(0.9,1.1)
# ------------------------------------------------X------------------------------------------------
ax1.set_xlim(4,6)
ax1.set_xlabel('ΔT in °C')
ax1.set_xticks(np.arange(4, 6.01, 0.25))
ax1.set_xticks(np.arange(4, 6.01, 0.05), minor = True)
# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Verteilung\ der\ Temperaturdifferenz\ auf\ der \ warmen \ Seite}$' + '\nZeitraum: 01.01.2022 - 31.01.2022')
ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
# ax1.legend(fontsize = 8)
plt.show()

#%%% Gütegrad
for i in df.index:
    df.loc[i, 'COP_the'] = (df.loc[i,'44:TempWP1Vorlauf']+273.15) / (df.loc[i,'44:TempWP1Vorlauf'] - df.loc[i,'43:TempWP1Quelleein'])

for i in df.index:
    df.loc[i, 'Guetegrad'] = df.loc[i, '48:TempWP1COP'] / df.loc[i, 'COP_the']*0.93

fig, ax1 = plt.subplots(figsize=(6.4, 4.8/3))
# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.boxplot(df['Guetegrad'], vert = False)
# labels
ax1.set_yticks([])
ax1.set_ylim(0.9,1.1)
# ------------------------------------------------X------------------------------------------------

# ax1.set_xlim(0.35,0.7)

ax1.set_xticks(np.arange(0.35, 0.575, 0.025),np.arange(35, 57.5, 2.5))
ax1.set_xticks(np.arange(0.35, 0.575, 0.005), minor = True)
ax1.set_xlim(0.35,0.575)
ax1.set_xlabel('\u03B7' + ' in %')
# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Verteilung\ des\ Gütegrades}$' + '\nZeitraum: 01.01.2022 - 31.01.2022')

ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
plt.show()

print('Gütegrad für diesen Monat: ' + str(round(df['Guetegrad'].mean(),2)))

#%%% prim. Volume
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Modelle/02_WP/VolumeFlow_HM1.csv', index_col = 0)
df = df[df['HM1_Vpkt']>300]

fig, ax1 = plt.subplots(figsize=(6.4, 4.8/3))
# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.boxplot(df['HM1_Vpkt'], vert = False)
# labels
ax1.set_yticks([])
# ax1.set_ylim(0.9,1.1)
# ------------------------------------------------X------------------------------------------------
# ax1.set_xlim(4,6)
# ax1.set_xlabel('ΔT in °C')
# ax1.set_xticks(np.arange(4, 6.01, 0.25))
# ax1.set_xticks(np.arange(4, 6.01, 0.05), minor = True)
# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Verteilung\ des\ Volumenstromes\ auf\ der \ kalten \ Seite}$' + '\nZeitraum: 01.01.2022 - 31.03.2022')
ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
# ax1.legend(fontsize = 8)
plt.show()

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
param['WP'] = {
    'etha'  :   0.47,
    }

#%% Model
def WP(WP_prim_theta_in, WP_sek_theta_in, WP_sek_delta_theta, WP_prim_Vpkt, WP_sek_Vpkt, WP_OI):
    #                        _______________
    # Prim_theta_in [°C]    |               | COP            [-]
    # --------------------->|               |--------------------->
    # Prim_theta_in [°C]    |               | Pel            [W]
    # --------------------->|               |--------------------->
    # ΔSek_theta    [°C]    |               | Qpkt_prim      [W]
    # --------------------->|               |--------------------->
    # Prim_Vpkt     [m³/h]  |     WP        | Prim_theta_out [°C]
    # --------------------->|               |--------------------->
    # Sek_theta_in  [°C]    |               | Qpkt_sek       [W] 
    # --------------------->|               |--------------------->
    # Sek_Vpkt      [°C]    |               | PSek_theta_out [°C]
    # --------------------->|_______________|--------------------->
    if WP_OI == 0:
        WP_COP = np.nan
        WP_Pel = np.nan
        WP_Qdot_prim = np.nan
        WP_Qdot_sec = np.nan
        WP_prim_theta_out = np.nan 
        WP_sek_theta_out = np.nan
        WP_prim_Vpkt = np.nan
        WP_sek_Vpkt = np.nan
    else:
        WP_prim_Vpkt = WP_prim_Vpkt
        WP_sek_Vpkt = WP_sek_Vpkt  
        WP_sek_theta_out = WP_sek_theta_in + WP_sek_delta_theta
        WP_COP = param['WP']['etha'] * (WP_sek_theta_out+273.15)/(WP_sek_theta_out - WP_prim_theta_in)
        WP_Qdot_sec = prop['wasser']['rho']*WP_sek_Vpkt/3600* prop['wasser']['cp']*(WP_sek_theta_out-WP_sek_theta_in)
        WP_Pel = (WP_Qdot_sec/WP_COP)
        WP_Qdot_prim = (WP_Pel-WP_Qdot_sec)
        WP_prim_theta_out = WP_prim_theta_in-(-WP_Qdot_prim/(prop['Antifrogen']['rho']*(WP_prim_Vpkt/3600)*prop['Antifrogen']['cp']))   
    WP_prim_theta_in = WP_prim_theta_in
    WP_sek_theta_in = WP_sek_theta_in
    WP_prim_Vpkt = WP_prim_Vpkt
    WP_sek_Vpkt = WP_sek_Vpkt
    return WP_COP, WP_Pel, WP_Qdot_prim, WP_Qdot_sec, WP_prim_theta_in, WP_prim_theta_out, WP_sek_theta_in, WP_sek_theta_out, WP_prim_Vpkt, WP_sek_Vpkt
#%% Single Step
wp = WP(WP_prim_theta_in = 10, 
        WP_sek_theta_in = 50, 
        WP_sek_delta_theta = 5,
        WP_prim_Vpkt = 0.3,
        WP_sek_Vpkt = 0.3,
        WP_OI = 1)

#%% Winter
#%%% Simulation
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_winter.csv', index_col = 0)
df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')
df.loc[df['WP_COP']<=3, 'WP_COP'] = 0
df.loc[df['WP_COP']==0, ['WP_COP','WP_sek_T_RET', 'WP_sek_T_SUP', 'WP_prim_T_RET', 'WP_prim_T_SUP', 'WP_Pel', 'WP_Qpkt']] = np.nan
a = df[['WP_sek_T_RET', 'WP_sek_T_SUP', 'WP_prim_T_RET', 'WP_prim_T_SUP', 'WP_COP', 'HM1_Vpkt','HM10_Vpkt']]

WP_COP = []
WP_Pel = []
WP_Qdot_prim = []
WP_Qdot_sec = []
WP_prim_theta_in = []
WP_prim_theta_out = []
WP_sek_theta_in = []
WP_sek_theta_out = []
WP_prim_Vpkt = []
WP_sek_Vpkt = []

for i in df.index:
    if df.loc[i, 'WP_COP'] == 0:
        WP_OI = 0
    else:
        WP_OI = 1
    wp = WP(WP_prim_theta_in = df.loc[i, 'WP_prim_T_SUP'], 
            WP_sek_theta_in = df.loc[i, 'WP_sek_T_RET'], 
            # WP_sek_delta_theta = df.loc[i, 'WP_sek_T_SUP']-df.loc[i, 'WP_sek_T_RET'], 
            WP_sek_delta_theta = 5, 
            WP_prim_Vpkt = df.loc[i, 'HM1_Vpkt']/1000, 
            WP_sek_Vpkt = df.loc[i, 'HM10_Vpkt']/1000, 
            WP_OI = WP_OI)
    
    WP_COP.append(wp[0])
    WP_Pel.append(wp[1])
    WP_Qdot_prim.append(wp[2])
    WP_Qdot_sec.append(wp[3])
    WP_prim_theta_in.append(wp[4])
    WP_prim_theta_out.append(wp[5])
    WP_sek_theta_in.append(wp[6])
    WP_sek_theta_out.append(wp[7])
    WP_prim_Vpkt.append(wp[8])
    WP_sek_Vpkt.append(wp[9])

df['s_WP_COP'] = WP_COP
df['s_WP_Pel'] = WP_Pel
df['s_WP_Qdot_prim'] = WP_Qdot_prim
df['s_WP_Qdot_sec'] = WP_Qdot_sec
df['s_WP_prim_theta_in'] = WP_prim_theta_in
df['s_WP_prim_theta_out'] = WP_prim_theta_out
df['s_WP_sek_theta_in'] = WP_sek_theta_in
df['s_WP_sek_theta_out'] = WP_sek_theta_out
df['s_WP_prim_Vpkt'] = WP_prim_Vpkt
df['s_WP_sek_Vpkt'] = WP_sek_Vpkt

#%%% Plotting TimeSeries
Tage = 7
for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:Tage], [df.index[0] + timedelta(days = i+1) for i in range(7)][:Tage]):
    # ----------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.6, 0.6]},figsize=(6.4, 4.8*1.5))
    # ------------------------------------------------Y1------------------------------------------------
    # Data
    ax1.scatter(df.index, df['WP_prim_T_SUP'], label = 'Messwert: Eintritt', c = 'tab:red', s = 12, alpha = 0.5)
    ax1.plot(df.index, df['WP_prim_T_SUP'], c = 'tab:red', linestyle = '--', alpha = 0.2)
    ax1.plot(df.index[0], -10, label = ' ', c = 'white')
    ax1.scatter(df.index, df['WP_prim_T_RET'], label = 'Messwert: Austritt', c = 'tab:blue', s = 12, alpha = 0.5)
    ax1.plot(df.index, df['WP_prim_T_RET'], c = 'tab:blue', linestyle = '--', alpha = 0.2)
    ax1.scatter(df.index, df['s_WP_prim_theta_out'], label = 'Simulation: Austritt', c = 'tab:blue', s = 12, marker = 'x')
    ax1.plot(df.index, df['s_WP_prim_theta_out'], c = 'tab:blue', linestyle = '--')
    
    ax1.set_ylabel('Temperatur in °C') 
    ax1.yaxis.set_major_formatter('{:.0f}'.format)
    ax1.set_ylim(-2,16)
    ax1.set_yticks(range(-2,16,2), minor = True)
    ax1.set_yticks(range(-2,16,4))
    ax1.set_title(r'$\bf{Kalte \ Seite}$' )
    ax1.grid(True)
    ax1.grid(which='minor', alpha = 0.3)
    ax1.legend(fontsize = 8, loc = 'upper right', ncol = 2)
    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], ['']*84,
                    rotation = 90)
    ax1.set_xlim(von, bis)
    ax1.set_title(r'$\bf{Temperaturen\ der\ Wärmepumpe}$' + '\nZeitraum: ' + str(von.strftime('%Y-%m-%d')) + ' bis ' + str(bis.strftime('%Y-%m-%d'))
                  + '\n' + '\n' + r'$\bf{Kalte \ Seite}$')
    
    
    
    ax2.scatter(df.index, df['WP_sek_T_RET'], label = 'Messwert: Eintritt', c = 'tab:blue', s = 12, alpha = 0.5)
    ax2.plot(df.index, df['WP_sek_T_RET'], c = 'tab:blue', linestyle = '--', alpha = 0.5)
    ax2.plot(df.index[0], -10, label = ' ', c = 'white')
    
    ax2.scatter(df.index, df['WP_sek_T_SUP'], label = 'Messwert: Austritt', c = 'tab:red', s = 12, alpha = 0.5)
    ax2.plot(df.index, df['WP_sek_T_SUP'], c = 'tab:red', linestyle = '--', alpha = 0.5)
            
    ax2.scatter(df.index, df['s_WP_sek_theta_out'], label = 'Simulation: Austritt', c = 'tab:red', s = 12, marker = 'x')
    ax2.plot(df.index, df['s_WP_sek_theta_out'], c = 'tab:red', linestyle = '--')   
    ax2.set_ylabel('Temperatur in °C') 
    ax2.yaxis.set_major_formatter('{:.0f}'.format)
    ax2.set_ylim(20,38)
    ax2.set_yticks(range(20,38,2), minor = True)
    ax2.set_yticks(range(20,38,4))
    ax2.set_title(r'$\bf{Warme \ Seite}$' )
    ax2.grid(True)
    ax2.grid(which='minor', alpha = 0.3)
    ax2.legend(fontsize = 8, loc = 'upper right', ncol = 2)
    ax2.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))], rotation = 90)
    ax2.set_xlim(von, bis)
    plt.show()

#%%% Leistungsdaten
for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:Tage], [df.index[0] + timedelta(days = i+1) for i in range(7)][:Tage]):
    # ----------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [1, 1, 1]},figsize=(6.4, 4.8*1.5))
    
    # ------------------------------------------------Y1------------------------------------------------
    # Data
    ax1.scatter(df.index, df['WP_COP'], label = 'Messwert', c = 'black', s = 12, alpha = 0.5)
    # ax1.plot(df.index, df['WP_COP'], c = 'black', linestyle = '--', alpha = 0.5)
    ax1.scatter(df.index, df['s_WP_COP'], label = 'Simulation', c = 'black', s = 12, marker = 'x')
    # ax1.plot(df.index, df['s_WP_COP'], c = 'black', linestyle = '--')
    # SetUp
    ax1.set_ylabel('COP in -') 
    ax1.yaxis.set_major_formatter('{:.1f}'.format)
    ax1.set_yticks(np.arange(0,16,1), minor = True)
    ax1.set_yticks(np.arange(0,16,3))
    ax1.set_ylim(0,11.9)
    ax1.set_title(r'$\bf{Leistungsdaten\ der\ Wärmepumpe}$' + '\nZeitraum: ' + str(von.strftime('%Y-%m-%d')) + ' bis ' + str(bis.strftime('%Y-%m-%d'))
                  + '\n' + '\n' + r'$\bf{COP}$')
    ax1.grid(True)
    ax1.grid(which='minor', alpha = 0.3)
    ax1.legend(fontsize = 8, loc = 'upper right')
    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], ['']*84,
                    rotation = 90)
    ax1.set_xlim(von, bis)
    
    
    
    ax2.scatter(df.index, df['WP_Pel'], label = 'Messwert', c = 'tab:red', s = 12, alpha = 0.5)
    # ax2.plot(df.index, df['WP_Pel'], c = 'tab:red', linestyle = '--', alpha = 0.5)
    ax2.scatter(df.index, [i/1000 for i in df['s_WP_Pel']], label = 'Simulation', c = 'tab:red', s = 12,marker = 'x')
    # ax2.plot(df.index, [i/1000 for i in df['s_WP_Pel']], c = 'tab:red', linestyle = '--')
    # SetUp
 
    ax2.set_ylabel('el. Leistung in kW') 
    ax2.yaxis.set_major_formatter('{:.1f}'.format)
    ax2.set_yticks(np.arange(0,1.25, 0.05), minor = True)
    ax2.set_yticks(np.arange(0, 1.25, 0.25))
    ax2.set_ylim(0,0.99)
    ax2.set_title(r'$\bf{Elektrische \ Leistungsaufnahme}$' )
    ax2.grid(True)
    ax2.grid(which='minor', alpha = 0.3)
    ax2.legend(fontsize = 8, loc = 'upper right')
    ax2.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], ['']*84,
                    rotation = 90)
    ax2.set_xlim(von, bis)


    ax3.scatter(df.index, df['WP_Qpkt'], label = 'Messwert', c = 'black', s = 12, alpha = 0.5)
    # ax3.plot(df.index, df['WP_Qpkt'], c = 'black', linestyle = '--', alpha = 0.5)
    ax3.scatter(df.index, [i/1000 for i in df['s_WP_Qdot_sec']], label = 'Simulation', c = 'black', s = 12, marker = 'x')
    # ax3.plot(df.index, df['s_WP_Qdot_sec'], c = 'black', linestyle = '--')
    ax3.set_ylabel('therm. Leistung in kW') 
    ax3.yaxis.set_major_formatter('{:.1f}'.format)
    ax3.set_yticks(np.arange(0,14,0.5), minor = True)
    ax3.set_yticks(range(0,14,2))
    ax3.set_ylim(0,5)
    ax3.set_title(r'$\bf{Thermische \ Leistung}$' )    
    ax3.grid(True)
    ax3.grid(which='minor', alpha = 0.3)
    ax3.legend(fontsize = 8, loc = 'upper right')
    ax3.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))], rotation = 90)

    ax3.set_xlim(von, bis)
    

    plt.show()


#%%% Plotting Energybalance
P_mess = []
P_sim = []

for i in df.index.day.unique():
    a = df[df.index.day == i]
    P_mess.append(a['WP_Pel'].sum()/3600*1000)
    P_sim.append(a['s_WP_Pel'].sum()/3600)

# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------SINGLE-PLOT-------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(6.4, 4.8/2))

# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.bar([i-0.2 for i in range(1,8,1)], P_mess, width = 0.4, label = 'Messwert')
ax1.bar([i+0.2 for i in range(1,8,1)], P_sim, width = 0.4, label = 'Simulationswert')
# labels
ax1.set_ylabel('el. Leistungsaufnahme in Wh', fontsize = 8)
# ax1.yaxis.set_major_formatter('{:.1f}'.format)
ax1.set_yticks(range(0,6,1), minor = True)
ax1.set_yticks(range(0,6,2))
ax1.set_ylim(0,6)

# ------------------------------------------------X------------------------------------------------
ax1.set_xticks(range(1,8,1), ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], rotation = 60,fontsize = 8)

# ax1.set_xticks([i / 2 for i in range(0,7*2,1)], minor = True)

# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Tagesbilanzen\ der\ elektrischen\ Leistungsaufnahme}$' + '\nZeitraum: 27.12.2021 - 02.01.2022')

ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
ax1.legend(fontsize = 8, loc = 'upper left')
ax1.text(6.05, 4.35,  r'$\bf{Gesamte\ Woche:}$' + '\n' + str(round(df['WP_Pel'].sum()/3600*1000,2)) + ' Wh gemessen'
          + '\n' + str(round(df['s_WP_Pel'].sum()/3600,2)) + ' Wh simuliert' + '\n' +
          str(round((df['s_WP_Pel'].sum()/3600-df['WP_Pel'].sum()/3600*1000)/(df['WP_Pel'].sum()/3600*1000)*100,2)) + ' % Abweichung',
          fontsize = 7, bbox=dict(facecolor='white', edgecolor='black'))
plt.show()

#%% Sommer
#%%% Simulation
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_summer.csv', index_col = 0)
df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')
df.loc[df['WP_COP']<=3, 'WP_COP'] = 0
df.loc[df['WP_COP']==0, ['WP_COP','WP_sek_T_RET', 'WP_sek_T_SUP', 'WP_prim_T_RET', 'WP_prim_T_SUP', 'WP_Pel', 'WP_Qpkt']] = np.nan
df.loc[df['WP_prim_T_SUP'] > 17, ['WP_prim_T_SUP', 'WP_prim_T_RET']] -= 2
a = df[['WP_sek_T_RET', 'WP_sek_T_SUP', 'WP_prim_T_RET', 'WP_prim_T_SUP', 'WP_COP', 'HM1_Vpkt', 'HM10_Vpkt', 'KS_T']]

WP_COP = []
WP_Pel = []
WP_Qdot_prim = []
WP_Qdot_sec = []
WP_prim_theta_in = []
WP_prim_theta_out = []
WP_sek_theta_in = []
WP_sek_theta_out = []
WP_prim_Vpkt = []
WP_sek_Vpkt = []

for i in df.index:
    print(i)
    if df.loc[i, 'WP_COP'] == 0:
        WP_OI = 0
    else:
        WP_OI = 1
    wp = WP(WP_prim_theta_in = df.loc[i, 'WP_prim_T_SUP'], 
            WP_sek_theta_in = df.loc[i, 'WP_sek_T_RET'], 
            WP_sek_delta_theta = 5, 
            WP_prim_Vpkt = 715/1000, 
            WP_sek_Vpkt = df.loc[i, 'HM10_Vpkt']/1000*1.1, 
            WP_OI = WP_OI)
    
    WP_COP.append(wp[0])
    WP_Pel.append(wp[1])
    WP_Qdot_prim.append(wp[2])
    WP_Qdot_sec.append(wp[3])
    WP_prim_theta_in.append(wp[4])
    WP_prim_theta_out.append(wp[5])
    WP_sek_theta_in.append(wp[6])
    WP_sek_theta_out.append(wp[7])
    WP_prim_Vpkt.append(wp[8])
    WP_sek_Vpkt.append(wp[9])

df['s_WP_COP'] = WP_COP
df['s_WP_Pel'] = WP_Pel
df['s_WP_Qdot_prim'] = WP_Qdot_prim
df['s_WP_Qdot_sec'] = WP_Qdot_sec
df['s_WP_prim_theta_in'] = WP_prim_theta_in
df['s_WP_prim_theta_out'] = WP_prim_theta_out
df['s_WP_sek_theta_in'] = WP_sek_theta_in
df['s_WP_sek_theta_out'] = WP_sek_theta_out
df['s_WP_prim_Vpkt'] = WP_prim_Vpkt
df['s_WP_sek_Vpkt'] = WP_sek_Vpkt

#%%% Plotting TimeSeries
Tage = 7
for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:Tage], [df.index[0] + timedelta(days = i+1) for i in range(7)][:Tage]):
    # ----------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.6, 0.6]},figsize=(6.4, 4.8*1.5))
    # ------------------------------------------------Y1------------------------------------------------
    # Data
    ax1.scatter(df.index, df['WP_prim_T_SUP'], label = 'Messwert: Eintritt', c = 'tab:red', s = 12, alpha = 0.5)
    ax1.plot(df.index, df['WP_prim_T_SUP'], c = 'tab:red', linestyle = '--', alpha = 0.2)
    ax1.plot(df.index[0], -10, label = ' ', c = 'white')
    ax1.scatter(df.index, df['WP_prim_T_RET'], label = 'Messwert: Austritt', c = 'tab:blue', s = 12, alpha = 0.5)
    ax1.plot(df.index, df['WP_prim_T_RET'], c = 'tab:blue', linestyle = '--', alpha = 0.2)
    ax1.scatter(df.index, df['s_WP_prim_theta_out'], label = 'Simulation: Austritt', c = 'tab:blue', s = 12, marker = 'x')
    ax1.plot(df.index, df['s_WP_prim_theta_out'], c = 'tab:blue', linestyle = '--')
    
    ax1.set_ylabel('Temperatur in °C') 
    ax1.yaxis.set_major_formatter('{:.0f}'.format)
    # ax1.set_ylim(-2,16)
    # ax1.set_yticks(range(-2,16,2), minor = True)
    # ax1.set_yticks(range(-2,16,4))
    ax1.set_title(r'$\bf{Kalte \ Seite}$' )
    ax1.grid(True)
    ax1.grid(which='minor', alpha = 0.3)
    ax1.legend(fontsize = 8, loc = 'upper right', ncol = 2)
    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], ['']*84,
                    rotation = 90)
    ax1.set_xlim(von, bis)
    ax1.set_title(r'$\bf{Temperaturen\ der\ Wärmepumpe}$' + '\nZeitraum: ' + str(von.strftime('%Y-%m-%d')) + ' bis ' + str(bis.strftime('%Y-%m-%d'))
                  + '\n' + '\n' + r'$\bf{Kalte \ Seite}$')
    
    
    
    ax2.scatter(df.index, df['WP_sek_T_RET'], label = 'Messwert: Eintritt', c = 'tab:blue', s = 12, alpha = 0.5)
    ax2.plot(df.index, df['WP_sek_T_RET'], c = 'tab:blue', linestyle = '--', alpha = 0.5)
    ax2.plot(df.index[0], -10, label = ' ', c = 'white')
    
    ax2.scatter(df.index, df['WP_sek_T_SUP'], label = 'Messwert: Austritt', c = 'tab:red', s = 12, alpha = 0.5)
    ax2.plot(df.index, df['WP_sek_T_SUP'], c = 'tab:red', linestyle = '--', alpha = 0.5)
            
    ax2.scatter(df.index, df['s_WP_sek_theta_out'], label = 'Simulation: Austritt', c = 'tab:red', s = 12, marker = 'x')
    ax2.plot(df.index, df['s_WP_sek_theta_out'], c = 'tab:red', linestyle = '--')   
    ax2.set_ylabel('Temperatur in °C') 
    ax2.yaxis.set_major_formatter('{:.0f}'.format)
    # ax2.set_ylim(20,38)
    # ax2.set_yticks(range(20,38,2), minor = True)
    # ax2.set_yticks(range(20,38,4))
    ax2.set_title(r'$\bf{Warme \ Seite}$' )
    ax2.grid(True)
    ax2.grid(which='minor', alpha = 0.3)
    ax2.legend(fontsize = 8, loc = 'upper right', ncol = 2)
    ax2.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))], rotation = 90)
    ax2.set_xlim(von, bis)
    plt.show()

#%%% Leistungsdaten
for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:Tage], [df.index[0] + timedelta(days = i+1) for i in range(7)][:Tage]):
    # ----------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [1, 1, 1]},figsize=(6.4, 4.8*1.5))
    
    # ------------------------------------------------Y1------------------------------------------------
    # Data
    ax1.scatter(df.index, df['WP_COP'], label = 'Messwert', c = 'black', s = 12, alpha = 0.5)
    # ax1.plot(df.index, df['WP_COP'], c = 'black', linestyle = '--', alpha = 0.5)
    ax1.scatter(df.index, df['s_WP_COP'], label = 'Simulation', c = 'black', s = 12, marker = 'x')
    # ax1.plot(df.index, df['s_WP_COP'], c = 'black', linestyle = '--')
    # SetUp
    ax1.set_ylabel('COP in -') 
    ax1.yaxis.set_major_formatter('{:.1f}'.format)
    ax1.set_yticks(np.arange(0,16,1), minor = True)
    ax1.set_yticks(np.arange(0,16,3))
    ax1.set_ylim(0,11.9)
    ax1.set_title(r'$\bf{Leistungsdaten\ der\ Wärmepumpe}$' + '\nZeitraum: ' + str(von.strftime('%Y-%m-%d')) + ' bis ' + str(bis.strftime('%Y-%m-%d'))
                  + '\n' + '\n' + r'$\bf{COP}$')
    ax1.grid(True)
    ax1.grid(which='minor', alpha = 0.3)
    ax1.legend(fontsize = 8, loc = 'upper right')
    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], ['']*84,
                    rotation = 90)
    ax1.set_xlim(von, bis)
    
    
    
    ax2.scatter(df.index, df['WP_Pel'], label = 'Messwert', c = 'tab:red', s = 12, alpha = 0.5)
    # ax2.plot(df.index, df['WP_Pel'], c = 'tab:red', linestyle = '--', alpha = 0.5)
    ax2.scatter(df.index, [i/1000 for i in df['s_WP_Pel']], label = 'Simulation', c = 'tab:red', s = 12,marker = 'x')
    # ax2.plot(df.index, [i/1000 for i in df['s_WP_Pel']], c = 'tab:red', linestyle = '--')
    # SetUp
 
    ax2.set_ylabel('el. Leistung in kW') 
    ax2.yaxis.set_major_formatter('{:.1f}'.format)
    ax2.set_yticks(np.arange(0,1.25, 0.05), minor = True)
    ax2.set_yticks(np.arange(0, 1.25, 0.25))
    ax2.set_ylim(0,0.99)
    ax2.set_title(r'$\bf{Elektrische \ Leistungsaufnahme}$' )
    ax2.grid(True)
    ax2.grid(which='minor', alpha = 0.3)
    ax2.legend(fontsize = 8, loc = 'upper right')
    ax2.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], ['']*84,
                    rotation = 90)
    ax2.set_xlim(von, bis)


    ax3.scatter(df.index, df['WP_Qpkt'], label = 'Messwert', c = 'black', s = 12, alpha = 0.5)
    # ax3.plot(df.index, df['WP_Qpkt'], c = 'black', linestyle = '--', alpha = 0.5)
    ax3.scatter(df.index, [i/1000 for i in df['s_WP_Qdot_sec']], label = 'Simulation', c = 'black', s = 12, marker = 'x')
    # ax3.plot(df.index, df['s_WP_Qdot_sec'], c = 'black', linestyle = '--')
    ax3.set_ylabel('therm. Leistung in kW') 
    ax3.yaxis.set_major_formatter('{:.1f}'.format)
    ax3.set_yticks(np.arange(0,14,0.5), minor = True)
    ax3.set_yticks(range(0,14,2))
    ax3.set_ylim(0,5)
    ax3.set_title(r'$\bf{Thermische \ Leistung}$' )    
    ax3.grid(True)
    ax3.grid(which='minor', alpha = 0.3)
    ax3.legend(fontsize = 8, loc = 'upper right')
    ax3.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))], rotation = 90)

    ax3.set_xlim(von, bis)
    plt.show()


#%%% Plotting Energybalance
P_mess = []
P_sim = []

for i in df.index.day.unique():
    a = df[df.index.day == i]
    P_mess.append(a['WP_Pel'].sum()/3600*1000)
    P_sim.append(a['s_WP_Pel'].sum()/3600)

# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------SINGLE-PLOT-------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(6.4, 4.8/2))

# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.bar([i-0.2 for i in range(1,8,1)], P_mess, width = 0.4, label = 'Messwert')
ax1.bar([i+0.2 for i in range(1,8,1)], P_sim, width = 0.4, label = 'Simulationswert')
# labels
ax1.set_ylabel('el. Leistungsaufnahme in Wh', fontsize = 8)
# ax1.yaxis.set_major_formatter('{:.1f}'.format)
ax1.set_yticks(range(0,9,1), minor = True)
ax1.set_yticks(range(0,9,2))
ax1.set_ylim(0,9)

# ------------------------------------------------X------------------------------------------------
ax1.set_xticks(range(1,8,1), ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], rotation = 60,fontsize = 8)

# ax1.set_xticks([i / 2 for i in range(0,7*2,1)], minor = True)

# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Tagesbilanzen\ der\ elektrischen\ Leistungsaufnahme}$' + '\nZeitraum: 27.12.2021 - 02.01.2022')

ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
ax1.legend(fontsize = 8, loc = 'upper left')
ax1.text(6.05, 6.5,  r'$\bf{Gesamte\ Woche:}$' + '\n' + str(round(df['WP_Pel'].sum()/3600*1000,2)) + ' Wh gemessen'
          + '\n' + str(round(df['s_WP_Pel'].sum()/3600,2)) + ' Wh simuliert' + '\n' +
          str(round((df['s_WP_Pel'].sum()/3600-df['WP_Pel'].sum()/3600*1000)/(df['WP_Pel'].sum()/3600*1000)*100,2)) + ' % Abweichung',
          fontsize = 7, bbox=dict(facecolor='white', edgecolor='black'))
plt.show()
    
    