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
param['PVT'] = {
    'quantity'  : 6,                                                            # -         number of modules    
    'A'         : 1.8,                                                          # m²        Area of one module        
    'etha_opt'  : 0.468*0.5,                                                        # -         optical efficiency
    'U'         : 40,                                                           # W/m²/K    Coupling between T_mass and T_Cell
    'etha_el'   : 0.175,                                                        # -         electrical efficiency
    'beta_el'   : 0.0039,                                                       # -         correction for el. efficiency
    'c1'        : 7.6*3,                                                          # W/m²/K    Heat loss coefficient 
    'c5'        : 26.05*1000,                                                   # J/m²K     Effective heat capacity
    }

#%% Model
def PVT(theta_mean_prev, PVT_theta_in, PVT_Vpkt, theta_e, BAT_SOC, G):
    #                        _______________
    # theta_mean[-1] [°C]   |               |
    # --------------------->|               |
    # theta_in       [°C]   |               | theta_out     [°C]
    # --------------------->|               |--------------------->
    # Vpkt           [m³/h] |               | Pel           [W]
    # --------------------->|     PVT       |--------------------->
    # theta_AMB      [°C]   |               | Qpkt          [W]
    # --------------------->|               |--------------------->
    # G_POA          [W/m²] |               | theta_cell    [°C]       
    # --------------------->|               |--------------------->
    #                       |_______________|  
                                                 
    A = param['PVT']['etha_opt']*G                                               
    D = 0
    F = (param['PVT']['c1']+D+param['PVT']['c5']/900)*0.5
    E = PVT_Vpkt/3600*prop['Antifrogen']['rho']*prop['Antifrogen']['cp'] / (param['PVT']['A'] * param['PVT']['quantity']) 
    if (PVT_Vpkt == 0 or str(PVT_theta_in) == 'nan'):
        PVT_theta_in = theta_mean_prev
    PVT_theta_out = (A+theta_e*(D+ param['PVT']['c1']) + param['PVT']['c5']/900 * theta_mean_prev - PVT_theta_in*(F-E))/(E+F) 
    PVT_Qpkt = PVT_Vpkt/3600 *prop['Antifrogen']['rho'] * prop['Antifrogen']['cp'] * (PVT_theta_out-PVT_theta_in)                    
    PVT_theta_cell = PVT_Qpkt / (param['PVT']['A']*param['PVT']['quantity'])/param['PVT']['U'] + (PVT_theta_in+PVT_theta_out)/2
    if BAT_SOC == 1:
        PVT_Pel = 0    
    else: 
        PVT_Pel = G*param['PVT']['etha_el']*(1-param['PVT']['beta_el']*(PVT_theta_cell-25))*(param['PVT']['A']*param['PVT']['quantity'])   
    if PVT_Pel >= 1040:
        PVT_Pel = 1040
    return PVT_Pel,PVT_theta_in,PVT_theta_out,PVT_Qpkt,PVT_theta_cell

#%% Single Step
# PVT_Pel = [0]
# PVT_theta_in = [20]
# PVT_theta_out = [20]
# PVT_Qpkt = [0]
# PVT_theta_cell = [0]

# pvt = PVT(theta_mean_prev = (PVT_theta_in[-1]+PVT_theta_out[-1])/2,
#           PVT_theta_in = 20,                                        
#           PVT_Vpkt= 0.5,                                         
#           theta_e= 12,                                     
#           i_dir= 100,                          
#           i_dif= 100,                         
#           sun_h=40,                       
#           sun_az=160,                     
#           u_wind= 0,                        
#           BAT_SOC=0)                                         
# print(pvt)

#%% Winter
#%%% Simulation
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_winter.csv', index_col = 0)
df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')
# df.loc[df['HM1_Vpkt']<=300, 'HM1_Vpkt'] = 0
df.loc[(df['HM1_Vpkt'] <= 300) | (df['HM1_T_RET'] > df['AMB_T']), 'HM1_Vpkt'] = 0
df.loc[df['HM1_Vpkt'] == 0, ['HM1_T_RET', 'HM1_T_SUP']] = np.nan
df['HM1_Qpkt'] = [vpkt/1000/3600 * prop['Antifrogen']['rho'] *prop['Antifrogen']['cp'] * (SUP-RET) for vpkt, SUP, RET in zip(df['HM1_Vpkt'], df['HM1_T_SUP'],df['HM1_T_RET'])]



a = df[['HM1_Vpkt', 'HM1_T_RET', 'HM1_T_SUP']]

PVT_Pel = [0]
PVT_theta_in = [20]
PVT_theta_out = [20]
PVT_Qpkt = [0]
PVT_theta_cell = [20]


for i in df.index:
    pvt = PVT(theta_mean_prev = (PVT_theta_in[-1]+PVT_theta_out[-1])/2, 
              PVT_theta_in = df.loc[i, 'HM1_T_RET'], 
              PVT_Vpkt = df.loc[i, 'HM1_Vpkt']/1000, 
              theta_e = df.loc[i, 'AMB_T'], 
              BAT_SOC = 0, 
              G = df.loc[i, 'Pyr_hor'])
    PVT_Pel.append(pvt[0])
    PVT_theta_in.append(pvt[1])
    PVT_theta_out.append(pvt[2])
    PVT_Qpkt.append(pvt[3])
    PVT_theta_cell.append(pvt[4])

df['s_PVT_Pel'] = PVT_Pel[1:]
df['PVT_theta_in'] = PVT_theta_in[1:]
df['PVT_theta_out'] = PVT_theta_out[1:]
df['PVT_Qpkt'] = PVT_Qpkt[1:]
df['PVT_theta_cell'] = PVT_theta_cell[1:]
  
#%%% Therm Timeseries
Tag = 7
for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:Tag], [df.index[0] + timedelta(days = i+1) for i in range(7)][:Tag]):
    # ----------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(6.4, 4.8/1.3))
    # fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]})
    ax2 = plt.twinx()
    # ------------------------------------------------Y1------------------------------------------------
    ax1.scatter(df.index, df['HM1_T_SUP'], label = 'Messwert: Outflow', c = 'tab:red', s = 12, alpha = 0.5)
    ax1.plot(df.index, df['HM1_T_SUP'], c = 'tab:red', linestyle = '--', alpha = 0.5)
    
    ax1.scatter(df.index, df['HM1_T_RET'], label = 'Messwert: Inflow', c = 'tab:blue', s = 12, alpha = 0.5)
    ax1.plot(df.index, df['HM1_T_RET'], c = 'tab:blue', linestyle = '--', alpha = 0.5)
    
    ax1.scatter(df.index, df['PVT_theta_out'], label = 'Simulation: Outflow', c = 'tab:red', s = 12, marker = 'x')
    ax1.plot(df.index, df['PVT_theta_out'], c = 'tab:red', linestyle = '--', alpha = 0.5)
    
    ax1.plot(df.index, df['AMB_T'], c = 'tab:green', alpha = 0.7, label = 'Umgebungstemperatur')

    # labels
    ax1.set_ylabel('Temperatur in °C')
    ax1.yaxis.set_major_formatter('{:.0f}'.format)
    ax1.set_ylim(-10,25)
    ax1.set_yticks(np.arange(-10,25,2.5), minor = True)
    ax1.set_yticks(range(-10,25,5))
    
    # ------------------------------------------------Y2------------------------------------------------
    # Data
    ax2.plot(df.index, [i/1 for i in df['Pyr_hor']], c = 'tab:orange', alpha = 0.7, label = 'Globalstrahlung')
    
    # labels
    ax2.set_ylabel('Globalstrahlung in W/m²')
    ax2.yaxis.set_major_formatter('{:.0f}'.format)
    ax2.set_ylim(0,1050)
    ax2.set_yticks(range(0,1050,150))
    # ------------------------------------------------X------------------------------------------------
    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))],
                    rotation = 90)
    ax1.set_xlim(von, bis)
    # ------------------------------------------------Global------------------------------------------------
    ax1.set_title(r'$\bf{Messwerte\ und\ Simulation\ des\ Wärmemengenzählers\ 4}$' + '\nZeitraum: ' + str(von.strftime('%Y-%m-%d')) + ' bis ' + str(bis.strftime('%Y-%m-%d')))
    ax1.grid(True)
    ax1.grid(which='minor', alpha = 0.3)
    ax1.legend(fontsize = 8, loc = 'upper left', ncol = 2)
    ax2.legend(fontsize = 8, loc = 'upper right', ncol = 2)
    
    plt.show()
#%%% Therm Energybalance
Q_mess = []
Q_sim = []

for i in df.index.day.unique():
    a = df[df.index.day == i]
    Q_mess.append(a['HM1_Qpkt'].sum()/3600)
    Q_sim.append(a['PVT_Qpkt'].sum()/3600)

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

ax1.text(6.1, 18.2,  r'$\bf{Gesamte\ Woche:}$' + '\n' + str(round(df['HM1_Qpkt'].sum()/3600,2)) + ' Wh gemessen'
          + '\n' + str(round(df['PVT_Qpkt'].sum()/3600,2)) + ' Wh simuliert' + '\n' +
          str(round((df['PVT_Qpkt'].sum()/3600-df['HM1_Qpkt'].sum()/3600)/(df['HM1_Qpkt'].sum()/3600)*100,2)) + ' % Abweichung',
          fontsize = 7, bbox=dict(facecolor='white', edgecolor='black'))
plt.show()    
    
#%%% El. Timeseries
Tag = 7
for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:Tag], [df.index[0] + timedelta(days = i+1) for i in range(7)][:Tag]):
    # ----------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(6.4, 4.8/1.3))
    # fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]})
    ax2 = plt.twinx()
    # ------------------------------------------------Y1------------------------------------------------
    ax1.scatter(df.index, [i*1000 for i in df['PVT_Pel']], label = 'Messwert: Outflow', c = 'tab:red', s = 12, alpha = 0.5)
    ax1.plot(df.index, [i*1000 for i in df['PVT_Pel']], c = 'tab:red', linestyle = '--', alpha = 0.5)
    
    
    ax1.scatter(df.index, df['s_PVT_Pel'], label = 'Simulation: Outflow', c = 'tab:red', s = 12, marker = 'x')
    ax1.plot(df.index, df['s_PVT_Pel'], c = 'tab:red', linestyle = '--')
    

    # labels
    ax1.set_ylabel('Temperatur in °C')
    ax1.yaxis.set_major_formatter('{:.0f}'.format)
    # ax1.set_ylim(-10,25)
    # ax1.set_yticks(np.arange(-10,25,2.5), minor = True)
    # ax1.set_yticks(range(-10,25,5))
    
    # ------------------------------------------------Y2------------------------------------------------
    # Data
    # ax2.plot(df.index, [i/1 for i in df['Pyr_hor']], c = 'tab:orange', alpha = 0.7, label = 'Globalstrahlung')
    
    # # labels
    # ax2.set_ylabel('Globalstrahlung in W/m²')
    # ax2.yaxis.set_major_formatter('{:.0f}'.format)
    # ax2.set_ylim(0,1050)
    # ax2.set_yticks(range(0,1050,150))
    # ------------------------------------------------X------------------------------------------------
    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))],
                    rotation = 90)
    ax1.set_xlim(von, bis)
    # ------------------------------------------------Global------------------------------------------------
    ax1.set_title(r'$\bf{Messwerte\ und\ Simulation\ des\ Wärmemengenzählers\ 4}$' + '\nZeitraum: ' + str(von.strftime('%Y-%m-%d')) + ' bis ' + str(bis.strftime('%Y-%m-%d')))
    ax1.grid(True)
    ax1.grid(which='minor', alpha = 0.3)
    ax1.legend(fontsize = 8, loc = 'upper left', ncol = 2)
    ax2.legend(fontsize = 8, loc = 'upper right', ncol = 2)
    
    plt.show()

#%%% El Energybalance
Q_mess = []
Q_sim = []

for i in df.index.day.unique():
    a = df[df.index.day == i]
    Q_mess.append(a['PVT_Pel'].sum()/3600*1000)
    Q_sim.append(a['s_PVT_Pel'].sum()/3600)

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
ax1.set_yticks(range(0,4,1), minor = True)
ax1.set_yticks(range(0,4,2))
ax1.set_ylim(0,4)

# ------------------------------------------------X------------------------------------------------
ax1.set_xticks(range(1,8,1), ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], rotation = 60,fontsize = 8)

# ax1.set_xticks([i / 2 for i in range(0,7*2,1)], minor = True)

# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Tagesbilanzen\ der\ abgegebenen\ Wärmemenge}$' + '\nZeitraum: 27.12.2021 - 02.01.2022')

ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
ax1.legend(fontsize = 8, loc = 'upper left')

ax1.text(6.1, 3.2,  r'$\bf{Gesamte\ Woche:}$' + '\n' + str(round(df['PVT_Pel'].sum()/3.600,2)) + ' Wh gemessen'
          + '\n' + str(round(df['s_PVT_Pel'].sum()/3600,2)) + ' Wh simuliert' + '\n' +
          str(round((df['s_PVT_Pel'].sum()/3600-df['PVT_Pel'].sum()/3.600)/(df['PVT_Pel'].sum()/3.600)*100,2)) + ' % Abweichung',
          fontsize = 7, bbox=dict(facecolor='white', edgecolor='black'))
plt.show()    


#%% Sommer
#%%% Simulation
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_summer.csv', index_col = 0)
df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')
df.loc[(df['HM1_Vpkt'] <= 300) | (df['HM1_T_RET'] <= 31), 'HM1_Vpkt'] = 0
df.loc[df['HM1_Vpkt'] == 0, ['HM1_T_RET', 'HM1_T_SUP']] = np.nan

df['HM1_Qpkt'] = [vpkt/1000/3600 * prop['Antifrogen']['rho'] *prop['Antifrogen']['cp'] * (SUP-RET) for vpkt, SUP, RET in zip(df['HM1_Vpkt'], df['HM1_T_SUP'],df['HM1_T_RET'])]

PVT_Pel = [0]
PVT_theta_in = [20]
PVT_theta_out = [20]
PVT_Qpkt = [0]
PVT_theta_cell = [20]


for i in df.index:
    pvt = PVT(theta_mean_prev = (PVT_theta_in[-1]+PVT_theta_out[-1])/2, 
              PVT_theta_in = df.loc[i, 'HM1_T_RET'], 
              PVT_Vpkt = df.loc[i, 'HM1_Vpkt']/1000, 
              theta_e = df.loc[i, 'AMB_T'], 
              BAT_SOC = 0, 
              G = df.loc[i, 'Pyr_ver'])
    PVT_Pel.append(pvt[0])
    PVT_theta_in.append(pvt[1])
    PVT_theta_out.append(pvt[2])
    PVT_Qpkt.append(pvt[3])
    PVT_theta_cell.append(pvt[4])

df['s_PVT_Pel'] = PVT_Pel[1:]
df['PVT_theta_in'] = PVT_theta_in[1:]
df['PVT_theta_out'] = PVT_theta_out[1:]
df['PVT_Qpkt'] = PVT_Qpkt[1:]
df['PVT_theta_cell'] = PVT_theta_cell[1:]
  

#%%% Therm Timeseries
Tag = 7
for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:Tag], [df.index[0] + timedelta(days = i+1) for i in range(7)][:Tag]):
    # ----------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(6.4, 4.8/1.3))
    # fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]})
    ax2 = plt.twinx()
    # ------------------------------------------------Y1------------------------------------------------
    ax1.scatter(df.index, df['HM1_T_SUP'], label = 'Messwert: Outflow', c = 'tab:red', s = 12, alpha = 0.5)
    ax1.plot(df.index, df['HM1_T_SUP'], c = 'tab:red', linestyle = '--', alpha = 0.5)
    
    ax1.scatter(df.index, df['HM1_T_RET'], label = 'Messwert: Inflow', c = 'tab:blue', s = 12, alpha = 0.5)
    ax1.plot(df.index, df['HM1_T_RET'], c = 'tab:blue', linestyle = '--', alpha = 0.5)
    
    ax1.scatter(df.index, df['PVT_theta_out'], label = 'Simulation: Outflow', c = 'tab:red', s = 12, marker = 'x')
    ax1.plot(df.index, df['PVT_theta_out'], c = 'tab:red', linestyle = '--', alpha = 0.5)
    
    ax1.plot(df.index, df['AMB_T'], c = 'tab:green', alpha = 0.7, label = 'Umgebungstemperatur')

    # labels
    ax1.set_ylabel('Temperatur in °C')
    ax1.yaxis.set_major_formatter('{:.0f}'.format)
    ax1.set_ylim(0,60)
    ax1.set_yticks(np.arange(0,60,2.5), minor = True)
    ax1.set_yticks(range(0,60,5))
    
    # ------------------------------------------------Y2------------------------------------------------
    # Data
    ax2.plot(df.index, [i/1 for i in df['Pyr_hor']], c = 'tab:orange', alpha = 0.7, label = 'Globalstrahlung')
    
    # labels
    ax2.set_ylabel('Globalstrahlung in W/m²')
    ax2.yaxis.set_major_formatter('{:.0f}'.format)
    ax2.set_ylim(0,2100)
    ax2.set_yticks(range(0,2100,300))
    # ------------------------------------------------X------------------------------------------------
    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))],
                    rotation = 90)
    ax1.set_xlim(von, bis)
    # ------------------------------------------------Global------------------------------------------------
    ax1.set_title(r'$\bf{Messwerte\ und\ Simulation\ des\ Wärmemengenzählers\ 4}$' + '\nZeitraum: ' + str(von.strftime('%Y-%m-%d')) + ' bis ' + str(bis.strftime('%Y-%m-%d')))
    ax1.grid(True)
    ax1.grid(which='minor', alpha = 0.3)
    ax1.legend(fontsize = 8, loc = 'upper left', ncol = 1)
    ax2.legend(fontsize = 8, loc = 'upper right', ncol = 2)
    
    plt.show()


#%%% Therm Energybalance
Q_mess = []
Q_sim = []

for i in df.index.day.unique():
    a = df[df.index.day == i]
    Q_mess.append(abs(a['HM1_Qpkt'].sum()/3600))
    Q_sim.append(abs(a['PVT_Qpkt'].sum()/3600))

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

ax1.text(6.1, 18.2,  r'$\bf{Gesamte\ Woche:}$' + '\n' + str(round(df['HM1_Qpkt'].sum()/3600,2)) + ' Wh gemessen'
          + '\n' + str(round(df['PVT_Qpkt'].sum()/3600,2)) + ' Wh simuliert' + '\n' +
          str(round((df['PVT_Qpkt'].sum()/3600-df['HM1_Qpkt'].sum()/3600)/(df['HM1_Qpkt'].sum()/3600)*100,2)) + ' % Abweichung',
          fontsize = 7, bbox=dict(facecolor='white', edgecolor='black'))
plt.show() 


#%%% El. Timeseries
Tag = 7
for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:Tag], [df.index[0] + timedelta(days = i+1) for i in range(7)][:Tag]):
    # ----------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(6.4, 4.8/1.3))
    # fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]})
    ax2 = plt.twinx()
    # ------------------------------------------------Y1------------------------------------------------
    ax1.scatter(df.index, [i*1000 for i in df['PVT_Pel']], label = 'Messwert: Outflow', c = 'tab:red', s = 12, alpha = 0.5)
    ax1.plot(df.index, [i*1000 for i in df['PVT_Pel']], c = 'tab:red', linestyle = '--', alpha = 0.5)
    
    
    ax1.scatter(df.index, df['s_PVT_Pel'], label = 'Simulation: Outflow', c = 'tab:red', s = 12, marker = 'x')
    ax1.plot(df.index, df['s_PVT_Pel'], c = 'tab:red', linestyle = '--')
    

    # labels
    ax1.set_ylabel('Temperatur in °C')
    ax1.yaxis.set_major_formatter('{:.0f}'.format)
    ax1.set_ylim(0,1750)
    # ax1.set_yticks(np.arange(-10,25,2.5), minor = True)
    # ax1.set_yticks(range(-10,25,5))
    
    # ------------------------------------------------Y2------------------------------------------------
    # Data
    ax2.plot(df.index, [i/1 for i in df['BAT_U']], c = 'tab:blue', alpha = 0.7, label = 'Batteriespannung')
    # # labels
    # ax2.set_ylabel('Globalstrahlung in W/m²')
    # ax2.yaxis.set_major_formatter('{:.0f}'.format)
    # ax2.set_ylim(0,1050)
    # ax2.set_yticks(range(0,1050,150))
    # ------------------------------------------------X------------------------------------------------
    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))],
                    rotation = 90)
    ax1.set_xlim(von, bis)
    # ------------------------------------------------Global------------------------------------------------
    ax1.set_title(r'$\bf{Messwerte\ und\ Simulation\ des\ Wärmemengenzählers\ 4}$' + '\nZeitraum: ' + str(von.strftime('%Y-%m-%d')) + ' bis ' + str(bis.strftime('%Y-%m-%d')))
    ax1.grid(True)
    ax1.grid(which='minor', alpha = 0.3)
    ax1.legend(fontsize = 8, loc = 'upper left', ncol = 2)
    ax2.legend(fontsize = 8, loc = 'upper right', ncol = 2)
    
    plt.show()

#%%% El Energybalance
Q_mess = []
Q_sim = []

for i in df.index.day.unique():
    a = df[df.index.day == i]
    Q_mess.append(a['PVT_Pel'].sum()/3600*1000)
    Q_sim.append(a['s_PVT_Pel'].sum()/3600)

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
ax1.set_yticks(range(0,18,1), minor = True)
ax1.set_yticks(range(0,18,2))
ax1.set_ylim(0,18)

# ------------------------------------------------X------------------------------------------------
ax1.set_xticks(range(1,8,1), ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], rotation = 60,fontsize = 8)

# ax1.set_xticks([i / 2 for i in range(0,7*2,1)], minor = True)

# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Tagesbilanzen\ der\ abgegebenen\ Wärmemenge}$' + '\nZeitraum: 27.12.2021 - 02.01.2022')

ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
ax1.legend(fontsize = 8, loc = 'upper left')

ax1.text(6.1, 3.2,  r'$\bf{Gesamte\ Woche:}$' + '\n' + str(round(df['PVT_Pel'].sum()/3.600,2)) + ' Wh gemessen'
          + '\n' + str(round(df['s_PVT_Pel'].sum()/3600,2)) + ' Wh simuliert' + '\n' +
          str(round((df['s_PVT_Pel'].sum()/3600-df['PVT_Pel'].sum()/3.600)/(df['PVT_Pel'].sum()/3.600)*100,2)) + ' % Abweichung',
          fontsize = 7, bbox=dict(facecolor='white', edgecolor='black'))
plt.show()



