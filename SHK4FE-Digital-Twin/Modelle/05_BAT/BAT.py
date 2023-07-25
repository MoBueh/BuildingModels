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

#%%
param = {}
param['BAT'] = {
    'C' : 8000                                                                  # Wh
    }

#%% Model
def BAT(BAT_Pel_in, BAT_Pel_out, BAT_SOC_prev):
    #                        _______________
    # Pel_in      [W]       |               |
    # --------------------->|               |
    # Pel_out     [W]       |               | SOC
    # --------------------->|     BAT       |--------------------->
    # SOC[-1]     [%]       |               | 
    # --------------------->|               |
    #                       |_______________|
    C_prev = BAT_SOC_prev*param['BAT']['C']*3600
    C = (BAT_Pel_in-BAT_Pel_out)*900 + C_prev
    BAT_SOC = C/(param['BAT']['C']*3600)
    if BAT_SOC > 1:
        BAT_SOC = 1
    if BAT_SOC < 0:
        BAT_SOC = 0
    return BAT_SOC
#%% Winter
# df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_summer.csv', index_col = 0)
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_winter.csv', index_col = 0)
df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')

BAT_SOC = [1]
for i in df.index:
    bat = BAT(BAT_Pel_in = (df.loc[i, 'PVT_Pel'] + df.loc[i, 'INV_P_IN'])*1000, 
              BAT_Pel_out = df.loc[i, 'Bat_Pel']*1000, 
              BAT_SOC_prev = BAT_SOC[-1])
    BAT_SOC.append(bat)
df['BAT_SOC'] = BAT_SOC[1:]

# %%% PLotting Timeseries
Tag = 7
for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:Tag], [df.index[0] + timedelta(days = i+1) for i in range(7)][:Tag]):
    # ----------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------
    # fig, ax1 = plt.subplots(figsize=(6.4, 4.8/1.3))
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]})

    # ------------------------------------------------Y1------------------------------------------------
    # Data
    ax1.scatter(df.index, df['BAT_U'], label = 'Messwert: Batteriespannung', c = 'grey', s = 12, alpha = 0.5)
    ax1.plot(df.index, df['BAT_U'], c = 'grey', linestyle = '--', alpha = 0.5)
    
    
    # labels
    ax1.set_ylabel('Batteriespannung in V', fontsize = 8)
    # ax1.yaxis.set_major_formatter('{:.0f}'.format)
    # ax1.set_ylim(55.06,57.76)
    ax1.set_ylim(55,58)
    ax1.set_yticks(np.arange(55,58.01,0.25), minor = True, fontsize = 8)
    ax1.set_yticks(np.arange(55,58.01,0.5), fontsize = 8)
    

    
    
    ax1_twin = ax1.twinx()
    ax1_twin.scatter(df.index, df['BAT_SOC'], label = 'Simulation: SOC', c = 'black', s = 12, marker = 'x')
    ax1_twin.plot(df.index, df['BAT_SOC'], c = 'black', linestyle = '--')
    ax1_twin.set_ylim(0,1)
    ax1_twin.legend(fontsize = 6, loc = 'lower right')
    ax1_twin.set_yticks(np.arange(0,1.01,0.25), minor = True, fontsize = 8)
    ax1_twin.set_yticks(np.arange(0,1.01,0.5), np.arange(0,101,50), fontsize = 8)

    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], ['']*84,
                    rotation = 90)
    ax1.legend(fontsize = 6, loc = 'lower left')
    ax1.grid(True)
    ax1.set_xlim(von, bis)
    
    
    # ------------------------------------------------Y2------------------------------------------------
    ax2.scatter(df.index, [i *1 for i in df['Bat_Pel']], label = 'Out', c = 'tab:red', s = 12, marker = 'x')
    # ax2.plot(df.index, df['TABS_Qpkt'], c = 'tab:green', linestyle = '--')
    ax2.scatter(df.index, [i *1 for i in df['INV_P_IN']], label = 'In Grid', c = 'tab:green', s = 12, alpha = 0.5)
    ax2.scatter(df.index, [i *1 for i in df['PVT_Pel']], label = 'In PVT', c = 'tab:green', s = 12, alpha = 0.5)
    ax2.grid(True)
    ax2.set_xlim(von, bis)

    # labels
    ax2.set_ylabel('el. Leistung in kW', fontsize = 8)
    ax2.yaxis.set_major_formatter('{:.0f}'.format)
    # ax2.set_ylim(0,7000)
    # ax2.set_yticks(range(0,7000,800))
    
    
    # ------------------------------------------------X------------------------------------------------
    ax2.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))],
                    rotation = 90)
    ax2.set_xlim(von, bis)
    ax2.grid(which='minor', alpha = 0.3)
    ax2.legend(fontsize = 8, loc = 'upper right') 
    plt.show()
#%%% whole Week
# ----------------------------------------------------------------------------------------------------------------
# ------------------------------------------------TWIN-PLOT-------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
# fig, ax1 = plt.subplots(figsize=(6.4, 4.8/1.3))
fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]})

# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.plot(df.index, df['BAT_U'], label = 'Messwert: Batteriespannung', c = 'grey', alpha = 0.5)


# labels
ax1.set_ylabel('Batteriespannung in V', fontsize = 8)
# ax1.yaxis.set_major_formatter('{:.0f}'.format)
# ax1.set_ylim(55.06,57.76)
ax1.set_ylim(55,58)
ax1.set_yticks(np.arange(55,58.01,0.25), minor = True, fontsize = 8)
ax1.set_yticks(np.arange(55,58.01,0.5), fontsize = 8)




ax1_twin = ax1.twinx()
ax1_twin.plot(df.index, df['BAT_SOC'], label = 'Simulation: SOC', c = 'black')

ax1_twin.set_ylim(0,1)
ax1_twin.legend(fontsize = 6, loc = 'lower right')
ax1_twin.set_yticks(np.arange(0,1.01,0.25), minor = True, fontsize = 8)
ax1_twin.set_yticks(np.arange(0,1.01,0.5), np.arange(0,101,50), fontsize = 8)

ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], ['']*84,
                rotation = 90)
ax1.legend(fontsize = 6, loc = 'lower left')
ax1.grid(True)


# ------------------------------------------------Y2------------------------------------------------
ax2.plot(df.index, [i *1 for i in df['Bat_Pel']], label = 'Out', c = 'tab:red')
# ax2.plot(df.index, df['TABS_Qpkt'], c = 'tab:green', linestyle = '--')
ax2.plot(df.index, [i *1 for i in df['INV_P_IN']], label = 'In Grid', c = 'tab:green')
ax2.plot(df.index, [i *1 for i in df['PVT_Pel']], label = 'In PVT', c = 'tab:green')
ax2.grid(True)

# labels
ax2.set_ylabel('el. Leistung in kW', fontsize = 8)
ax2.yaxis.set_major_formatter('{:.0f}'.format)
# ax2.set_ylim(0,7000)
# ax2.set_yticks(range(0,7000,800))


# ------------------------------------------------X------------------------------------------------
ax2.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))],
                rotation = 90)
ax2.grid(which='minor', alpha = 0.3)
ax2.legend(fontsize = 8, loc = 'upper right') 
plt.show()


#%% Summer
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_summer.csv', index_col = 0)
# df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_winter.csv', index_col = 0)
df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')

BAT_SOC = [0.7]

for i in df.index:
    bat = BAT(BAT_Pel_in = (df.loc[i, 'PVT_Pel'] + df.loc[i, 'INV_P_IN'])*1000, 
              BAT_Pel_out = df.loc[i, 'Bat_Pel']*1000, 
              BAT_SOC_prev = BAT_SOC[-1])
    BAT_SOC.append(bat)

df['BAT_SOC'] = BAT_SOC[1:]


# %%% PLotting Timeseries
Tag = 7
for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:Tag], [df.index[0] + timedelta(days = i+1) for i in range(7)][:Tag]):
    # ----------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------
    # fig, ax1 = plt.subplots(figsize=(6.4, 4.8/1.3))
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]})

    # ------------------------------------------------Y1------------------------------------------------
    # Data
    ax1.scatter(df.index, df['BAT_U'], label = 'Messwert: Batteriespannung', c = 'grey', s = 12, alpha = 0.5)
    ax1.plot(df.index, df['BAT_U'], c = 'grey', linestyle = '--', alpha = 0.5)
    
    
    # labels
    ax1.set_ylabel('Batteriespannung in V', fontsize = 8)
    # ax1.yaxis.set_major_formatter('{:.0f}'.format)
    # ax1.set_ylim(55.06,57.76)
    ax1.set_ylim(55,58)
    ax1.set_yticks(np.arange(55,58.01,0.25), minor = True, fontsize = 8)
    ax1.set_yticks(np.arange(55,58.01,0.5), fontsize = 8)
    

    
    
    ax1_twin = ax1.twinx()
    ax1_twin.scatter(df.index, df['BAT_SOC'], label = 'Simulation: SOC', c = 'black', s = 12, marker = 'x')
    ax1_twin.plot(df.index, df['BAT_SOC'], c = 'black', linestyle = '--')
    ax1_twin.set_ylim(0,1)
    ax1_twin.legend(fontsize = 6, loc = 'lower right')
    ax1_twin.set_yticks(np.arange(0,1.01,0.25), minor = True, fontsize = 8)
    ax1_twin.set_yticks(np.arange(0,1.01,0.5), np.arange(0,101,50), fontsize = 8)

    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], ['']*84,
                    rotation = 90)
    ax1.legend(fontsize = 6, loc = 'lower left')
    ax1.grid(True)
    ax1.set_xlim(von, bis)
    
    
    # ------------------------------------------------Y2------------------------------------------------
    ax2.scatter(df.index, [i *1 for i in df['Bat_Pel']], label = 'Out', c = 'tab:red', s = 12, marker = 'x')
    # ax2.plot(df.index, df['TABS_Qpkt'], c = 'tab:green', linestyle = '--')
    ax2.scatter(df.index, [i *1 for i in df['INV_P_IN']], label = 'In Grid', c = 'tab:green', s = 12, alpha = 0.5)
    ax2.scatter(df.index, [i *1 for i in df['PVT_Pel']], label = 'In PVT', c = 'tab:green', s = 12, alpha = 0.5)
    ax2.grid(True)
    ax2.set_xlim(von, bis)

    # labels
    ax2.set_ylabel('el. Leistung in kW', fontsize = 8)
    ax2.yaxis.set_major_formatter('{:.0f}'.format)
    # ax2.set_ylim(0,7000)
    # ax2.set_yticks(range(0,7000,800))
    
    
    # ------------------------------------------------X------------------------------------------------
    ax2.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))],
                    rotation = 90)
    ax2.set_xlim(von, bis)
    ax2.grid(which='minor', alpha = 0.3)
    ax2.legend(fontsize = 8, loc = 'upper right') 
    plt.show()

#%%% whole Week
# ----------------------------------------------------------------------------------------------------------------
# ------------------------------------------------TWIN-PLOT-------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
# fig, ax1 = plt.subplots(figsize=(6.4, 4.8/1.3))
fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]})

# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.plot(df.index, df['BAT_U'], label = 'Messwert: Batteriespannung', c = 'grey', alpha = 0.5)

ax1_twin = ax1.twinx()
ax1_twin.plot(df.index, df['BAT_SOC'], label = 'Simulation: SOC', c = 'black')
ax1_twin.set_ylim(0,1)
ax1_twin.legend(fontsize = 6, loc = 'lower right')
ax1_twin.set_yticks(np.arange(0,1.01,0.166666), minor = True)
ax1_twin.set_yticks(np.arange(0,1.01,0.5), np.arange(0,101,50))

# labels
ax1.set_ylabel('Batteriespannung in V', fontsize = 8)
ax1.set_ylim(55,58)
ax1.set_yticks(np.arange(55,58.01,0.25), minor = True, fontsize = 8)
ax1.set_yticks(np.arange(55,58.01,0.5), fontsize = 8)

ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*4) for i in range(int(672/4/4))], minor = True)
ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*24) for i in range(int(672/4/24))], 
                ['']*7, rotation = 60)

ax1.legend(fontsize = 6, loc = 'lower left')
ax1.grid(True)
ax1.set_xlim(df.index[0], df.index[-1])

# ------------------------------------------------Y2------------------------------------------------
ax2.plot(df.index, [i *1 for i in df['Bat_Pel']], label = 'Out', c = 'tab:red')
# ax2.plot(df.index, df['TABS_Qpkt'], c = 'tab:green', linestyle = '--')
# ax2.plot(df.index, [i *1 for i in df['INV_P_IN']], label = 'In Grid', c = 'lightgreen')
ax2.plot(df.index, [i *1 for i in df['PVT_Pel']], label = 'In PVT', c = 'tab:green')
ax2.grid(True)
# labels
ax2.set_ylabel('el. Leistung in kW', fontsize = 8)
# ax2.yaxis.set_major_formatter('{:.0f}'.format)
ax2.set_ylim(0,1.2)
ax2.set_yticks(np.arange(0, 1.21, 0.1),fontsize = 8, minor = True)
ax2.set_yticks(np.arange(0, 1.21, 0.5),fontsize = 8)
ax2.legend(fontsize = 6, loc = 'upper left', ncol = 2) 


# ------------------------------------------------X------------------------------------------------
# ax2.set_xticks([df.index[0] + timedelta(minutes=i*60*4) for i in range(int(672/4/4))], minor = True)

ax2.set_xticks([df.index[0] + timedelta(minutes=i*60*24) for i in range(int(672/4/24))], 
                [(df.index[0] + timedelta(minutes=i*60*24)).strftime('%d-%m') for i in range(int(672/4/24))],
                rotation = 60)
ax2.grid(which='minor', alpha = 0.3)
ax2.set_xlim(df.index[0], df.index[-1])

plt.show()
    
    
    
    
    
    
    
    
    