# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 12:45:48 2023

@author: MoBueh
"""

#%% Directory
import os
os.chdir('C:/Users/49157/Documents/Python/')

#%% Libraries
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np
plt.rcParams["figure.dpi"] = 1000
plt.rc('axes', axisbelow=True)

#%% Sommer
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_summer.csv', index_col = 0)
df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')
df.loc[df['WP_COP']<  2, 'WP_COP'] = 0
df.loc[datetime.strptime('2022-08-11 15:00:00', '%Y-%m-%d %H:%M:%S'),'HM1_T_RET'] = 47

von = datetime.strptime('2022-08-11 13:15:00', '%Y-%m-%d %H:%M:%S')
bis = von +timedelta(hours = 6)

# ----------------------------------------------------------------------------------------------------------
# ------------------------------------------------COP-------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(6.4, 4.8/3))
# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.scatter(df.index, df['WP_COP'], c = 'black', s = 20)
ax1.plot(df.index, df['WP_COP'], linestyle = '--', c = 'black', alpha = 0.5)
# labels
ax1.set_ylabel('COP')
ax1.yaxis.set_major_formatter('{:.0f}'.format)
ax1.set_yticks(range(0,7,2))
ax1.set_yticks(range(0,7,1), minor = True)
# ------------------------------------------------X------------------------------------------------
ax1.set_xlim(von - timedelta(minutes = 10) , bis - timedelta(minutes = 5))
ax1.set_xticks([von + timedelta(minutes=i*15) for i in range(6*4)], minor = True)
ax1.set_xticks([von + timedelta(minutes=i*30) for i in range(6*2)], [(von + timedelta(minutes=i*30)).strftime('%H:%M') for i in range(6*2)], rotation = 90)
# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Coefficiant\ of\ Performance}$')
ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
# ax1.legend()
plt.show()

# ---------------------------------------------------------------------------------------------------------------
# ------------------------------------------------Speicher-------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(6.4, 4.8/2))
# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.scatter(df.index, df['WS_T'], label = 'Wärmespeicher', c = 'tab:red', s = 20)
ax1.plot(df.index, df['WS_T'], linestyle = '--', c = 'tab:red', alpha = 0.5)
ax1.scatter(df.index, df['KS_T'], label = 'Kältespeicher', c = 'tab:blue', s = 20)
ax1.plot(df.index, df['KS_T'], linestyle = '--', c = 'tab:blue', alpha = 0.5)
# labels
ax1.set_ylabel('Temperatur in °C')
ax1.yaxis.set_major_formatter('{:.0f}'.format)
ax1.set_yticks(range(0,60,10))
ax1.set_yticks(range(0,60,5), minor = True)
# ax1.set_ylim(0,6)
# ------------------------------------------------X------------------------------------------------
ax1.set_xlim(von - timedelta(minutes = 10) , bis - timedelta(minutes = 5))
ax1.set_xticks([von + timedelta(minutes=i*15) for i in range(6*4)], minor = True)
ax1.set_xticks([von + timedelta(minutes=i*30) for i in range(6*2)], [(von + timedelta(minutes=i*30)).strftime('%H:%M') for i in range(6*2)], rotation = 90)
# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Speichertemperaturen}$')
ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
ax1.legend(fontsize = 8)
plt.show()



# ---------------------------------------------------------------------------------------------------------------
# ------------------------------------------------Recooling------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(6.4, 4.8/2))
# fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]})
ax2 = plt.twinx()
# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.scatter(df.index, df['HM1_T_RET'], label = 'Rücklauf', c = 'tab:red', s = 20)
ax1.plot(df.index, df['HM1_T_RET'], linestyle = '--', c = 'tab:red', alpha = 0.5)
ax1.scatter(df.index, df['HM1_T_SUP'], label = 'Vorlauf', c = 'tab:blue', s = 20)
ax1.plot(df.index, df['HM1_T_SUP'], linestyle = '--', c = 'tab:blue', alpha = 0.5)
# labels
ax1.set_ylabel('Temperatur in °C')
ax1.yaxis.set_major_formatter('{:.0f}'.format)
ax1.set_yticks(range(0,60,10))
ax1.set_yticks(range(0,60,5), minor = True)
ax1.set_ylim(0,55)
# ------------------------------------------------Y2------------------------------------------------
# Data
ax2.scatter(df.index, df['HM1_Vpkt'], label = 'Volumenstrom', c = 'tab:green', s = 20)
ax2.plot(df.index, df['HM1_Vpkt'], linestyle = '--', c = 'tab:green', alpha = 0.5)
# labels☻
ax2.set_ylabel('Volumenstrom in l/h')
# ax2.yaxis.set_major_formatter('{:.1f}'.format)
ax2.set_ylim(0,1375)
# ------------------------------------------------X------------------------------------------------
ax1.set_xlim(von - timedelta(minutes = 10) , bis - timedelta(minutes = 5))
ax1.set_xticks([von + timedelta(minutes=i*15) for i in range(6*4)], minor = True)
ax1.set_xticks([von + timedelta(minutes=i*30) for i in range(6*2)], [(von + timedelta(minutes=i*30)).strftime('%H:%M') for i in range(6*2)], rotation = 60)
# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{PVT-Kollektoren}$')
ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
ax1.legend(fontsize = 8)
ax2.legend(fontsize = 8, loc = 'lower right')
plt.show()






#%% Winter
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_winter.csv', index_col = 0)
df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')
df.loc[df['HM10_Vpkt']<200, 'HM10_Vpkt'] = 0
df.loc[df['HM10_Vpkt']== 0, ['HM10_T_RET', 'HM10_T_SUP']] = np.nan
df.loc[df['WP_COP']<  2, 'WP_COP'] = 0
von = datetime.strptime('2021-12-29 07:00:00', '%Y-%m-%d %H:%M:%S')
bis = von +timedelta(hours = 6)

# ----------------------------------------------------------------------------------------------------------
# ------------------------------------------------COP-------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# fig, ax1 = plt.subplots(figsize=(6.4, 4.8/3))
# # ------------------------------------------------Y1------------------------------------------------
# # Data
# ax1.scatter(df.index, df['WP_COP'], c = 'black', s = 20)
# ax1.plot(df.index, df['WP_COP'], linestyle = '--', c = 'black', alpha = 0.5)
# # labels
# ax1.set_ylabel('COP')
# ax1.yaxis.set_major_formatter('{:.0f}'.format)
# ax1.set_yticks(range(0,7,2))
# ax1.set_yticks(range(0,7,1), minor = True)
# # ------------------------------------------------X------------------------------------------------
# ax1.set_xlim(von - timedelta(minutes = 10) , bis - timedelta(minutes = 5))
# ax1.set_xticks([von + timedelta(minutes=i*15) for i in range(6*4)], minor = True)
# ax1.set_xticks([von + timedelta(minutes=i*30) for i in range(6*2)], [(von + timedelta(minutes=i*30)).strftime('%H:%M') for i in range(6*2)], rotation = 90)
# # ------------------------------------------------Global------------------------------------------------
# ax1.set_title(r'$\bf{Coefficiant\ of\ Performance}$')
# ax1.grid(True)
# ax1.grid(which='minor', alpha = 0.3)
# # ax1.legend()
# plt.show()

# ---------------------------------------------------------------------------------------------------------------
# -------------------------------------------Speicher-------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(6.4, 4.8/2))
# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.scatter(df.index, df['WS_T'], label = 'Wärmespeicher', c = 'tab:red', s = 20)
ax1.plot(df.index, df['WS_T'], linestyle = '--', c = 'tab:red', alpha = 0.5)
# ax1.scatter(df.index, df['KS_T'], label = 'Kältespeicher', c = 'tab:blue', s = 20)
# ax1.plot(df.index, df['KS_T'], linestyle = '--', c = 'tab:blue', alpha = 0.5)
# labels
ax1.set_ylabel('Temperatur in °C')
ax1.yaxis.set_major_formatter('{:.0f}'.format)
ax1.set_yticks(range(20,30,2))
ax1.set_yticks(range(20,30,1), minor = True)
# ax1.set_ylim(0,6)
# ------------------------------------------------X------------------------------------------------
ax1.set_xlim(von - timedelta(minutes = 10) , bis - timedelta(minutes = 5))
ax1.set_xticks([von + timedelta(minutes=i*15) for i in range(6*4)], minor = True)
ax1.set_xticks([von + timedelta(minutes=i*30) for i in range(6*2)], [(von + timedelta(minutes=i*30)).strftime('%H:%M') for i in range(6*2)], rotation = 90)
# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Speichertemperatur}$')
ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
ax1.legend(fontsize = 8)
plt.show()
#%%

# ---------------------------------------------------------------------------------------------------------------
# ------------------------------------------------Zone------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(6.4, 4.8/2))
# fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]})
ax2 = plt.twinx()
# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.scatter(df.index, df['Zone_T'], label = 'Zone', c = 'black', s = 20)
ax1.plot(df.index, df['Zone_T'], linestyle = '--', c = 'black', alpha = 0.5)
ax1.scatter(df.index, df['HM4_T_SUP'], label = 'Vorlauf', c = 'tab:red', s = 20)
ax1.plot(df.index, df['HM4_T_SUP'], linestyle = '--', c = 'tab:red', alpha = 0.5)
ax1.scatter(df.index, df['HM4_T_RET'], label = 'Rücklauf', c = 'tab:blue', s = 20)
ax1.plot(df.index, df['HM4_T_RET'], linestyle = '--', c = 'tab:blue', alpha = 0.5)
# labels
ax1.set_ylabel('Temperatur in °C')
ax1.yaxis.set_major_formatter('{:.0f}'.format)
ax1.set_yticks(range(5,40,10))
ax1.set_yticks(range(5,40,5), minor = True)
# ax1.set_ylim(0,55)
# ------------------------------------------------Y2------------------------------------------------
# Data
ax2.scatter(df.index, df['HM4_Vpkt'], label = 'Volumenstrom', c = 'tab:green', s = 20)
ax2.plot(df.index, df['HM4_Vpkt'], linestyle = '--', c = 'tab:green', alpha = 0.5)
# labels
ax2.set_ylabel('Volumenstrom in l/h')
# ax2.yaxis.set_major_formatter('{:.1f}'.format)
ax2.set_ylim(0,1500)
# ------------------------------------------------X------------------------------------------------
ax1.set_xlim(von - timedelta(minutes = 10) , bis - timedelta(minutes = 5))
ax1.set_xticks([von + timedelta(minutes=i*15) for i in range(6*4)], minor = True)
ax1.set_xticks([von + timedelta(minutes=i*30) for i in range(6*2)], [(von + timedelta(minutes=i*30)).strftime('%H:%M') for i in range(6*2)], rotation = 60)
# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Innentemperatur\ und\ Fußbodenheizung}$')
ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
ax1.legend(fontsize = 8, loc = 'upper left')
ax2.legend(fontsize = 8, loc = 'upper right')
plt.show()
#%%
# ---------------------------------------------------------------------------------------------------------------
# ------------------------------------------------WP------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(6.4, 4.8/2))
# fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]})
ax2 = plt.twinx()
# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.plot(df.index[0], -99, label = r'$\bf{Temperaturen}$', c = 'white', alpha = 0)
ax1.scatter(df.index, df['HM10_T_SUP'], label = 'Vorlauf Senke', c = 'tab:red', s = 20)
ax1.plot(df.index, df['HM10_T_SUP'], linestyle = '--', c = 'tab:red', alpha = 0.5)
ax1.scatter(df.index, df['HM10_T_RET'], label = 'Rücklauf Senke', c = 'tab:blue', s = 20)
ax1.plot(df.index, df['HM10_T_RET'], linestyle = '--', c = 'tab:blue', alpha = 0.5)

ax1.scatter(df.index, df['HM1_T_SUP'], label = 'Vorlauf Quelle (PVT)', c = 'tab:red', s = 20, alpha = 0.5)
ax1.plot(df.index, df['HM1_T_SUP'], linestyle = '--', c = 'tab:red', alpha = 0.25)
ax1.scatter(df.index, df['HM1_T_RET'], label = 'Rücklauf Quelle (PVT)', c = 'tab:blue', s = 20, alpha = 0.5)
ax1.plot(df.index, df['HM1_T_RET'], linestyle = '--', c = 'tab:blue', alpha = 0.25)
# labels
ax1.set_ylabel('Temperatur in °C')
ax1.yaxis.set_major_formatter('{:.0f}'.format)
ax1.set_yticks(range(-50,40,10))
ax1.set_yticks(range(-50,40,5), minor = True)
ax1.set_ylim(-50,35)
# ------------------------------------------------Y2------------------------------------------------
# Data
ax2.plot(df.index[0], -99, label = r'$\bf{Volumenströme}$', c = 'white', alpha = 0)
ax2.scatter(df.index, df['HM1_Vpkt'], label = 'Quelle', c = 'tab:green', s = 20, alpha = 0.5)
ax2.plot(df.index, df['HM1_Vpkt'], linestyle = '--', c = 'tab:green', alpha = 0.25)

ax2.scatter(df.index, df['HM10_Vpkt'], label = 'Senke', c = 'tab:green', s = 20)
ax2.plot(df.index, df['HM10_Vpkt'], linestyle = '--', c = 'tab:green', alpha = 0.5)
# # labels
ax2.set_ylabel('Volumenstrom in l/h')
# # ax2.yaxis.set_major_formatter('{:.1f}'.format)
ax2.set_ylim(0,2125)
# ax2.set_yticks(range(0,2040,204))
# ------------------------------------------------X------------------------------------------------
ax1.set_xlim(von - timedelta(minutes = 10) , bis - timedelta(minutes = 5))
ax1.set_xticks([von + timedelta(minutes=i*15) for i in range(6*4)], minor = True)
ax1.set_xticks([von + timedelta(minutes=i*30) for i in range(6*2)], [(von + timedelta(minutes=i*30)).strftime('%H:%M') for i in range(6*2)], rotation = 60)
# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Wärmepumpe}$')
ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
ax1.legend(fontsize = 8, loc = 'upper left')
ax2.legend(fontsize = 8, loc = 'upper right')
plt.show()
