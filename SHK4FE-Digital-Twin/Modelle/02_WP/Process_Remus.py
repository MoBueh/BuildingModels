# -*- coding: utf-8 -*-
"""
Created on Mon May 23 10:01:50 2022

@author: 49157
"""

#%% Directory
import os
os.chdir('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Modelle/02_WP')

#%% Libraries
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np
import math
plt.rcParams["figure.dpi"] = 400
plt.rc('axes', axisbelow=True)


#%% Merge Data
#%%% ExWirkungsgrad
columns = {}
# columns['HM1'] = ['42:TempWP1Quelleaus', '43:TempWP1Quelleein', '44:TempWP1Vorlauf', '45:TempWP1Ruecklauf','46:TempWP1PTherm', '47:TempWP1PElek', '48:TempWP1COP','98:1_Vdot(lh)', '97:1_T_RL(C)', '96:1_T_VL(C)']
columns['HM1'] = ['98:1_Vdot(lh)', '97:1_T_RL(C)', '96:1_T_VL(C)']
data = {'HM1': pd.DataFrame()}

j = 1
for file in os.listdir('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Modelle/02_WP/Remus')[:31]:
    print('adding: ' + file + ' ' + str(j) + '/' + str(len([name for name in os.listdir('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Modelle/02_WP/Remus')])))
    df = pd.read_table('Remus/' + file, sep = ' ')
    index = [datetime.strptime(df.iloc[i,0] + ' ' + df.iloc[i,1], '%Y.%m.%d %H:%M:%S') for i in df.index]
    for i in list(data):
        data[i] = data[i].append(pd.DataFrame(index = index, data = df[columns[i]].to_numpy(), columns=columns[i]))
    j = j+1
del i, file, index, df, j, columns

df = data['HM1']
del data

df1 = df.rename(columns={'98:1_Vdot(lh)': 'HM1_Vpkt', '97:1_T_RL(C)' : 'HM1_T_RET', '96:1_T_VL(C)' : 'HM1_T_SUP'})
df2 = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Modelle/02_WP/Remus/ExWirkungsgrad.csv', index_col = 0)
df2.index = pd.to_datetime(df2.index)

df = pd.merge(df1, df2, left_index=True, right_index=True, how='outer')
df.to_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Modelle/02_WP/ExWirkungsgrad.csv')

#%%% VolumeFlow
columns = {}
# columns['HM1'] = ['42:TempWP1Quelleaus', '43:TempWP1Quelleein', '44:TempWP1Vorlauf', '45:TempWP1Ruecklauf','46:TempWP1PTherm', '47:TempWP1PElek', '48:TempWP1COP','98:1_Vdot(lh)', '97:1_T_RL(C)', '96:1_T_VL(C)']
columns['HM1'] = ['98:1_Vdot(lh)', '97:1_T_RL(C)', '96:1_T_VL(C)']
data = {'HM1': pd.DataFrame()}

j = 1
for file in os.listdir('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Modelle/02_WP/Remus')[:-1]:
    print('adding: ' + file + ' ' + str(j) + '/' + str(len([name for name in os.listdir('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Modelle/02_WP/Remus')])))
    df = pd.read_table('Remus/' + file, sep = ' ')
    index = [datetime.strptime(df.iloc[i,0] + ' ' + df.iloc[i,1], '%Y.%m.%d %H:%M:%S') for i in df.index]
    for i in list(data):
        data[i] = data[i].append(pd.DataFrame(index = index, data = df[columns[i]].to_numpy(), columns=columns[i]))
    j = j+1
del i, file, index, df, j, columns

df = data['HM1']
df = df.rename(columns={'98:1_Vdot(lh)': 'HM1_Vpkt', '97:1_T_RL(C)' : 'HM1_T_RET', '96:1_T_VL(C)' : 'HM1_T_SUP'})
df.to_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Modelle/02_WP/VolumeFlow_HM1.csv')
del data

#%%
df = df[df['HM1_Vpkt']>300]
#%%
for i in [0, 100, 200, 300, 400, 500, 600]:
    
    df1 = df[df['HM1_Vpkt']>i]
    fig, ax1 = plt.subplots(figsize=(6.4, 4.8/3))
    
    # ------------------------------------------------Y1------------------------------------------------
    # Data
    ax1.boxplot(df1['HM1_Vpkt'], vert = False)
    # labels
    ax1.set_yticks([])
    # ax1.set_ylim(0.9,1.1)
    # ------------------------------------------------X------------------------------------------------
    # ax1.set_xlim(4,6)
    # ax1.set_xlabel('ΔT in °C')
    # ax1.set_xticks(np.arange(4, 6.01, 0.25))
    # ax1.set_xticks(np.arange(4, 6.01, 0.05), minor = True)
    # ------------------------------------------------Global------------------------------------------------
    ax1.set_title(r'$\bf{Verteilung\ des\ Volumesntromes\ auf\ der \ kalten \ Seite}$' + '\nZeitraum: 01.01.2022 - 31.01.2022')
    ax1.grid(True)
    ax1.grid(which='minor', alpha = 0.3)
    ax1.set_xlim(0,800)
    ax1.set_title(str(i))
    # ax1.legend(fontsize = 8)
    plt.show()






