# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 10:39:42 2023

@author: MoBueh
"""

#%% Directory
import os
os.chdir('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin')

#%% Libraries
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
plt.rcParams["figure.dpi"] = 1000
plt.rc('axes', axisbelow=True)

#%%
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_Winter.csv')
df.loc[df['HM10_Vpkt'] < 30, 'HM10_Vpkt'] = 0


#%%
plt.plot(df.index, df['HM4_Vpkt'], label = 'HM4_Vpkt')
plt.plot(df.index, df['HM10_Vpkt'], label = 'HM10_Vpkt')
plt.plot(df.index, df['WP_COP'], label = 'WP_COP')
plt.xlim(450, 479)
plt.grid()
plt.legend()
plt.show()

#%%
plt.plot(df.index, df['HM4_T_SUP'], label = 'HM4_T_SUP')
plt.plot(df.index, df['HM4_T_RET'], label = 'HM4_T_RET')
plt.plot(df.index, df['HM10_T_SUP'], label = 'HM10_T_SUP')
plt.plot(df.index, df['WS_T'], label = 'WS_T')
plt.xlim(450, 479)
plt.grid()
plt.legend()
plt.show()


