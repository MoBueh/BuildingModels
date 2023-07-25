# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 11:03:57 2023

@author: MoBueh
"""

#%% Directory
import os
os.chdir('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin')

#%% Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
plt.rcParams["figure.dpi"] = 1000
plt.rc('axes', axisbelow=True)

#%% Same Columns?
# columns = {}

# for file in os.listdir(os.getcwd()+'/Data/Monitoring/Raw_5')[:-1]:
#     # print(file)
#     columns[file.split('.')[0]] = list(pd.read_csv('Data/Monitoring/Raw_5/'+file, index_col=0, delimiter = ';',skiprows = 1, nrows = 0).columns)

# amount_per_column = []
# for i in range(0,36):
#     elements = []
#     for key in columns:
#         elements.append(columns[key][i])  
#     amount_per_column.append(len(list(set(elements))))
#     # print(list(set(elements)))
# print(amount_per_column)
    
# del columns, elements, file, i, key
# # del amount_per_column


     
    
#%% Merge Raw Data
# -> 104 days/df

df = pd.DataFrame()
# for file in os.listdir(os.getcwd()+'/Data/Monitoring/Raw_5')[1:]:
for file in os.listdir(os.getcwd()+'/Data/Monitoring/Raw_15')[1:]:
    print(file)
    # append = pd.read_csv('Data/Monitoring/Raw_5/'+file, index_col=0, sep = ';', skiprows = 1)
    append = pd.read_csv('Data/Monitoring/Raw_15/'+file, index_col=0, sep = ';', skiprows = 1)
    append.index = pd.to_datetime(append.index, format='%Y/%m/%d %H:%M:%S')
    # append.index = pd.to_datetime(append.index, format='%d.%m.%Y %H:%M')
    df = pd.concat([df, append], axis = 0)
del append, file
df.sort_index(inplace=True)
df.replace(-99, np.nan, inplace=True)
df = df.drop('GEN.EL-PV----CALC.POW.EL-', axis=1)
# df = df.resample('15min').asfreq()
# df.to_csv('Data/Monitoring/Processed_5/total.csv')
df.to_csv('Data/Monitoring/Processed_15/total.csv')

#%% DataFrame per Day
# for year in df.index.year.unique():
#     df_y = df[df.index.year == year] 
#     for month in range(1,13):
#         df_m = df_y[df_y.index.month == month]    
#         for day in df_m.index.day.unique():
#             df_store = df_m[df_m.index.day == day]
#             # df_store.to_csv('Data/Monitoring/Processed_5/'+df_store.index[0].strftime('%y-%m-%d')+'.csv')
#             df_store.to_csv('Data/Monitoring/Processed_15/'+df_store.index[0].strftime('%y-%m-%d')+'.csv')
#             print(df_store.index[0].strftime('%y-%m-%d'))
# del day, month, year, df_y, df_m, df_store

#%% Change Columns
from Data.Monitoring.Rename_Dictionary import rename
df = df.rename(columns=rename)
del rename

#%%
summer = df[df.index >= datetime.strptime('2022-08-08 00:00:00', '%Y-%m-%d %H:%M:%S')]
summer = summer[summer.index < datetime.strptime('2022-08-15 00:00:00', '%Y-%m-%d %H:%M:%S')]
summer.loc[summer.index.hour < 6, 'Pyr_hor'] = 0
summer.loc[summer.index.hour > 20, 'Pyr_hor'] = 0

plt.plot(summer.index, summer['AMB_T'], c = 'TAB:GREEN')
plt.plot(summer.index, summer['Zone_T'], c = 'black')
plt.plot(summer.index, [y/50 for y in summer['HM4_Vpkt']], c = 'TAB:BLUE')
plt.xticks(rotation = 90)
plt.title('2022: Week ' + str(summer.index[0].isocalendar().week))
plt.ylim(-10, 40)
plt.grid()
plt.show()
summer.to_csv('Data/Monitoring/TestSet/TestSet_summer.csv')



#%%
winter = df[df.index >= datetime.strptime('2021-12-27 00:00:00', '%Y-%m-%d %H:%M:%S')]
winter = winter[winter.index < datetime.strptime('2022-01-03 00:00:00', '%Y-%m-%d %H:%M:%S')]
winter.loc[winter.index.hour < 8, 'Pyr_hor'] = 0
winter.loc[winter.index.hour > 16, 'Pyr_hor'] = 0

plt.plot(winter.index, winter['AMB_T'], c = 'TAB:GREEN')
plt.plot(winter.index, winter['Zone_T'], c = 'black')
plt.plot(winter.index, [y/50 for y in winter['HM4_Vpkt']], c = 'TAB:BLUE')
# plt.plot(y23_week.index, y23_week['Pyr_hor'], c = 'y')
plt.xticks(rotation = 90)
plt.title('2021: Week ' + str(winter.index[0].isocalendar().week))
plt.ylim(-10, 30)
plt.grid()
plt.show()
winter.to_csv('Data/Monitoring/TestSet/TestSet_winter.csv')
    
    
    
    
    
    
    
    

