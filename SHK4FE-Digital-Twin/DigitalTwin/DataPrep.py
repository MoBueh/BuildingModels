"""
@author: MoBueh
"""

#%% Libraries
import os
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np
import copy
import math
from scipy.integrate import odeint
import matplotlib.image as mpimg
plt.rcParams["figure.dpi"] = 100
plt.rc('axes', axisbelow=True)


#%% Global
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Wetter_HSO/global_irr.csv', index_col = 0, sep = ';')
df.index = pd.to_datetime(df.index, format='%Y-%m-%d %H:%M:%S')  
df = df.loc[df.index.year == 2021, :] 
df = df.rename(columns={'global_irr [kWh/m2]' : 'G'})
df1 = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Wetter_HSO/air_temperature.csv', index_col = 0, sep = ';')
df1.index = pd.to_datetime(df1.index, format='%Y-%m-%d %H:%M:%S')  
df1 = df1.loc[df1.index.year == 2021, :] 
df1 = df1.rename(columns={'air_temperature' : 't'})
df['t'] = df1['t']

#%%
df.to_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/DigitalTwin/21_wetter_HSO.csv')