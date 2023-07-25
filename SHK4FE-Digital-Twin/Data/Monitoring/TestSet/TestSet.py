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

#%% Libraries
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
plt.rcParams["figure.dpi"] = 400
plt.rc('axes', axisbelow=True)

def SunPos(lon, lat, year, month, day, hour, minute, shift):
    # only wintertime (UTC +1) is valid; sun hight and azimith in degree, lon = 7.9498017 and lat = 48.473451 for Offenburg,
    # year, month day, hour as integer, shift for different UTC
    time = datetime(year, month, day, hour, minute)
    time = time + timedelta(hours = shift)
    cal_day = time.timetuple().tm_yday
    oz = time.hour+time.minute/60+time.second/60/60
    moz = oz - 4*(15-lon)/60
    zgl = 0.0066+7.325*math.cos(math.radians(360/365*cal_day+85.9))+9.9359*math.cos(math.radians(2*360/365*cal_day+108.9))+0.3387*math.cos(math.radians(3*360/365*cal_day+105.2))
    woz = moz + zgl/60
    sun_dek = 0.3948-23.2559*math.cos(math.radians(360/365*cal_day+9.1))-0.3915*math.cos(math.radians(2*360/365*cal_day+5.4))-0.1764*math.cos(math.radians(3*360/365*cal_day+26.0))    
    sun_h = math.degrees(math.asin(math.cos(math.radians((12-woz)*15))*math.cos(math.radians(lat))*math.cos(math.radians(sun_dek))+math.sin(math.radians(lat))*math.sin(math.radians(sun_dek))))
    if sun_h<0:
        sun_h = 0
    if woz <= 12:
        sun_az = 180 - math.degrees(math.acos((math.sin(math.radians(sun_h))*math.sin((math.radians(lat)))-math.sin(math.radians(sun_dek)))/(math.cos(math.radians(sun_h))*math.cos(math.radians(lat)))))
    else:
        sun_az = 180 + math.degrees(math.acos((math.sin(math.radians(sun_h))*math.sin((math.radians(lat)))-math.sin(math.radians(sun_dek)))/(math.cos(math.radians(sun_h))*math.cos(math.radians(lat)))))
    return sun_az, sun_h

param = {}
param['Position'] = {
    'lon' : 7.9498017,
    'lat' : 48.473451,
    }

#%% Winter
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_winter.csv', index_col = 0)
df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')

i_dif = []
i_dir = []

for i in df.index:
    i_dif.append((1/(1 + 2.71828182846**(7.997*(df.loc[i, 'Pyr_hor']/1000 - 0.586))))*df.loc[i, 'Pyr_hor'])
    i_dir.append(df.loc[i, 'Pyr_hor']-i_dif[-1])
    
df['i_dif'] = i_dif
df['i_dir'] = i_dir

#%%% DataPrep
Tmin = []
Tmax = []
Tmean = []
Tschwankung = []

G_tot = []
G_mean = []

for i in df.index.day.unique():
    day = df[df.index.day == i]
    Tmin.append(day['AMB_T'].min())
    Tmax.append(day['AMB_T'].max())
    Tmean.append(day['AMB_T'].mean())
    G_tot.append(day['Pyr_hor'].sum()*900/3600) #Wh
    G_mean.append(day[day['Pyr_hor'] != 0]['Pyr_hor'].mean())
    Tschwankung.append((Tmax[-1] - Tmin[-1])/2)

#%%% Temperatur
# ----------------------------------------------------------------------------------------------------------------
# ------------------------------------------------Zeitreihe-------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 1]},figsize=(6.4, 4.8/1.25))

# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.plot(df.index, df['AMB_T'], color = 'tab:green', label = 'Umgebungstemperatur')
# labels
ax1.set_ylabel('Temperatur in °C', fontsize = 9)
ax1.set_ylim(0,20)

# ------------------------------------------------X------------------------------------------------
ax1.set_xlim(df.index[0],df.index[-1])
ax1.set_xticks([df.index[0] + timedelta(days=i) for i in range(7)], [' '] * 7, rotation = 60)

# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Zeitreihe}$', pad=-10).set_size(10) 
ax1.grid(True)

# --------------------------------------------------------------------------------------------------------------
# ------------------------------------------------Analyse-------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------

# ------------------------------------------------Y1------------------------------------------------
# Data
ax2.scatter([i+0.5 for i in range(0,7,1)], Tmin, color = 'tab:blue', label = 'Minimalwert', s = 12)
ax2.plot([i+0.5 for i in range(0,7,1)], Tmin, color = 'tab:blue', alpha = 0.4, linestyle = '--')
ax2.scatter([i+0.5 for i in range(0,7,1)], Tmax, color = 'tab:red', label = 'Maximalwert', s = 12)
ax2.plot([i+0.5 for i in range(0,7,1)], Tmax, color = 'tab:red', alpha = 0.4, linestyle = '--')
ax2.scatter([i+0.5 for i in range(0,7,1)], Tmean, color = 'tab:green', label = 'Mittelwert', s = 12)
ax2.plot([i+0.5 for i in range(0,7,1)], Tmean, color = 'tab:green', alpha = 0.4, linestyle = '--')
# labels
ax2.set_ylabel('Temperatur in °C', fontsize = 9)
ax2.yaxis.set_major_formatter('{:.1f}'.format)
ax2.set_ylim(0,20)
ax2_twin = ax2.twinx()
# ------------------------------------------------Y2------------------------------------------------
# Data
ax2_twin.bar([i+0.5 for i in range(0,7,1)], Tschwankung, color = 'tab:green', label = 'Schwankung', alpha = 0.2)
# labels
ax2_twin.set_ylabel('Schwankung in K', fontsize = 9)
ax2_twin.yaxis.set_major_formatter('{:.1f}'.format)
ax2_twin.set_ylim(0,10)
ax2_twin.set_yticks([0, 2.5, 5, 7.5, 10])
# ------------------------------------------------X------------------------------------------------
ax2.set_xlim(0,7)
ax2.set_xticks(range(0,7,1), ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], rotation = 60, fontsize = 8)
# ------------------------------------------------Global------------------------------------------------
ax2.set_title(r'$\bf{Tageswerte}$', pad=-10).set_size(10) 
# ax.set_title(r'$\bf{First\ Row}$' + '\nSecond Row')
ax2.grid(True)
ax2.legend(fontsize = 6, loc = 'upper left')
ax2_twin.legend(fontsize = 6, loc = 'upper right')
plt.show()



#%%% Plot Strahlung
# ----------------------------------------------------------------------------------------------------------------
# ------------------------------------------------Zeitreihe-------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(6.4, 4.8/2.5))
# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.fill_between(df.index, df['Pyr_hor'], df['i_dir'], color='tab:orange')
ax1.fill_between(df.index, df['i_dir'], [0]*672, color='tab:red')
ax1.plot(df.index, df['Pyr_hor'], color = 'black', linewidth = 1, label = 'Global')
ax1.plot(df.index, df['i_dir'], color = 'black', linewidth = 0.2)
plt.bar(df.index[0], -200, color = 'tab:red', label = 'Direkt')
plt.bar(df.index[0], -200, color = 'tab:orange', label = 'Diffus')
# labels
ax1.set_ylabel('Globalstrahlung in W/m²', fontsize = 7)
ax1.set_ylim(0,400)
# ------------------------------------------------X------------------------------------------------
ax1.set_xlim(df.index[0],df.index[-1])
ax1.set_xticks([df.index[0] + timedelta(days=i) for i in range(7)],['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], rotation = 60, fontsize = 7)

# ------------------------------------------------Global------------------------------------------------
ax1.grid(True)
ax1.legend(fontsize = 7, ncol = 3, loc = 'upper right')
plt.show()


#%% Sommer
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_summer.csv', index_col = 0)
df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')

i_dif = []
i_dir = []

for i in df.index:
    i_dif.append((1/(1 + 2.71828182846**(7.997*(df.loc[i, 'Pyr_hor']/1000 - 0.586))))*df.loc[i, 'Pyr_hor'])
    i_dir.append(df.loc[i, 'Pyr_hor']-i_dif[-1])
df['i_dif'] = i_dif
df['i_dir'] = i_dir

#%%% DataPrep
Tmin = []
Tmax = []
Tmean = []
Tschwankung = []

G_tot = []
G_mean = []

for i in df.index.day.unique():
    day = df[df.index.day == i]
    Tmin.append(day['AMB_T'].min())
    Tmax.append(day['AMB_T'].max())
    Tmean.append(day['AMB_T'].mean())
    G_tot.append(day['Pyr_hor'].sum()*900/3600) #Wh
    G_mean.append(day[day['Pyr_hor'] != 0]['Pyr_hor'].mean())
    Tschwankung.append((Tmax[-1] - Tmin[-1])/2)

#%%% Temperatur
# ----------------------------------------------------------------------------------------------------------------
# ------------------------------------------------Zeitreihe-------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 1]},figsize=(6.4, 4.8/1.25))

# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.plot(df.index, df['AMB_T'], color = 'tab:green', label = 'Umgebungstemperatur')
# labels
ax1.set_ylabel('Temperatur in °C', fontsize = 9)
ax1.set_ylim(10,40)

# ------------------------------------------------X------------------------------------------------
ax1.set_xlim(df.index[0],df.index[-1])
ax1.set_xticks([df.index[0] + timedelta(days=i) for i in range(7)], [' '] * 7, rotation = 60)

# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Zeitreihe}$', pad=-10).set_size(10) 
ax1.grid(True)

# --------------------------------------------------------------------------------------------------------------
# ------------------------------------------------Analyse-------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------

# ------------------------------------------------Y1------------------------------------------------
# Data
ax2.scatter([i+0.5 for i in range(0,7,1)], Tmin, color = 'tab:blue', label = 'Minimalwert', s = 12)
ax2.plot([i+0.5 for i in range(0,7,1)], Tmin, color = 'tab:blue', alpha = 0.4, linestyle = '--')
ax2.scatter([i+0.5 for i in range(0,7,1)], Tmax, color = 'tab:red', label = 'Maximalwert', s = 12)
ax2.plot([i+0.5 for i in range(0,7,1)], Tmax, color = 'tab:red', alpha = 0.4, linestyle = '--')
ax2.scatter([i+0.5 for i in range(0,7,1)], Tmean, color = 'tab:green', label = 'Mittelwert', s = 12)
ax2.plot([i+0.5 for i in range(0,7,1)], Tmean, color = 'tab:green', alpha = 0.4, linestyle = '--')
# labels
ax2.set_ylabel('Temperatur in °C', fontsize = 9)
ax2.yaxis.set_major_formatter('{:.1f}'.format)
ax2.set_ylim(10,40)

ax2_twin = ax2.twinx()
# ------------------------------------------------Y2------------------------------------------------
# Data
ax2_twin.bar([i+0.5 for i in range(0,7,1)], Tschwankung, color = 'tab:green', label = 'Schwankung', alpha = 0.2)
# labels
ax2_twin.set_ylabel('Schwankung in K', fontsize = 9)
ax2_twin.yaxis.set_major_formatter('{:.1f}'.format)
ax2_twin.set_ylim(0,12)
# ax2_twin.set_yticks([0, 2.5, 5, 7.5, 10])

# ------------------------------------------------X------------------------------------------------
ax2.set_xlim(0,7)
ax2.set_xticks(range(0,7,1), ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], rotation = 60, fontsize = 8)

# ------------------------------------------------Global------------------------------------------------
ax2.set_title(r'$\bf{Tageswerte}$', pad=-10).set_size(10) 
# ax.set_title(r'$\bf{First\ Row}$' + '\nSecond Row')
ax2.grid(True)
ax2.legend(fontsize = 6, loc = 'upper left')
ax2_twin.legend(fontsize = 6, loc = 'upper right')
plt.show()

#%%% Plot Strahlung
fig, ax1 = plt.subplots(figsize=(6.4, 4.8/2.5))
# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.fill_between(df.index, df['Pyr_hor'], df['i_dir'], color='tab:orange')
ax1.fill_between(df.index, df['i_dir'], [0]*672, color='tab:red')
ax1.plot(df.index, df['Pyr_hor'], color = 'black', linewidth = 1, label = 'Global')
ax1.plot(df.index, df['i_dir'], color = 'black', linewidth = 0.2)
plt.bar(df.index[0], -200, color = 'tab:red', label = 'Direkt')
plt.bar(df.index[0], -200, color = 'tab:orange', label = 'Diffus')
# labels
ax1.set_ylabel('Globalstrahlung in W/m²', fontsize = 7)
ax1.set_ylim(0,1050)
# ------------------------------------------------X------------------------------------------------
ax1.set_xlim(df.index[0],df.index[-1])
ax1.set_xticks([df.index[0] + timedelta(days=i) for i in range(7)],['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], rotation = 60, fontsize = 7)

# ------------------------------------------------Global------------------------------------------------
ax1.grid(True)
ax1.legend(fontsize = 7, ncol = 3, loc = 'upper right')
plt.show()