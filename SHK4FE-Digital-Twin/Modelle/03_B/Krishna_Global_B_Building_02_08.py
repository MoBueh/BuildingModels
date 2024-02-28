"""
@author: MoBueh & Krishna
"""

#%% Libraries
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np
import math
plt.rcParams["figure.dpi"] = 400
plt.rc('axes', axisbelow=True)
def prRed(skk): print("\033[31;1;m {}\033[00m" .format(skk))
def prYellow(skk): print("\033[33;1;m {}\033[00m" .format(skk))

from influxdb_client import InfluxDBClient, Point, WriteOptions, WriteApi
from influxdb_client.client.write_api import SYNCHRONOUS

import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import pickle
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
param['Position'] = {
    'lon' : 7.942275907924158,
    'lat' : 48.457691216034874,
    }
l = 3.60*2.0; b = 3.60 ;
h = 3.00

param['B'] = {
    'A_f': l*b,                                                                 # m²        conditioned area
    'Coe_Am': 3.5,                                                              # -         coefficient for the determination of mass-related area according to table 12
    'Coe_Cm': 370000,                                                           # J/K*m²    coefficiant for the determination of Internal heat storage capacity according to table 12
    'transparent components': {
        # each column is one Component 
        'shading'    : [0,      0,      0,         1,      0,       0],           # -         1 for no shading, 0 for absolutely shadowed
        'g_tot'      : [0.4,    0.4,    0.4,     0.4,    0.4,       0],           # -         DIN 4108-2 8.3.2 (3) und DIN EN ISO 13790 Tabelle 7
        'F_c'        : [1,      1,      1,         1,      1,       1],           # -         DIN 4108-2 8.3.2 (3) und DIN EN ISO 13790 Tabelle 7
        'F_f'        : [0,      0,      0,      0.05,      0,       0],           # -         Frame part of the component; 1 -> 100 % frame, 0 -> 0 % frame
        'A'          : [0,      0,      0,      7.35,      0,       0],           # m²        Area of the component
        'direction'  : [41,     131,    221,     311,    270,     270],           # °         0° -> north, 90° -> east, 180° ->south, 270° -> west
        'tilt'       : [90,     90,     90,       90,      0,     180],           # °         0° horizontal top, 90° -> vertical, 180° -> horizontal down
        'U'          : [5.0]*6,                                                   # W/m²/K    Heat transfer coefficient determined according to ISO 6946
        },
    'opaque components': {
        # each column is one Component
        'shading'    : [0,      0,      0,      0,      0,    0],               # -         1 for no shading, 0 for absolutely shadowed
        'R'          : [0.05]*6,                                                #           Product of surface heat transfer resistance(ISO 6946) and absorption coefficient for solar radiation on opaque component, "Baehr, H.D.; Stephan, K. ; Wärme- und Stoffübertragung; Auszug Kap. 5.5 Strahlungsaustausch, aus Tab. 5.8,S.633" https://www.baunetzwissen.de/bauphysik/fachwissen/waermeschutz/absorption-und-waerme-auf-oberflaechen-4733315
        'U'          : [0.2]*6,                                                 # W/m²/K    Heat transfer coefficient determined according to ISO 6946
        'A'          : [b*h,    l*h,    b*h,    l*h,    l*b,    l*b],           # m²        Area of the component
        'direction'  : [41,     131,    221,    311,    270,    270],           # °         0° -> north, 90° -> east, 180° ->south, 270° -> west
        'tilt'       : [90,     90,     90,     90,     0,      180]            # °         0° horizontal top, 90° -> vertical, 180° -> horizontal down
        },
    }

#%% Model
def B(B_phi_int, B_phi_HC, B_theta_m_prev, theta_e, sun_h, sun_az, I_dir, I_diff, V_Vdot, V_theta_sup):
    #                        _______________
    # Int Qpkt        [W]   |               |
    # --------------------->|               |
    # HC_Qpkt         [W]   |               |
    # --------------------->|               |
    # T_mass[-1]      [°C]  |               | T_air         [°C]
    # --------------------->|               |--------------------->
    # T_AMB           [°C]  |               | T_surface     [°C]
    # --------------------->|               |--------------------->
    # Sonnenstand     [°]   |     5R1C      | T_mass        [°C]
    # --------------------->|               |--------------------->
    # Solarstrahlung  [W/m2]|               | T_operativ    [°C]
    # --------------------->|               |--------------------->
    # Vpkt_Lüftung  [m3/hr] |               |
    # --------------------->|               |
    # T_Lüftung       [°C]  |               |
    # --------------------->|               |
    #                       |_______________|
    if str(B_phi_HC) == 'nan':
        B_phi_HC = 0
    else:
        B_phi_HC = B_phi_HC
    A_m = param['B']['A_f']*param['B']['Coe_Am'] # 12.3.1 Moritz 5R1C; effective mass area, expressed in m2.
    A_tot = 4.5 * param['B']['A_f'] # 7.2.2.2 Moritz 5R1C; area of all surfaces facing the building zone in m2
    c_m = param['B']['Coe_Cm']*param['B']['A_f'] # internal heat capacity of the building or building zone, calculated in accordance with 12.3.1
    H_tr_w = sum([U*A for U,A in zip(param['B']['transparent components']['U'], param['B']['transparent components']['A'])]) #overall transmission heat transfer coefficient W/K
    if V_Vdot == 0:
        H_ve = 0.00000001 # overall ventilation heat transfer coefficient W/K , else Eq C.6 becomes infinity
    else:
        H_ve = 1200*V_Vdot/3600                                                 # keine Wärmerückgewinnung
        # H_ve = (1-0.85)*1200*V_Vdot/3600                                      # mit Wärmerückgewinnung
    H_tr_ms = 9.1*A_m #coupling conductance between nodes m and s, expressed in watts per kelvin 12.2.2
    H_op = sum([U*A for U,A in zip(param['B']['opaque components']['U'], param['B']['opaque components']['A'])]) #overall ventilation heat transfer coefficient W/K
    H_tr_em = 1/(1/H_op+1/H_tr_ms) # heat transfer coefficient W/K between e and m
    H_tr_is = 3.45*A_tot # heat transfer coefficient W/K between air and surface
# Phi_sol: transparent components
    phi_sol_trans = []
    for shading, g_tot, F_c, F_f, A, direction, tilt in zip(param['B']['transparent components']['shading'],
                                                                       param['B']['transparent components']['g_tot'],
                                                                       param['B']['transparent components']['F_c'],
                                                                       param['B']['transparent components']['F_f'],
                                                                       param['B']['transparent components']['A'],
                                                                       param['B']['transparent components']['direction'],
                                                                       param['B']['transparent components']['tilt']):
        if sun_h < 5:
            I_f = 0
        else:
            I_f = (math.sin(math.radians(sun_h))*math.cos(math.radians(tilt)) + math.cos(math.radians(sun_h))*math.sin(math.radians(tilt)) * math.cos(math.radians(abs(direction-sun_az)))) * I_dir / math.sin(math.radians(sun_h))
            # I_f = (math.sin(math.radians(sun_h))*math.cos(math.radians(direction)) + math.cos(math.radians(sun_h))*math.sin(math.radians(direction)) * math.cos(math.radians(abs(tilt-sun_az)))) * I_dir / math.sin(math.radians(sun_h))
        if I_f <0:
            I_f = 0
        if I_f > 1200:
            I_f = 1200
        I_umg = (I_diff + I_dir) * 0.5 * 0.2 * (1-math.cos(math.radians(tilt))) #Ground reflexion onto the surface area with any orientation and inclination VDI6007
        I_d = I_diff * 0.5 * (1 + math.cos(math.radians(tilt)))
        I_sol = I_f + I_umg + I_d
        phi_sol_trans.append(shading*g_tot*F_c*0.9*(1-F_f)*A*I_sol)
# Phi_sol: opaque components
    phi_sol_op = []
    for shading, R, U, A, direction, tilt in zip(param['B']['opaque components']['shading'], # alpha and R = Rse
                                                                   param['B']['opaque components']['R'],
                                                                   param['B']['opaque components']['U'],
                                                                   param['B']['opaque components']['A'],
                                                                   param['B']['opaque components']['direction'],
                                                                   param['B']['opaque components']['tilt']):
        # print(R)
        if sun_h < 5:
            I_f = 0
        else:
            I_f = (math.sin(math.radians(sun_h))*math.cos(math.radians(tilt)) + math.cos(math.radians(sun_h))*math.sin(math.radians(tilt)) * math.cos(math.radians(abs(direction-sun_az)))) * I_dir / math.sin(math.radians(sun_h))
        if I_f <0:
            I_f = 0
        if I_f > 1200:
            I_f = 1200
        # print('I_f:  ' + str(I_f))
        I_umg = (I_diff + I_dir) * 0.5 * 0.2 * (1-math.cos(math.radians(tilt)))
        # print('I_umg: ' + str(I_umg))
        I_d = I_diff * 0.5 * (1 + math.cos(math.radians(tilt)))
        # print('I_d:   ' + str(I_d))
        # print('')
        I_sol = I_f +I_umg+I_d
        phi_sol_op.append(shading*R*U*A*I_sol)
# Phi_sol: sum of opaque and transparent
    phi_sol = sum(phi_sol_trans) + sum(phi_sol_op)
# Equation for Model regarding Appendix C
    H_tr_1 = 1 / (1/H_ve + 1/H_tr_is)                                           # C.6
    H_tr_2 = H_tr_1 + H_tr_w                                                    # C.7
    H_tr_3 = 1 / (1/H_tr_2 + 1/ H_tr_ms)                                        # C.8
    phi_ia = 0.5*B_phi_int+0*B_phi_HC                                   + 0*B_phi_HC                                                    # C.1
    phi_m = A_m/A_tot * (0.5 * B_phi_int + phi_sol)                     + 1*B_phi_HC                          # C.2
    phi_st = (1-A_m/A_tot-H_tr_w/9.1/A_tot)*(0.5 * B_phi_int + phi_sol) + 0*B_phi_HC          # C.3
    phi_mtot = phi_m + H_tr_em * theta_e + H_tr_3 *(phi_st +H_tr_w * theta_e + H_tr_1 * ((phi_ia)/H_ve + V_theta_sup))/ H_tr_2   # C.5
    B_theta_m_t = (B_theta_m_prev * (c_m/900-0.5 * (H_tr_3 +H_tr_em))+phi_mtot)/(c_m/900+0.5*(H_tr_3+H_tr_em))  # C.4
        # Output
    B_theta_m = (B_theta_m_t + B_theta_m_prev)/2                                      # C.9
    # B_theta_m = B_theta_m_t                                     # C.9
    B_theta_s = (H_tr_ms * B_theta_m + phi_st + H_tr_w * theta_e + H_tr_1 *(V_theta_sup + (phi_ia)/H_ve))/(H_tr_ms + H_tr_w + H_tr_1)     # C.10
    B_theta_air = (H_tr_is * B_theta_s + H_ve * V_theta_sup + phi_ia)/(H_tr_is + H_ve)          # C.11
    B_theta_operativ = 0.3 * B_theta_air + 0.7 * B_theta_s                            # C.12
    return B_theta_m, B_theta_s, B_theta_air, B_theta_operativ, phi_sol_op, phi_sol_trans, phi_sol

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

#%% Single Step
b = B(B_phi_int = 0, #  heat flow rate from internal heat source
      B_phi_HC = 0,     # heat flow rate from heating, cooling and ventilation systems,
      B_theta_m_prev = 20, #mass temperature previous step
      theta_e = 30, # Outdoor temperature
      sun_h = 10, # sun position
      sun_az = 66, # Azimuth
      I_dir = 100, # direct radiation
      I_diff = 0,  # diffuse radation
      V_Vdot = 30, # m3/hr
      V_theta_sup = 30) #temperature of supply air
#%% Heating

#%%%% Data Acquisition

token = "clqtOMwwAIMuVlMnBySs4Am5py_afCRB2K0yXGuu-ORykw8kANEGEeoyd-t7eUcbUSd2skYnAmj6a1-duEmJLw=="
my_org = "HS-Offenburg"
buckety = "heatGUIde"
measurement = "70B3D57ED005E201"
url="https://heatguide-ines.hs-offenburg.de:8086/"

#%%%%% Weather
query = '''from(bucket: "Mosmix")
  |> range(start: 2023-07-31T00:00:00Z, stop: 2023-08-06T00:00:00Z)
  |> filter(fn: (r) => r["_measurement"] == "K2301")
  |> filter(fn: (r) => r["_field"] == "temperature_air_mean_005")
  |> filter(fn: (r) => r["hour"] == "0")
  |> pivot(rowKey:["_time"], columnKey: ["_field"], valueColumn: "_value")'''


client = InfluxDBClient(url=url, token=token, org=my_org, debug=False, verify_ssl=False ) # connects to the Influx DB

weather_data = client.query_api().query_data_frame(org=my_org, query=query) # gets the data in the form of a dataframe

prYellow(weather_data.head()) # displays only the first 10 rows

weather_data = weather_data.loc[:,["_time", "temperature_air_mean_005"]]
weather_data["temperature_air_mean_005"] = weather_data["temperature_air_mean_005"] - 273.15
weather_data['_time'] = weather_data['_time'].dt.tz_convert('Europe/Berlin').dt.tz_localize(None)
weather_data.set_index('_time', inplace=True)
df_15min_avg = weather_data.resample('15min').mean().fillna(method="pad")
df_15min_avg.columns = ["AMB_T"] # Outdoor Temperature from DWD

#%%%%% Radiation
query = '''from(bucket: "Mosmix")
  |> range(start: 2023-07-31T00:00:00Z, stop: 2023-08-06T00:00:00Z)
  |> filter(fn: (r) => r["_measurement"] == "K2301")
  |> filter(fn: (r) => r["_field"] == "radiation_global")
  |> filter(fn: (r) => r["hour"] == "0")
  |> pivot(rowKey:["_time"], columnKey: ["_field"], valueColumn: "_value")'''


client = InfluxDBClient(url=url, token=token, org=my_org, debug=False, verify_ssl=False ) # connects to the Influx DB

radiation_data = client.query_api().query_data_frame(org=my_org, query=query) # gets the data in the form of a dataframe

prYellow(radiation_data.head()) # displays only the first 10 rows

radiation_data = radiation_data.loc[:,["_time", "radiation_global"]]
radiation_data['_time'] = radiation_data['_time'].dt.tz_convert('Europe/Berlin').dt.tz_localize(None)
radiation_data.set_index('_time', inplace=True)
radiation_data = (radiation_data*1000)/3600 # Convert Mosmix data from KJ/m2 to W/m2
df_rad_15min_avg = radiation_data.resample('15min').mean().fillna(method="pad")
df_rad_15min_avg.columns = ["Pyr_hor"] # Global Radiation Temperature from DWD
#%%%%% Indoor Temperature

query = '''from(bucket: "Campus-OG")
  |> range(start: 2023-07-31T00:00:00Z, stop: 2023-08-06T00:00:00Z)
  |> filter(fn: (r) => r["_measurement"] == "70B3D57ED005E201")
  |> filter(fn: (r) => r["_field"] == "temp_htu21")
  |> pivot(rowKey:["_time"], columnKey: ["_field"], valueColumn: "_value")'''



client = InfluxDBClient(url=url, token=token, org=my_org, debug=False, verify_ssl=False ) # connects to the Influx DB

node_data = client.query_api().query_data_frame(org=my_org, query=query) # gets the data in the form of a dataframe

prYellow(node_data.head()) # displays only the first 10 rows
node_data = node_data.loc[:,["_time", "temp_htu21"]]
node_data['_time'] = node_data['_time'].dt.tz_convert('Europe/Berlin').dt.tz_localize(None)
node_data.set_index('_time', inplace=True)
df_hourly_avg = node_data.resample('15min').mean().fillna(method="pad")
df_hourly_avg.columns = ["Zone_T"] # Indoor Temperature from Node Measurement
df_hourly_avg = df_hourly_avg.truncate(before=str(df_15min_avg.index[0]), after=str(df_15min_avg.index[-1]), axis=0, copy=True)


#%%% Actual Simulation


# =============================================================================
# df = pd.read_csv('C:/Users/rdevinen/Music/GitHub/BuildingModels/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_winter.csv', index_col = 0)
# df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')
# df['HM4_Qpkt'] = [vpkt/1000/3600 * prop['wasser']['rho'] *prop['wasser']['cp'] * (SUP-RET) for vpkt, SUP, RET in zip(df['HM4_Vpkt'], df['HM4_T_SUP'],df['HM4_T_RET'])]
# df['HM4_Qpkt_delay'] = df['HM4_Qpkt'].shift(2)
# =============================================================================

df = pd.concat([df_15min_avg, df_hourly_avg, df_rad_15min_avg], axis = 1)
df['HM4_Qpkt'] = 0



B_theta_m = [df.loc[df.index[0], 'Zone_T']]
B_theta_s = [20]
B_theta_air = [20]
B_theta_operativ = [20]
sun_h = []
sun_az = []
phi_sol_op = [0]
phi_sol_trans = [0]
i_dif = []
i_dir = []

for i in df.index:
    sp = SunPos(param['Position']['lon'],
                param['Position']['lat'],
                year = i.year,
                month = i.month,
                day = i.day,
                hour = i.hour,
                minute = i.minute,
                shift = 0)
    sun_az.append(sp[0])
    sun_h.append(sp[1])

    i_dif.append((1/(1 + 2.71828182846**(7.997*(df.loc[i, 'Pyr_hor']/1000 - 0.586))))*df.loc[i, 'Pyr_hor'])
    i_dir.append(df.loc[i, 'Pyr_hor']-i_dif[-1])


    b = B(B_phi_int = 50, #  heat flow rate from internal heat source in Watts
          B_phi_HC = 0, # heat flow rate from heating, cooling and ventilation systems,
          B_theta_m_prev = B_theta_m[-1], # mass temperature previous step
          theta_e = df.loc[i, 'AMB_T'], # Outdoor temperature
          sun_h = sun_h[-1], # sun position
          sun_az = sun_az[-1], # Azimuth
          # I_dir = 0.5 *df.loc[i, 'Pyr_hor'],
          I_dir = i_dir[-1], # direct radiation
          # I_diff = 0.5*df.loc[i, 'Pyr_hor'],
          I_diff = i_dif[-1], # diffuse radation
          V_Vdot = 30, # m3/hr volumne flow
          V_theta_sup = df.loc[i, 'AMB_T'] ) # temperature of supply air(ambient in our case)

    if i.hour in [6, 12, 0] and i.minute == 0 and i.second == 0:
        B_theta_m.append(df.loc[i, 'Zone_T'])
    else:
        B_theta_m.append(b[0])

    B_theta_s.append(b[1])
    B_theta_air.append(b[2])
    B_theta_operativ.append(b[3])
    phi_sol_op.append(b[4])
    phi_sol_trans.append(b[5])


df['i_dif'] = i_dif
df['i_dir'] = i_dir

df['s_B_theta_m'] = B_theta_m[1:]
df['s_B_theta_s'] = B_theta_s[1:]
df['s_B_theta_air'] = B_theta_air[1:]
df['s_B_theta_operativ'] = B_theta_operativ[1:]
df['phi_sol_op'] = phi_sol_op[1:]
df['sun_h'] = sun_h
df['sun_az'] = sun_az

AE = [abs(mess-sim) for mess, sim in zip(df['Zone_T'], df['s_B_theta_air'])]
print('MAE: ' + str(round(sum(AE) / len(AE),2)))
print('Max. Error: ' + str(round(max(AE),2)))

#%%% Plotting TimeSeries

# ----------------------------------------------------------------------------------------------------------------
# ------------------------------------------------TWIN-PLOT-------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(6.4, 4.8/1.2))
ax2 = plt.twinx()
# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.plot(df.index, df['Zone_T'], c = 'tab:green', linestyle = '-', alpha = 1, label = 'Messwert Innentemperatur')
# Simulation
ax1.plot(df.index, df['s_B_theta_air'], c = 'tab:blue', linestyle = '-', alpha = 1, label = 'Simulation: Innenraumtemperatur')

# labels
ax1.legend(fontsize = 8, loc = 'lower left')
ax1.set_ylabel('Temperatur in °C')
ax1.yaxis.set_major_formatter('{:.0f}'.format)
ax1.set_ylim(0,40)
ax1.set_yticks(range(0,40,2), minor = True)
ax1.set_yticks(range(0,40,4))
# ------------------------------------------------Y2------------------------------------------------
# Data
ax2.plot(df.index, df['HM4_Qpkt'], c = 'tab:red', linestyle = '-', alpha = 1, label = 'Messwert: Wärmeleistung')
ax2.plot(df.index, df['Pyr_hor'], c = 'tab:orange', linestyle = '-', alpha = 1, label = 'Messwert: Globalstrahlung')
# labels
ax2.legend(fontsize = 8, loc = 'lower right')

ax2.set_ylabel('Wärmeleistung in W \n Globalstrahlung in W/m²')
ax2.yaxis.set_major_formatter('{:.0f}'.format)
ax2.set_ylim(0,2800)
ax2.set_yticks(range(0,2801,400))
# ------------------------------------------------X------------------------------------------------
ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*6) for i in range(int(672/4/6))], minor = True)
ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*24) for i in range(int(672/4/24))],
                # [(df.index[0] + timedelta(minutes=i*60*24)).strftime('%y-%m-%d') for i in range(int(672/4/24))],
                ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'],
                rotation = 0)
# ax1.set_xticks([i / 2 for i in range(0,7*2,1)], minor = True)
ax1.set_xlim(df.index[0], df.index[-1])
# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Raumtemperatur\ und\ Umgebung}$' + '\nZeitraum: ' + str(df.index[0].strftime('%Y-%m-%d')) + ' bis ' + str(df.index[-1].strftime('%Y-%m-%d')))
ax1.grid(True)
ax1.grid(which='minor', alpha = 0.3)
plt.show()

#%%% Plotting OneDay
for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:], [df.index[0] + timedelta(days = i+1) for i in range(7)][:]):

    # ----------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(6.4, 4.8/1.2))
    ax2 = plt.twinx()
    # ------------------------------------------------Y1------------------------------------------------
    # Data
    ax1.plot(df.index, df['AMB_T'], c = 'grey', linestyle = '-', alpha = 0.8, label = 'Messwert: Außentemperatur')
    ax1.plot(df.index, df['Zone_T'], c = 'tab:green', linestyle = '-', alpha = 1, label = 'Messwert Innentemperatur')

    # Simulation
    ax1.plot(df.index, df['s_B_theta_air'], c = 'tab:blue', linestyle = '-', alpha = 1, label = 'Simulation: Innenraumtemperatur')

    # labels
    ax1.legend(fontsize = 8, loc = 'upper left')
    ax1.set_ylabel('Temperatur in °C')
    ax1.yaxis.set_major_formatter('{:.0f}'.format)
    ax1.set_ylim(0,34)
    ax1.set_yticks(range(0,34,2), minor = True)
    ax1.set_yticks(range(0,36,4))
    # ------------------------------------------------Y2------------------------------------------------
    # Data
    ax2.plot(df.index, [i*1 for i in df['HM4_Qpkt']], c = 'tab:red', linestyle = '-', alpha = 1, label = 'Messwert: Wärmeleistung')
    ax2.plot(df.index, df['Pyr_hor'], c = 'tab:orange', linestyle = '-', alpha = 1, label = 'Messwert: Globalstrahlung')
    # labels
    ax2.legend(fontsize = 8, loc = 'upper right')

    ax2.set_ylabel('Wärmeleistung in W \n Globalstrahlung in W/m²')
    ax2.yaxis.set_major_formatter('{:.0f}'.format)
    ax2.set_ylim(0,2800)
    ax2.set_yticks(range(0,2801,400))
    # ------------------------------------------------X------------------------------------------------
    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*6) for i in range(int(672/4/6))], minor = True)
    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*24) for i in range(int(672/4/24))],
                    # [(df.index[0] + timedelta(minutes=i*60*24)).strftime('%y-%m-%d') for i in range(int(672/4/24))],
                    ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'],
                    rotation = 60)
    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))],
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))],
                    rotation = 90)
    # ax1.set_xticks([i / 2 for i in range(0,7*2,1)], minor = True)
    ax1.set_xlim(von, bis)
    # ------------------------------------------------Global------------------------------------------------
    ax1.set_title(r'$\bf{Raumtemperatur\ und\ Umgebung}$' + '\nZeitraum: ' + str(von.strftime('%Y-%m-%d')) + ' bis ' + str(bis.strftime('%Y-%m-%d')))
    ax1.grid(True)
    ax1.grid(which='minor', alpha = 0.3)
    plt.show()

#%%% Plotting Temperatures per Day
Mess_Tmin = []
Mess_Tmax = []
Mess_Tmean = []
Mess_Tschwankung = []
Sim_Tmin = []
Sim_Tmax = []
Sim_Tmean = []
Sim_Tschwankung = []

for i in df.index.day.unique():
    day = df[df.index.day == i]
    Mess_Tmin.append(day['Zone_T'].min())
    Mess_Tmax.append(day['Zone_T'].max())
    Mess_Tmean.append(day['Zone_T'].mean())
    Mess_Tschwankung.append((Mess_Tmax[-1] - Mess_Tmin[-1])/2)
    Sim_Tmin.append(day['s_B_theta_air'].min())
    Sim_Tmax.append(day['s_B_theta_air'].max())
    Sim_Tmean.append(day['s_B_theta_air'].mean())
    Sim_Tschwankung.append((Sim_Tmax[-1] - Sim_Tmin[-1])/2)

# --------------------------------------------------------------------------------------------------------------
# ------------------------------------------------Analyse-------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(6.4, 4.8))
ax2 = plt.twinx()

# ------------------------------------------------Y1------------------------------------------------
# Data
# ax1.scatter(0, 0, label = '', color = 'white')
ax1.scatter([i+0.5 for i in range(0,7,1)], Mess_Tmin, color = 'tab:blue', label = '   ', alpha = 0.5)
ax1.plot([i+0.5 for i in range(0,7,1)], Mess_Tmin, color = 'tab:blue', alpha = 0.4, linestyle = '--')
ax1.scatter([i+0.5 for i in range(0,7,1)], Mess_Tmax, color = 'tab:red', label = '   ', alpha = 0.5)
ax1.plot([i+0.5 for i in range(0,7,1)], Mess_Tmax, color = 'tab:red', alpha = 0.4, linestyle = '--')
ax1.scatter([i+0.5 for i in range(0,7,1)], Mess_Tmean, color = 'tab:green', label = '   ', alpha = 0.5)
ax1.plot([i+0.5 for i in range(0,7,1)], Mess_Tmean, color = 'tab:green', alpha = 0.4, linestyle = '--')

# ax1.scatter(0, 0, label = ' ', color = 'white')
ax1.scatter([i+0.5 for i in range(0,7,1)], Sim_Tmin, color = 'tab:blue', label = '     Minimalwert', marker = 'x')
ax1.plot([i+0.5 for i in range(0,7,1)], Sim_Tmin, color = 'tab:blue', alpha = 1, linestyle = '--')
ax1.scatter([i+0.5 for i in range(0,7,1)], Sim_Tmax, color = 'tab:red', label = '     Maximalwert', marker = 'x')
ax1.plot([i+0.5 for i in range(0,7,1)], Sim_Tmax, color = 'tab:red', alpha = 1, linestyle = '--')
ax1.scatter([i+0.5 for i in range(0,7,1)], Sim_Tmean, color = 'tab:green', label = '     Mittelwert', marker = 'x')
ax1.plot([i+0.5 for i in range(0,7,1)], Sim_Tmean, color = 'tab:green', alpha = 1, linestyle = '--')
# labels
ax1.set_ylabel('Temperatur in °C')
ax1.yaxis.set_major_formatter('{:.1f}'.format)
ax1.set_ylim(20,40)
ax1.legend(fontsize = 8, loc = 'upper left', ncol = 2, title = 'Messwert    Simulation              ', title_fontsize = 8)
# ------------------------------------------------Y2------------------------------------------------
# Data
ax2.bar([i+0.35 for i in range(0,7,1)], Mess_Tschwankung, width = 0.3, color = 'tab:green', label = 'Messwert:   Schwankung', alpha = 0.3)
ax2.bar([i+0.65 for i in range(0,7,1)], Sim_Tschwankung, width = 0.3, color = 'tab:green', label =  'Simulation: Schwankung', alpha = 0.6)
# labels
ax2.set_ylabel('Schwankung in K')
ax2.yaxis.set_major_formatter('{:.1f}'.format)
ax2.set_ylim(0,10)
ax2.set_yticks([0, 2.5, 5, 7.5, 10])
ax2.legend(fontsize = 8, loc = 'upper right')
# ------------------------------------------------X------------------------------------------------
ax1.set_xlim(0,7)
ax1.set_xticks(range(0,7,1), ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], rotation = 0)
# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Tagesbilanz\ der\ Innenraumtemperatur}$' + '\nZeitraum:' + str(df.index[0].strftime('%Y-%m-%d')) + ' bis ' + str(df.index[-1].strftime('%Y-%m-%d')))
# ax.set_title(r'$\bf{First\ Row}$' + '\nSecond Row')
ax1.grid(True)


plt.show()
# plt.close("all")

#%%

with open("b_building_sim.pkl", 'wb') as file:
    pickle.dump(df, file)

