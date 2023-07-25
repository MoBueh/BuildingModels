"""
@author: MoBueh
"""

#%% Libraries
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np
import copy
import math
from scipy.integrate import odeint
plt.rcParams["figure.dpi"] = 100
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
param['WS'] = {
    'H'             : 1.2,      # m         Speicherhöhe
    'N'             : 8,        # -         Anzahl Schichten
    'D'             : 0.6,      # m         Speicherdurchmesser
    'th'            : 0.02,      # m²        Wanddicke
    'k'             : 0.02*100,    # W/m²K     Wärmetransferkoeffizient
    'lambda_eff'    : 0.0015*1000    # W/mK      vertikale Wärmeleitfähigkeit von Wasser   
    }
print('Layer: ' + str(param['WS']['N']))

prop_water_rho = 980
prop_water_cp  = 4180
prop_brine_rho = 1026
prop_brine_cp = 3500

param_HS_V = 0.29556


#%% Modell
def WS(Vpkt_prim, T_prim_in, Vpkt_sek, T_sek_in, T_vektor_prev, Tamb):  
    #                        _______________
    # Vpkt_prim, T_prim_in  |     T_n+2     | T_sek_out
    # --------------------->|_______________|--------------------->
    #                       |     T_n+1     |
    #                       |_______________|
    #                       |     T_n       |
    #                       |_______________|
    #            T_prim_out |     T_n-1     | Vpkt_sek, T_sek_in
    # <---------------------|_______________|<---------------------   
    if str(T_prim_in) == 'nan':
        T_prim_in = 20
    if str(T_sek_in) == 'nan':
        T_sek_in = 20
    # -------------------Parameter-Massenströme-Randbedingungen----------------------------------------------------
    z_i = param['WS']['H']/param['WS']['N']                                     # (25) Höhe einer Schicht
    A_ext = np.pi*param['WS']['D']*z_i                                          # (26) Außenfläche einer Schicht
    A_i = np.pi*(param['WS']['D']-2*param['WS']['th'])**2/4                     # (27) Fläche zwischen zwei Layern
    m_i = A_i*z_i*prop['wasser']['rho']                                         # (28) Masse pro Schicht
    mpkt_prim = Vpkt_prim*prop['wasser']['rho']/3600
    mpkt_sek = Vpkt_sek*prop['wasser']['rho']/3600
    mpkt_i = mpkt_prim-mpkt_sek    
    Ti = copy.deepcopy(T_vektor_prev)
    # -------------------------------------------------------------------------------------------------------------
    # Die folgenden Gleichungen entsprechen Gleichung (29) 
    
    # -------------------------------------------------------------------------------------------------------------
    # ------------------------------BOTTOM-Layer-------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------
    def EB_BOTTOM(y, t):
    # -----------------------------------DGL-----------------------------------------------------------------------
        Ti[0] = y
        if mpkt_i>=0:                                                           # Resultierender Massenstrom ↓           
            dTdt = (0
                    -mpkt_sek*prop['wasser']['cp']*(Ti[0] - T_sek_in)           # Sek: mpkt-bedingte Änderung
                    -(param['WS']['k'] * A_ext * (Ti[0] - Tamb))                # Wärmeverluste an Umgebung
                    + (mpkt_i * prop['wasser']['cp'] * (Ti[0+1] - Ti[0]))       # Stofftransport
                    +((A_i*param['WS']['lambda_eff']/z_i) * (Ti[0+1] - Ti[0]))  # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt

        else:                                                                   # Resultierender Massenstrom ↑    
            dTdt = (0
                    -(mpkt_sek * prop['wasser']['cp'] * (Ti[0] - T_sek_in))     # Sek: mpkt-bedingte Änderung
                    -(param['WS']['k']*A_ext*(Ti[0] - Tamb))                    # Wärmeverluste an Umgebung
                    +0  # keine Schicht unterhalb vorhanden                     # Stofftransport                    
                    +(A_i*param['WS']['lambda_eff']/z_i*(Ti[0+1] - Ti[0]))      # Wärmeleitung zw. Layer 
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt 
        return dTdt
    # -------------------------------Solving-DGL-------------------------------------------------------------------
    Ti[0] = float(odeint(EB_BOTTOM, Ti[0], np.linspace(0, 900))[-1]) # solve ODE
    # ------------------------------------------------------------------------------------------------------------- 
    
    
    # -------------------------------------------------------------------------------------------------------------
    # --------------------------------TOP-Layer--------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------
    def EB_TOP(y, t):
    # -----------------------------------DGL-----------------------------------------------------------------------
        Ti[-1] = y
        if mpkt_i>=0:                                                           # Resultierender Massenstrom ↓            
            dTdt = (0
                    +(mpkt_prim * prop['wasser']['cp'] * (T_prim_in - Ti[-1]))  # Prim. mpkt-bedingte Änderung
                    -(param['WS']['k']*A_ext*(Ti[-1] - Tamb))                   # Wärmeverluste an Umgebung 
                    +0 # keine Schicht oberhalb vorhanden                       # Stofftransport
                    +((A_i*param['WS']['lambda_eff']/z_i)*(Ti[-1-1] - Ti[-1]))  # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt  
        
        else:                                                                   # Resultierender Massenstrom ↑   
            dTdt = (0
                    +(mpkt_prim * prop['wasser']['cp'] * (T_prim_in - Ti[-1]))  # Prim. mpkt-bedingte Änderung
                    -(param['WS']['k']*A_ext*(Ti[-1] - Tamb))                   # Wärmeverluste an Umgebung 
                    +(mpkt_i * prop['wasser']['cp'] * (Ti[-1] - Ti[-1-1]))      # Stofftransport
                    +(A_i*param['WS']['lambda_eff']/z_i*(Ti[-1-1] - Ti[-1]))     # Wärmeleitung zw. Layer
                    )/ (m_i*prop['wasser']['cp'])                               # aus dQ/dt 
        return dTdt
    # -------------------------------Solving-DGL-------------------------------------------------------------------
    Ti[-1] = float(odeint(EB_TOP, Ti[-1], np.linspace(0, 900))[-1]) # solve ODE
    # -------------------------------------------------------------------------------------------------------------
    
    
    # -------------------------------------------------------------------------------------------------------------
    # ------------------------------MIDDLE-Layer-------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------
    def EB_MIDDLE(y, t): # Function for the middle layer
        Ti[i] = y
    # -----------------------------------DGL-----------------------------------------------------------------------
        if mpkt_i>=0:                                                           # Resultierender Massenstrom ↓ 
            dTdt = (0
                    +0 # kein Anschluss in diesen Layern                        # Sek/prim mpkt-bedingte Änderung
                    -(param['WS']['k']*A_ext*(Ti[i] - Tamb))                    # Wärmeverluste an Umgebung
                    +(mpkt_i*prop['wasser']['cp']*(Ti[i+1] - Ti[i]))            # Stofftransport
                    +((A_i*param['WS']['lambda_eff'])/z_i * (Ti[i+1] - 2*Ti[i] + Ti[i-1])) # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt
   
        else:                                                                   # Resultierender Massenstrom ↑ 
            dTdt = (0
                    +0 # kein Anschluss in diesen Layern                        # Sek/prim mpkt-bedingte Änderung
                    -(param['WS']['k']*A_ext*(Ti[i] - Tamb))                    # Wärmeverluste an Umgebung
                    +(mpkt_i*prop['wasser']['cp']*(Ti[i] - Ti[i-1]))            # Stofftransport
                    +((A_i*param['WS']['lambda_eff'])/z_i * (Ti[i+1] - 2*Ti[i] + Ti[i-1])) # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt 
        return dTdt
    # -------------------------------Solving-DGL-------------------------------------------------------------------
    # for i in range(1, param['WS']['N'] - 1):
    #     Ti[i] = float(odeint(EB_MIDDLE, Ti[i], np.linspace(0, 900))[-1])
    if mpkt_i>=0:                                                               # Lösen der Schichtung in Richtung mpkt_i ↓
        for i in list(reversed(range(1, param['WS']['N'] - 1))):
            Ti[i] = float(odeint(EB_MIDDLE, Ti[i], np.linspace(0, 900))[-1])
    else:                                                                       # Lösen der Schichtung in Richtung mpkt_i ↑
        for i in list(range(1, param['WS']['N'] - 1)):
            Ti[i] = float(odeint(EB_MIDDLE, Ti[i], np.linspace(0, 900))[-1])
    # -------------------------------------------------------------------------------------------------------------    
     
    
    # ---------------------------Return-Temperatures---------------------------------------------------------------
    return Ti

def HS(HS_prim_theta_in, HS_prim_Vdot, HS_sec_theta_in, HS_sec_Vdot, HS_ter_theta_in, HS_ter_Vdot, HS_theta_prev):      # Speichertemperatur im Zeitschritt zuvor
    # print('nan')    
    if str(HS_prim_theta_in) == 'nan':
        HS_prim_theta_in = 20
        HS_prim_Vdot = 0
    if str(HS_sec_theta_in) == 'nan':
        HS_sec_theta_in = 20
        HS_sec_Vdot = 0
    
    A = HS_prim_Vdot/3600 * prop_water_rho*prop_water_cp
    B = HS_sec_Vdot/3600 * prop_water_rho*prop_water_cp
    C = HS_ter_Vdot/3600 * prop_brine_rho*prop_brine_cp
    D = param_HS_V*prop_water_rho*prop_water_cp/900
    HS_theta = (A*HS_prim_theta_in + B*HS_sec_theta_in + C*HS_ter_theta_in + D * HS_theta_prev)/(A + B + C + D)
    return HS_theta


#%% Testing
# WS_theta = [[10]*param['WS']['N']]

# dauer = 45
# for i in range(0,dauer):
#     # print('Loop: ' + str(i) )
#     ws = WS(Vpkt_prim       = 800/1000, 
#             T_prim_in       = 25, 
#             Vpkt_sek        = 500/1000, 
#             T_sek_in        = 12, 
#             T_vektor_prev   = WS_theta[-1],
#             Tamb            = 20)
#     WS_theta.append(ws)

# for i in range(param['WS']['N']):
#     plt.plot(range(0, dauer+1), [p[i] for p in WS_theta], label=[str(a) for a in list(range(param['WS']['N']))][i], color = plt.cm.RdYlBu(np.linspace(0, 1, param['WS']['N']))[::-1][i])
# plt.legend()
# plt.grid()



#%% Winter
#%%% Simulation
df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_winter.csv', index_col = 0)
df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')
df.loc[df['HM10_Vpkt'] < 80, 'HM10_Vpkt'] = 0
df.loc[df['WP_COP'] < 2, 'HM10_Vpkt'] = 0
df.loc[df['HM10_Vpkt'] == 0, ['HM10_T_RET', 'HM10_T_SUP']] = np.nan
df.loc[df['HM4_Vpkt'] < 20, 'HM4_Vpkt'] = 0
df.loc[df['HM4_Vpkt'] == 0, ['HM4_T_RET', 'HM4_T_SUP']] = np.nan

a = df[['HM4_Vpkt', 'HM4_T_RET', 'HM4_T_SUP', 'HM10_Vpkt', 'HM10_T_RET', 'HM10_T_SUP', 'Warm Unten', 'Warm Oben', 'Kalt Oben', 'Kalt Unten']]


WS_theta = [np.linspace(df['Warm Unten'][0], df['Warm Oben'][0], param['WS']['N']).tolist()]
Mixed = [(df['Warm Unten'][0] +  df['Warm Oben'][0])/2]

for i in df.index:
    ws = WS(Vpkt_prim       = df.loc[i, 'HM10_Vpkt']/1000, 
            T_prim_in       = df.loc[i, 'HM10_T_SUP'], 
            Vpkt_sek        = df.loc[i, 'HM4_Vpkt']/1000, 
            T_sek_in        = df.loc[i, 'HM4_T_RET'], 
            T_vektor_prev   = WS_theta[-1],
            Tamb            = df.loc[i, 'Zone_T'])
    WS_theta.append(ws)
    
    hs = HS(HS_prim_theta_in = df.loc[i, 'HM10_T_SUP'], 
            HS_prim_Vdot =df.loc[i, 'HM10_Vpkt']/1000, 
            HS_sec_theta_in = df.loc[i, 'HM4_T_RET'], 
            HS_sec_Vdot = df.loc[i, 'HM4_Vpkt']/1000, 
            HS_ter_theta_in = 0,
            HS_ter_Vdot = 0, 
            HS_theta_prev = Mixed[-1])
    Mixed.append(hs)
    
df['WS_theta'] = WS_theta[1:]
df['Mixed'] = Mixed[1:]


# %%% Plotting Timeseries
ymax = 35
ymin = 15

Tag = 7

for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:Tag], [df.index[0] + timedelta(days = i+1) for i in range(7)][:Tag]):
#     for i in range(param['WS']['N']):
#         plt.plot(df.index, [p[i] for p in df['WS_theta']], label=[str(a) for a in list(range(param['WS']['N']))][i], color = plt.cm.RdYlBu(np.linspace(0, 1, param['WS']['N']))[::-1][i])
#     plt.xlim(von, bis)
#     plt.legend()
#     plt.grid()
#     plt.show()

    # ----------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(6.4, 4.8/1.3))
    # fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.9, 0.1]})
    
    # ------------------------------------------------Y1------------------------------------------------
    ax1.plot(df.index, df['Warm Oben'], label = 'Messwert: Top', c = 'tab:red', alpha = 0.3, linestyle = '--', marker = 'o', markersize = 4)
    ax1.plot(df.index[::4], df['Warm Unten'][::4], label = 'Messwert: Bottom', c = 'tab:blue', alpha = 0.3, linestyle = '--', marker = 'x', markersize = 4)
    
    ax1.plot(df.index, df['Mixed'], label = '1Layer', c = 'grey', alpha = 0.6)
    
    ax1.plot(df.index, [posEl[0] for posEl in df['WS_theta']], label = 'Simulation: Bottom', c = 'tab:blue')
    ax1.plot(df.index, [posEl[-1] for posEl in df['WS_theta']], label = 'Simulation: Top', c = 'tab:red')


    # labels
    ax1.set_ylabel('Temperatur in °C')
    ax1.set_ylim(ymin,ymax)
    # ax1.set_yticks(np.arange(-10,25,2.5), minor = True)
    # ax1.set_yticks(range(-10,25,5))
    ax1.set_title(r'$\bf{Speichertemperatur}$' + '\nZeitraum: ' + str(von.strftime('%Y-%m-%d')) + ' bis ' + str(bis.strftime('%Y-%m-%d')))

    
    # ------------------------------------------------X1------------------------------------------------
    ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    ['']*84,
                    rotation = 90)
    ax1.set_xlim(von, bis)
    ax1.grid(True)
    ax1.grid(which='minor', alpha = 0.3)
    ax1.legend(fontsize = 6, loc = 'upper left', ncol = 2)
    
    
    ax1_twin = ax1.twinx()
    
    ax1_twin.plot(df.index, [1 if value > 0 else value for value in df['HM4_Vpkt']],label = 'Entladung (Tabs)', c = 'tab:green')
    ax1_twin.fill_between(df.index, [1 if value > 0 else value for value in df['HM4_Vpkt']],color = 'tab:green', alpha = 0.3)
    
    ax1_twin.plot(df.index, [1 if value > 0 else value for value in df['HM10_Vpkt']], label = 'Beladung (WP)', c = 'tab:orange')
    ax1_twin.fill_between(df.index, [1 if value > 0 else value for value in df['HM10_Vpkt']],color = 'tab:orange', alpha = 0.3)

    ax1_twin.legend(fontsize = 6, loc = 'upper right', ncol = 2)
    ax1_twin.set_ylabel('AN / AUS')
    ax1_twin.set_yticks(range(0,2,1))
    ax1_twin.set_ylim(0,8)
    plt.show()
    
    # abweichung_top = [abs(a-b) for a,b in zip(df['Warm Oben'][von:bis], [posEl[-1] for posEl in df['WS_theta'][von:bis]])] 
    # MAE_top = sum(abweichung_top) / len(abweichung_top)
    # print('MAE_top: ' +str(MAE_top))

    # abweichung_bottom = [abs(a-b) for a,b in zip(df['Warm Unten'][von:bis][::4], [posEl[0] for posEl in df['WS_theta'][von:bis][::4]])] 
    # MAE_bottom = sum(abweichung_bottom) / len(abweichung_bottom)
    # print('MAE_bottom: ' +str(MAE_bottom))
    # print(' ')
    
# #%%% 7 Tage, 1 Plot    
# ymax = 35
# ymin = 15

# Tag = 7
# fig, axes = plt.subplots(4, 2, figsize=(12, 15))
# fig.suptitle(r'$\bf{Wärmespeicher\ während\ der\ Sommerwoche}$', fontsize = 16)
# for von,bis,z in zip([df.index[0] + timedelta(days = i) for i in range(7)][:Tag], [df.index[0] + timedelta(days = i+1) for i in range(7)][:Tag], range(0,7)):
#     row = z // 2  # Bestimmen der aktuellen Zeile
#     col = z % 2   # Bestimmen der aktuellen Spalte
#     ymax = 35
#     ymin = 15
#     # Erstellen des Plots in der entsprechenden Achse
#     ax1 = axes[row, col]
#     # ax.plot(df['x'], df['y'])
#     ax1.plot(df.index, df['Warm Oben'], label = 'Messwert: Top', c = 'tab:red', alpha = 0.3, linestyle = '--', marker = 'o', markersize = 4)
#     ax1.plot(df.index[::4], df['Warm Unten'][::4], label = 'Messwert: Bottom', c = 'tab:blue', alpha = 0.3, linestyle = '--', marker = 'x', markersize = 4)
    
#     ax1.plot(df.index, df['Mixed'], label = '1Layer', c = 'grey', alpha = 0.6)
    
#     ax1.plot(df.index, [posEl[0] for posEl in df['WS_theta']], label = 'Simulation: Bottom', c = 'tab:blue')
#     ax1.plot(df.index, [posEl[-1] for posEl in df['WS_theta']], label = 'Simulation: Top', c = 'tab:red')
       
       
#     # labels
#     ax1.set_ylabel('Temperatur in °C')
#     ax1.set_ylim(ymin,ymax)
#     # ax1.set_yticks(np.arange(-10,25,2.5), minor = True)
#     # ax1.set_yticks(range(-10,25,5))
#     ax1.set_title('Datum: ' + str(von.strftime('%Y-%m-%d')))
    
#     # ------------------------------------------------X1------------------------------------------------
#     ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
#                     ['']*84,
#                     rotation = 90)
#     ax1.set_xlim(von, bis)
#     ax1.grid(True)
#     ax1.grid(which='minor', alpha = 0.3)
#     ax1.legend(fontsize = 6, loc = 'upper left', ncol = 2)
    
    
#     ax1_twin = ax1.twinx()
    
#     ax1_twin.plot(df.index, [1 if value > 0 else value for value in df['HM4_Vpkt']],label = 'Entladung (Tabs)', c = 'tab:green')
#     ax1_twin.fill_between(df.index, [1 if value > 0 else value for value in df['HM4_Vpkt']],color = 'tab:green', alpha = 0.3)
    
#     ax1_twin.plot(df.index, [1 if value > 0 else value for value in df['HM10_Vpkt']], label = 'Beladung (WP)', c = 'tab:orange')
#     ax1_twin.fill_between(df.index, [1 if value > 0 else value for value in df['HM10_Vpkt']],color = 'tab:orange', alpha = 0.3)
       
#     ax1_twin.legend(fontsize = 6, loc = 'upper right', ncol = 2)
#     ax1_twin.set_ylabel('AN / AUS')
#     ax1_twin.set_yticks(range(0,2,1))
#     ax1_twin.set_ylim(0,8)
# axes[3, 1].axis('off')     
# # Das Layout anpassen
# plt.tight_layout()
# # Anzeigen oder Speichern des Plots
# plt.show()

#%%% Plotting Temperatures per Day
Top_Mess_Tmin = []
Top_Mess_Tmax = []
Top_Mess_Tmean = []
Top_Mess_Tschwankung = []
Top_Sim_Tmin = []
Top_Sim_Tmax = []
Top_Sim_Tmean = []
Top_Sim_Tschwankung = []
Top_MAE = []

Bot_Mess_Tmin = []
Bot_Mess_Tmax = []
Bot_Mess_Tmean = []
Bot_Mess_Tschwankung = []
Bot_Sim_Tmin = []
Bot_Sim_Tmax = []
Bot_Sim_Tmean = []
Bot_Sim_Tschwankung = []
Bot_MAE = []

for i in df.index.day.unique():
    day = df[df.index.day == i]
    Top_Mess_Tmin.append(day['Warm Oben'].min())
    Top_Mess_Tmax.append(day['Warm Oben'].max())
    Top_Mess_Tmean.append(day['Warm Oben'].mean())
    Top_Mess_Tschwankung.append((Top_Mess_Tmax[-1] - Top_Mess_Tmin[-1])/2)
    Top_Sim_Tmin.append(min([posEl[-1] for posEl in day['WS_theta']]))
    Top_Sim_Tmax.append(max([posEl[-1] for posEl in day['WS_theta']]))
    Top_Sim_Tmean.append(sum([posEl[-1] for posEl in day['WS_theta']])/len([posEl[-1] for posEl in day['WS_theta']]))
    Top_Sim_Tschwankung.append((Top_Sim_Tmax[-1] - Top_Sim_Tmin[-1])/2)
    abweichung = [abs(a-b) for a,b in zip(day['Warm Oben'], [posEl[-1] for posEl in day['WS_theta']])]
    Top_MAE.append((sum(abweichung) / len(abweichung)))

    Bot_Mess_Tmin.append(day['Warm Unten'].min())
    Bot_Mess_Tmax.append(day['Warm Unten'].max())
    Bot_Mess_Tmean.append(day['Warm Unten'].mean())
    Bot_Mess_Tschwankung.append((Bot_Mess_Tmax[-1] - Bot_Mess_Tmin[-1])/2)
    Bot_Sim_Tmin.append(min([posEl[0] for posEl in day['WS_theta']]))
    Bot_Sim_Tmax.append(max([posEl[0] for posEl in day['WS_theta']]))
    Bot_Sim_Tmean.append(sum([posEl[0] for posEl in day['WS_theta']])/len([posEl[0] for posEl in day['WS_theta']]))
    Bot_Sim_Tschwankung.append((Bot_Sim_Tmax[-1] - Bot_Sim_Tmin[-1])/2)
    abweichung_bot = [abs(a-b) for a,b in zip(day['Warm Unten'][::4], [posEl[0] for posEl in day['WS_theta'][::4]])]
    Bot_MAE.append((sum(abweichung_bot) / len(abweichung_bot)))

# --------------------------------------------------------------------------------------------------------------
# ------------------------------------------------Analyse-------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]}, figsize=(6.4, 4.8*2))
ax1_twin = ax1.twinx()
ax2_twin = ax2.twinx()
# ------------------------------------------------Y1------------------------------------------------
# Data
ax1.plot([i+0.5 for i in range(0,7,1)], Top_Mess_Tmin, color = 'tab:blue', label = '   ', marker = 'o', markersize = 5 ,alpha = 0.5, linestyle = '--')
ax1.plot([i+0.5 for i in range(0,7,1)], Top_Mess_Tmax, color = 'tab:red', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )
ax1.plot([i+0.5 for i in range(0,7,1)], Top_Mess_Tmean, color = 'tab:green', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )

ax1.plot([i+0.5 for i in range(0,7,1)], Top_Sim_Tmin, color = 'tab:blue', alpha = 1, linestyle = '--', label = 'Minimalwert', marker = 'x', markersize = 5)
ax1.plot([i+0.5 for i in range(0,7,1)], Top_Sim_Tmax, color = 'tab:red', alpha = 1, linestyle = '--', label = 'Maximalwert', marker = 'x', markersize = 5)
ax1.plot([i+0.5 for i in range(0,7,1)], Top_Sim_Tmean, color = 'tab:green', alpha = 1, linestyle = '--', label = 'Mittelwert', marker = 'x', markersize = 5)

# labels
ax1.set_ylabel('Temperatur in °C')
ax1.set_yticks(np.arange(15,36,2.5))
ax1.set_yticks(np.arange(15,36,0.5), minor = True)
ax1_twin.yaxis.set_major_formatter('{:.1f}'.format)
ax1.set_ylim(15,35)
ax1.legend(fontsize = 8, loc = 'upper left', ncol = 2, title = 'Messwert    Simulation              ', title_fontsize = 8)
# ------------------------------------------------Y1_Twin------------------------------------------------
# Data
ax1_twin.bar([i+0.35 for i in range(0,7,1)], Top_Mess_Tschwankung, width = 0.3, color = 'tab:green', label = 'Messwert:   Schwankung', alpha = 0.3)
ax1_twin.bar([i+0.65 for i in range(0,7,1)], Top_Sim_Tschwankung, width = 0.3, color = 'tab:green', label =  'Simulation: Schwankung', alpha = 0.6)
ax1_twin.plot([i+0.5 for i in range(0,7,1)], Top_MAE, color = 'grey', alpha = 1, linestyle = '--', label = 'MAE', marker = 'x', markersize = 5)

# labels
ax1_twin.set_ylabel('Schwankung und MAE in K')
ax1_twin.yaxis.set_major_formatter('{:.1f}'.format)
ax1_twin.set_ylim(0,16)
ax1_twin.set_yticks(range(0,16,2))
ax1_twin.legend(fontsize = 8, loc = 'upper right')
# ------------------------------------------------X------------------------------------------------
ax1.set_xlim(0,7)
ax1.set_xticks(range(0,7,1), ['']*7, rotation = 60)
# ------------------------------------------------Global------------------------------------------------
ax1.set_title(r'$\bf{Tagesbilanz\ der\ obersten\ Schicht}$')
ax1.grid(True)


# ------------------------------------------------Y2------------------------------------------------
# Data
ax2.plot([i+0.5 for i in range(0,7,1)], Bot_Mess_Tmin, color = 'tab:blue', label = '   ', marker = 'o', markersize = 5 ,alpha = 0.5, linestyle = '--')
ax2.plot([i+0.5 for i in range(0,7,1)], Bot_Mess_Tmax, color = 'tab:red', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )
ax2.plot([i+0.5 for i in range(0,7,1)], Bot_Mess_Tmean, color = 'tab:green', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )

ax2.plot([i+0.5 for i in range(0,7,1)], Bot_Sim_Tmin, color = 'tab:blue', alpha = 1, linestyle = '--', label = 'Minimalwert', marker = 'x', markersize = 5)
ax2.plot([i+0.5 for i in range(0,7,1)], Bot_Sim_Tmax, color = 'tab:red', alpha = 1, linestyle = '--', label = 'Maximalwert', marker = 'x', markersize = 5)
ax2.plot([i+0.5 for i in range(0,7,1)], Bot_Sim_Tmean, color = 'tab:green', alpha = 1, linestyle = '--', label = 'Mittelwert', marker = 'x', markersize = 5)

# labels
ax2.set_ylabel('Temperatur in °C')
ax2.set_yticks(np.arange(15,36,2.5))
ax2.set_yticks(np.arange(15,36,0.5), minor = True)
ax2_twin.yaxis.set_major_formatter('{:.1f}'.format)
ax2.set_ylim(15,35)
ax2.legend(fontsize = 8, loc = 'upper left', ncol = 2, title = 'Messwert    Simulation              ', title_fontsize = 8)
# ------------------------------------------------Y2_Twin------------------------------------------------
# Data
ax2_twin.bar([i+0.35 for i in range(0,7,1)], Bot_Mess_Tschwankung, width = 0.3, color = 'tab:green', label = 'Messwert:   Schwankung', alpha = 0.3)
ax2_twin.bar([i+0.65 for i in range(0,7,1)], Bot_Sim_Tschwankung, width = 0.3, color = 'tab:green', label =  'Simulation: Schwankung', alpha = 0.6)
ax2_twin.plot([i+0.5 for i in range(0,7,1)], Bot_MAE, color = 'grey', alpha = 1, linestyle = '--', label = 'MAE', marker = 'x', markersize = 5)

# labels
ax2_twin.set_ylabel('Schwankung und MAE in K')
ax2_twin.yaxis.set_major_formatter('{:.1f}'.format)
ax2_twin.set_ylim(0,16)
ax2_twin.set_yticks(range(0,16,2))
ax2_twin.legend(fontsize = 8, loc = 'upper right')
# ------------------------------------------------X------------------------------------------------
ax2.set_xlim(0,7)
ax2.set_xticks(range(0,7,1), ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], rotation = 60)
# ------------------------------------------------Global------------------------------------------------
ax2.set_title(r'$\bf{Tagesbilanz\ der\ untersten\ Schicht}$')

ax2.grid(True)

E_Top = [(b-a) for a,b in zip(df['Warm Oben'], [posEl[-1] for posEl in df['WS_theta']])]
E_Bot = [(b-a) for a,b in zip(df['Warm Unten'][::4], [posEl[0] for posEl in df['WS_theta'][::4]])]


print('')
print('Winterfall:')
print('MAE Top: ' + str(round(sum(Top_MAE)/len(Top_MAE),2)))
print('ME Top: ' + str(round(sum(E_Top)/len(E_Top),2)))
print('Top_Max: ' + str(round(max(E_Top),2)))
print('Top_Min: ' + str(round(min(E_Top),2)))
print('-')

print('MAE Bot: ' + str(round(sum(Bot_MAE)/len(Bot_MAE),2)))
print('ME Bot: ' + str(round(sum(E_Bot)/len(E_Bot),2)))
print('Bot_Max: ' + str(round(max(E_Bot),2)))
print('Bot_Min: ' + str(round(min(E_Bot),2)))

plt.show()





#%% Sommer
# #%%% Simulation
# df = pd.read_csv('C:/Users/49157/Documents/Python/SHK4FE-Digital-Twin/Data/Monitoring/TestSet/TestSet_summer.csv', index_col = 0)
# df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')
# df.loc[(df['HM10_Vpkt'] < 80) | (df['WP_COP'] < 2), 'HM10_Vpkt'] = 0
# df.loc[df['HM10_Vpkt'] == 0, ['HM10_T_RET', 'HM10_T_SUP']] = np.nan
# df.loc[df['HM1_Vpkt'] < 400, 'HM1_Vpkt'] = 0
# df.loc[df['HM1_Vpkt'] == 0, ['HM1_T_RET', 'HM1_T_SUP']] = np.nan


# WS_theta = [np.linspace(df['Warm Unten'][0], df['Warm Oben'][0], param['WS']['N']).tolist()]
# Mixed = [(df['Warm Unten'][0] +  df['Warm Oben'][0])/2]

# for i in df.index:
#     # print(i)
#     # print(b)
#     ws = WS(Vpkt_prim       = df.loc[i, 'HM10_Vpkt']/1000, 
#             T_prim_in       = df.loc[i, 'HM10_T_SUP']-2, 
#             Vpkt_sek        = df.loc[i, 'HM1_Vpkt']/1000, 
#             T_sek_in        = df.loc[i, 'HM1_T_RET'], 
#             T_vektor_prev   = WS_theta[-1],
#             Tamb            = df.loc[i, 'Zone_T'])
#     WS_theta.append(ws)
    
#     if i.hour == 0 and i.minute == 0 and i.second == 0:
#         WS_theta[-1] = [df.loc[i, 'Warm Oben']]*param['WS']['N']    
    
#     hs = HS(HS_prim_theta_in = df.loc[i, 'HM10_T_SUP'], 
#             HS_prim_Vdot =df.loc[i, 'HM10_Vpkt']/1000, 
#             HS_sec_theta_in = df.loc[i, 'HM1_T_RET'], 
#             HS_sec_Vdot = df.loc[i, 'HM1_Vpkt']/1000, 
#             HS_ter_theta_in = 0,
#             HS_ter_Vdot = 0, 
#             HS_theta_prev = Mixed[-1])
#     Mixed.append(hs)
#     if i.hour == 0 and i.minute == 0 and i.second == 0:
#         Mixed[-1] = (df.loc[i, 'Warm Oben'] + df.loc[i, 'Warm Unten'])/2
    
# df['WS_theta'] = WS_theta[1:]
# df['Mixed'] = Mixed[1:]

# #%%% Plotting Timeseries

# Tag = 7

# for von,bis in zip([df.index[0] + timedelta(days = i) for i in range(7)][:Tag], [df.index[0] + timedelta(days = i+1) for i in range(7)][:Tag]):
#     # for i in range(param['WS']['N']):
#     #     plt.plot(df.index, [p[i] for p in df['WS_theta']], label=[str(a) for a in list(range(param['WS']['N']))][i], color = plt.cm.RdYlBu(np.linspace(0, 1, param['WS']['N']))[::-1][i])
#     # plt.xlim(von, bis)
#     # plt.legend()
#     # plt.grid()
#     # plt.show()

#     # ----------------------------------------------------------------------------------------------------------------
#     # ------------------------------------------------TWIN-PLOT-------------------------------------------------------
#     # ----------------------------------------------------------------------------------------------------------------
#     fig, ax1 = plt.subplots(figsize=(6.4, 4.8/1.3))
#     # fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.9, 0.1]})
    
#     # ------------------------------------------------Y1------------------------------------------------
#     ax1.plot(df.index, df['Warm Oben'], label = 'Messwert: Top', c = 'tab:red', alpha = 0.3, linestyle = '--', marker = 'o', markersize = 4)
#     ax1.plot(df.index[::4], df['Warm Unten'][::4], label = 'Messwert: Bottom', c = 'tab:blue', alpha = 0.3, linestyle = '--', marker = 'x', markersize = 4)
    
#     # ax1.plot(df.index, df['Mixed'], label = '1Layer', c = 'grey', alpha = 0.6)
    
#     ax1.plot(df.index, [posEl[0] for posEl in df['WS_theta']], label = 'Simulation: Bottom', c = 'tab:blue')
#     ax1.plot(df.index, [posEl[-1] for posEl in df['WS_theta']], label = 'Simulation: Top', c = 'tab:red')


#     # labels
#     ax1.set_ylabel('Temperatur in °C')
#     ax1.set_ylim(25,55)
#     ax1.set_yticks(np.arange(25,55,2.5), minor = True)
#     ax1.set_yticks(range(25,55,5))
#     ax1.set_title(r'$\bf{Speichertemperatur}$' + '\nZeitraum: ' + str(von.strftime('%Y-%m-%d')) + ' bis ' + str(bis.strftime('%Y-%m-%d')))

    
#     # ------------------------------------------------X1------------------------------------------------
#     ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
#                     ['']*84,
#                     rotation = 90)
#     ax1.set_xlim(von, bis)
#     ax1.grid(True)
#     ax1.grid(which='minor', alpha = 0.3)
#     ax1.legend(fontsize = 6, loc = 'upper left', ncol = 2)
    
#     ax1_twin = ax1.twinx()
#     ax1_twin.plot(df.index, [2 if value > 0 else value for value in df['HM1_Vpkt']],label = 'Rückkühlung (PVT)', c = 'tab:green')
#     ax1_twin.fill_between(df.index, [2 if value > 0 else value for value in df['HM1_Vpkt']],color = 'tab:green', alpha = 0.3)
    
#     ax1_twin.plot(df.index, [2 if value > 0 else value for value in df['HM10_Vpkt']], label = 'Beladung (WP)', c = 'tab:orange')
#     ax1_twin.fill_between(df.index, [2 if value > 0 else value for value in df['HM10_Vpkt']],color = 'tab:orange', alpha = 0.3)

#     ax1_twin.legend(fontsize = 6, loc = 'upper right', ncol = 2)
#     ax1_twin.set_ylabel('AN / AUS')
#     ax1_twin.set_yticks(range(0,3,2), [0,1])
#     ax1_twin.set_ylim(0,12)

#     plt.show()
    
#     # abweichung_top = [abs(a-b) for a,b in zip(df['Warm Oben'][von:bis], [posEl[-1] for posEl in df['WS_theta'][von:bis]])] 
#     # MAE_top = sum(abweichung_top) / len(abweichung_top)
#     # print('MAE_top: ' +str(MAE_top))

#     # abweichung_bottom = [abs(a-b) for a,b in zip(df['Warm Unten'][von:bis][::4], [posEl[0] for posEl in df['WS_theta'][von:bis][::4]])] 
#     # MAE_bottom = sum(abweichung_bottom) / len(abweichung_bottom)
#     # print('MAE_bottom: ' +str(MAE_bottom))
#     # print(' ')

# #%%% 7 Tage, 1 Plot    
# ymax = 35
# ymin = 15

# Tag = 7
# fig, axes = plt.subplots(4, 2, figsize=(12, 15))
# fig.suptitle(r'$\bf{Wärmespeicher\ während\ der\ Winterwoche}$', fontsize = 16)
# for von,bis,z in zip([df.index[0] + timedelta(days = i) for i in range(7)][:Tag], [df.index[0] + timedelta(days = i+1) for i in range(7)][:Tag], range(0,7)):
#     row = z // 2  # Bestimmen der aktuellen Zeile
#     col = z % 2   # Bestimmen der aktuellen Spalte
#     ax1 = axes[row, col]
#     ymax = 35
#     ymin = 15
#     ax1.plot(df.index, df['Warm Oben'], label = 'Messwert: Top', c = 'tab:red', alpha = 0.3, linestyle = '--', marker = 'o', markersize = 4)
#     ax1.plot(df.index[::4], df['Warm Unten'][::4], label = 'Messwert: Bottom', c = 'tab:blue', alpha = 0.3, linestyle = '--', marker = 'x', markersize = 4)
    
#     # ax1.plot(df.index, df['Mixed'], label = '1Layer', c = 'grey', alpha = 0.6)
    
#     ax1.plot(df.index, [posEl[0] for posEl in df['WS_theta']], label = 'Simulation: Bottom', c = 'tab:blue')
#     ax1.plot(df.index, [posEl[-1] for posEl in df['WS_theta']], label = 'Simulation: Top', c = 'tab:red')


#     # labels
#     ax1.set_ylabel('Temperatur in °C')
#     ax1.set_ylim(25,55)
#     ax1.set_yticks(np.arange(25,55,2.5), minor = True)
#     ax1.set_yticks(range(25,55,5))
#     ax1.set_title('\nDatum: ' + str(von.strftime('%Y-%m-%d')))

    
#     # ------------------------------------------------X1------------------------------------------------
#     ax1.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
#                     ['']*84,
#                     rotation = 90)
#     ax1.set_xlim(von, bis)
#     ax1.grid(True)
#     ax1.grid(which='minor', alpha = 0.3)
#     ax1.legend(fontsize = 6, loc = 'upper left', ncol = 2)
    
#     ax1_twin = ax1.twinx()
#     ax1_twin.plot(df.index, [2 if value > 0 else value for value in df['HM1_Vpkt']],label = 'Rückkühlung (PVT)', c = 'tab:green')
#     ax1_twin.fill_between(df.index, [2 if value > 0 else value for value in df['HM1_Vpkt']],color = 'tab:green', alpha = 0.3)
    
#     ax1_twin.plot(df.index, [2 if value > 0 else value for value in df['HM10_Vpkt']], label = 'Beladung (WP)', c = 'tab:orange')
#     ax1_twin.fill_between(df.index, [2 if value > 0 else value for value in df['HM10_Vpkt']],color = 'tab:orange', alpha = 0.3)

#     ax1_twin.legend(fontsize = 6, loc = 'upper right', ncol = 2)
#     ax1_twin.set_ylabel('AN / AUS')
#     ax1_twin.set_yticks(range(0,3,2), [0,1])
#     ax1_twin.set_ylim(0,12)
# axes[3, 1].axis('off') 
# # Das Layout anpassen
# plt.tight_layout()
# # Anzeigen oder Speichern des Plots
# plt.show()
  


# #%%% Plotting Temperatures per Day
# Top_Mess_Tmin = []
# Top_Mess_Tmax = []
# Top_Mess_Tmean = []
# Top_Mess_Tschwankung = []
# Top_Sim_Tmin = []
# Top_Sim_Tmax = []
# Top_Sim_Tmean = []
# Top_Sim_Tschwankung = []
# Top_MAE = []

# Bot_Mess_Tmin = []
# Bot_Mess_Tmax = []
# Bot_Mess_Tmean = []
# Bot_Mess_Tschwankung = []
# Bot_Sim_Tmin = []
# Bot_Sim_Tmax = []
# Bot_Sim_Tmean = []
# Bot_Sim_Tschwankung = []
# Bot_MAE = []

# for i in df.index.day.unique():
#     day = df[df.index.day == i]
#     Top_Mess_Tmin.append(day['Warm Oben'].min())
#     Top_Mess_Tmax.append(day['Warm Oben'].max())
#     Top_Mess_Tmean.append(day['Warm Oben'].mean())
#     Top_Mess_Tschwankung.append((Top_Mess_Tmax[-1] - Top_Mess_Tmin[-1])/2)
#     Top_Sim_Tmin.append(min([posEl[-1] for posEl in day['WS_theta']]))
#     Top_Sim_Tmax.append(max([posEl[-1] for posEl in day['WS_theta']]))
#     Top_Sim_Tmean.append(sum([posEl[-1] for posEl in day['WS_theta']])/len([posEl[-1] for posEl in day['WS_theta']]))
#     Top_Sim_Tschwankung.append((Top_Sim_Tmax[-1] - Top_Sim_Tmin[-1])/2)
#     abweichung = [abs(a-b) for a,b in zip(day['Warm Oben'], [posEl[-1] for posEl in day['WS_theta']])]
#     Top_MAE.append((sum(abweichung) / len(abweichung)))

#     Bot_Mess_Tmin.append(day['Warm Unten'].min())
#     Bot_Mess_Tmax.append(day['Warm Unten'].max())
#     Bot_Mess_Tmean.append(day['Warm Unten'].mean())
#     Bot_Mess_Tschwankung.append((Bot_Mess_Tmax[-1] - Bot_Mess_Tmin[-1])/2)
#     Bot_Sim_Tmin.append(min([posEl[0] for posEl in day['WS_theta']]))
#     Bot_Sim_Tmax.append(max([posEl[0] for posEl in day['WS_theta']]))
#     Bot_Sim_Tmean.append(sum([posEl[0] for posEl in day['WS_theta']])/len([posEl[0] for posEl in day['WS_theta']]))
#     Bot_Sim_Tschwankung.append((Bot_Sim_Tmax[-1] - Bot_Sim_Tmin[-1])/2)
#     abweichung_bot = [abs(a-b) for a,b in zip(day['Warm Unten'][::4], [posEl[0] for posEl in day['WS_theta'][::4]])]
#     Bot_MAE.append((sum(abweichung_bot) / len(abweichung_bot)))

# # --------------------------------------------------------------------------------------------------------------
# # ------------------------------------------------Analyse-------------------------------------------------------
# # --------------------------------------------------------------------------------------------------------------
# fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]}, figsize=(6.4, 4.8*2))
# ax1_twin = ax1.twinx()
# ax2_twin = ax2.twinx()
# # ------------------------------------------------Y1------------------------------------------------
# # Data
# ax1.plot([i+0.5 for i in range(0,7,1)], Top_Mess_Tmin, color = 'tab:blue', label = '   ', marker = 'o', markersize = 5 ,alpha = 0.5, linestyle = '--')
# ax1.plot([i+0.5 for i in range(0,7,1)], Top_Mess_Tmax, color = 'tab:red', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )
# ax1.plot([i+0.5 for i in range(0,7,1)], Top_Mess_Tmean, color = 'tab:green', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )

# ax1.plot([i+0.5 for i in range(0,7,1)], Top_Sim_Tmin, color = 'tab:blue', alpha = 1, linestyle = '--', label = 'Minimalwert', marker = 'x', markersize = 5)
# ax1.plot([i+0.5 for i in range(0,7,1)], Top_Sim_Tmax, color = 'tab:red', alpha = 1, linestyle = '--', label = 'Maximalwert', marker = 'x', markersize = 5)
# ax1.plot([i+0.5 for i in range(0,7,1)], Top_Sim_Tmean, color = 'tab:green', alpha = 1, linestyle = '--', label = 'Mittelwert', marker = 'x', markersize = 5)

# # labels
# ax1.set_ylabel('Temperatur in °C')
# ax1.set_yticks(np.arange(15,65,2.5))
# ax1.set_yticks(np.arange(15,65,0.5), minor = True)
# ax1.set_ylim(15,65)
# ax1_twin.yaxis.set_major_formatter('{:.1f}'.format)
# ax1.legend(fontsize = 8, loc = 'upper left', ncol = 2, title = 'Messwert    Simulation              ', title_fontsize = 8)
# # ------------------------------------------------Y1_Twin------------------------------------------------
# # Data
# ax1_twin.bar([i+0.35 for i in range(0,7,1)], Top_Mess_Tschwankung, width = 0.3, color = 'tab:green', label = 'Messwert:   Schwankung', alpha = 0.3)
# ax1_twin.bar([i+0.65 for i in range(0,7,1)], Top_Sim_Tschwankung, width = 0.3, color = 'tab:green', label =  'Simulation: Schwankung', alpha = 0.6)
# ax1_twin.plot([i+0.5 for i in range(0,7,1)], Top_MAE, color = 'grey', alpha = 1, linestyle = '--', label = 'MAE', marker = 'x', markersize = 5)

# # labels
# ax1_twin.set_ylabel('Schwankung und MAE in K')
# ax1_twin.yaxis.set_major_formatter('{:.1f}'.format)
# ax1_twin.set_ylim(0,30)
# ax1_twin.set_yticks(range(0,30,2))
# ax1_twin.legend(fontsize = 8, loc = 'upper right')
# # ------------------------------------------------X------------------------------------------------
# ax1.set_xlim(0,7)
# ax1.set_xticks(range(0,7,1), ['']*7, rotation = 60)
# # ------------------------------------------------Global------------------------------------------------
# ax1.set_title(r'$\bf{Tagesbilanz\ der\ obersten\ Schicht}$')
# ax1.grid(True)


# # ------------------------------------------------Y2------------------------------------------------
# # Data
# ax2.plot([i+0.5 for i in range(0,7,1)], Bot_Mess_Tmin, color = 'tab:blue', label = '   ', marker = 'o', markersize = 5 ,alpha = 0.5, linestyle = '--')
# ax2.plot([i+0.5 for i in range(0,7,1)], Bot_Mess_Tmax, color = 'tab:red', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )
# ax2.plot([i+0.5 for i in range(0,7,1)], Bot_Mess_Tmean, color = 'tab:green', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )

# ax2.plot([i+0.5 for i in range(0,7,1)], Bot_Sim_Tmin, color = 'tab:blue', alpha = 1, linestyle = '--', label = 'Minimalwert', marker = 'x', markersize = 5)
# ax2.plot([i+0.5 for i in range(0,7,1)], Bot_Sim_Tmax, color = 'tab:red', alpha = 1, linestyle = '--', label = 'Maximalwert', marker = 'x', markersize = 5)
# ax2.plot([i+0.5 for i in range(0,7,1)], Bot_Sim_Tmean, color = 'tab:green', alpha = 1, linestyle = '--', label = 'Mittelwert', marker = 'x', markersize = 5)

# # labels
# ax2.set_ylabel('Temperatur in °C')
# ax2.set_yticks(np.arange(15,60,2.5))
# ax2.set_yticks(np.arange(15,60,0.5), minor = True)
# ax2.set_ylim(15,60)
# ax2_twin.yaxis.set_major_formatter('{:.1f}'.format)
# ax2.legend(fontsize = 8, loc = 'upper left', ncol = 2, title = 'Messwert    Simulation              ', title_fontsize = 8)
# # ------------------------------------------------Y2_Twin------------------------------------------------
# # Data
# ax2_twin.bar([i+0.35 for i in range(0,7,1)], Bot_Mess_Tschwankung, width = 0.3, color = 'tab:green', label = 'Messwert:   Schwankung', alpha = 0.3)
# ax2_twin.bar([i+0.65 for i in range(0,7,1)], Bot_Sim_Tschwankung, width = 0.3, color = 'tab:green', label =  'Simulation: Schwankung', alpha = 0.6)
# ax2_twin.plot([i+0.5 for i in range(0,7,1)], Bot_MAE, color = 'grey', alpha = 1, linestyle = '--', label = 'MAE', marker = 'x', markersize = 5)

# # labels
# ax2_twin.set_ylabel('Schwankung und MAE in K')
# ax2_twin.yaxis.set_major_formatter('{:.1f}'.format)
# ax2_twin.set_ylim(0,30)
# ax2_twin.set_yticks(range(0,30,2))
# ax2_twin.legend(fontsize = 8, loc = 'upper right')
# # ------------------------------------------------X------------------------------------------------
# ax2.set_xlim(0,7)
# ax2.set_xticks(range(0,7,1), ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], rotation = 60)
# # ------------------------------------------------Global------------------------------------------------
# ax2.set_title(r'$\bf{Tagesbilanz\ der\ untersten\ Schicht}$')

# ax2.grid(True)

# E_Top = [b-a for a,b in zip(df['Warm Oben'], [posEl[-1] for posEl in df['WS_theta']])]
# E_Bot = [b-a for a,b in zip(df['Warm Unten'][::4], [posEl[0] for posEl in df['WS_theta'][::4]])]

# plt.show()

# print('')
# print('Sommerfall:')
# print('MAE Top: ' + str(round(sum(Top_MAE)/len(Top_MAE),2)))
# print('ME Top: ' + str(round(sum(E_Top)/len(E_Top),2)))
# print('Top_Max: ' + str(round(max(E_Top),2)))
# print('Top_Min: ' + str(round(min(E_Top),2)))
# print('-')
# print('MAE Bot: ' + str(round(sum(Bot_MAE)/len(Bot_MAE),2)))
# print('ME Bot: ' + str(round(sum(E_Bot)/len(E_Bot),2)))
# print('Bot_Max: ' + str(round(max(E_Bot),2)))
# print('Bot_Min: ' + str(round(min(E_Bot),2)))
# print('--------------------------------------------')














