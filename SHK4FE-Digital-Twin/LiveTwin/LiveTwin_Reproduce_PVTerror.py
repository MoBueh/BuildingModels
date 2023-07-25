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
plt.rcParams["figure.dpi"] = 400
plt.rc('axes', axisbelow=True)

#%% Directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))                            #GoTo aktuelles Arbeitsverzeichnis

#%% Modelle

def WS(Vpkt_prim, T_prim_in, Vpkt_sek, T_sek_in, T_vektor_prev, Tamb):
    # "Development and experimental evaluation of grey-box models of a microscale polygeneration system for application in optimal control" 
    # https://www.sciencedirect.com/science/article/abs/pii/S0378778819331342
    
    #                        _______________
    # Vpkt_prim, T_prim_in  |     T_n+2     | T_sek_out
    # --------------------->|_______________|--------------------->
    #                       |     T_n+1     |
    #                       |_______________|
    #                       |     T_n       |
    #                       |_______________|
    #            T_prim_out |     T_n-1     | Vpkt_sek, T_sek_in
    # <---------------------|_______________|<---------------------   
    
    # -------------------------------Fehlende-Messwerte------------------------------------------------------------
    if str(T_prim_in) == 'nan':
        T_prim_in = 20
    if str(T_sek_in) == 'nan':
        T_sek_in = 20
    if str(Vpkt_prim) == 'nan':
        Ti = T_vektor_prev
        return Ti
    if str(Vpkt_sek) == 'nan':
        Ti = T_vektor_prev
        return Ti
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

def KS(Vpkt_prim, T_prim_in, Vpkt_sek, T_sek_in, T_vektor_prev, Tamb):
    # "Development and experimental evaluation of grey-box models of a microscale polygeneration system for application in optimal control" 
    # https://www.sciencedirect.com/science/article/abs/pii/S0378778819331342
    
    #                        _______________
    # Vpkt_prim, T_prim_in  |     T_n+2     | T_sek_out
    # --------------------->|_______________|--------------------->
    #                       |     T_n+1     |
    #                       |_______________|
    #                       |     T_n       |
    #                       |_______________|
    #            T_prim_out |     T_n-1     | Vpkt_sek, T_sek_in
    # <---------------------|_______________|<---------------------  
    
    # -------------------------------Fehlende-Messwerte------------------------------------------------------------
    if str(T_prim_in) == 'nan':
        T_prim_in = 20
    if str(T_sek_in) == 'nan':
        T_sek_in = 20
    if str(Vpkt_prim) == 'nan':
        Ti = T_vektor_prev
        return Ti
    if str(Vpkt_sek) == 'nan':
        Ti = T_vektor_prev
        return Ti
    # -------------------Parameter-Massenströme-Randbedingungen----------------------------------------------------
    z_i = param['KS']['H']/param['KS']['N']                                     # (25) Höhe einer Schicht
    A_ext = np.pi*param['KS']['D']*z_i                                          # (26) Außenfläche einer Schicht
    A_i = np.pi*(param['KS']['D']-2*param['KS']['th'])**2/4                     # (27) Fläche zwischen zwei Layern
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
                    -(param['KS']['k'] * A_ext * (Ti[0] - Tamb))                # Wärmeverluste an Umgebung
                    + (mpkt_i * prop['wasser']['cp'] * (Ti[0+1] - Ti[0]))       # Stofftransport
                    +((A_i*param['KS']['lambda_eff']/z_i) * (Ti[0+1] - Ti[0]))  # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt

        else:                                                                   # Resultierender Massenstrom ↑    
            dTdt = (0
                    -(mpkt_sek * prop['wasser']['cp'] * (Ti[0] - T_sek_in))     # Sek: mpkt-bedingte Änderung
                    -(param['KS']['k']*A_ext*(Ti[0] - Tamb))                    # Wärmeverluste an Umgebung
                    +0  # keine Schicht unterhalb vorhanden                     # Stofftransport                    
                    +(A_i*param['KS']['lambda_eff']/z_i*(Ti[0+1] - Ti[0]))      # Wärmeleitung zw. Layer 
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
                    -(param['KS']['k']*A_ext*(Ti[-1] - Tamb))                   # Wärmeverluste an Umgebung 
                    +0 # keine Schicht oberhalb vorhanden                       # Stofftransport
                    +((A_i*param['KS']['lambda_eff']/z_i)*(Ti[-1-1] - Ti[-1]))  # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt  
        
        else:                                                                   # Resultierender Massenstrom ↑   
            dTdt = (0
                    +(mpkt_prim * prop['wasser']['cp'] * (T_prim_in - Ti[-1]))  # Prim. mpkt-bedingte Änderung
                    -(param['KS']['k']*A_ext*(Ti[-1] - Tamb))                   # Wärmeverluste an Umgebung 
                    +(mpkt_i * prop['wasser']['cp'] * (Ti[-1] - Ti[-1-1]))      # Stofftransport
                    +(A_i*param['KS']['lambda_eff']/z_i*(Ti[-1-1] - Ti[-1]))     # Wärmeleitung zw. Layer
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
                    -(param['KS']['k']*A_ext*(Ti[i] - Tamb))                    # Wärmeverluste an Umgebung
                    +(mpkt_i*prop['wasser']['cp']*(Ti[i+1] - Ti[i]))            # Stofftransport
                    +((A_i*param['KS']['lambda_eff'])/z_i * (Ti[i+1] - 2*Ti[i] + Ti[i-1])) # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt
   
        else:                                                                   # Resultierender Massenstrom ↑ 
            dTdt = (0
                    +0 # kein Anschluss in diesen Layern                        # Sek/prim mpkt-bedingte Änderung
                    -(param['KS']['k']*A_ext*(Ti[i] - Tamb))                    # Wärmeverluste an Umgebung
                    +(mpkt_i*prop['wasser']['cp']*(Ti[i] - Ti[i-1]))            # Stofftransport
                    +((A_i*param['KS']['lambda_eff'])/z_i * (Ti[i+1] - 2*Ti[i] + Ti[i-1])) # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt 
        return dTdt
    # -------------------------------Solving-DGL-------------------------------------------------------------------
    # for i in range(1, param['KS']['N'] - 1):
    #     Ti[i] = float(odeint(EB_MIDDLE, Ti[i], np.linspace(0, 900))[-1])
    if mpkt_i>=0:                                                               # Lösen der Schichtung in Richtung mpkt_i ↓
        for i in list(reversed(range(1, param['KS']['N'] - 1))):
            Ti[i] = float(odeint(EB_MIDDLE, Ti[i], np.linspace(0, 900))[-1])
    else:                                                                       # Lösen der Schichtung in Richtung mpkt_i ↑
        for i in list(range(1, param['KS']['N'] - 1)):
            Ti[i] = float(odeint(EB_MIDDLE, Ti[i], np.linspace(0, 900))[-1])
    # -------------------------------------------------------------------------------------------------------------    

    # ---------------------------Return-Temperatures---------------------------------------------------------------
    return Ti

def TABS(TABS_theta_in, TABS_Vpkt, B_theta_air):
    #                        _______________
    # theta_in      [°C]    |               | theta_out [°C]
    # --------------------->|               |--------------------->
    # Vpkt          [m³/h]  |               | 
    # --------------------->|     TABS      |
    # B_theta_air   [°C]    |               | Qpkt      [W]
    # --------------------->|               |--------------------->
    #                       |_______________|
    if TABS_Vpkt > 0:
        R_w = (param['TABS']['dx']**0.13*(param['TABS']['da']-2*param['TABS']['dp'])/(TABS_Vpkt*prop['wasser']['rho']/3600*param['TABS']['A']/param['TABS']['dx'])**0.87)/(8*math.pi)
        R_p = (param['TABS']['dx']*math.log(param['TABS']['da']/(param['TABS']['da']-2*param['TABS']['dp'])))/(2*math.pi*param['TABS']['lambda_p'])
        R_x = 1**param['TABS']['dx'] * math.log(param['TABS']['dx']/(math.pi*param['TABS']['da']))/(2*math.pi*param['TABS']['lambda_s'])
        R_z = 1/(2*(TABS_Vpkt*prop['wasser']['rho']/3600)*prop['wasser']['cp'])
        R_d = 1/(1/(1/8 + param['TABS']['d']/param['TABS']['lambda_s']))
        k = 1/(R_w + R_p + R_x + R_z + R_d)
        TABS_theta_out= ((TABS_theta_in - B_theta_air)/math.exp(k*param['TABS']['A']/prop['wasser']['rho']/(TABS_Vpkt/3600)/prop['wasser']['cp']))+B_theta_air
        TABS_Qpkt=prop['wasser']['rho']*TABS_Vpkt/3600*prop['wasser']['cp']* (TABS_theta_in-TABS_theta_out)
    else:
        TABS_Qpkt=np.nan
        TABS_theta_out=np.nan
    TABS_theta_in = TABS_theta_in
    TABS_Vpkt = TABS_Vpkt
    return TABS_theta_in, TABS_theta_out, TABS_Vpkt, TABS_Qpkt

def BAT(BAT_Pel_in, BAT_Pel_out, BAT_SOC_prev):
    #                        _______________
    # Pel_in      [W]       |               |
    # --------------------->|               |
    # Pel_out     [W]       |               | SOC
    # --------------------->|      BAT      |--------------------->
    # SOC[-1]     [%]       |               | 
    # --------------------->|               |
    #                       |_______________|
    BAT_Pel_in, BAT_Pel_out = (0 if math.isnan(val) else val for val in (BAT_Pel_in, BAT_Pel_out))
    C_prev = BAT_SOC_prev*param['BAT']['C']*3600
    C = (BAT_Pel_in-BAT_Pel_out)*900 + C_prev
    BAT_SOC = C/(param['BAT']['C']*3600)
    if BAT_SOC > 1:
        BAT_SOC = 1
    if BAT_SOC < 0:
        BAT_SOC = 0
    return BAT_SOC

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
    if str(PVT_Vpkt) == 'nan':
        PVT_Pel, PVT_theta_in, PVT_theta_out, PVT_Qpkt, PVT_theta_cell = 0, theta_mean_prev, theta_mean_prev, 0, theta_mean_prev
        return PVT_Pel,PVT_theta_in,PVT_theta_out,PVT_Qpkt,PVT_theta_cell                                        
    A = param['PVT']['etha_opt']*G                                               
    D = 0
    F = (param['PVT']['c1']+D+param['PVT']['c5']/900)*0.5
    E = PVT_Vpkt/3600*prop['Antifrogen']['rho']*prop['Antifrogen']['cp'] / (param['PVT']['A'] * param['PVT']['quantity'])
    if (PVT_Vpkt == 0 or str(PVT_theta_in) == 'nan'):
        PVT_theta_in = theta_mean_prev
    PVT_theta_out = max([(A+theta_e*(D+ param['PVT']['c1']) + param['PVT']['c5']/900 * theta_mean_prev - PVT_theta_in*(F-E))/(E+F), [(PVT_theta_in-2000/(prop['Antifrogen']['rho']*prop['Antifrogen']['cp']*PVT_Vpkt/3600)) if PVT_Vpkt != 0 else np.nan][0]])
    PVT_Qpkt = PVT_Vpkt/3600 *prop['Antifrogen']['rho'] * prop['Antifrogen']['cp'] * (PVT_theta_out-PVT_theta_in)                    
    PVT_theta_cell = PVT_Qpkt / (param['PVT']['A']*param['PVT']['quantity'])/param['PVT']['U'] + (PVT_theta_in+PVT_theta_out)/2
    if BAT_SOC == 1:
        PVT_Pel = 0    
    else: 
        PVT_Pel = G*param['PVT']['etha_el']*(1-param['PVT']['beta_el']*(PVT_theta_cell-25))*(param['PVT']['A']*param['PVT']['quantity'])   
    if PVT_Pel >= 1040:
        PVT_Pel = 1040
    return PVT_Pel,PVT_theta_in,PVT_theta_out,PVT_Qpkt,PVT_theta_cell

def WP(WP_prim_theta_in, WP_sek_theta_in, WP_sek_delta_theta, WP_prim_Vpkt, WP_sek_Vpkt, WP_OI):
    #                        _______________
    # Prim_theta_in [°C]    |               | COP            [-]
    # --------------------->|               |--------------------->
    # Prim_theta_in [°C]    |               | Pel            [W]
    # --------------------->|               |--------------------->
    # ΔSek_theta    [°C]    |               | Qpkt_prim      [W]
    # --------------------->|               |--------------------->
    # Prim_Vpkt     [m³/h]  |     WP        | Prim_theta_out [°C]
    # --------------------->|               |--------------------->
    # Sek_theta_in  [°C]    |               | Qpkt_sek       [W] 
    # --------------------->|               |--------------------->
    # Sek_Vpkt      [°C]    |               | PSek_theta_out [°C]
    # --------------------->|_______________|--------------------->
    if str(WP_sek_Vpkt) == 'nan':
        WP_OI = 0
    
    if WP_OI == 0:
        WP_COP = np.nan
        WP_Pel = np.nan
        WP_Qdot_prim = np.nan
        WP_Qdot_sec = np.nan
        WP_prim_theta_out = np.nan 
        WP_sek_theta_out = np.nan
        WP_prim_Vpkt = np.nan
        WP_sek_Vpkt = np.nan
    else:
        WP_prim_Vpkt = WP_prim_Vpkt
        WP_sek_Vpkt = WP_sek_Vpkt  
        WP_sek_theta_out = WP_sek_theta_in + WP_sek_delta_theta
        WP_COP = param['WP']['etha'] * (WP_sek_theta_out+273.15)/(WP_sek_theta_out - WP_prim_theta_in)
        WP_Qdot_sec = prop['wasser']['rho']*WP_sek_Vpkt/3600* prop['wasser']['cp']*(WP_sek_theta_out-WP_sek_theta_in)
        WP_Pel = (WP_Qdot_sec/WP_COP)
        WP_Qdot_prim = (WP_Pel-WP_Qdot_sec)
        WP_prim_theta_out = WP_prim_theta_in-(-WP_Qdot_prim/(prop['Antifrogen']['rho']*(WP_prim_Vpkt/3600)*prop['Antifrogen']['cp']))   
    WP_prim_theta_in = WP_prim_theta_in
    WP_sek_theta_in = WP_sek_theta_in
    WP_prim_Vpkt = WP_prim_Vpkt
    WP_sek_Vpkt = WP_sek_Vpkt
    return WP_COP, WP_Pel, WP_Qdot_prim, WP_Qdot_sec, WP_prim_theta_in, WP_prim_theta_out, WP_sek_theta_in, WP_sek_theta_out, WP_prim_Vpkt, WP_sek_Vpkt

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
    # Sonnenstand     [°]   |     5R1C       | T_mass        [°C]
    # --------------------->|               |--------------------->
    # Solarstrahlung  [°C]  |               | T_operativ    [°C]
    # --------------------->|               |--------------------->
    # Vpkt_Lüftung    [°C]  |               |
    # --------------------->|               |
    # T_Lüftung       [°C]  |               | 
    # --------------------->|               |
    #                       |_______________|  
    if str(I_dir) == 'nan':
        B_theta_m, B_theta_s, B_theta_air, B_theta_operativ = [B_theta_m_prev]*4
        return B_theta_m, B_theta_s, B_theta_air, B_theta_operativ
    if str(B_phi_HC) == 'nan':
        B_phi_HC = 0
    else:
        B_phi_HC = B_phi_HC
    A_m = param['B']['A_f']*param['B']['Coe_Am']
    A_tot = 4.5 * param['B']['A_f']
    c_m = param['B']['Coe_Cm']*param['B']['A_f']*2.5
    H_tr_w = sum([U*A for U,A in zip(param['B']['transparent components']['U'], param['B']['transparent components']['A'])])
    if V_Vdot == 0:
        H_ve = 0.00000001
    else:
        # H_ve = 1200*V_Vdot/3600                                                 # keine Wärmerückgewinnung
        H_ve = (1-0.85)*1200*V_Vdot/3600                                        # mit Wärmerückgewinnung
    H_tr_ms = 9.1*A_m
    H_op = sum([U*A for U,A in zip(param['B']['opaque components']['U'], param['B']['opaque components']['A'])])
    H_tr_em = 1/(1/H_op+1/H_tr_ms)
    H_tr_is = 3.45*A_tot
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
        I_umg = (I_diff + I_dir) * 0.5 * 0.2 * (1-math.cos(math.radians(tilt)))
        I_d = I_diff * 0.5 * (1 + math.cos(math.radians(tilt)))
        I_sol = I_f +I_umg+I_d
        phi_sol_trans.append(shading*g_tot*F_c*0.9*(1-F_f)*A*I_sol)
# Phi_sol: opaque components
    phi_sol_op = []
    for shading, R, U, A, direction, tilt in zip(param['B']['opaque components']['shading'],
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
    phi_ia = 0.5*B_phi_int+0*B_phi_HC                                                      # C.1
    phi_m = A_m/A_tot * (0.5 * B_phi_int + phi_sol)                     + 1*B_phi_HC                          # C.2
    phi_st = (1-A_m/A_tot-H_tr_w/9.1/A_tot)*(0.5 * B_phi_int + phi_sol) + 0*B_phi_HC          # C.3
    phi_mtot = phi_m + H_tr_em * theta_e + H_tr_3 *(phi_st +H_tr_w * theta_e + H_tr_1 * ((phi_ia)/H_ve + V_theta_sup))/ H_tr_2   # C.5
    B_theta_m_t = (B_theta_m_prev * (c_m/3600-0.5 * (H_tr_3 +H_tr_em))+phi_mtot)/(c_m/3600+0.5*(H_tr_3+H_tr_em))  # C.4
        # Output
    B_theta_m = (B_theta_m_t + B_theta_m_prev)/2                                      # C.9
    # B_theta_m = B_theta_m_t                                     # C.9
    B_theta_s = (H_tr_ms * B_theta_m + phi_st + H_tr_w * theta_e + H_tr_1 *(V_theta_sup + (phi_ia)/H_ve))/(H_tr_ms + H_tr_w + H_tr_1)     # C.10
    B_theta_air = (H_tr_is * B_theta_s + H_ve * V_theta_sup + phi_ia)/(H_tr_is + H_ve)          # C.11
    B_theta_operativ = 0.3 * B_theta_air + 0.7 * B_theta_s                            # C.12
    return B_theta_m, B_theta_s, B_theta_air, B_theta_operativ, phi_sol_op, phi_sol_trans, phi_sol

def SunPos(lon, lat, year, month, day, hour, minute, shift):
    #                        _______________
    # Zeit                  |               | Sonnenhöhe
    # --------------------->|               |--------------------->
    # Position              |               | Sonnenazimuth
    # --------------------->|_______________|--------------------->
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
    'lon' : 7.9498017,
    'lat' : 48.473451,
    }

param['WS'] = {
    'H'             : 1.2,                                                      # m         Speicherhöhe
    'N'             : 4,                                                        # -         Anzahl Schichten
    'D'             : 0.6,                                                      # m         Speicherdurchmesser
    'th'            : 0.02,                                                     # m²        Wanddicke
    'k'             : 0.02*100,                                                 # W/m²K     Wärmetransferkoeffizient
    'lambda_eff'    : 0.0015*1000                                               # W/mK      vertikale Wärmeleitfähigkeit von Wasser   
    }

param['KS'] = {
    'H'             : 1.2,                                                      # m         Speicherhöhe
    'N'             : 8,                                                        # -         Anzahl Schichten
    'D'             : 0.6,                                                      # m         Speicherdurchmesser
    'th'            : 0.02,                                                     # m²        Wanddicke
    'k'             : 0.02*100,                                                 # W/m²K     Wärmetransferkoeffizient
    'lambda_eff'    : 0.0015*1000                                               # W/mK      vertikale Wärmeleitfähigkeit von Wasser     
    }

param['TABS'] = {
    'A'         : 15,                                                           # Area of floor
    'da'        : 0.03,                                                         # outside diameter of pipe
    'dp'        : 0.001,                                                        # thickness of pipe
    'dx'        : 0.125,                                                        # Laying distance of pipes
    'lambda_p'  : 340,                                                          # Heat transfer pipe (copper)
    'lambda_s'  : 2,                                                            # Heat transfer floor
    'd'         : 0.01,                                                         # Laying depth
    }

param['BAT'] = {
    'C' : 8000                                                                  # Wh
    }

param['PVT'] = {
    'quantity'  : 3,                                                            # -         number of modules    
    'A'         : 1.8,                                                          # m²        Area of one module        
    'etha_opt'  : 0.234,                                                        # -         optical efficiency
    'U'         : 40,                                                           # W/m²/K    Coupling between T_mass and T_Cell
    'etha_el'   : 0.175,                                                        # -         electrical efficiency
    'beta_el'   : 0.0039,                                                       # -         correction for el. efficiency
    'c1'        : 22.8,                                                         # W/m²/K    Heat loss coefficient 
    'c5'        : 26050,                                                        # J/m²K     Effective heat capacity
    }

param['WP'] = {
    'etha'  :   0.5,                                                           # #-        Exergetischer Wirkungsgrad
    }

param['B'] = {
    'A_f': 14.8,                                                                # m²        conditioned area
    'Coe_Am': 2.5,                                                              # -         coefficient for the determination of mass-related area according to table 12
    'Coe_Cm': 120000,                                                           # J/K*m²    coefficiant for the determination of Internal heat storage capacity according to table 12
    'transparent components': {
        # each column is one Component
        'shading'    : [1,      1,      1,      1,      0.1,    0],             # -         1 for no shading, 0 for absolutely shadowed
        'g_tot'      : [0.4,    0.4,    0.4,    0.4,    0.4,    0],             # -         DIN 4108-2 8.3.2 (3) und DIN EN ISO 13790 Tabelle 7   
        'F_c'        : [1,      1,      1,      1,      1,      1],             # -         DIN 4108-2 8.3.2 (3) und DIN EN ISO 13790 Tabelle 7
        'F_f'        : [0,      0,      0,      0,      0,      0],             # -         Frame part of the component; 1 -> 100 % frame, 0 -> 0 % frame
        'A'          : [0,      0,      0,      0,      0,      0],             # m²        Area of the component
        'direction'  : [66,     156,    246,    336,    270,    270],           # °         0° -> north, 90° -> east, 180° ->south, 270° -> west
        'tilt'       : [90,     90,     90,     90,     0,      180],           # °         0° horizontal top, 90° -> vertical, 180° -> horizontal down
        'U'          : [0.8,    0.8,    0.8,    0.8,    0.8,    0.8]            # W/m²/K    Heat transfer coefficient determined according to ISO 6946
        },
    'opaque components': {
        # each column is one Component
        'shading'    : [1,      1,      1,      1,      0.1,    1],             # -         1 for no shading, 0 for absolutely shadowed
        'R'          : [0.05,   0.05,   0.05,   0.05,   0.05,   0.05],          #           product of surface heat transfer resistance(ISO 6946) and absorption coefficient for solar radiation on opaque component, "Baehr, H.D.; Stephan, K. ; Wärme- und Stoffübertragung; Auszug Kap. 5.5 Strahlungsaustausch, aus Tab. 5.8,S.633" https://www.baunetzwissen.de/bauphysik/fachwissen/waermeschutz/absorption-und-waerme-auf-oberflaechen-4733315
        'U'          : [0.6,    0.6,    0.6,    0.6,    0.6,    0.6],           # W/m²/K    Heat transfer coefficient determined according to ISO 6946
        'A'          : [5.3,    14.6,   6.3,    15.7,   14.8,   14.8],          # m²        Area of the component
        'direction'  : [66,     156,    246,    336,    270,    270],           # °         0° -> north, 90° -> east, 180° ->south, 270° -> west
        'tilt'       : [90,     90,     90,     90,     0,      180]            # °         0° horizontal top, 90° -> vertical, 180° -> horizontal down
        },
    }

#%% DataPrep
a = 0   # Welche Datei ist Startdatei?
for file in os.listdir(os.getcwd()+'\Data')[a:7]:
    a = a+1
    df = pd.read_csv('Data/' + file, index_col = 0, delimiter = ';', header = 1)     # Messwerte von Mondas
    df.index = pd.to_datetime(df.index, format='%Y/%m/%d %H:%M:%S')                         # index als DateTimeObject
    df.replace(-99, np.nan, inplace=True)                                                   # -99 (Mondas) zu nan
    df = df.drop('GEN.EL-PV----CALC.POW.EL-', axis=1)                                       # Dopplung!!
    from Rename_Dictionary import rename
    df = df.rename(columns=rename)                                                          # rename Columns
    del rename
    #-----Pyranometer-----
    df.loc[(df.index.hour < 7) | (df.index.hour > 23), ['Pyr_hor', 'Pyr_ver']] = 0 

    #-----Wärmepumpe-----
    df.loc[df['WP_COP'] < 3, ['WP_COP', 'WP_Pel', 'WP_Qpkt']] = 0
    df.loc[df['WP_COP'] == 0, ['WP_prim_T_RET', 'WP_prim_T_SUP', 'WP_sek_T_RET', 'WP_sek_T_SUP', 'WP_Pel', 'WP_Qpkt']] = np.nan
    df.loc[df['WP_prim_T_SUP'] > 17, ['WP_prim_T_SUP', 'WP_prim_T_RET']] -= 2
    
    #-----HM10:WP-----
    df.loc[(df['HM10_Vpkt'] < 80) | (df['WP_COP'] < 3), 'HM10_Vpkt'] = 0                        
    df.loc[df['HM10_Vpkt'] == 0, ['HM10_T_RET', 'HM10_T_SUP']] = np.nan
    
    #-----HM1:PVT-----
    df.loc[(df['HM1_Vpkt'] < 200) | (df['HM1_T_RET'] <= 31), 'HM1_Vpkt'] = 0
    df.loc[df['HM1_Vpkt'] == 0, ['HM1_T_RET', 'HM1_T_SUP']] = np.nan
    df['HM1_Qpkt'] = [math.nan if x > 0 else x for x in [vpkt/1000/3600 * prop['Antifrogen']['rho'] *prop['Antifrogen']['cp'] * (SUP-RET) for vpkt, SUP, RET in zip(df['HM1_Vpkt'], df['HM1_T_SUP'],df['HM1_T_RET'])]]
    
    #-----HM4:TABS-----    
    df.loc[(df['HM4_Vpkt'] < 40) | (df['HM4_T_RET'] > 20), 'HM4_Vpkt'] = 0
    df.loc[df['HM4_Vpkt'] == 0, ['HM4_T_RET', 'HM4_T_SUP']] = np.nan
    df['HM4_Qpkt'] = [vpkt/1000/3600 * prop['wasser']['rho'] *prop['wasser']['cp'] * (SUP-RET) for vpkt, SUP, RET in zip(df['HM4_Vpkt'], df['HM4_T_SUP'],df['HM4_T_RET'])]
    
    shift = 3
    df['HM4_B'] = [vpkt/1000/3600 * prop['wasser']['rho'] *prop['wasser']['cp'] * (SUP-RET) for vpkt, SUP, RET in zip(df['HM4_Vpkt'].shift(shift), df['HM4_T_SUP'].shift(shift),df['HM4_T_RET'].shift(shift))]
    df.loc[df['HM4_B'] < -900, 'HM4_B'] = -900


    #----------------------------------------Startwerte------------------------------------------------
    WS_theta = [np.linspace(df['Warm Unten'][0], df['Warm Oben'][0], param['WS']['N']).tolist()]
    KS_theta = [np.linspace(df['Kalt Unten'][0], df['Kalt Oben'][0], param['KS']['N']).tolist()]
    TABS_theta_in = []
    TABS_theta_out = []
    TABS_Vpkt = []
    TABS_Qpkt = []
    BAT_SOC = [0.4348*df.loc[df['WP_COP'] == 0,  'BAT_U'][0] -24.088]
    PVT_Pel = [0]
    PVT_theta_in = [20]
    PVT_theta_out = [20]
    PVT_Qpkt = [0]
    PVT_theta_cell = [20]
    WP_COP = []
    WP_Pel = []
    WP_Qdot_prim = []
    WP_Qdot_sec = []
    WP_prim_theta_in = []
    WP_prim_theta_out = []
    WP_sek_theta_in = []
    WP_sek_theta_out = []
    WP_prim_Vpkt = []
    WP_sek_Vpkt = []
    B_theta_m = [df.loc[df.index[0], 'Zone_T']]
    B_theta_s = [20] 
    B_theta_air = [20] 
    B_theta_operativ = [20] 
    sun_h = []
    sun_az = []
    i_dif = []
    i_dir = []
    
    for i in df.index:
        # print(i.strftime("%H:%M"))
        #----------------------------------------Sonnenposition----------------------------------------
        sp = SunPos(param['Position']['lon'], 
                    param['Position']['lat'], 
                    year = i.year, 
                    month = i.month, 
                    day = i.day, 
                    hour = i.hour,
                    minute = i.minute,
                    shift = -1)
        sun_az.append(sp[0])
        sun_h.append(sp[1])
        
        #----------------------------------------Solarstrahlung----------------------------------------
        i_dif.append((1/(1 + 2.71828182846**(7.997*(df.loc[i, 'Pyr_hor']/1000 - 0.586))))*df.loc[i, 'Pyr_hor'])
        i_dir.append(df.loc[i, 'Pyr_hor']-i_dif[-1])
        
        #----------------------------------------Wärmespeicher----------------------------------------
        ws = WS(Vpkt_prim       = df.loc[i, 'HM10_Vpkt']/1000*0.8, 
                T_prim_in       = df.loc[i, 'WP_sek_T_SUP']-1, 
                Vpkt_sek        = df.loc[i, 'HM1_Vpkt']/1000*1.5, 
                T_sek_in        = df.loc[i, 'HM1_T_RET'], 
                T_vektor_prev   = WS_theta[-1],
                Tamb            = df.loc[i, 'Zone_T'])
        WS_theta.append(ws)
        if i.hour == 0 and i.minute == 0 and i.second == 0:
            WS_theta[-1] = np.linspace(df['Warm Unten'][i], df['Warm Oben'][i], param['WS']['N']).tolist() 
            
        #----------------------------------------Kältespeicher----------------------------------------
        ks = KS(Vpkt_prim       = df.loc[i, 'HM4_Vpkt']/1000, 
                T_prim_in       = df.loc[i, 'HM4_T_RET'], 
                Vpkt_sek        = 715/1000*[1 if df.loc[i, 'WP_COP']>0 else 0][0], 
                T_sek_in        = df.loc[i, 'WP_prim_T_RET']+3,                     # Anpassung durch IWT
                T_vektor_prev   = KS_theta[-1],
                Tamb            = df.loc[i, 'Zone_T'])
        KS_theta.append(ks)
        if i.hour == 0 and i.minute == 0 and i.second == 0:
            KS_theta[-1] = np.linspace(df['Kalt Unten'][i], df['Kalt Oben'][i], param['KS']['N']).tolist()
            
        #-------------------------------------------TABS----------------------------------------------
        tabs = TABS(TABS_theta_in = df.loc[i, 'HM4_T_SUP'], 
                    TABS_Vpkt = df.loc[i, 'HM4_Vpkt']/1000,
                    B_theta_air = df.loc[i, 'Zone_T'])
        TABS_theta_in.append(tabs[0])
        TABS_theta_out.append(tabs[1])
        TABS_Vpkt.append(tabs[2])
        TABS_Qpkt.append(tabs[3])
        
        #-------------------------------------------BAT-----------------------------------------------
        bat = BAT(BAT_Pel_in = (df.loc[i, 'PVT_Pel'] + df.loc[i, 'INV_P_IN'])*1000, 
                  BAT_Pel_out = df.loc[i, 'Bat_Pel']*1000, 
                  BAT_SOC_prev = BAT_SOC[-1])
        BAT_SOC.append(bat)
        
        #-------------------------------------------PVT-----------------------------------------------
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
        
        #-------------------------------------------WP-----------------------------------------------
        wp = WP(WP_prim_theta_in = df.loc[i, 'WP_prim_T_SUP'], 
                WP_sek_theta_in = df.loc[i, 'WP_sek_T_RET'], 
                WP_sek_delta_theta = 6, 
                WP_prim_Vpkt = 715/1000, 
                WP_sek_Vpkt = df.loc[i, 'HM10_Vpkt']/1000, 
                WP_OI = [1 if df.loc[i, 'WP_COP']>0 else 0][0])
        WP_COP.append(wp[0])
        WP_Pel.append(wp[1])
        WP_Qdot_prim.append(wp[2])
        WP_Qdot_sec.append(wp[3])
        WP_prim_theta_in.append(wp[4])
        WP_prim_theta_out.append(wp[5])
        WP_sek_theta_in.append(wp[6])
        WP_sek_theta_out.append(wp[7])
        WP_prim_Vpkt.append(wp[8])
        WP_sek_Vpkt.append(wp[9])

        #-------------------------------------------B-----------------------------------------------
        b = B(B_phi_int = 90,
              B_phi_HC = df.loc[i, 'HM4_B'], 
              B_theta_m_prev = B_theta_m[-1], 
              theta_e = df.loc[i, 'AMB_T'], 
              sun_h = sun_h[-1], 
              sun_az = sun_az[-1], 
              I_dir = i_dir[-1], 
              I_diff = i_dif[-1], 
              V_Vdot = 30/4, 
              V_theta_sup = df.loc[i, 'AMB_T'] )

        # if i.hour in [12, 0] and i.minute == 0 and i.second == 0:
        if i.hour in [0] and i.minute == 0 and i.second == 0:
            B_theta_m.append(df.loc[i, 'Zone_T'])
        else:
            B_theta_m.append(b[0])
        B_theta_s.append(b[1])
        B_theta_air.append(b[2])
        B_theta_operativ.append(b[3])
               
    df['WS_theta'] = WS_theta[1:]
    df['KS_theta'] = KS_theta[1:]
    df['TABS_theta_in'] = TABS_theta_in
    df['TABS_theta_out'] = TABS_theta_out
    df['TABS_Vpkt'] = TABS_Vpkt
    df['TABS_Qpkt'] = TABS_Qpkt
    df['BAT_SOC'] = BAT_SOC[1:]
    df['s_PVT_Pel'] = PVT_Pel[1:]
    df['PVT_theta_in'] = PVT_theta_in[1:]
    df['PVT_theta_out'] = PVT_theta_out[1:]
    df['PVT_Qpkt'] = PVT_Qpkt[1:]
    df['PVT_theta_cell'] = PVT_theta_cell[1:]
    df['s_WP_COP'] = WP_COP
    df['s_WP_Pel'] = WP_Pel
    df['s_WP_Qdot_prim'] = WP_Qdot_prim
    df['s_WP_Qdot_sec'] = WP_Qdot_sec
    df['s_WP_prim_theta_in'] = WP_prim_theta_in
    df['s_WP_prim_theta_out'] = WP_prim_theta_out
    df['s_WP_sek_theta_in'] = WP_sek_theta_in
    df['s_WP_sek_theta_out'] = WP_sek_theta_out
    df['s_WP_prim_Vpkt'] = WP_prim_Vpkt
    df['s_WP_sek_Vpkt'] = WP_sek_Vpkt
    df['sun_h'] = sun_h
    df['sun_az'] = sun_az  
    df['i_dif'] = i_dif
    df['i_dir'] = i_dir
    df['s_B_theta_m'] = B_theta_m[1:]
    df['s_B_theta_s'] = B_theta_s[1:]
    df['s_B_theta_air'] = B_theta_air[1:]
    df['s_B_theta_operativ'] = B_theta_operativ[1:]
    
    df.to_csv(os.path.dirname(os.path.abspath(__file__))+'/Simulationsergebnisse/Sim_' + file.split('.')[0] + '.csv')
   
    font = 14
    font_title = 12
    labelfont = 9
    legendfont = 9
    tickfont = 9
    
    fig, ax32 = plt.subplots(figsize=(6.4, 4.8/1.3))
    ax32_twin = ax32.twinx()
    # Data
    ax32.scatter(df.index, [i*1000 for i in df['PVT_Pel']], label = 'Mess', c = 'tab:red', s = 12, alpha = 0.5)
    ax32.plot(df.index, [i*1000 for i in df['PVT_Pel']], c = 'tab:red', linestyle = '--', alpha = 0.5)
    ax32.scatter(df.index, df['s_PVT_Pel'], label = 'Sim', c = 'tab:red', s = 12, marker = 'x')
    ax32.plot(df.index, df['s_PVT_Pel'], c = 'tab:red', linestyle = '--')
    # labels
    ax32.set_ylabel('elektrische Leistung in kW', fontsize = labelfont)
    ax32.set_ylim(0,1200)
    ax32.set_yticks(np.arange(0,1200,100), minor = True)
    ax32.set_yticks(np.arange(0,1201,200), [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    # ------------------------------------------------Y2------------------------------------------------
    # Data
    ax32_twin.plot(df.index, [i/1000 for i in df['Pyr_hor']], c = 'tab:orange', alpha = 0.7, label = 'Globalstrahlung')
    # # labels
    ax32_twin.set_ylabel('Globalstrahlung in kW/m²', fontsize = labelfont)
    # ax2.yaxis.set_major_formatter('{:.0f}'.format)
    ax32_twin.set_yticks(np.arange(0,1.200, 0.100), minor = True)
    ax32_twin.set_yticks(np.arange(0,1.300,0.200))
    ax32_twin.set_ylim(0,1.200)
    # ------------------------------------------------X------------------------------------------------
    ax32.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], ['']*84, rotation = 90)
    ax32.set_xlim(df.index[0], df.index[0]+timedelta(hours = 24))
    # ------------------------------------------------Global------------------------------------------------
    # ax32.set_title(r'$\bf{PVT-Kollektoren:\ elektrische\ Seite}$', fontsize = font_title)    
    ax32.grid(True)
    ax32.grid(which='minor', alpha = 0.3)        
    ax32.legend(fontsize = legendfont, loc = 'upper left', ncol = 2)
    ax32_twin.legend(fontsize = legendfont, loc = 'upper right', ncol = 2)
    plt.show()
    
    
    if df.index[0].weekday() == 6:
        print('week')
        df = pd.DataFrame()                     # index als DateTimeObject
        for weeks in os.listdir(os.getcwd()+'\Data')[a-7:a]:
            print(weeks)
            zwischen = pd.read_csv('Data/' + weeks, index_col = 0, delimiter = ';', header = 1)
            zwischen.index = pd.to_datetime(zwischen.index, format='%Y/%m/%d %H:%M:%S')
            df = pd.concat([df, zwischen], axis = 0)
        df.replace(-99, np.nan, inplace=True)                                                   # -99 (Mondas) zu nan
        df = df.drop('GEN.EL-PV----CALC.POW.EL-', axis=1)                                       # Dopplung!!
        from Rename_Dictionary import rename
        df = df.rename(columns=rename)                                                          # rename Columns
        del rename
        #-----Pyranometer-----
        df.loc[(df.index.hour < 7) | (df.index.hour > 23), ['Pyr_hor', 'Pyr_ver']] = 0 
    
        #-----Wärmepumpe-----
        df.loc[df['WP_COP'] < 3, ['WP_COP', 'WP_Pel', 'WP_Qpkt']] = 0
        df.loc[df['WP_COP'] == 0, ['WP_prim_T_RET', 'WP_prim_T_SUP', 'WP_sek_T_RET', 'WP_sek_T_SUP', 'WP_Pel', 'WP_Qpkt']] = np.nan
        df.loc[df['WP_prim_T_SUP'] > 17, ['WP_prim_T_SUP', 'WP_prim_T_RET']] -= 2
        
        #-----HM10:WP-----
        df.loc[(df['HM10_Vpkt'] < 80) | (df['WP_COP'] < 3), 'HM10_Vpkt'] = 0                        
        df.loc[df['HM10_Vpkt'] == 0, ['HM10_T_RET', 'HM10_T_SUP']] = np.nan
        
        #-----HM1:PVT-----
        df.loc[(df['HM1_Vpkt'] < 200) | (df['HM1_T_RET'] <= 31), 'HM1_Vpkt'] = 0
        df.loc[df['HM1_Vpkt'] == 0, ['HM1_T_RET', 'HM1_T_SUP']] = np.nan
        df['HM1_Qpkt'] = [math.nan if x > 0 else x for x in [vpkt/1000/3600 * prop['Antifrogen']['rho'] *prop['Antifrogen']['cp'] * (SUP-RET) for vpkt, SUP, RET in zip(df['HM1_Vpkt'], df['HM1_T_SUP'],df['HM1_T_RET'])]]
        
        #-----HM4:TABS-----    
        df.loc[(df['HM4_Vpkt'] < 40) | (df['HM4_T_RET'] > 20), 'HM4_Vpkt'] = 0
        df.loc[df['HM4_Vpkt'] == 0, ['HM4_T_RET', 'HM4_T_SUP']] = np.nan
        df['HM4_Qpkt'] = [vpkt/1000/3600 * prop['wasser']['rho'] *prop['wasser']['cp'] * (SUP-RET) for vpkt, SUP, RET in zip(df['HM4_Vpkt'], df['HM4_T_SUP'],df['HM4_T_RET'])]
        
        shift = 3
        df['HM4_B'] = [vpkt/1000/3600 * prop['wasser']['rho'] *prop['wasser']['cp'] * (SUP-RET) for vpkt, SUP, RET in zip(df['HM4_Vpkt'].shift(shift), df['HM4_T_SUP'].shift(shift),df['HM4_T_RET'].shift(shift))]
        df.loc[df['HM4_B'] < -900, 'HM4_B'] = -900
        
        #----------------------------------------Startwerte------------------------------------------------
        WS_theta = [np.linspace(df['Warm Unten'][0], df['Warm Oben'][0], param['WS']['N']).tolist()]
        KS_theta = [np.linspace(df['Kalt Unten'][0], df['Kalt Oben'][0], param['KS']['N']).tolist()]
        TABS_theta_in = []
        TABS_theta_out = []
        TABS_Vpkt = []
        TABS_Qpkt = []
        BAT_SOC = [0.4348*df.loc[df['WP_COP'] == 0,  'BAT_U'][0] -24.088]
        PVT_Pel = [0]
        PVT_theta_in = [20]
        PVT_theta_out = [20]
        PVT_Qpkt = [0]
        PVT_theta_cell = [20]
        WP_COP = []
        WP_Pel = []
        WP_Qdot_prim = []
        WP_Qdot_sec = []
        WP_prim_theta_in = []
        WP_prim_theta_out = []
        WP_sek_theta_in = []
        WP_sek_theta_out = []
        WP_prim_Vpkt = []
        WP_sek_Vpkt = []
        B_theta_m = [df.loc[df.index[0], 'Zone_T']]
        B_theta_s = [20] 
        B_theta_air = [20] 
        B_theta_operativ = [20] 
        sun_h = []
        sun_az = []
        i_dif = []
        i_dir = []
        
        for i in df.index:
            # print(i.strftime("%H:%M"))
            #----------------------------------------Sonnenposition----------------------------------------
            sp = SunPos(param['Position']['lon'], 
                        param['Position']['lat'], 
                        year = i.year, 
                        month = i.month, 
                        day = i.day, 
                        hour = i.hour,
                        minute = i.minute,
                        shift = -1)
            sun_az.append(sp[0])
            sun_h.append(sp[1])
            
            #----------------------------------------Solarstrahlung----------------------------------------
            i_dif.append((1/(1 + 2.71828182846**(7.997*(df.loc[i, 'Pyr_hor']/1000 - 0.586))))*df.loc[i, 'Pyr_hor'])
            i_dir.append(df.loc[i, 'Pyr_hor']-i_dif[-1])
            
            #----------------------------------------Wärmespeicher----------------------------------------
            ws = WS(Vpkt_prim       = df.loc[i, 'HM10_Vpkt']/1000*0.8, 
                    T_prim_in       = df.loc[i, 'WP_sek_T_SUP']-1, 
                    Vpkt_sek        = df.loc[i, 'HM1_Vpkt']/1000*1.5, 
                    T_sek_in        = df.loc[i, 'HM1_T_RET'], 
                    T_vektor_prev   = WS_theta[-1],
                    Tamb            = df.loc[i, 'Zone_T'])
            WS_theta.append(ws)
            if i.hour == 0 and i.minute == 0 and i.second == 0:
                WS_theta[-1] = np.linspace(df['Warm Unten'][i], df['Warm Oben'][i], param['WS']['N']).tolist() 
                
            #----------------------------------------Kältespeicher----------------------------------------
            ks = KS(Vpkt_prim       = df.loc[i, 'HM4_Vpkt']/1000, 
                    T_prim_in       = df.loc[i, 'HM4_T_RET'], 
                    Vpkt_sek        = 715/1000*[1 if df.loc[i, 'WP_COP']>0 else 0][0], 
                    T_sek_in        = df.loc[i, 'WP_prim_T_RET']+3,                     # Anpassung durch IWT
                    T_vektor_prev   = KS_theta[-1],
                    Tamb            = df.loc[i, 'Zone_T'])
            KS_theta.append(ks)
            if i.hour == 0 and i.minute == 0 and i.second == 0:
                KS_theta[-1] = np.linspace(df['Kalt Unten'][i], df['Kalt Oben'][i], param['KS']['N']).tolist()
                
            #-------------------------------------------TABS----------------------------------------------
            tabs = TABS(TABS_theta_in = df.loc[i, 'HM4_T_SUP'], 
                        TABS_Vpkt = df.loc[i, 'HM4_Vpkt']/1000,
                        B_theta_air = df.loc[i, 'Zone_T'])
            TABS_theta_in.append(tabs[0])
            TABS_theta_out.append(tabs[1])
            TABS_Vpkt.append(tabs[2])
            TABS_Qpkt.append(tabs[3])
            
            #-------------------------------------------BAT-----------------------------------------------
            bat = BAT(BAT_Pel_in = (df.loc[i, 'PVT_Pel'] + df.loc[i, 'INV_P_IN'])*1000, 
                      BAT_Pel_out = df.loc[i, 'Bat_Pel']*1000, 
                      BAT_SOC_prev = BAT_SOC[-1])
            BAT_SOC.append(bat)
            if i.hour == 0 and i.minute == 0 and i.second == 0:
                BAT_SOC[-1] = 0.4348*df.loc[df['WP_COP'] == 0,  'BAT_U'][0] -24.088
            
            #-------------------------------------------PVT-----------------------------------------------
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
            
            #-------------------------------------------WP-----------------------------------------------
            wp = WP(WP_prim_theta_in = df.loc[i, 'WP_prim_T_SUP'], 
                    WP_sek_theta_in = df.loc[i, 'WP_sek_T_RET'], 
                    WP_sek_delta_theta = 6, 
                    WP_prim_Vpkt = 715/1000, 
                    WP_sek_Vpkt = df.loc[i, 'HM10_Vpkt']/1000, 
                    WP_OI = [1 if df.loc[i, 'WP_COP']>0 else 0][0])
            WP_COP.append(wp[0])
            WP_Pel.append(wp[1])
            WP_Qdot_prim.append(wp[2])
            WP_Qdot_sec.append(wp[3])
            WP_prim_theta_in.append(wp[4])
            WP_prim_theta_out.append(wp[5])
            WP_sek_theta_in.append(wp[6])
            WP_sek_theta_out.append(wp[7])
            WP_prim_Vpkt.append(wp[8])
            WP_sek_Vpkt.append(wp[9])

            #-------------------------------------------B-----------------------------------------------
            b = B(B_phi_int = 90,
                  B_phi_HC = df.loc[i, 'HM4_B'], 
                  B_theta_m_prev = B_theta_m[-1], 
                  theta_e = df.loc[i, 'AMB_T'], 
                  sun_h = sun_h[-1], 
                  sun_az = sun_az[-1], 
                  I_dir = i_dir[-1], 
                  I_diff = i_dif[-1], 
                  V_Vdot = 30/4, 
                  V_theta_sup = df.loc[i, 'AMB_T'] )

            # if i.hour in [12, 0] and i.minute == 0 and i.second == 0:
            if i.hour in [0] and i.minute == 0 and i.second == 0:
                B_theta_m.append(df.loc[i, 'Zone_T'])
            else:
                B_theta_m.append(b[0])
            B_theta_s.append(b[1])
            B_theta_air.append(b[2])
            B_theta_operativ.append(b[3])
                   
        df['WS_theta'] = WS_theta[1:]
        df['KS_theta'] = KS_theta[1:]
        df['TABS_theta_in'] = TABS_theta_in
        df['TABS_theta_out'] = TABS_theta_out
        df['TABS_Vpkt'] = TABS_Vpkt
        df['TABS_Qpkt'] = TABS_Qpkt
        df['BAT_SOC'] = BAT_SOC[1:]
        df['s_PVT_Pel'] = PVT_Pel[1:]
        df['PVT_theta_in'] = PVT_theta_in[1:]
        df['PVT_theta_out'] = PVT_theta_out[1:]
        df['PVT_Qpkt'] = PVT_Qpkt[1:]
        df['PVT_theta_cell'] = PVT_theta_cell[1:]
        df['s_WP_COP'] = WP_COP
        df['s_WP_Pel'] = WP_Pel
        df['s_WP_Qdot_prim'] = WP_Qdot_prim
        df['s_WP_Qdot_sec'] = WP_Qdot_sec
        df['s_WP_prim_theta_in'] = WP_prim_theta_in
        df['s_WP_prim_theta_out'] = WP_prim_theta_out
        df['s_WP_sek_theta_in'] = WP_sek_theta_in
        df['s_WP_sek_theta_out'] = WP_sek_theta_out
        df['s_WP_prim_Vpkt'] = WP_prim_Vpkt
        df['s_WP_sek_Vpkt'] = WP_sek_Vpkt
        df['sun_h'] = sun_h
        df['sun_az'] = sun_az  
        df['i_dif'] = i_dif
        df['i_dir'] = i_dir
        df['s_B_theta_m'] = B_theta_m[1:]
        df['s_B_theta_s'] = B_theta_s[1:]
        df['s_B_theta_air'] = B_theta_air[1:]
        df['s_B_theta_operativ'] = B_theta_operativ[1:]
        
        # -----------------------------------------------------------------------------------------------------------------
        # -------------------------------------PVT-el. Seite------------------------------------------------------------
        # -----------------------------------------------------------------------------------------------------------------
        Q_mess = []
        Q_sim = []
        
        for i in df.index.day.unique():
            aa = df[df.index.day == i]
            Q_mess.append(aa['PVT_Pel'].sum()*0.25)
            Q_sim.append(aa['s_PVT_Pel'].sum()*0.25/1000)
            
        fig, ax42 = plt.subplots(figsize=(6.4, 4.8/1.3))
        # ------------------------------------------------Y1------------------------------------------------
        # Data
        ax42.bar([i-0.2 for i in range(1,8,1)], Q_mess, width = 0.4, label = 'Messwert')
        ax42.bar([i+0.2 for i in range(1,8,1)], Q_sim, width = 0.4, label = 'Simulationswert')
        # labels
        ax42.set_ylabel('Wärmemenge in kWh', fontsize = labelfont)
        # ax1.yaxis.set_major_formatter('{:.1f}'.format)
        ax42.set_yticks(range(0,30,1), minor = True)
        ax42.set_yticks(range(0,30,2))
        ax42.set_ylim(0,10)
        
        # ------------------------------------------------X------------------------------------------------
        ax42.set_xticks(range(1,8,1), ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], rotation = 60, fontsize = labelfont)
        
        # ax1.set_xticks([i / 2 for i in range(0,7*2,1)], minor = True)
        
        # ------------------------------------------------Global------------------------------------------------
        
        ax42.grid(True)
        ax42.grid(which='minor', alpha = 0.3)
        ax42.legend(fontsize = 8, loc = 'upper left')

        ax42.text(5.6, 7.2,  r'$\bf{Gesamte\ Woche:}$' + '\n' + str(round(df['PVT_Pel'].sum()*0.25,2)) + ' kWh gemessen'
                  + '\n' + str(round(df['s_PVT_Pel'].sum()*0.25/1000,2)) + ' Wh simuliert' + '\n' +
                  str(round((df['s_PVT_Pel'].sum()*0.25/1000-df['PVT_Pel'].sum()*0.25)/(df['PVT_Pel'].sum()*0.25)*100,2)) + ' % Abweichung',
                  fontsize = labelfont, bbox=dict(facecolor='white', edgecolor='black'))
        
        
        
       