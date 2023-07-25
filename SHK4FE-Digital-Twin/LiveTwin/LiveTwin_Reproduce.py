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
        PVT_Pel = G*param['PVT']['etha_el']*(1-param['PVT']['beta_el']*(PVT_theta_cell-25))*(param['PVT']['A']*param['PVT']['quantity']/2)   
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
    'quantity'  : 6,                                                            # -         number of modules    
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
for file in os.listdir(os.getcwd()+'\Data')[a:]:
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
    
    fig, ((ax01, ax02), (ax11, ax12), (ax21, ax22), (ax31, ax32), (ax41, ax42)) = plt.subplots(5, 2, figsize = (6.4*2/0.9, 4.8/1.5*5/0.842))
    ax02_twin = ax02.twinx()
    ax11_twin = ax11.twinx()
    ax12_twin = ax12.twinx()
    ax21_twin = ax21.twinx()
    ax22_twin = ax22.twinx()
    ax31_twin = ax31.twinx()
    ax32_twin = ax32.twinx()
    ax02.tick_params(axis='both', labelsize=tickfont)
    ax02_twin.tick_params(axis='both', labelsize=tickfont)
    ax11.tick_params(axis='both', labelsize=tickfont)
    ax11_twin.tick_params(axis='both', labelsize=tickfont)
    ax12.tick_params(axis='both', labelsize=tickfont)
    ax12_twin.tick_params(axis='both', labelsize=tickfont)
    ax21.tick_params(axis='both', labelsize=tickfont)
    ax21_twin.tick_params(axis='both', labelsize=tickfont)
    ax22.tick_params(axis='both', labelsize=tickfont)
    ax22_twin.tick_params(axis='both', labelsize=tickfont)
    ax31.tick_params(axis='both', labelsize=tickfont)
    ax31_twin.tick_params(axis='both', labelsize=tickfont)
    ax32.tick_params(axis='both', labelsize=tickfont)
    ax32_twin.tick_params(axis='both', labelsize=tickfont)
    ax41.tick_params(axis='both', labelsize=tickfont)
    ax42.tick_params(axis='both', labelsize=tickfont)
    fig.suptitle(r'$\bf{Verknüpfung\ mit\ Temperatur-\ und\ Volumenstrommesswerten}$' + '\nDatum: ' + str(df.index[0].strftime('%d.%m.%Y')), 
                  fontsize=font)
    plt.subplots_adjust(top=0.94)
    
    #--------------------------------------------------------------------------------------------------
    #-----------------------------------------Bilanz---------------------------------------------------
    #--------------------------------------------------------------------------------------------------
    ax01.set_title(r'$\bf{Energiebilanz\ für\ den\ Tag}$', fontsize = font_title)
    a_Grid = round(df['INV_P_IN'].sum()*0.25,1)
    a_PVT = round(df['PVT_Pel'].sum()*0.25,1)
    a_PVT_s = round(df['s_PVT_Pel'].sum()*0.25/1000,1)
    a_WP = round(df.loc[df['WP_COP']>0, 'Bat_Pel'].sum()*0.25,1)
    a_WP_s = round(df['s_WP_Pel'].sum()*0.25/1000,1)
    a_user = round(df.loc[df['WP_COP']==0, 'Bat_Pel'].sum()*0.25,1)
    a_room = round(df['HM4_Qpkt'].sum()*0.25/-1000,1)
    a_room_s = round(df['TABS_Qpkt'].sum()*0.25/-1000,1)
    a_PVT_th = round(df['HM1_Qpkt'].sum()*0.25/-1000,1)
    a_PVT_th_s = round(df['PVT_Qpkt'].sum()*0.25/-1000,1)
    
    ax01.imshow(mpimg.imread(os.path.dirname(os.path.abspath(__file__)) + '/Energiebilanz.JPG'), extent=[0, 7.21, 0, 3.8])
    ax01.set_yticks([])
    ax01.set_xticks([])
    # Grid
    ax01.text(0.7, 1.7, str(a_Grid) + ' kWh',c = 'tab:blue', fontsize=legendfont)
    # WP
    ax01.text(2.8, 1.7, str(a_WP) + ' kWh',c = 'tab:blue', fontsize=legendfont)
    ax01.text(2.8, 1.45, str(a_WP_s) + ' kWh',c = 'tab:red', fontsize=legendfont)
    # User
    ax01.text(2.8+1, 1.7-1.17, str(a_user) + ' kWh',c = 'tab:blue', fontsize=legendfont)
    ax01.text(2.8+1, 1.45-1.17, str(round(0.15*len(df.loc[df['WP_COP']==0, :])*0.25,1)) + ' kWh',c = 'tab:red', fontsize=legendfont)
    # PVT_el
    ax01.text(1.2, 1.7+1, str(a_PVT) + ' kWh',c = 'tab:blue', fontsize=legendfont)
    ax01.text(1.2, 2.7-0.25, str(a_PVT_s) + ' kWh',c = 'tab:red', fontsize=legendfont)
    # PVT_th
    ax01.text(4.6, 3, str(a_PVT_th) + ' kWh',c = 'tab:blue', fontsize=legendfont)
    ax01.text(4.6, 3-0.25, str(a_PVT_th_s) + ' kWh',c = 'tab:red', fontsize=legendfont)
    # Room
    ax01.text(5.3, 2.05, str(a_room) + ' kWh',c = 'tab:blue', fontsize=legendfont)
    ax01.text(5.3, 2.05-0.25, str(a_room_s) + ' kWh',c = 'tab:red', fontsize=legendfont)
    # Room
    ax01.text(5.3, 2.05-0.745, str(0.0) + ' kWh',c = 'tab:blue', fontsize=legendfont)
    ax01.text(5.3, 2.05-0.745-0.25, str(0.0) + ' kWh',c = 'tab:red', fontsize=legendfont)
    # Legend
    ax01.text(0.1, 3.6, 'Messwert',c = 'tab:blue', fontsize=legendfont)
    ax01.text(0.1, 3.6-0.25, 'Simulationswert' ,c = 'tab:red', fontsize=legendfont)
    
    #--------------------------------------------------------------------------------------------------
    #---------------------------------------Gebäudemodell----------------------------------------------
    #--------------------------------------------------------------------------------------------------
    ax02.plot(df.index, df['AMB_T'], c = 'grey', linestyle = '-', alpha = 0.5, label = 'Mess: Außentemperatur', marker = 'o', markersize = 4)
    ax02.plot(df.index, df['Zone_T'], c = 'tab:green', linestyle = '-', alpha = 0.5, label = 'Mess Innentemperatur', marker = 'o', markersize = 4)
    # Simulation
    ax02.plot(df.index, df['s_B_theta_operativ'], c = 'tab:green', linestyle = '-', alpha = 1, label = 'Sim: Innenraumtemperatur', marker = 'x', markersize = 4)
    # labels
    ax02.legend(fontsize=legendfont, loc='upper left')
    ax02.set_ylabel('Temperatur in °C', fontsize = labelfont)
    ax02.yaxis.set_major_formatter('{:.0f}'.format)
    ax02.set_yticks(np.arange(10,40,2.5), minor = True)
    ax02.set_yticks(range(10,40,5))
    ax02.set_ylim(10,40)
    # ------------------------------------------------Y2------------------------------------------------
    # Data
    ax02_twin.plot(df.index, [i*-1 for i in df['HM4_Qpkt']], c = 'tab:blue', linestyle = '-', alpha = 0.5, label = 'Mess: Kälteleistung', marker = 'o', markersize = 4)
    ax02_twin.plot(df.index, df['Pyr_hor'], c = 'tab:orange', linestyle = '-', alpha = 0.5, label = 'Mess: Globalstrahlung', marker = 'o', markersize = 4)
    # labels
    ax02_twin.legend(loc = 'upper right', fontsize = legendfont) 
    ax02_twin.set_ylabel('Kälteleistung in -W \n Globalstrahlung in W/m²', fontsize = labelfont)
    ax02_twin.yaxis.set_major_formatter('{:.0f}'.format)
    ax02_twin.set_yticks(range(0,3601,600))
    ax02_twin.set_ylim(0,3600)
    # ------------------------------------------------X------------------------------------------------
    ax02.set_xticks([df.index[0] + timedelta(minutes=i*60*6) for i in range(int(672/4/6))], minor = True)
    ax02.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], ['']*84)
    ax02.set_xlim(df.index[0], df.index[0]+timedelta(hours = 24))
    ax02.set_title(r'$\bf{Innentemperatur}$', fontsize = font_title)
    # ------------------------------------------------Global------------------------------------------------
    ax02.grid(True)
    ax02.grid(which='minor', alpha = 0.3)
    MAE = (sum([abs(a - b) for a, b in zip(df['AMB_T'], df['s_B_theta_operativ'])]) /
            len([abs(a - b) for a, b in zip(df['AMB_T'], df['s_B_theta_operativ'])]))
    # ax02.text(df.index[0]+timedelta(hours = 0.31), 15, r'$\bf{Mean\ Absolute\ Error:}$' + str(round(MAE, 2)) + ' K\n' , bbox={'boxstyle': 'round', 'facecolor': 'white', 'edgecolor': 'grey'}, fontsize=legendfont)
    data = [
    ['Mean', str(round(df['Zone_T'].mean(),1)) + ' °C', str(round(df['s_B_theta_operativ'].mean(),1)) + ' °C', str(round(abs(df['Zone_T'].mean() - df['s_B_theta_operativ'].mean())/df['Zone_T'].mean()*100,1)) + ' %'],
    ['Max' , str(round(df['Zone_T'].max(),1)) + ' °C' , str(round(df['s_B_theta_operativ'].max(),1)) + ' °C' , str(round(abs(df['Zone_T'].max() - df['s_B_theta_operativ'].max())/df['Zone_T'].max()*100,1)) + ' %'],
    ['Min' , str(round(df['Zone_T'].min(),1)) + ' °C' , str(round(df['s_B_theta_operativ'].min(),1)) + ' °C' , str(round(abs(df['Zone_T'].min() - df['s_B_theta_operativ'].min())/df['Zone_T'].min()*100,1)) + ' %']
    ]
    table = ax02.table(cellText=data, colLabels=['', 'Mess', 'Sim', 'DEV'], loc='center', cellLoc='center', colColours=['lightgray']*4,
                      bbox=[0.01, 0.61, 0.45, 0.2]) # Tabelle formatieren
    table.auto_set_font_size(False)
    table.set_fontsize(legendfont)
    table.set_zorder(100)
    
    #--------------------------------------------------------------------------------------------------
    #---------------------------------------Wärmespeicher----------------------------------------------
    #--------------------------------------------------------------------------------------------------
    abweichung_top = [x for x in [abs(a-b) for a,b in zip(df['Warm Oben'][df.index[0]:df.index[-1]], [posEl[-1] for posEl in df['WS_theta'][df.index[0]:df.index[-1]]])] if not math.isnan(x)]
    MAE_top = sum(abweichung_top) / len(abweichung_top)
    abweichung_bottom = [x for x in [abs(a-b) for a,b in zip(df['Warm Unten'][df.index[0]:df.index[-1]], [posEl[0] for posEl in df['WS_theta'][df.index[0]:df.index[-1]]])]  if not math.isnan(x)]
    MAE_bottom = sum(abweichung_bottom) / len(abweichung_bottom)
    # ------------------------------------------------Y1------------------------------------------------
    ax11.plot(df.index, df['Warm Oben'], label = 'Mess: Top', c = 'tab:red', alpha = 0.5, linestyle = '--', marker = 'o', markersize = 4)
    ax11.plot(df.index, df['Warm Unten'], label = 'Mess: Bottom', c = 'tab:blue', alpha = 0.5, linestyle = '--', marker = 'x', markersize = 4)
    ax11.plot(df.index, [posEl[-1] for posEl in df['WS_theta']], label = 'Sim: Top', c = 'tab:red')
    ax11.plot(df.index, [posEl[0] for posEl in df['WS_theta']], label = 'Sim: Bottom', c = 'tab:blue')
    # labels
    ax11.set_ylabel('Temperatur in °C', fontsize = labelfont)
    ax11.yaxis.set_major_formatter('{:.1f}'.format)
    ax11.set_ylim(25,60)
    ax11.set_yticks(np.arange(25,61,2.5), minor = True)
    ax11.set_yticks(range(25,61,5))
    ax11.set_title(r'$\bf{Wärmespeicher}$', fontsize = font_title)
    # ------------------------------------------------X1------------------------------------------------
    ax11.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    ['']*84,
                    rotation = 90)
    ax11.set_xlim(df.index[0], df.index[0]+timedelta(hours = 24))
    ax11.grid(True)
    ax11.grid(which='minor', alpha = 0.3)
    ax11.legend(fontsize = legendfont, loc = 'upper left', ncol = 2)
    ax11_twin.plot(df.index, [2 if value > 0 else value for value in df['HM1_Vpkt']],label = 'Rückkühlung (PVT)', c = 'tab:green')
    ax11_twin.fill_between(df.index, [2 if value > 0 else value for value in df['HM1_Vpkt']],color = 'tab:green', alpha = 0.3)
    ax11_twin.plot(df.index, [2 if value > 0 else value for value in df['HM10_Vpkt']], label = 'Beladung (WP)', c = 'tab:orange')
    ax11_twin.fill_between(df.index, [2 if value > 0 else value for value in df['HM10_Vpkt']],color = 'tab:orange', alpha = 0.3)
    ax11_twin.legend(fontsize = legendfont, loc = 'upper right', ncol = 1)
    ax11_twin.set_ylabel('AN / AUS', fontsize = labelfont)
    ax11_twin.set_yticks(range(0,3,2), [0,1])
    ax11_twin.set_ylim(0,14)
    ax11.text(df.index[0]+timedelta(hours = 0.5), 49., 'MAE @ Top: ' + str(round(MAE_top, 2)) + ' K' + '\nMAE @ Bot: ' + str(round(MAE_bottom, 2))+ ' K' , bbox={'boxstyle': 'round', 'facecolor': 'white', 'edgecolor': 'grey'}, fontsize=legendfont)

    
    #--------------------------------------------------------------------------------------------------
    #---------------------------------------Kältespeicher----------------------------------------------
    #--------------------------------------------------------------------------------------------------
    abweichung_top = [x for x in [abs(a-b) for a,b in zip(df['Kalt Oben'][df.index[0]:df.index[-1]], [posEl[-1] for posEl in df['KS_theta'][df.index[0]:df.index[-1]]])] if not math.isnan(x)]
    MAE_top = sum(abweichung_top) / len(abweichung_top)
    abweichung_bottom = [x for x in [abs(a-b) for a,b in zip(df['Kalt Unten'][df.index[0]:df.index[-1]], [posEl[0] for posEl in df['KS_theta'][df.index[0]:df.index[-1]]])]  if not math.isnan(x)]
    MAE_bottom = sum(abweichung_bottom) / len(abweichung_bottom)
    # ------------------------------------------------Y1------------------------------------------------
    ax12.plot(df.index, df['Kalt Oben'], label = 'Mess: Top', c = 'tab:red', alpha = 0.5, linestyle = '--', marker = 'o', markersize = 4)
    ax12.plot(df.index, df['Kalt Unten'], label = 'Mess: Bottom', c = 'tab:blue', alpha = 0.5, linestyle = '--', marker = 'x', markersize = 4)
    ax12.plot(df.index, [posEl[-1] for posEl in df['KS_theta']], label = 'Sim: Top', c = 'tab:red')
    ax12.plot(df.index, [posEl[0] for posEl in df['KS_theta']], label = 'Sim: Bottom', c = 'tab:blue')
    # labels
    ax12.set_ylabel('Temperatur in °C', fontsize = labelfont)
    ax12.set_ylim(5,25)
    # ax1.set_yticks(np.arange(25,55,2.5), minor = True)
    # ax1.set_yticks(range(25,55,5))
    ax12.set_title(r'$\bf{Kältespeicher}$', fontsize = font_title)
    # ------------------------------------------------X1------------------------------------------------
    ax12.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    ['']*84,
                    rotation = 90)
    ax12.set_xlim(df.index[0], df.index[0]+timedelta(hours = 24))
    ax12.grid(True)
    ax12.grid(which='minor', alpha = 0.3)
    ax12.legend(fontsize = legendfont, loc = 'upper left', ncol = 2)
    
    ax12_twin.plot(df.index, df['WP_COP'].apply(lambda x: 2 if x > 0 else x),label = 'Kühlung (WP)', c = 'tab:green')
    ax12_twin.fill_between(df.index, df['WP_COP'].apply(lambda x: 2 if x > 0 else x), color = 'tab:green', alpha = 0.3)
    ax12_twin.plot(df.index, [2 if value > 0 else value for value in df['HM4_Vpkt']], label = 'Entladung (TABS)', c = 'tab:orange')
    ax12_twin.fill_between(df.index, [2 if value > 0 else value for value in df['HM4_Vpkt']],color = 'tab:orange', alpha = 0.3)
    
    ax12_twin.legend(fontsize = legendfont, loc = 'upper right', ncol = 1)
    ax12_twin.set_ylabel('AN / AUS', fontsize = labelfont)
    ax12_twin.set_yticks(range(0,3,2), [0,1])
    ax12_twin.set_ylim(0,16)
    ax12.text(df.index[0]+timedelta(hours = 0.5), ((49.-25)/(60-25)*20)+5, 'MAE @ Top: ' + str(round(MAE_top, 2)) + ' K' + '\nMAE @ Bot: ' + str(round(MAE_bottom, 2))+ ' K' , bbox={'boxstyle': 'round', 'facecolor': 'white', 'edgecolor': 'grey'}, fontsize=legendfont) 
    
    #--------------------------------------------------------------------------------------------------
    #---------------------------------------------TABS-------------------------------------------------
    #--------------------------------------------------------------------------------------------------
    # Data
    ax21.scatter(df.index, df['HM4_T_SUP'], label = 'Mess: Vorlauftemperatur', c = 'tab:red', s = 12, alpha = 0.5)
    ax21.plot(df.index, df['HM4_T_SUP'], c = 'tab:red', linestyle = '--', alpha = 0.5)
    ax21.scatter(df.index, df['HM4_T_RET'], label = 'Mess: Rücklauftemperatur', c = 'tab:blue', s = 12, alpha = 0.5)
    ax21.plot(df.index, df['HM4_T_RET'], c = 'tab:blue', linestyle = '--', alpha = 0.5)
    # Simulation
    ax21.scatter(df.index, df['TABS_theta_out'], label = 'Sim: Rücklauftemperatur', c = 'tab:blue', s = 12, marker = 'x')
    ax21.plot(df.index, df['TABS_theta_out'], c = 'tab:blue', linestyle = '--')
    # labels
    ax21.set_ylabel('Temperatur in °C', fontsize = labelfont)
    ax21.yaxis.set_major_formatter('{:.0f}'.format)
    ax21.set_ylim(0,30)
    ax21.set_yticks(range(0,30,2), minor = True)
    ax21.set_yticks(range(0,30,4))
    # ------------------------------------------------Y2------------------------------------------------
    # Data
    ax21_twin.scatter(df.index, [i/1000 for i in abs(df['TABS_Qpkt'])], label = 'Sim: Kälteleistung', c = 'tab:green', s = 12, marker = 'x')
    # ax21_twin.plot(df.index, abs(df['TABS_Qpkt']), c = 'tab:green', linestyle = '--')
    ax21_twin.scatter(df.index, [i/1000 for i in abs(df['HM4_Qpkt'])], label = 'Mess: Kälteleistung', c = 'tab:green', s = 12, alpha = 0.5)
    # ax21_twin.plot(df.index, abs(df['HM4_Qpkt']), c = 'tab:green', linestyle = '--', alpha = 0.5)
    # labels
    ax21_twin.set_ylabel('Kälteleistung in -kW', fontsize = labelfont)
    ax21_twin.yaxis.set_major_formatter('{:.1f}'.format)
    ax21_twin.set_yticks(np.arange(0,7,2))
    ax21_twin.set_ylim(0,6)
    # ------------------------------------------------X------------------------------------------------
    ax21.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], ['']*84, rotation = 90)
    ax21.set_xlim(df.index[0], df.index[0]+timedelta(hours = 24))
    # ------------------------------------------------Global------------------------------------------------
    ax21.set_title(r'$\bf{Thermoaktives\ Bauteilsystem}$', fontsize = font_title)
    ax21.grid(True)
    ax21.grid(which='minor', alpha = 0.3)
    ax21.legend(fontsize = legendfont, loc = 'upper left')
    ax21_twin.legend(fontsize = legendfont, loc = 'upper right') 
    
    ax21.text(df.index[0] + timedelta(hours = 0.5), 4,  r'$\bf{Tagesbilanz:}$' + '\n' + str(round(df['HM4_Qpkt'].sum()*0.25,1)) + ' Wh gemessen'
          + '\n' + str(round(df['TABS_Qpkt'].sum()*0.25,1)) + ' Wh simuliert' + '\n' +
          str(round((df['TABS_Qpkt'].sum()*0.25-df['HM4_Qpkt'].sum()*0.25)/(df['HM4_Qpkt'].sum()*0.25)*100,2)) + ' % Abweichung',
          bbox={'boxstyle': 'round', 'facecolor': 'white', 'edgecolor': 'grey'}, fontsize=legendfont)
    
    #--------------------------------------------------------------------------------------------------
    #---------------------------------------------BAT--------------------------------------------------
    #--------------------------------------------------------------------------------------------------
    # Data
    ax22.scatter(df.index, df['BAT_U'], label = 'Mess: Batteriespannung', c = 'grey', s = 12, alpha = 0.5)
    ax22.plot(df.index, df['BAT_U'], c = 'grey', linestyle = '--', alpha = 0.5)
    # labels
    ax22.set_ylabel('Batteriespannung in V', fontsize = labelfont)
    # ax1.yaxis.set_major_formatter('{:.0f}'.format)
    # ax1.set_ylim(55.06,57.76)
    ax22.set_ylim(55,58)
    # ax22.set_yticks(np.arange(55,58.01,0.33), minor = True)
    ax22.set_yticks(np.arange(55,58.01,0.3))

    ax22_twin.scatter(df.index, df['BAT_SOC'], label = 'Sim: SOC', c = 'black', s = 12, marker = 'x')
    ax22_twin.plot(df.index, df['BAT_SOC'], c = 'black', linestyle = '--')
    
    ax22_twin.set_yticks(np.arange(0,1.01,0.1), minor = True)
    ax22_twin.set_yticks(np.arange(0,1.01,0.2), np.arange(0,101,20))
    ax22_twin.set_ylim(0,1)
    ax22_twin.legend(fontsize = legendfont, loc = 'upper right')
    ax22.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], ['']*84,
                    rotation = 90)
    ax22.legend(fontsize = legendfont, loc = 'upper left')
    ax22.grid(True)
    ax22.set_title(r'$\bf{Batteriespeicher}$', fontsize = font_title)
    ax22.set_xlim(df.index[0], df.index[0]+timedelta(hours = 24))
    

    #--------------------------------------------------------------------------------------------------
    #------------------------------------------PVT-therm-----------------------------------------------
    #--------------------------------------------------------------------------------------------------
    # Data
    ax31.scatter(df.index, df['HM1_T_SUP'], label = 'Mess: Outflow', c = 'tab:red', s = 12, alpha = 0.5)
    ax31.plot(df.index, df['HM1_T_SUP'], c = 'tab:red', linestyle = '--', alpha = 0.5)
    ax31.scatter(df.index, df['HM1_T_RET'], label = 'Mess: Inflow', c = 'tab:blue', s = 12, alpha = 0.5)
    ax31.plot(df.index, df['HM1_T_RET'], c = 'tab:blue', linestyle = '--', alpha = 0.5)
    ax31.scatter(df.index, df['PVT_theta_out'], label = 'Sim: Outflow', c = 'tab:red', s = 12, marker = 'x')
    ax31.plot(df.index, df['PVT_theta_out'], c = 'tab:red', linestyle = '--', alpha = 0.5)
    ax31.plot(df.index, df['AMB_T'], c = 'tab:green', alpha = 0.7, label = 'Umgebungstemperatur')
    # labels
    ax31.set_ylabel('Temperatur in °C', fontsize = labelfont)
    ax31.yaxis.set_major_formatter('{:.0f}'.format)
    ax31.set_ylim(0,60)
    ax31.set_yticks(np.arange(0,60,2.5), minor = True)
    ax31.set_yticks(range(0,60,5))
    # ------------------------------------------------Y2------------------------------------------------
    # Data
    ax31_twin.plot(df.index, [i/1000 for i in df['Pyr_hor']], c = 'tab:orange', alpha = 0.7, label = 'Globalstrahlung')
    # labels
    ax31_twin.set_ylabel('Globalstrahlung in kW/m²', fontsize = labelfont)
    # ax31_twin.yaxis.set_major_formatter('{:.0f}'.format)
    ax31_twin.set_yticks(np.arange(0,2.1,0.3))
    ax31_twin.set_ylim(0,2.1)
    # ------------------------------------------------X------------------------------------------------
    ax31.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                ['']*84, rotation = 90)
    ax31.set_xlim(df.index[0], df.index[0]+timedelta(hours = 24))
    # ------------------------------------------------Global------------------------------------------------
    ax31.set_title(r'$\bf{PVT-Kollektoren:\ thermische\ Seite}$', fontsize = font_title)
    ax31.grid(True)
    ax31.grid(which='minor', alpha = 0.3)
    ax31.legend(fontsize = 8, loc = 'upper left', ncol = 1)
    ax31_twin.legend(fontsize = 8, loc = 'upper right', ncol = 2)
    ax31.text(df.index[0] + timedelta(hours = 0.5), 3,  r'$\bf{Tagesbilanz:}$' + '\n' + str(round(df['HM1_Qpkt'].sum()*0.25,1)) + ' Wh gemessen'
              + '\n' + str(round(df['PVT_Qpkt'].sum()*0.25,1)) + ' Wh simuliert' + '\n' +
              str(round((df['PVT_Qpkt'].sum()*0.25-df['HM1_Qpkt'].sum()*0.25)/(df['HM1_Qpkt'].sum()*0.25)*100,2)) + ' % Abweichung',
              bbox={'boxstyle': 'round', 'facecolor': 'white', 'edgecolor': 'grey'}, fontsize=legendfont)
    
    #--------------------------------------------------------------------------------------------------
    #------------------------------------------PVT-el--------------------------------------------------
    #--------------------------------------------------------------------------------------------------
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
    ax32.set_title(r'$\bf{PVT-Kollektoren:\ elektrische\ Seite}$', fontsize = font_title)    
    ax32.grid(True)
    ax32.grid(which='minor', alpha = 0.3)
    ax32.text(df.index[0] + timedelta(hours = 0.5), 750,  r'$\bf{Tagesbilanz:}$' + '\n' + str(round(df['PVT_Pel'].sum()*0.25,1)) + ' kWh gemessen'
              + '\n' + str(round(df['s_PVT_Pel'].sum()*0.25/1000,1)) + ' kWh simuliert' + '\n' +
              str(round((df['s_PVT_Pel'].sum()*0.25/1000-df['PVT_Pel'].sum()*0.25)/(df['PVT_Pel'].sum()*0.25)*100,2)) + ' % Abweichung',
              bbox={'boxstyle': 'round', 'facecolor': 'white', 'edgecolor': 'grey'}, fontsize=legendfont)
        
    ax32.legend(fontsize = legendfont, loc = 'upper left', ncol = 2)
    ax32_twin.legend(fontsize = legendfont, loc = 'upper right', ncol = 2)
    
    
    
    #--------------------------------------------------------------------------------------------------
    #------------------------------------------WP-leistung---------------------------------------------
    #--------------------------------------------------------------------------------------------------
    ax41.scatter(df.index[0], -99, label = 'COP', alpha = 0)    
    ax41.scatter(df.index, df['WP_COP'], label = 'Mess', c = 'black', s = 12, alpha = 0.5)
    ax41.scatter(df.index, df['s_WP_COP'], label = 'Sim', c = 'black', s = 12, marker = 'x')

    ax41.scatter(df.index[0], -99, label = 'therm. Leistung', alpha = 0)    
    ax41.scatter(df.index, df['WP_Qpkt'], label = 'Mess', c = 'tab:red', s = 12, alpha = 0.5)
    ax41.scatter(df.index, [i/1000 for i in df['s_WP_Qdot_sec']], label = 'Sim', c = 'tab:red', s = 12, marker = 'x')
    ax41.scatter(df.index[0], -99, label = 'el. Leistung', alpha = 0)  
    ax41.scatter(df.index, df['WP_Pel'], label = 'Mess', c = 'tab:orange', s = 12, alpha = 0.5)
    ax41.scatter(df.index, [i/1000 for i in df['s_WP_Pel']], label = 'Sim', c = 'tab:orange', s = 12,marker = 'x')
    # SetUp
    ax41.set_ylabel('COP in -, el. und therm. Leistung in kW', fontsize = labelfont) 
    # ax41.yaxis.set_major_formatter('{:.1f}'.format)
    ax41.set_yticks(np.arange(0,9,0.5), minor = True)
    ax41.set_yticks(np.arange(0,9,1))
    ax41.set_ylim(0,9)
    ax41.set_title(r'$\bf{Wärmepumpe: Leistungsdaten}$', fontsize = font_title)
    ax41.grid(True)
    ax41.grid(which='minor', alpha = 0.3)
    legend = ax41.legend(fontsize = legendfont, loc = 'upper left', ncol = 3)
    legend_texts = legend.get_texts()
    legend_texts[0].set_weight('bold')  # 1. Element fett (Index 0)
    legend_texts[3].set_weight('bold')  # 3. Element fett (Index 2)
    legend_texts[6].set_weight('bold')  # 3. Element fett (Index 2)
    
    ax41.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))],
                    rotation = 60)
    ax41.set_xlim(df.index[0], df.index[0]+timedelta(hours = 24))
    
    MAE_COP = (sum([abs(a-b) for a, b in zip(df.loc[df['WP_COP'] > 0,'s_WP_COP'], df.loc[df['WP_COP'] > 0,'WP_COP'])])/
                len([abs(a-b) for a, b in zip(df.loc[df['WP_COP'] > 0,'s_WP_COP'], df.loc[df['WP_COP'] > 0,'WP_COP'])]))
    MAE_Pel = (sum([abs(a-b*1000) for a, b in zip(df.loc[df['WP_COP'] > 0,'s_WP_Pel'], df.loc[df['WP_COP'] > 0,'WP_Pel'])])/
                len([abs(a-b*1000) for a, b in zip(df.loc[df['WP_COP'] > 0,'s_WP_Pel'], df.loc[df['WP_COP'] > 0,'WP_Pel'])]))
    MAE_Qpkt = (sum([abs(a-b*1000) for a, b in zip(df.loc[df['WP_COP'] > 0,'s_WP_Qdot_sec'], df.loc[df['WP_COP'] > 0,'WP_Qpkt'])])/
                len([abs(a-b*1000) for a, b in zip(df.loc[df['WP_COP'] > 0,'s_WP_Qdot_sec'], df.loc[df['WP_COP'] > 0,'WP_Qpkt'])]))
    ax41.text(df.index[0]+timedelta(hours = 0.5), 0.5, r'$\bf{Mean\ Absolute\ Error:}$' + 
              '\nCOP:                        ' + str(round(MAE_COP, 2)) + 
              '\nel. Leistung:          ' + str(round(MAE_Pel, 2)) + ' W' +
              '\ntherm. Leistung: ' + str(round(MAE_Qpkt, 2)) + ' W' , bbox={'boxstyle': 'round', 'facecolor': 'white', 'edgecolor': 'grey'}, fontsize=legendfont)
    
    #--------------------------------------------------------------------------------------------------
    #------------------------------------------WP-Temperaturen-----------------------------------------
    #--------------------------------------------------------------------------------------------------
    ax42_twin = ax42.twinx()
    # ------------------------------------------------Y1------------------------------------------------
    # Data
    ax42.scatter(df.index, df['WP_prim_T_SUP'], label = 'Mess: Eintritt', c = 'tab:red', s = 12, alpha = 0.5)
    ax42.plot(df.index, df['WP_prim_T_SUP'], c = 'tab:red', linestyle = '--', alpha = 0.2)

    ax42.scatter(df.index, df['WP_prim_T_RET'], label = 'Mess: Austritt', c = 'tab:blue', s = 12, alpha = 0.5)
    ax42.plot(df.index, df['WP_prim_T_RET'], c = 'tab:blue', linestyle = '--', alpha = 0.2)
    
    ax42.plot(df.index[0], -99, alpha = 0, label = ' ')
    
    ax42.scatter(df.index, df['s_WP_prim_theta_out'], label = 'Sim: Austritt', c = 'tab:blue', s = 12, marker = 'x')
    ax42.plot(df.index, df['s_WP_prim_theta_out'], c = 'tab:blue', linestyle = '--')
    
    
    offset = 0
    ax42_twin.scatter(df.index, [i-offset for i in df['WP_sek_T_RET']], label = 'Mess: Eintritt', c = 'tab:blue', s = 12, alpha = 0.5)
    ax42_twin.plot(df.index, [i-offset for i in df['WP_sek_T_RET']], c = 'tab:blue', linestyle = '--', alpha = 0.5)
    
    
    ax42_twin.scatter(df.index, [i-offset for i in df['WP_sek_T_SUP']], label = 'Mess: Austritt', c = 'tab:red', s = 12, alpha = 0.5)
    ax42_twin.plot(df.index, [i-offset for i in df['WP_sek_T_SUP']], c = 'tab:red', linestyle = '--', alpha = 0.5)

    ax42_twin.plot(df.index[0], -99, alpha = 0, label = ' ')
      
    ax42_twin.scatter(df.index, [i-offset for i in df['s_WP_sek_theta_out']], label = 'Sim: Austritt', c = 'tab:red', s = 12, marker = 'x')
    ax42_twin.plot(df.index, [i-offset for i in df['s_WP_sek_theta_out']], c = 'tab:red', linestyle = '--')   
    
    ax42.set_ylabel('Temperatur in °C', fontsize = labelfont) 
    # ax42_twin.set_yticks(np.arange(5,65,100), minor = True)
    ax42_twin.set_yticks(np.arange(0,66,100))
    ax42_twin.set_ylim(5,65)
    
    ax42.set_yticks(np.arange(5,65,2.5), minor = True)
    ax42.set_yticks(np.arange(5,66,5))
    ax42.set_ylim(5,65)
    
    ax42.plot([df.index[0],df.index[0]+timedelta(hours = 24)], [25,25], c = 'grey')
    ax42.text(df.index[0]+timedelta(hours = 18.5), 25+1.5, '↑ Wärmesenke', fontsize=legendfont+1)
    ax42.text(df.index[0]+timedelta(hours = 18.5), 20+1.5, '↓ Wärmequelle', fontsize=legendfont+1)
    
    ax42.set_title(r'$\bf{Wärmepumpe: Temperaturen}$' )
    ax42.grid(True)
    ax42.grid(which='minor', alpha = 0.3)
    ax42.legend(fontsize = legendfont, loc = 'lower left', ncol = 2)
    ax42_twin.legend(fontsize = legendfont, loc = 'upper left', ncol = 2)
    ax42.set_xticks([df.index[0] + timedelta(minutes=i*60*2) for i in range(int(672/4/2))], 
                    [(df.index[0] + timedelta(minutes=i*60*2)).strftime('%H:%M') for i in range(int(672/4/2))],
                    rotation = 60)    
    ax42.set_xlim(df.index[0], df.index[0]+timedelta(hours = 24))
    
    
    MAE_Quelle = (sum([abs(a-b) for a, b in zip(df.loc[df['WP_COP'] > 0,'s_WP_prim_theta_out'], df.loc[df['WP_COP'] > 0,'WP_prim_T_RET'])])/
                len([abs(a-b) for a, b in zip(df.loc[df['WP_COP'] > 0,'s_WP_prim_theta_out'], df.loc[df['WP_COP'] > 0,'WP_prim_T_RET'])]))
    MAE_Senke = (sum([abs(a-b) for a, b in zip(df.loc[df['WP_COP'] > 0,'s_WP_sek_theta_out'], df.loc[df['WP_COP'] > 0,'WP_sek_T_SUP'])])/
                len([abs(a-b) for a, b in zip(df.loc[df['WP_COP'] > 0,'s_WP_sek_theta_out'], df.loc[df['WP_COP'] > 0,'WP_sek_T_SUP'])]))
    ax42.text(df.index[0]+timedelta(hours = 0.5), 18.5, r'$\bf{Mean\ Absolute\ Error:}$' + str(round(MAE_Quelle, 2)) + ' K' , bbox={'boxstyle': 'round', 'facecolor': 'white', 'edgecolor': 'grey'}, fontsize=legendfont)
    ax42.text(df.index[0]+timedelta(hours = 0.5), 50.3, r'$\bf{Mean\ Absolute\ Error:}$' + str(round(MAE_Senke, 2)) + ' K' , bbox={'boxstyle': 'round', 'facecolor': 'white', 'edgecolor': 'grey'}, fontsize=legendfont)
    plt.subplots_adjust(wspace=0.22)
    plt.show()
    
    
    fig.savefig(os.path.dirname(os.path.abspath(__file__)) + '/Grafik/' + file.split('.')[0] + '.png')
    
    print('day')
    
    
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
        

        font = 14
        font_title = 12
        labelfont = 9
        legendfont = 9
        tickfont = 9
        
        fig, ((ax01, ax02), (ax11, ax12), (ax21, ax22), (ax31, ax32), (ax41, ax42)) = plt.subplots(5, 2, figsize = (6.4*2/0.9, 4.8/1.5*5/0.842))
        ax02_twin = ax02.twinx()
        ax11_twin = ax11.twinx()
        ax12_twin = ax12.twinx()
        ax21_twin = ax21.twinx()
        ax22_twin = ax22.twinx()

        ax02.tick_params(axis='both', labelsize=tickfont)
        ax02_twin.tick_params(axis='both', labelsize=tickfont)
        ax11.tick_params(axis='both', labelsize=tickfont)
        ax11_twin.tick_params(axis='both', labelsize=tickfont)
        ax12.tick_params(axis='both', labelsize=tickfont)
        ax12_twin.tick_params(axis='both', labelsize=tickfont)
        ax21.tick_params(axis='both', labelsize=tickfont)
        ax21_twin.tick_params(axis='both', labelsize=tickfont)
        ax22.tick_params(axis='both', labelsize=tickfont)
        ax22_twin.tick_params(axis='both', labelsize=tickfont)
        ax31.tick_params(axis='both', labelsize=tickfont)
        ax32.tick_params(axis='both', labelsize=tickfont)
        ax41.tick_params(axis='both', labelsize=tickfont)
        ax42.tick_params(axis='both', labelsize=tickfont)



        fig.suptitle(r'$\bf{Verknüpfung\ mit\ Temperatur-\ und\ Volumenstrommesswerten}$' + '\nDatum: ' + str(df.index[0].strftime('%d.%m.%Y')) + ' bis ' + str(df.index[-1].strftime('%d.%m.%Y')), fontsize=font)
        plt.subplots_adjust(top=0.94)
        
        # -------------------------------------Energiebilanz-------------------------------------
        ax01.set_title(r'$\bf{Energiebilanz\ für\ den\ Tag}$', fontsize = font_title)
        a_Grid = round(df['INV_P_IN'].sum()*0.25,1)
        a_PVT = round(df['PVT_Pel'].sum()*0.25,1)
        a_PVT_s = round(df['s_PVT_Pel'].sum()*0.25/1000,1)
        a_WP = round(df.loc[df['WP_COP']>-0.2, 'WP_Pel'].sum()*0.25,1)
        a_WP_s = round(df['s_WP_Pel'].sum()*0.25/1000,1)
        a_user = round(df.loc[df['WP_COP']==0, 'Bat_Pel'].sum()*0.25,1)
        a_room = round(df['HM4_Qpkt'].sum()*0.25/-1000,1)
        a_room_s = round(df['TABS_Qpkt'].sum()*0.25/-1000,1)
        a_PVT_th = round(df['HM1_Qpkt'].sum()*0.25/-1000,1)
        a_PVT_th_s = round(df['PVT_Qpkt'].sum()*0.25/-1000,1)
        
        ax01.imshow(mpimg.imread(os.path.dirname(os.path.abspath(__file__)) + '/Energiebilanz.JPG'), extent=[0, 7.21, 0, 3.8])
        ax01.set_yticks([])
        ax01.set_xticks([])
        # Grid
        ax01.text(0.7, 1.7, str(a_Grid) + ' kWh',c = 'tab:blue', fontsize=legendfont)
        #WP
        a_WP = round(df.loc[df['WP_COP']>0, 'WP_Pel'].sum()*0.25,1)
        a_WP_s = round(df['s_WP_Pel'].sum()*0.25/1000,1)
        # User
        ax01.text(2.8+1, 1.7-1.17, str(a_user) + ' kWh',c = 'tab:blue', fontsize=legendfont)
        ax01.text(2.8+1, 1.45-1.17, str(round(0.15*len(df.loc[df['WP_COP']==0, :])*0.25,1)) + ' kWh',c = 'tab:red', fontsize=legendfont)
        # PVT_el
        ax01.text(1.2, 1.7+1, str(a_PVT) + ' kWh',c = 'tab:blue', fontsize=legendfont)
        ax01.text(1.2, 2.7-0.25, str(a_PVT_s) + ' kWh',c = 'tab:red', fontsize=legendfont)
        # PVT_th
        ax01.text(4.6, 3, str(a_PVT_th) + ' kWh',c = 'tab:blue', fontsize=legendfont)
        ax01.text(4.6, 3-0.25, str(a_PVT_th_s) + ' kWh',c = 'tab:red', fontsize=legendfont)
        # Room
        ax01.text(5.3, 2.05, str(a_room) + ' kWh',c = 'tab:blue', fontsize=legendfont)
        ax01.text(5.3, 2.05-0.25, str(a_room_s) + ' kWh',c = 'tab:red', fontsize=legendfont)
        # Room
        ax01.text(5.3, 2.05-0.745, str(0.0) + ' kWh',c = 'tab:blue', fontsize=legendfont)
        ax01.text(5.3, 2.05-0.745-0.25, str(0.0) + ' kWh',c = 'tab:red', fontsize=legendfont)
        # Legend
        ax01.text(0.1, 3.6, 'Messwert',c = 'tab:blue', fontsize=legendfont)
        ax01.text(0.1, 3.6-0.25, 'Simulationswert' ,c = 'tab:red', fontsize=legendfont)
            
        # -----------------------------------------------------------------------------------------------------------------
        # -------------------------------------Innenraum-------------------------------------------------------------------
        # -----------------------------------------------------------------------------------------------------------------
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
        
        # ------------------------------------------------Y1------------------------------------------------
        # Data
        # ax1.scatter(0, 0, label = '', color = 'white')
        ax02.scatter([i+0.5 for i in range(0,7,1)], Mess_Tmin, color = 'tab:blue', label = '   ', alpha = 0.5)
        ax02.plot([i+0.5 for i in range(0,7,1)], Mess_Tmin, color = 'tab:blue', alpha = 0.4, linestyle = '--')
        ax02.scatter([i+0.5 for i in range(0,7,1)], Mess_Tmax, color = 'tab:red', label = '   ', alpha = 0.5)
        ax02.plot([i+0.5 for i in range(0,7,1)], Mess_Tmax, color = 'tab:red', alpha = 0.4, linestyle = '--')
        ax02.scatter([i+0.5 for i in range(0,7,1)], Mess_Tmean, color = 'tab:green', label = '   ', alpha = 0.5)
        ax02.plot([i+0.5 for i in range(0,7,1)], Mess_Tmean, color = 'tab:green', alpha = 0.4, linestyle = '--')
        
        # ax1.scatter(0, 0, label = ' ', color = 'white')
        ax02.scatter([i+0.5 for i in range(0,7,1)], Sim_Tmin, color = 'tab:blue', label = '     Minimalwert', marker = 'x')
        ax02.plot([i+0.5 for i in range(0,7,1)], Sim_Tmin, color = 'tab:blue', alpha = 1, linestyle = '--')
        ax02.scatter([i+0.5 for i in range(0,7,1)], Sim_Tmax, color = 'tab:red', label = '     Maximalwert', marker = 'x')
        ax02.plot([i+0.5 for i in range(0,7,1)], Sim_Tmax, color = 'tab:red', alpha = 1, linestyle = '--')
        ax02.scatter([i+0.5 for i in range(0,7,1)], Sim_Tmean, color = 'tab:green', label = '     Mittelwert', marker = 'x')
        ax02.plot([i+0.5 for i in range(0,7,1)], Sim_Tmean, color = 'tab:green', alpha = 1, linestyle = '--')
        # labels
        ax02.set_ylabel('Temperatur in °C', fontsize = labelfont)
        # ax02.yaxis.set_major_formatter('{:.1f}'.format)
        ax02.set_yticks(np.arange(10,50, 5))
        ax02.set_yticks(np.arange(10,50, 2.5), minor = True)
        ax02.set_ylim(10,40)
        ax02.legend(fontsize = legendfont, loc = 'upper left', ncol = 2, title = 'Messwert    Simulation              ', title_fontsize = legendfont)
        # ------------------------------------------------Y2------------------------------------------------
        # Data
        ax02_twin.bar([i+0.35 for i in range(0,7,1)], Mess_Tschwankung, width = 0.3, color = 'tab:green', label = 'Messwert:   Schwankung', alpha = 0.3)
        ax02_twin.bar([i+0.65 for i in range(0,7,1)], Sim_Tschwankung, width = 0.3, color = 'tab:green', label =  'Simulation: Schwankung', alpha = 0.6)
        # labels
        ax02_twin.set_ylabel('Schwankung in K', fontsize = labelfont)
        # ax02_twin.yaxis.set_major_formatter('{:.1f}'.format)
        ax02_twin.set_yticks(range(0,12,2))
        ax02_twin.set_ylim(0,12)

        ax02_twin.legend(fontsize = legendfont, loc = 'upper right')
        # ------------------------------------------------X------------------------------------------------
        ax02.set_xlim(0,7)
        ax02.set_xticks(range(0,7,1), ['']*7)
        # ------------------------------------------------Global------------------------------------------------
        ax02.set_title(r'$\bf{Innentemeratur}$', fontsize = font_title)
        # ax.set_title(r'$\bf{First\ Row}$' + '\nSecond Row')
        ax02.grid(True)
        ax02.grid(which='minor', alpha = 0.3)
     
        
        
        
        
        # -----------------------------------------------------------------------------------------------------------------
        # -------------------------------------Kaltspeicher------------------------------------------------------------
        # -----------------------------------------------------------------------------------------------------------------
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
            Top_Mess_Tmin.append(day['Kalt Oben'].min())
            Top_Mess_Tmax.append(day['Kalt Oben'].max())
            Top_Mess_Tmean.append(day['Kalt Oben'].mean())
            Top_Mess_Tschwankung.append((Top_Mess_Tmax[-1] - Top_Mess_Tmin[-1])/2)
            Top_Sim_Tmin.append(min([posEl[-1] for posEl in day['KS_theta']]))
            Top_Sim_Tmax.append(max([posEl[-1] for posEl in day['KS_theta']]))
            Top_Sim_Tmean.append(sum([posEl[-1] for posEl in day['KS_theta']])/len([posEl[-1] for posEl in day['KS_theta']]))
            Top_Sim_Tschwankung.append((Top_Sim_Tmax[-1] - Top_Sim_Tmin[-1])/2)
            abweichung = [abs(a-b) for a,b in zip(day['Kalt Oben'], [posEl[-1] for posEl in day['KS_theta']])]
            Top_MAE.append((sum(abweichung) / len(abweichung)))
        
            Bot_Mess_Tmin.append(day['Kalt Unten'].min())
            Bot_Mess_Tmax.append(day['Kalt Unten'].max())
            Bot_Mess_Tmean.append(day['Kalt Unten'].mean())
            Bot_Mess_Tschwankung.append((Bot_Mess_Tmax[-1] - Bot_Mess_Tmin[-1])/2)
            Bot_Sim_Tmin.append(min([posEl[0] for posEl in day['KS_theta']]))
            Bot_Sim_Tmax.append(max([posEl[0] for posEl in day['KS_theta']]))
            Bot_Sim_Tmean.append(sum([posEl[0] for posEl in day['KS_theta']])/len([posEl[0] for posEl in day['KS_theta']]))
            Bot_Sim_Tschwankung.append((Bot_Sim_Tmax[-1] - Bot_Sim_Tmin[-1])/2)
            abweichung_bot = [abs(a-b) for a,b in zip(day['Kalt Unten'], [posEl[0] for posEl in day['KS_theta']])]
            Bot_MAE.append((sum(abweichung_bot) / len(abweichung_bot)))
        

        # ------------------------------------------------Y1------------------------------------------------
        # Data
        ax12.plot([i+0.5 for i in range(0,7,1)], Top_Mess_Tmin, color = 'tab:blue', label = '   ', marker = 'o', markersize = 5 ,alpha = 0.5, linestyle = '--')
        ax12.plot([i+0.5 for i in range(0,7,1)], Top_Mess_Tmax, color = 'tab:red', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )
        ax12.plot([i+0.5 for i in range(0,7,1)], Top_Mess_Tmean, color = 'tab:green', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )
        
        ax12.plot([i+0.5 for i in range(0,7,1)], Top_Sim_Tmin, color = 'tab:blue', alpha = 1, linestyle = '--', label = 'Minimalwert', marker = 'x', markersize = 5)
        ax12.plot([i+0.5 for i in range(0,7,1)], Top_Sim_Tmax, color = 'tab:red', alpha = 1, linestyle = '--', label = 'Maximalwert', marker = 'x', markersize = 5)
        ax12.plot([i+0.5 for i in range(0,7,1)], Top_Sim_Tmean, color = 'tab:green', alpha = 1, linestyle = '--', label = 'Mittelwert', marker = 'x', markersize = 5)
        
        # labels
        ax12.set_ylabel('Temperatur in °C', fontsize = labelfont)
        ax12.set_yticks(np.arange(0,30,2.5))
        ax12.set_yticks(np.arange(0,30,0.5), minor = True)
        ax12.set_ylim(7.5,27.5)
        ax12_twin.yaxis.set_major_formatter('{:.1f}'.format)
        ax12.legend(fontsize = legendfont, loc = 'upper left', ncol = 2, title = 'Messwert    Simulation              ', title_fontsize = legendfont)
        # ------------------------------------------------Y1_Twin------------------------------------------------
        # Data
        ax12_twin.bar([i+0.35 for i in range(0,7,1)], Top_Mess_Tschwankung, width = 0.3, color = 'tab:green', label = 'Messwert:   Schwankung', alpha = 0.3)
        ax12_twin.bar([i+0.65 for i in range(0,7,1)], Top_Sim_Tschwankung, width = 0.3, color = 'tab:green', label =  'Simulation: Schwankung', alpha = 0.6)
        ax12_twin.plot([i+0.5 for i in range(0,7,1)], Top_MAE, color = 'grey', alpha = 1, linestyle = '--', label = 'MAE', marker = 'x', markersize = 5)
        
        # labels
        ax12_twin.set_ylabel('Schwankung und MAE in K', fontsize = labelfont)
        ax12_twin.yaxis.set_major_formatter('{:.1f}'.format)
        ax12_twin.set_ylim(0,20)
        # ax1_twin.set_yticks(range(0,30,2))
        ax12_twin.legend(fontsize = legendfont, loc = 'upper right')
        # ------------------------------------------------X------------------------------------------------
        ax12.set_xlim(0,7)
        ax12.set_xticks(range(0,7,1), ['']*7, rotation = 60)
        # ------------------------------------------------Global------------------------------------------------
        ax12.set_title(r'$\bf{Kältespeicher:\ Oberste\ Schicht}$', fontsize = font_title)
        ax12.grid(True)
        

        # ------------------------------------------------Unterste Schicht------------------------------------------------
        # Data
        ax22.plot([i+0.5 for i in range(0,7,1)], Bot_Mess_Tmin, color = 'tab:blue', label = '   ', marker = 'o', markersize = 5 ,alpha = 0.5, linestyle = '--')
        ax22.plot([i+0.5 for i in range(0,7,1)], Bot_Mess_Tmax, color = 'tab:red', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )
        ax22.plot([i+0.5 for i in range(0,7,1)], Bot_Mess_Tmean, color = 'tab:green', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )
        
        ax22.plot([i+0.5 for i in range(0,7,1)], Bot_Sim_Tmin, color = 'tab:blue', alpha = 1, linestyle = '--', label = 'Minimalwert', marker = 'x', markersize = 5)
        ax22.plot([i+0.5 for i in range(0,7,1)], Bot_Sim_Tmax, color = 'tab:red', alpha = 1, linestyle = '--', label = 'Maximalwert', marker = 'x', markersize = 5)
        ax22.plot([i+0.5 for i in range(0,7,1)], Bot_Sim_Tmean, color = 'tab:green', alpha = 1, linestyle = '--', label = 'Mittelwert', marker = 'x', markersize = 5)
        
        # labels
        ax22.set_ylabel('Temperatur in °C', fontsize = labelfont)
        ax22.set_yticks(np.arange(0,35,2.5))
        ax22.set_yticks(np.arange(0,35,0.5), minor = True)
        ax22.set_ylim(5,30)
        # ax22.yaxis.set_major_formatter('{:.1f}'.format)
        ax22.legend(fontsize = legendfont, loc = 'upper left', ncol = 2, title = 'Messwert    Simulation              ', title_fontsize = legendfont)
        # ------------------------------------------------Y2_Twin------------------------------------------------
        # Data
        ax22_twin.bar([i+0.35 for i in range(0,7,1)], Bot_Mess_Tschwankung, width = 0.3, color = 'tab:green', label = 'Messwert:   Schwankung', alpha = 0.3)
        ax22_twin.bar([i+0.65 for i in range(0,7,1)], Bot_Sim_Tschwankung, width = 0.3, color = 'tab:green', label =  'Simulation: Schwankung', alpha = 0.6)
        ax22_twin.plot([i+0.5 for i in range(0,7,1)], Bot_MAE, color = 'grey', alpha = 1, linestyle = '--', label = 'MAE', marker = 'x', markersize = 5)
        
        # labels
        ax22_twin.set_ylabel('Schwankung und MAE in K', fontsize = labelfont)
        # ax22_twin.yaxis.set_major_formatter('{:.1f}'.format)
        ax22_twin.set_yticks(np.arange(0,20,3))
        ax22_twin.set_ylim(0,15)
        ax22_twin.legend(fontsize = legendfont, loc = 'upper right')
        # ------------------------------------------------X------------------------------------------------
        ax22.set_xlim(0,7)
        ax22.set_xticks(range(0,7,1), ['']*7, rotation = 60)
        # ------------------------------------------------Global------------------------------------------------
        ax22.set_title(r'$\bf{Kältespeicher: \ Unterste\ Schicht}$', fontsize = font_title)
        
        ax22.grid(True)
        
        E_Top = [b-a for a,b in zip(df['Kalt Oben'], [posEl[-1] for posEl in df['KS_theta']])]
        E_Bot = [b-a for a,b in zip(df['Kalt Unten'], [posEl[0] for posEl in df['KS_theta']])]
        
        
        # -----------------------------------------------------------------------------------------------------------------
        # -------------------------------------Wärmespeicher------------------------------------------------------------
        # -----------------------------------------------------------------------------------------------------------------
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
        

        # fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [0.5, 0.5]}, figsize=(6.4, 4.8*2))
        # ax1_twin = ax1.twinx()
        # ax2_twin = ax2.twinx()
        # ------------------------------------------------Y1------------------------------------------------
        # Data
        ax11.plot([i+0.5 for i in range(0,7,1)], Top_Mess_Tmin, color = 'tab:blue', label = '   ', marker = 'o', markersize = 5 ,alpha = 0.5, linestyle = '--')
        ax11.plot([i+0.5 for i in range(0,7,1)], Top_Mess_Tmax, color = 'tab:red', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )
        ax11.plot([i+0.5 for i in range(0,7,1)], Top_Mess_Tmean, color = 'tab:green', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )
        
        ax11.plot([i+0.5 for i in range(0,7,1)], Top_Sim_Tmin, color = 'tab:blue', alpha = 1, linestyle = '--', label = 'Minimalwert', marker = 'x', markersize = 5)
        ax11.plot([i+0.5 for i in range(0,7,1)], Top_Sim_Tmax, color = 'tab:red', alpha = 1, linestyle = '--', label = 'Maximalwert', marker = 'x', markersize = 5)
        ax11.plot([i+0.5 for i in range(0,7,1)], Top_Sim_Tmean, color = 'tab:green', alpha = 1, linestyle = '--', label = 'Mittelwert', marker = 'x', markersize = 5)
        
        # labels
        ax11.set_ylabel('Temperatur in °C', fontsize = labelfont)
        ax11.set_yticks(np.arange(10,80,5))
        ax11.set_yticks(np.arange(10,80,2.5), minor = True)
        ax11.set_ylim(15,70)
        # ax1_twin.yaxis.set_major_formatter('{:.1f}'.format)
        ax11.legend(fontsize = legendfont, loc = 'upper left', ncol = 2, title = 'Messwert    Simulation              ', title_fontsize = legendfont)
        # ------------------------------------------------Y1_Twin------------------------------------------------
        # Data
        ax11_twin.bar([i+0.35 for i in range(0,7,1)], Top_Mess_Tschwankung, width = 0.3, color = 'tab:green', label = 'Messwert:   Schwankung', alpha = 0.3)
        ax11_twin.bar([i+0.65 for i in range(0,7,1)], Top_Sim_Tschwankung, width = 0.3, color = 'tab:green', label =  'Simulation: Schwankung', alpha = 0.6)
        ax11_twin.plot([i+0.5 for i in range(0,7,1)], Top_MAE, color = 'grey', alpha = 1, linestyle = '--', label = 'MAE', marker = 'x', markersize = 5)
        
        # labels
        ax11_twin.set_ylabel('Schwankung und MAE in K', fontsize = labelfont)
        # ax11_twin.yaxis.set_major_formatter('{:.1f}'.format)
        ax11_twin.set_yticks(range(0,44,4))
        ax11_twin.set_ylim(0,44)
        ax11_twin.legend(fontsize = legendfont, loc = 'upper right')
        # ------------------------------------------------X------------------------------------------------
        ax11.set_xlim(0,7)
        ax11.set_xticks(range(0,7,1), ['']*7, rotation = 60)
        # ------------------------------------------------Global------------------------------------------------
        ax11.set_title(r'$\bf{Wärmespeicher: \ Oberste\ Schicht}$', fontsize = font_title)
        ax11.grid(True)


        # ------------------------------------------------Y2------------------------------------------------
        # Data
        ax21.plot([i+0.5 for i in range(0,7,1)], Bot_Mess_Tmin, color = 'tab:blue', label = '   ', marker = 'o', markersize = 5 ,alpha = 0.5, linestyle = '--')
        ax21.plot([i+0.5 for i in range(0,7,1)], Bot_Mess_Tmax, color = 'tab:red', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )
        ax21.plot([i+0.5 for i in range(0,7,1)], Bot_Mess_Tmean, color = 'tab:green', label = '   ', alpha = 0.4, linestyle = '--', marker = 'o', markersize = 5 )
        
        ax21.plot([i+0.5 for i in range(0,7,1)], Bot_Sim_Tmin, color = 'tab:blue', alpha = 1, linestyle = '--', label = 'Minimalwert', marker = 'x', markersize = 5)
        ax21.plot([i+0.5 for i in range(0,7,1)], Bot_Sim_Tmax, color = 'tab:red', alpha = 1, linestyle = '--', label = 'Maximalwert', marker = 'x', markersize = 5)
        ax21.plot([i+0.5 for i in range(0,7,1)], Bot_Sim_Tmean, color = 'tab:green', alpha = 1, linestyle = '--', label = 'Mittelwert', marker = 'x', markersize = 5)
        
        # labels
        ax21.set_ylabel('Temperatur in °C', fontsize = labelfont)
        ax21.set_yticks(np.arange(10,80,5))
        ax21.set_yticks(np.arange(10,80,2.5), minor = True)
        ax21.set_ylim(15,70)
        # ax2_twin.yaxis.set_major_formatter('{:.1f}'.format)
        ax21.legend(fontsize = legendfont, loc = 'upper left', ncol = 2, title = 'Messwert    Simulation              ', title_fontsize = legendfont)
        # ------------------------------------------------Y2_Twin------------------------------------------------
        # Data
        ax21_twin.bar([i+0.35 for i in range(0,7,1)], Bot_Mess_Tschwankung, width = 0.3, color = 'tab:green', label = 'Messwert:   Schwankung', alpha = 0.3)
        ax21_twin.bar([i+0.65 for i in range(0,7,1)], Bot_Sim_Tschwankung, width = 0.3, color = 'tab:green', label =  'Simulation: Schwankung', alpha = 0.6)
        ax21_twin.plot([i+0.5 for i in range(0,7,1)], Bot_MAE, color = 'grey', alpha = 1, linestyle = '--', label = 'MAE', marker = 'x', markersize = 5)
        
        # labels
        ax21_twin.set_ylabel('Schwankung und MAE in K', fontsize = labelfont)
        # a12_twin.yaxis.set_major_formatter('{:.1f}'.format)

        ax21_twin.set_yticks(range(0,33,3))
        ax21_twin.set_ylim(0,33)
        ax21_twin.legend(fontsize = legendfont, loc = 'upper right')
        # ------------------------------------------------X------------------------------------------------
        ax21.set_xlim(0,7)
        ax21.set_xticks(range(0,7,1), ['']*7, rotation = 60)
        # ------------------------------------------------Global------------------------------------------------
        ax21.set_title(r'$\bf{Wärmespeicher\ Unterste\ Schicht}$', fontsize = font_title)
        ax21.grid(True)
        E_Top = [b-a for a,b in zip(df['Warm Oben'], [posEl[-1] for posEl in df['WS_theta']])]
        E_Bot = [b-a for a,b in zip(df['Warm Unten'][::4], [posEl[0] for posEl in df['WS_theta'][::4]])]



        
        # -----------------------------------------------------------------------------------------------------------------
        # -------------------------------------TABS------------------------------------------------------------
        # -----------------------------------------------------------------------------------------------------------------
        Q_mess = []
        Q_sim = []
        for i in df.index.day.unique():
            aa = df[df.index.day == i]
            Q_mess.append(abs(aa['HM4_Qpkt'].sum()*0.25/-1000))
            Q_sim.append(abs(aa['TABS_Qpkt'].sum()*0.25/-1000))

        # ------------------------------------------------Y1------------------------------------------------
        # Data
        ax31.bar([i-0.2 for i in range(1,8,1)], Q_mess, width = 0.4, label = 'Messwert')
        ax31.bar([i+0.2 for i in range(1,8,1)], Q_sim, width = 0.4, label = 'Simulationswert')
        # labels
        ax31.set_ylabel('Kältemenge in -kWh', fontsize = labelfont)
        # ax1.yaxis.set_major_formatter('{:.1f}'.format)
        ax31.set_yticks(range(0,30,2), minor = True)
        ax31.set_yticks(range(0,30,4))
        ax31.set_ylim(0,28)
        # ------------------------------------------------X------------------------------------------------
        ax31.set_xticks(range(1,8,1), ['']*7, rotation = 60)
        # ------------------------------------------------Global------------------------------------------------
        ax31.set_title(r'$\bf{Thermoaktives\ Bauteilsystem}$',fontsize = font_title)
        ax31.grid(True)
        ax31.grid(which='minor', alpha = 0.3)
        ax31.legend(fontsize = 8, loc = 'upper left')
        ax31.text(5.45, 20.5,  r'$\bf{Gesamte\ Woche:}$' + '\n' + str(round(df['HM4_Qpkt'].sum()*0.25/1000,2)) + ' kWh gemessen'
                  + '\n' + str(round(df['TABS_Qpkt'].sum()*0.25/1000,2)) + ' kWh simuliert' + '\n' +
                  str(round((df['TABS_Qpkt'].sum()*0.25-df['HM4_Qpkt'].sum()*0.25)/(df['HM4_Qpkt'].sum()*0.25)*100,2)) + ' % Abweichung',
                  fontsize = labelfont, bbox=dict(facecolor='white', edgecolor='black'))

        
        
        # -----------------------------------------------------------------------------------------------------------------
        # -------------------------------------WP------------------------------------------------------------
        # -----------------------------------------------------------------------------------------------------------------
        P_mess = []
        P_sim = []
        
        for i in df.index.day.unique():
            aa = df[df.index.day == i]
            P_mess.append(aa.loc[aa['WP_COP']>0, 'Bat_Pel'].sum()*0.25)
            P_sim.append(aa['s_WP_Pel'].sum()*0.25/1000)     
            

            
            # WP
        ax01.text(2.8, 1.7, str(a_WP) + ' kWh',c = 'tab:blue', fontsize=legendfont)
        ax01.text(2.8, 1.45, str(a_WP_s) + ' kWh',c = 'tab:red', fontsize=legendfont)
        # ------------------------------------------------Y1------------------------------------------------
        # Data
        ax32.bar([i-0.2 for i in range(1,8,1)], P_mess, width = 0.4, label = 'Messwert')
        ax32.bar([i+0.2 for i in range(1,8,1)], P_sim, width = 0.4, label = 'Simulationswert')
        # labels
        ax32.set_ylabel('el. Leistungsaufnahme in kWh', fontsize = labelfont)
        # ax1.yaxis.set_major_formatter('{:.1f}'.format)
        ax32.set_yticks(range(0,15,1), minor = True)
        ax32.set_yticks(range(0,15,2))
        ax32.set_ylim(0,15)
        
        # ------------------------------------------------X------------------------------------------------
        ax32.set_xticks(range(1,8,1), ['']*7)
        
        # ax1.set_xticks([i / 2 for i in range(0,7*2,1)], minor = True)
        
        # ------------------------------------------------Global------------------------------------------------
        ax32.set_title(r'$\bf{Wärmepumpe:\ Leistungsaufnahme}$', fontsize = font_title)
        
        ax32.grid(True)
        ax32.grid(which='minor', alpha = 0.3)
        ax32.legend(fontsize = legendfont, loc = 'upper left')
        ax32.text(5.57, 11,  r'$\bf{Gesamte\ Woche:}$' + '\n' + str(round(df['WP_Pel'].sum()*0.25,2)) + ' kWh gemessen'
                  + '\n' + str(round(df['s_WP_Pel'].sum()*0.25/1000,2)) + ' kWh simuliert' + '\n' +
                  str(round((df['s_WP_Pel'].sum()*0.25-df['WP_Pel'].sum()*0.25*1000)/(df['WP_Pel'].sum()*0.25*1000)*100,2)) + ' % Abweichung',
                  fontsize = labelfont, bbox=dict(facecolor='white', edgecolor='black'))

        
        
        # -----------------------------------------------------------------------------------------------------------------
        # -------------------------------------PVT-Thermisch------------------------------------------------------------
        # -----------------------------------------------------------------------------------------------------------------
        Q_mess = []
        Q_sim = []
        for i in df.index.day.unique():
            aa = df[df.index.day == i]
            Q_mess.append(abs(aa['HM1_Qpkt'].sum()*0.25/-1000))
            Q_sim.append(abs(aa['PVT_Qpkt'].sum()*0.25/-1000))
        # ------------------------------------------------Y1------------------------------------------------
        # Data
        ax41.bar([i-0.2 for i in range(1,8,1)], Q_mess, width = 0.4, label = 'Messwert')
        ax41.bar([i+0.2 for i in range(1,8,1)], Q_sim, width = 0.4, label = 'Simulationswert')
        # labels
        ax41.set_ylabel('Wärmemenge in kWh', fontsize = labelfont)
        # ax1.yaxis.set_major_formatter('{:.1f}'.format)
        ax41.set_yticks(np.arange(0,45,2.5), minor = True)
        ax41.set_yticks(range(0,45,5))
        ax41.set_ylim(0,35)
        
        # ------------------------------------------------X------------------------------------------------
        ax41.set_xticks(range(1,8,1), ['Montag', 'Dienstag', 'Mittwoch', 'Donnerstag', 'Freitag', 'Samstag', 'Sonntag'], fontsize = labelfont, rotation = 60)
        
        # ax1.set_xticks([i / 2 for i in range(0,7*2,1)], minor = True)
        
        # ------------------------------------------------Global------------------------------------------------
        ax41.set_title(r'$\bf{PVT-Kollektoren:\ thermische\ Seite}$', fontsize = font_title)
        
        ax41.grid(True)
        ax41.grid(which='minor', alpha = 0.3)
        ax41.legend(fontsize = 8, loc = 'upper left')
        
        ax41.text(5.4, 25.4,  r'$\bf{Gesamte\ Woche:}$' + '\n' + str(round(df['HM1_Qpkt'].sum()*0.25/1000,2)) + ' kWh gemessen'
                  + '\n' + str(round(df['PVT_Qpkt'].sum()*0.25/1000,2)) + ' kWh simuliert' + '\n' +
                  str(round((df['PVT_Qpkt'].sum()*0.25-df['HM1_Qpkt'].sum()*0.25)/(df['HM1_Qpkt'].sum()*0.25)*100,2)) + ' % Abweichung',
                  fontsize = labelfont, bbox=dict(facecolor='white', edgecolor='black'))

        
        # -----------------------------------------------------------------------------------------------------------------
        # -------------------------------------PVT-el. Seite------------------------------------------------------------
        # -----------------------------------------------------------------------------------------------------------------
        Q_mess = []
        Q_sim = []
        
        for i in df.index.day.unique():
            aa = df[df.index.day == i]
            Q_mess.append(aa['PVT_Pel'].sum()*0.25)
            Q_sim.append(aa['s_PVT_Pel'].sum()*0.25/1000)

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
        ax42.set_title(r'$\bf{PVT-Kollektoren:\ elektrische\ Seite}$', fontsize = font_title)
        
        ax42.grid(True)
        ax42.grid(which='minor', alpha = 0.3)
        ax42.legend(fontsize = 8, loc = 'upper left')

        ax42.text(5.6, 7.2,  r'$\bf{Gesamte\ Woche:}$' + '\n' + str(round(df['PVT_Pel'].sum()*0.25,2)) + ' kWh gemessen'
                  + '\n' + str(round(df['s_PVT_Pel'].sum()*0.25/1000,2)) + ' Wh simuliert' + '\n' +
                  str(round((df['s_PVT_Pel'].sum()*0.25/1000-df['PVT_Pel'].sum()*0.25)/(df['PVT_Pel'].sum()*0.25)*100,2)) + ' % Abweichung',
                  fontsize = labelfont, bbox=dict(facecolor='white', edgecolor='black'))
        
        plt.subplots_adjust(wspace=0.22)
        plt.show()
        
        fig.savefig(os.path.dirname(os.path.abspath(__file__)) + '/Grafik/KW_' + str(df.index[0].isocalendar()[1]) + '.png')

    
    
    
    
    
    
    
