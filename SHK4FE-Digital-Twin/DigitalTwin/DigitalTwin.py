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

#%%
os.chdir(os.path.dirname(os.path.abspath(__file__)))
#%% Modelle
def BAT(BAT_Pel_in, BAT_Pel_out, BAT_SOC_prev):
    #                        _______________
    # Pel_in      [W]       |               |
    # --------------------->|               |
    # Pel_out     [W]       |               | SOC
    # --------------------->|     BAT       |--------------------->
    # SOC[-1]     [%]       |               | 
    # --------------------->|               |
    #                       |_______________|
    C_prev = BAT_SOC_prev*param['BAT']['C']*3600
    C = (BAT_Pel_in-BAT_Pel_out)*900 + C_prev
    if C<0:
        P_Grid = -C/15/60
        C = (BAT_Pel_in-BAT_Pel_out + P_Grid)*900 + C_prev
    else:
        P_Grid = 0
    BAT_SOC = C/(param['BAT']['C']*3600)
    return BAT_SOC, P_Grid
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
    if str(T_prim_in) == 'nan':
        T_prim_in = 20
    if str(T_sek_in) == 'nan':
        T_sek_in = 20
    if str(Vpkt_prim) == 'nan':
        Vpkt_prim = 0
    if str(Vpkt_sek) == 'nan':
        Vpkt_sek = 0

    
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
                    -(param['KS']['k']*100 * A_ext * (Ti[0] - Tamb))                # Wärmeverluste an Umgebung
                    + (mpkt_i * prop['wasser']['cp'] * (Ti[0+1] - Ti[0]))       # Stofftransport
                    +((A_i*param['KS']['lambda_eff']/z_i) * (Ti[0+1] - Ti[0]))  # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt

        else:                                                                   # Resultierender Massenstrom ↑    
            dTdt = (0
                    -(mpkt_sek * prop['wasser']['cp'] * (Ti[0] - T_sek_in))     # Sek: mpkt-bedingte Änderung
                    -(param['KS']['k']*100*A_ext*(Ti[0] - Tamb))                    # Wärmeverluste an Umgebung
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
                    -(param['KS']['k']*100*A_ext*(Ti[-1] - Tamb))                   # Wärmeverluste an Umgebung 
                    +0 # keine Schicht oberhalb vorhanden                       # Stofftransport
                    +((A_i*param['KS']['lambda_eff']/z_i)*(Ti[-1-1] - Ti[-1]))  # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt  
        
        else:                                                                   # Resultierender Massenstrom ↑   
            dTdt = (0
                    +(mpkt_prim * prop['wasser']['cp'] * (T_prim_in - Ti[-1]))  # Prim. mpkt-bedingte Änderung
                    -(param['KS']['k']*100*A_ext*(Ti[-1] - Tamb))                   # Wärmeverluste an Umgebung 
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
                    -(param['KS']['k']*100*A_ext*(Ti[i] - Tamb))                    # Wärmeverluste an Umgebung
                    +(mpkt_i*prop['wasser']['cp']*(Ti[i+1] - Ti[i]))            # Stofftransport
                    +((A_i*param['KS']['lambda_eff'])/z_i * (Ti[i+1] - 2*Ti[i] + Ti[i-1])) # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt
   
        else:                                                                   # Resultierender Massenstrom ↑ 
            dTdt = (0
                    +0 # kein Anschluss in diesen Layern                        # Sek/prim mpkt-bedingte Änderung
                    -(param['KS']['k']*100*A_ext*(Ti[i] - Tamb))                    # Wärmeverluste an Umgebung
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
    if str(Vpkt_prim) == 'nan':
        Vpkt_prim = 0
    # -------------------Parameter-Massenströme-Randbedingungen----------------------------------------------------
    z_i = param['WS']['H']/param['WS']['N']                                     # (25) Höhe einer Schicht
    A_ext = np.pi*param['WS']['D']*z_i                                          # (26) Außenfläche einer Schicht
    A_i = np.pi*(param['WS']['D']-2*param['WS']['th']*10)**2/4                     # (27) Fläche zwischen zwei Layern
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
                    -(param['WS']['k']*100 * A_ext * (Ti[0] - Tamb))                # Wärmeverluste an Umgebung
                    + (mpkt_i * prop['wasser']['cp'] * (Ti[0+1] - Ti[0]))       # Stofftransport
                    +((A_i*param['WS']['lambda_eff']/z_i) * (Ti[0+1] - Ti[0]))  # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt

        else:                                                                   # Resultierender Massenstrom ↑    
            dTdt = (0
                    -(mpkt_sek * prop['wasser']['cp'] * (Ti[0] - T_sek_in))     # Sek: mpkt-bedingte Änderung
                    -(param['WS']['k']*100*A_ext*(Ti[0] - Tamb))                    # Wärmeverluste an Umgebung
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
                    -(param['WS']['k']*100*A_ext*(Ti[-1] - Tamb))                   # Wärmeverluste an Umgebung 
                    +0 # keine Schicht oberhalb vorhanden                       # Stofftransport
                    +((A_i*param['WS']['lambda_eff']/z_i)*(Ti[-1-1] - Ti[-1]))  # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt  
        
        else:                                                                   # Resultierender Massenstrom ↑   
            dTdt = (0
                    +(mpkt_prim * prop['wasser']['cp'] * (T_prim_in - Ti[-1]))  # Prim. mpkt-bedingte Änderung
                    -(param['WS']['k']*100*A_ext*(Ti[-1] - Tamb))                   # Wärmeverluste an Umgebung 
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
                    -(param['WS']['k']*100*A_ext*(Ti[i] - Tamb))                    # Wärmeverluste an Umgebung
                    +(mpkt_i*prop['wasser']['cp']*(Ti[i+1] - Ti[i]))            # Stofftransport
                    +((A_i*param['WS']['lambda_eff'])/z_i * (Ti[i+1] - 2*Ti[i] + Ti[i-1])) # Wärmeleitung zw. Layer
                    ) / (m_i*prop['wasser']['cp'])                              # aus dQ/dt

        else:                                                                   # Resultierender Massenstrom ↑ 
            dTdt = (0
                    +0 # kein Anschluss in diesen Layern                        # Sek/prim mpkt-bedingte Änderung
                    -(param['WS']['k']*100*A_ext*(Ti[i] - Tamb))                    # Wärmeverluste an Umgebung
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
        TABS_theta_in = np.nan
    TABS_theta_in = TABS_theta_in
    TABS_Vpkt = TABS_Vpkt
    return TABS_theta_in, TABS_theta_out, TABS_Vpkt, TABS_Qpkt
def PVT(theta_mean_prev, PVT_theta_in, PVT_Vpkt, theta_e, BAT_SOC, G, WP_Pel):
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
        PVT_Vpkt = 0
    A = param['PVT']['etha_opt']*G                                               
    F = (param['PVT']['c1']+param['PVT']['c5']/900)*0.5
    E = PVT_Vpkt/3600*prop['Antifrogen']['rho']*prop['Antifrogen']['cp'] / (param['PVT']['A'] * param['PVT']['quantity']) 
    if (PVT_Vpkt == 0 or str(PVT_theta_in) == 'nan'):
        PVT_theta_in = theta_mean_prev
    PVT_theta_out = (A+theta_e*(param['PVT']['c1']) + param['PVT']['c5']/900 * theta_mean_prev - PVT_theta_in*(F-E))/(E+F) 
    PVT_Qpkt = PVT_Vpkt/3600 *prop['Antifrogen']['rho'] * prop['Antifrogen']['cp'] * (PVT_theta_out-PVT_theta_in)                    
    PVT_theta_cell = PVT_Qpkt / (param['PVT']['A']*param['PVT']['quantity'])/param['PVT']['U'] + (PVT_theta_in+PVT_theta_out)/2
    if BAT_SOC == 1:
        PVT_Pel = 0   
    else: 
        PVT_Pel = G*param['PVT']['etha_el']*(1-param['PVT']['beta_el']*(PVT_theta_cell-25))*(param['PVT']['A']*param['PVT']['quantity'])   
    if PVT_Pel >= 1040:
        PVT_Pel = 1040
    PVT_max = (0 if str(WP_Pel) == 'nan' else WP_Pel) + (param['BAT']['C']*3600 - BAT_SOC*param['BAT']['C']*3600)/900
    if PVT_Pel > PVT_max:
        PVT_Pel = PVT_max
    
    
    
    return PVT_Pel,PVT_theta_in,PVT_theta_out,PVT_Qpkt,PVT_theta_cell, PVT_Vpkt
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
    # Sek_Vpkt      [m³/h]  |               | PSek_theta_out [°C]
    # --------------------->|_______________|--------------------->
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
        if WP_COP > 8:
            WP_COP = 8
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
    if str(B_phi_HC) == 'nan':
        B_phi_HC = 0
    else:
        B_phi_HC = B_phi_HC
    A_m = param['B']['A_f']*param['B']['Coe_Am']
    A_tot = 4.5 * param['B']['A_f']
    c_m = param['B']['Coe_Cm']*param['B']['A_f']*2
    H_tr_w = sum([U*A for U,A in zip(param['B']['transparent components']['U'], param['B']['transparent components']['A'])])
    if V_Vdot == 0:
        H_ve = 0.00000001
    else:
        # H_ve = 1200*V_Vdot/3600
        H_ve = (1-0.85)*1200*V_Vdot/3600   
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
        I_umg = (I_diff + I_dir) * 0.5 * 0.2 * (1-math.cos(math.radians(tilt)))
        I_d = I_diff * 0.5 * (1 + math.cos(math.radians(tilt)))
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
def radprocessor(i_dir, i_dif, sun_h, sun_az, pos_h, pos_az):
# i_dir   -->|            |
# i_dif   -->|            | 
# sun_h   -->| Radiation  |--> i_sol
# sun_az  -->| Processor  |
# pos_h   -->|            | 
# pos_az  -->|            |
    # Required packages: math
    # Based on the height of the sun and the positioning of a component, indicates the resulting radiation on the component 
    # sun_h -> Sun hight in °
    # sun_az -> Sun azimuth in °
    # i_dir -> Direct irradiation on horizontal plane 
    # i_dif -> Diffuse irradiation on horizontal plane
    # pos_az -> azimuthal orientation of the component, 0° -> north, 90° -> east, 180° -> south, 270° -> west
    # pos_h -> alignment in axis of sun elevation, 0° -> horizontal up, 90°-> vertical, 180° -> horizontal down    
    
    # Direct radiation
    # Angle of Incidence cos(zeta)
    aoi = math.sin(math.radians(sun_h))*math.cos(math.radians(pos_h)) + math.cos(math.radians(sun_h))*math.sin(math.radians(pos_h)) * math.cos(math.radians(abs(pos_az-sun_az)))
    if aoi < 0:
        aoi = 0
    if sun_h<=5:
        i_f = 0
    else:
        i_f = i_dir/math.sin(math.radians(sun_h))*aoi
    # Strahlung durch Bodenreflexion
    i_umg = (i_dir+i_dif)*0.5*0.2*(1-math.cos(math.radians(pos_h)))
    # Diffuse Strahlung
    i_d = i_dif * 0.5 *(1+math.cos(math.radians(pos_h)))
    # Summe
    i_sol = i_f + i_umg + i_d
    return i_sol
def CON_winter(t_MA, WP_OI_prev, WS_theta_top, WS_theta_bot, TABS_OI_prev, B_theta_air):
    # Set Point Hot Tank
    if t_MA >= param['CON']['WS_theta_e']:
        WS_theta_set = param['CON']['WS_theta_min']
    else:
        WS_theta_set =  param['CON']['WS_slope']*t_MA + (param['CON']['WS_theta_min']-param['CON']['WS_slope']*param['CON']['WS_theta_e'])
    if WS_theta_set > param['CON']['WS_theta_max']:
        WS_theta_set = param['CON']['WS_theta_max']
        
    # Hysterese der Wärmepumpe
    if WP_OI_prev == 0 and WS_theta_top <= WS_theta_set - param['CON']['WP_Hyst']:
    # if WP_OI_prev == 0 and WS_theta_top <= WS_theta_set:
        WP_OI = 1
    elif WP_OI_prev == 1 and WS_theta_bot >= WS_theta_set: 
    # elif WP_OI_prev == 1 and WS_theta_bot >= WS_theta_set + param['CON']['WP_Hyst']: 
        WP_OI = 0
    else:
        WP_OI = WP_OI_prev
        
    # Hysterese TABS
    if TABS_OI_prev == 0 and B_theta_air <= 20:
        TABS_OI = 1
    elif TABS_OI_prev == 1 and B_theta_air > 23: 
        TABS_OI = 0
    else:
        TABS_OI = TABS_OI_prev
        
    return WS_theta_set, WP_OI, TABS_OI
def CON_sommer(t_MA, WP_OI_prev, KS_theta_top, KS_theta_bot, TABS_OI_prev, B_theta_air, PVT_ReC_prev, HS_theta_top, HS_theta_bot):
    # set pont cold tank
    if t_MA <= param['CON']['KS_theta_e']:
        KS_theta_set = param['CON']['KS_theta_max']
    else:
        KS_theta_set =  param['CON']['KS_slope']*t_MA + (param['CON']['KS_theta_max']-param['CON']['KS_slope']*param['CON']['KS_theta_e'])
    if KS_theta_set < param['CON']['KS_theta_min']:
        KS_theta_set = param['CON']['KS_theta_min']
    
    # Hysterese der Wärmepumpe
    if WP_OI_prev == 0 and KS_theta_bot > KS_theta_set + param['CON']['WP_Hyst']:
        WP_OI = 1
    elif WP_OI_prev == 1 and KS_theta_top <= KS_theta_set: 
        WP_OI = 0
    else:
        WP_OI = WP_OI_prev    
        
    # Hysterese TABS
    if TABS_OI_prev == 0 and B_theta_air >= 26:
        TABS_OI = 1
    elif TABS_OI_prev == 1 and B_theta_air < 23: 
        TABS_OI = 0
    else:
        TABS_OI = TABS_OI_prev
        
    # Recooling PVT
    if PVT_ReC_prev == 0 and  HS_theta_bot >= param['CON']['WS_theta_max']+5:
        PVT_ReC = 1
    elif PVT_ReC_prev == 1 and HS_theta_top < param['CON']['WS_theta_max']:
        PVT_ReC = 0
    else:
        PVT_ReC = PVT_ReC_prev
    
    return KS_theta_set, WP_OI, TABS_OI, PVT_ReC
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
param['KS'] = {
    'H'             : 1.2,                                                      # m         Speicherhöhe
    'N'             : 6,                                                        # -         Anzahl Schichten
    'D'             : 0.6,                                                      # m         Speicherdurchmesser
    'th'            : 0.02,                                                     # m²        Wanddicke
    'k'             : 0.002    ,                                                 # W/m²K     Wärmetransferkoeffizient
    'lambda_eff'    : 1.5                                                       # W/mK      vertikale Wärmeleitfähigkeit von Wasser     
    }
param['WS'] = {
    'H'             : 1.2,                                                      # m         Speicherhöhe
    'N'             : 4,                                                        # -         Anzahl Schichten
    'D'             : 0.6,                                                      # m         Speicherdurchmesser
    'th'            : 0.002,                                                     # m²        Wanddicke
    'k'             : 0.002,                                                     # W/m²K     Wärmetransferkoeffizient
    'lambda_eff'    : 1.5                                                       # W/mK      vertikale Wärmeleitfähigkeit von Wasser   
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
param['TABS'] = {
    'A'         : 15,                                                           # Area of floor
    'da'        : 0.03,                                                         # outside diameter of pipe
    'dp'        : 0.001,                                                        # thickness of pipe
    'dx'        : 0.125,                                                        # Laying distance of pipes
    'lambda_p'  : 340,                                                          # Heat transfer pipe (copper)
    'lambda_s'  : 2,                                                            # Heat transfer floor
    'd'         : 0.025,                                                         # Laying depth
    }
param['WP'] = {
    'etha'  :   0.47,
    }
param['B'] = {
    'A_f': 14.8,                                                                # m²        conditioned area
    'Coe_Am': 2.5,                                                              # -         coefficient for the determination of mass-related area according to table 12
    'Coe_Cm': 120000,                                                           # J/K*m²    coefficiant for the determination of Internal heat storage capacity according to table 12
    'transparent components': {
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
        'shading'    : [1,      1,      1,      1,      0.1,    1],             # -         1 for no shading, 0 for absolutely shadowed
        'R'          : [0.05,   0.05,   0.05,   0.05,   0.05,   0.05],          #           product of surface heat transfer resistance(ISO 6946) and absorption coefficient for solar radiation on opaque component, "Baehr, H.D.; Stephan, K. ; Wärme- und Stoffübertragung; Auszug Kap. 5.5 Strahlungsaustausch, aus Tab. 5.8,S.633" https://www.baunetzwissen.de/bauphysik/fachwissen/waermeschutz/absorption-und-waerme-auf-oberflaechen-4733315
        'U'          : [0.6,    0.6,    0.6,    0.6,    0.6,    0.6],           # W/m²/K    Heat transfer coefficient determined according to ISO 6946
        'A'          : [5.3,    14.6,   6.3,    15.7,   14.8,   14.8],          # m²        Area of the component
        'direction'  : [66,     156,    246,    336,    270,    270],           # °         0° -> north, 90° -> east, 180° ->south, 270° -> west
        'tilt'       : [90,     90,     90,     90,     0,      180]            # °         0° horizontal top, 90° -> vertical, 180° -> horizontal down
        },
    }
param['BAT'] = {
    'C' : 8000                                                                  # Wh
    }

param['CON'] = {
    'WP_Hyst':5,                # deviation between set and current state of the storage at which HP starts
    # Winter
    'WS_theta_max' : 45,        # see characteristic curve of the storage temperatures 
    'WS_slope': -1,           # see characteristic curve of the storage temperatures 
    'WS_theta_min': 25,         # see characteristic curve of the storage temperatures 
    'WS_theta_e' : 16,          # see characteristic curve of the storage temperatures 

    # Summer
    'KS_theta_min' : 7,         # see characteristic curve of the storage temperatures 
    'KS_slope': -0.6,           # see characteristic curve of the storage temperatures 
    'KS_theta_max': 10,         # see characteristic curve of the storage temperatures 
    'KS_theta_e' : 16,          # see characteristic curve of the storage temperatures 
    }

#%% Data
df = pd.read_csv('21_wetter_HSO.csv', index_col = 0)
df.index = pd.to_datetime(df.index, format='%Y-%m-%d %H:%M:%S')  

# Direkt und Diffus
df['I_dif'] = [(1/(1 + 2.71828182846**(7.997*(G/1000 - 0.586))))*G for G in df['G']]
df['I_dir'] = [G-dif for G,dif in zip(df['G'],df['I_dif'])]

# Moving average der Außentemperatur
MA_fenster = 24                                                                #h
df['t_MA'] = df['t'].rolling(window=MA_fenster*4).mean()
df['t_MA'] = [next(x for x in df['t_MA'] if not pd.isna(x)) if pd.isna(x) else x for x in df['t_MA']]

# Winter/Sommer
df['case'] = ['Winter' if x<=param['CON']['WS_theta_e'] else ('Sommer' if x >= param['CON']['KS_theta_e'] else np.nan) for x in df['t_MA']]

# Sonnenposition
df['sun_h'] = [SunPos(lon = 7.9498017, lat = 48.473451, year = i.year, month = i.month, 
                      day = i.day, hour = i.hour, minute = i.minute, shift = 0)[1] for i in df.index]
df['sun_az'] = [SunPos(lon = 7.9498017, lat = 48.473451, year = i.year, month = i.month, 
                      day = i.day, hour = i.hour, minute = i.minute, shift = 0)[0] for i in df.index]
    
del MA_fenster

#%% Simulation

# Startwerte
CON_WS_set, CON_KS_set, CON_WP_OI, CON_TABS_OI, CON_PVT_ReC = [],[], [0], [0], [0]

WP_COP, WP_Pel, WP_Qdot_prim,WP_Qdot_sec, WP_prim_theta_in, WP_prim_theta_out, WP_prim_Vpkt, WP_sek_theta_in, WP_sek_theta_out, WP_sek_Vpkt = [],[],[],[],[],[],[],[],[],[]

WS_theta = [np.linspace(22, 22, param['WS']['N']).tolist()]

KS_theta = [np.linspace(20, 20, param['KS']['N']).tolist()]

PVT_Pel, PVT_theta_in, PVT_theta_out, PVT_Qpkt, PVT_theta_cell, PVT_G, PVT_Vpkt = [0],[df.loc[df.index[0], 't']],[df.loc[df.index[0], 't']],[0], [df.loc[df.index[0], 't']], [0], [0]

B_theta_m, B_theta_s, B_theta_air, B_theta_operativ = [25], [25], [25], [25]

TABS_theta_in, TABS_theta_out, TABS_Vpkt, TABS_Qpkt = [], [], [], []

BAT_SOC = [1]
P_Grid = [0]

von = 0
bis = von + (24*4)*365
# von = (24*4)*30*7
# bis = von + (24*4)
# von = 0
# bis = von + (24*4)*7
for i in df.index[von:bis]:
    print(str(i.month)+ '.' + str(i.day))
    if df.loc[i, 'case'] == 'Winter':
        # Controller
        conW = CON_winter(t_MA = df.loc[i, 't_MA'],
                          WP_OI_prev = CON_WP_OI[-1], 
                          WS_theta_top = WS_theta[-1][-1],
                          WS_theta_bot = WS_theta[-1][1],
                          TABS_OI_prev = CON_TABS_OI[-1],
                          B_theta_air = B_theta_air[-1])
        CON_WS_set.append(conW[0])
        CON_WP_OI.append(conW[1])
        CON_TABS_OI.append(conW[2])
        CON_KS_set.append(np.nan)
        CON_PVT_ReC.append(0)
        
        # Wärmepumpe
        wp = WP(WP_prim_theta_in = PVT_theta_out[-1], 
                WP_sek_theta_in = WS_theta[-1][0], 
                WP_sek_delta_theta = 5, 
                WP_prim_Vpkt = 0.5, 
                WP_sek_Vpkt = 0.5, 
                WP_OI = CON_WP_OI[-1])
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
        
        # PVT
        G_POA = radprocessor(i_dir = df.loc[i, 'I_dir'], i_dif  = df.loc[i, 'I_dif'], 
                             sun_h = df.loc[i, 'sun_h'], sun_az  = df.loc[i, 'sun_az'], 
                             pos_h = 35, pos_az = 156)
        PVT_G.append(G_POA)
        pvt = PVT(theta_mean_prev = (PVT_theta_in[-1]+PVT_theta_out[-1])/2, 
                  PVT_theta_in = WP_prim_theta_out[-1], 
                  PVT_Vpkt = WP_prim_Vpkt[-1], 
                  theta_e = df.loc[i, 't'], 
                  BAT_SOC = BAT_SOC[-1], 
                  G = G_POA,
                  WP_Pel = WP_Pel[-1])
        PVT_Pel.append(pvt[0])
        PVT_theta_in.append(pvt[1])
        PVT_theta_out.append(pvt[2])
        PVT_Qpkt.append(pvt[3])
        PVT_theta_cell.append(pvt[4])
        PVT_Vpkt.append(pvt[5])
        
        # TABS
        tabs = TABS(TABS_theta_in = WS_theta[-1][-1], 
                    TABS_Vpkt = CON_TABS_OI[-1] * 0.4,
                    B_theta_air = B_theta_air[-1])
        TABS_theta_in.append(tabs[0])
        TABS_theta_out.append(tabs[1])
        TABS_Vpkt.append(tabs[2])
        TABS_Qpkt.append(tabs[3])
        
        # Wärmespeicher
        ws = WS(Vpkt_prim = WP_sek_Vpkt[-1], 
                T_prim_in = WP_sek_theta_out[-1], 
                Vpkt_sek = TABS_Vpkt[-1],
                T_sek_in = TABS_theta_out[-1], 
                T_vektor_prev = WS_theta[-1], 
                Tamb = B_theta_air[-1])
        WS_theta.append(ws)
        
        # Kältespeicher
        ks = KS(Vpkt_prim       = 0, 
                T_prim_in       = 0, 
                Vpkt_sek        = 0, 
                T_sek_in        = 0,
                T_vektor_prev   = KS_theta[-1],
                Tamb            = B_theta_air[-1])
        KS_theta.append(ks)
        
        # Gebäude
        b = B(B_phi_int = 90, 
              B_phi_HC = TABS_Qpkt[-1], 
              B_theta_m_prev = B_theta_m[-1], 
              theta_e = df.loc[i, 't'], 
              sun_h = df.loc[i, 'sun_h'], 
              sun_az = df.loc[i, 'sun_az'], 
              I_dir = df.loc[i, 'I_dir'], 
              I_diff = df.loc[i, 'I_dif'], 
              V_Vdot = 30, 
              V_theta_sup = df.loc[i, 't'])
        B_theta_m.append(b[0])
        B_theta_s.append(b[1])
        B_theta_air.append(b[2])
        B_theta_operativ.append(b[3])
        
        # Batterie
        bat = BAT(BAT_Pel_in = (0 if str(PVT_Pel[-1]) == 'nan' else PVT_Pel[-1]), 
              BAT_Pel_out = (0 if str(WP_Pel[-1]) == 'nan' else WP_Pel[-1]),  
              BAT_SOC_prev = BAT_SOC[-1])
        BAT_SOC.append(bat[0])
        P_Grid.append(bat[1])
        
    elif df.loc[i, 'case'] == 'Sommer':       
        conS = CON_sommer(t_MA = df.loc[i, 't_MA'],
                          WP_OI_prev = CON_WP_OI[-1], 
                          KS_theta_top = KS_theta[-1][-1],
                          KS_theta_bot = KS_theta[-1][0],
                          TABS_OI_prev =  CON_TABS_OI[-1],
                          B_theta_air = B_theta_air[-1],
                          PVT_ReC_prev = CON_PVT_ReC[-1], 
                          HS_theta_top = WS_theta[-1][-1],
                          HS_theta_bot = WS_theta[-1][0])
        
        CON_WS_set.append(np.nan)
        CON_KS_set.append(conS[0])
        CON_WP_OI.append(conS[1])
        CON_TABS_OI.append(conS[2])
        CON_PVT_ReC.append(conS[3])
        
        wp = WP(WP_prim_theta_in = KS_theta[-1][-1], 
                WP_sek_theta_in = WS_theta[-1][0], 
                WP_sek_delta_theta = 5, 
                WP_prim_Vpkt = 0.5, 
                WP_sek_Vpkt = 0.5, 
                WP_OI = CON_WP_OI[-1])
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
        
        # PVT
        G_POA = radprocessor(i_dir = df.loc[i, 'I_dir'], i_dif  = df.loc[i, 'I_dif'], 
                             sun_h = df.loc[i, 'sun_h'], sun_az  = df.loc[i, 'sun_az'], 
                             pos_h = 35, pos_az = 156)
        PVT_G.append(G_POA)
        pvt = PVT(theta_mean_prev = (PVT_theta_in[-1]+PVT_theta_out[-1])/2, 
                  PVT_theta_in = WS_theta[-1][-1], 
                  PVT_Vpkt = CON_PVT_ReC[-1]*0.5, 
                  theta_e = df.loc[i, 't'], 
                  BAT_SOC = BAT_SOC[-1], 
                  G = G_POA,
                  WP_Pel = WP_Pel[-1])
        PVT_Pel.append(pvt[0])
        PVT_theta_in.append(pvt[1])
        PVT_theta_out.append(pvt[2])
        PVT_Qpkt.append(pvt[3])
        PVT_theta_cell.append(pvt[4])
        PVT_Vpkt.append(pvt[5])

        # TABS
        tabs = TABS(TABS_theta_in = KS_theta[-1][0], 
                    TABS_Vpkt = CON_TABS_OI[-1] * 0.4,
                    B_theta_air = B_theta_air[-1])
        TABS_theta_in.append(tabs[0])
        TABS_theta_out.append(tabs[1])
        TABS_Vpkt.append(tabs[2])
        TABS_Qpkt.append(tabs[3])

        # Kältespeicher
        ks = KS(Vpkt_prim       = TABS_Vpkt[-1], 
                T_prim_in       = TABS_theta_out[-1], 
                Vpkt_sek        = WP_prim_Vpkt[-1], 
                T_sek_in        = WP_prim_theta_out[-1],
                T_vektor_prev   = KS_theta[-1],
                Tamb            = B_theta_air[-1])
        KS_theta.append(ks)

        # Wärmespeicher
        ws = WS(Vpkt_prim = WP_sek_Vpkt[-1], 
                T_prim_in = WP_sek_theta_out[-1], 
                Vpkt_sek = PVT_Vpkt[-1],
                T_sek_in = PVT_theta_out[-1], 
                T_vektor_prev = WS_theta[-1], 
                Tamb = B_theta_air[-1])
        WS_theta.append(ws)

        b = B(B_phi_int = 90, 
              B_phi_HC = TABS_Qpkt[-1], 
              B_theta_m_prev = B_theta_m[-1], 
              theta_e = df.loc[i, 't'], 
              sun_h = df.loc[i, 'sun_h'], 
              sun_az = df.loc[i, 'sun_az'], 
              I_dir = df.loc[i, 'I_dir'], 
              I_diff = df.loc[i, 'I_dif'], 
              V_Vdot = 30, 
              V_theta_sup = df.loc[i, 't'])
        B_theta_m.append(b[0])
        B_theta_s.append(b[1])
        B_theta_air.append(b[2])
        B_theta_operativ.append(b[3])
        
        # Batterie
        bat = BAT(BAT_Pel_in = (0 if str(PVT_Pel[-1]) == 'nan' else PVT_Pel[-1]), 
              BAT_Pel_out = (0 if str(WP_Pel[-1]) == 'nan' else WP_Pel[-1]),  
              BAT_SOC_prev = BAT_SOC[-1])
        BAT_SOC.append(bat[0])
        P_Grid.append(bat[1])

a = pd.DataFrame({'G' : df['G'][von: bis],
                  'I_dif' : df['I_dif'][von: bis],
                  'I_dir' : df['I_dir'][von: bis],
                  't' : df['t'][von: bis],
                  't_MA' : df['t_MA'][von: bis],
                  'case' : df['case'][von: bis],
                  'sun_h' : df['sun_h'][von: bis],
                  'sun_az' : df['sun_az'][von: bis],
                  'CON_WS_set': CON_WS_set,
                   'CON_KS_set': CON_KS_set,
                   'CON_WP_OI': CON_WP_OI[1:],
                   'CON_PVT_ReC': CON_PVT_ReC[1:],
                   'WP_COP': WP_COP,
                   'WP_Pel': WP_Pel,
                   'WP_Qdot_prim': WP_Qdot_prim,
                   'WP_Qdot_sec': WP_Qdot_sec,
                   'WP_prim_theta_in': WP_prim_theta_in,
                   'WP_prim_theta_out': WP_prim_theta_out,
                   'WP_prim_Vpkt': WP_prim_Vpkt,
                   'WP_sek_theta_in': WP_sek_theta_in,
                   'WP_sek_theta_out': WP_sek_theta_out,
                   'WP_sek_Vpkt': WP_sek_Vpkt,
                   'WS_theta': WS_theta[1:],
                   'KS_theta': KS_theta[1:],
                    'PVT_G' : PVT_G[1:],
                    'PVT_Pel' : PVT_Pel[1:],
                    'PVT_theta_in' : PVT_theta_in[1:],
                    'PVT_theta_out' : PVT_theta_out[1:],
                    'PVT_Qpkt' : PVT_Qpkt[1:],
                   'PVT_theta_cell' : PVT_theta_cell[1:],
                   'PVT_Vpkt' : PVT_Vpkt[1:],
                   'CON_TABS_OI': CON_TABS_OI[1:],
                   'TABS_theta_in': TABS_theta_in[:],
                   'TABS_theta_out': TABS_theta_out[:],
                   'TABS_Vpkt': TABS_Vpkt[:],
                   'TABS_Qpkt': TABS_Qpkt[:],
                   'B_theta_m': B_theta_m[1:],
                   'B_theta_s': B_theta_s[1:],
                   'B_theta_air': B_theta_air[1:],
                   'B_theta_operativ': B_theta_operativ[1:],
                   'BAT_SOC': BAT_SOC[1:],
                   'P_Grid': P_Grid[1:],
                  }, index=df.index[von: bis])

del CON_WS_set, CON_WP_OI, CON_TABS_OI, WP_COP, WP_Pel, WP_Qdot_prim, WP_Qdot_sec, WP_prim_theta_in, WP_prim_theta_out
del WP_sek_theta_in, WP_sek_theta_out, WP_prim_Vpkt, WP_sek_Vpkt, PVT_Pel, PVT_theta_in, PVT_theta_out, PVT_Qpkt, PVT_theta_cell
del TABS_theta_in, TABS_theta_out, TABS_Vpkt, TABS_Qpkt, WS_theta, B_theta_m, B_theta_s, B_theta_air, B_theta_operativ
del conW, G_POA, b, bis, i, pvt, PVT_G, von, wp, ws, tabs, bat, BAT_SOC, P_Grid, PVT_Vpkt, ks, conS, CON_PVT_ReC, CON_KS_set
#%% Gebäude und TABS
colors_map = {'Sommer': 'tab:blue', 'Winter': 'tab:red', np.nan: 'gray'}
plt.scatter(a.index, a['B_theta_air'], label = 'Air', s = 1, c = [colors_map[val] for val in a['case']])
# plt.scatter(a.index, a['t'], label = 'Amb', s = 1, c = 'tab:green')
# plt.plot(a.index, a['TABS_theta_in'], label = 'Tabs_in')
# plt.plot(a.index, a['TABS_theta_out'], label = 'Tabs_out')
# plt.plot(a.index, a['TABS_theta_out'], label = 'Air')
plt.ylim(0,30)
plt.legend()
plt.grid()
plt.show()
#%%
too_low = sum([20-i for i in a[a['B_theta_air']<20]['B_theta_air']])/4
too_high = sum([i-26 for i in a[a['B_theta_air']>26]['B_theta_air']])/4
print(too_low+too_high)
#%% Batterie
plt.plot(a.index, a['BAT_SOC'], label = 'WP_Pel')
plt.plot(a.index, [i/1000 for i in a['WP_Pel']], label = 'WP_Pel')
plt.plot(a.index, [i/1000 for i in a['PVT_Pel']], label = 'PVT_Pel')
plt.plot(a.index, [i/1000 for i in a['P_Grid']], label = 'P_Grid', alpha = 0.5)
plt.legend()
plt.ylim(0,2)
plt.xticks(rotation = 90)
plt.grid()
#%% PVT

plt.plot(a.index, a['PVT_theta_in'],label = 'In')
plt.plot(a.index, a['PVT_theta_out'],label = 'Out')
plt.plot(a.index, a['PVT_Vpkt'],label = 'Vpkt')
plt.legend()
plt.grid()
plt.show()

#%% Recooling WS
for i in range(param['WS']['N']):
        plt.plot(a.index, [p[i] for p in a['WS_theta']], label=[str(a) for a in list(range(param['WS']['N']))][i], color = plt.cm.RdYlBu(np.linspace(0, 1, param['WS']['N']))[::-1][i])

plt.plot(a.index, a['CON_PVT_ReC'])

plt.plot([a.index[0], a.index[-1]], [param['CON']['WS_theta_max']+5, param['CON']['WS_theta_max']+5], c = 'black')
plt.plot([a.index[0], a.index[-1]], [param['CON']['WS_theta_max'], param['CON']['WS_theta_max']], c = 'black')

# plt.ylim(0,30)
plt.legend()
plt.grid()
plt.show()


#%% Speicher
for i in range(param['KS']['N']):
        plt.plot(a.index, [p[i] for p in a['KS_theta']], label=[str(a) for a in list(range(param['KS']['N']))][i], color = plt.cm.RdYlBu(np.linspace(0, 1, param['KS']['N']))[::-1][i])
for i in range(param['WS']['N']):
        plt.plot(a.index, [p[i] for p in a['WS_theta']], label=[str(a) for a in list(range(param['WS']['N']))][i], color = plt.cm.RdYlBu(np.linspace(0, 1, param['WS']['N']))[::-1][i])


plt.plot(a.index, a['CON_KS_set'], c = 'black', linestyle = '--')
plt.plot(a.index, [i + 5 for i in a['CON_KS_set']], c = 'black', linestyle = '--')
plt.plot(a.index, a['CON_WS_set'], c = 'black', linestyle = '--')
plt.plot(a.index, [i + 5 for i in a['CON_WS_set']], c = 'black', linestyle = '--')
plt.plot(a.index, a['CON_WP_OI'], c = 'black')

# plt.ylim(0,30)
# plt.legend()
plt.grid()
plt.show()





#%%
colors_map = {'Sommer': 'tab:blue', 'Winter': 'tab:red', np.nan: 'gray'}
plt.scatter(a['t_MA'], a['B_theta_air'], label = 'BAT_SOC', s = 1, c = [colors_map[val] for val in a['case']])
# plt.plot(a.index, [i/1000 for i in a['WP_Pel']], label = 'WP_Pel')
# plt.plot(a.index, [i/1000 for i in a['PVT_Pel']], label = 'PVT_Pel')
# plt.plot(a.index, [i/1000 for i in a['P_Grid']], label = 'P_Grid', alpha = 0.5)
plt.legend()
plt.xticks(rotation = 90)
plt.grid()

#%%
a['Uhrzeit'] = [datetime.strptime(c, '%H:%M') for c in [i.strftime("%H:%M") for i in a.index]]
plt.scatter(a['Uhrzeit'], a['WP_COP'], s = 2)
plt.xticks(rotation = 90)


#%% Durchschnittlicher COP je Uhrzeit
COP_perTime = []
for i in a['Uhrzeit'].unique():
    COP_perTime.append(a.loc[(a['Uhrzeit'] == i) & (a['WP_COP'] > 0), :]['WP_COP'].mean())
plt.ylim(3,6)
plt.plot(range(0,96), COP_perTime)

#%% Laufzeit der Wärmepumpe je Uhrzeit (15 Minuten)
RunTime = []
for i in a['Uhrzeit'].unique():
    RunTime.append(len(a.loc[(a['Uhrzeit'] == i) & (a['WP_COP'] > 0), :]['WP_COP'])*0.25)
plt.ylim(0,25)
plt.plot(range(0,96), RunTime)


















