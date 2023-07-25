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
        C = 0
    else:
        P_Grid = 0
    if C >=param['BAT']['C']*3600:
        C = param['BAT']['C']*3600
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
    # if str(Vpkt_prim) == 'nan':
    #     Vpkt_prim = 20
    if str(T_sek_in) == 'nan':
        T_sek_in = 20
    # if str(T_sek_in) == 'nan':
    #     T_sek_in = 20

    
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
    'k'             : 0.02    ,                                                 # W/m²K     Wärmetransferkoeffizient
    'lambda_eff'    : 1.5                                                       # W/mK      vertikale Wärmeleitfähigkeit von Wasser     
    }
param['WS'] = {
    'H'             : 1.2,                                                      # m         Speicherhöhe
    'N'             : 4,                                                        # -         Anzahl Schichten
    'D'             : 0.6,                                                      # m         Speicherdurchmesser
    'th'            : 0.002,                                                     # m²        Wanddicke
    'k'             : 0.02,                                                     # W/m²K     Wärmetransferkoeffizient
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
    'WS_slope': -0.7,           # see characteristic curve of the storage temperatures 
    'WS_theta_min': 20,         # see characteristic curve of the storage temperatures 
    'WS_theta_e' : 16,          # see characteristic curve of the storage temperatures 

    # Summer
    'KS_theta_min' : 7,         # see characteristic curve of the storage temperatures 
    'KS_slope': -0.6,           # see characteristic curve of the storage temperatures 
    'KS_theta_max': 15,         # see characteristic curve of the storage temperatures 
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
def CON_winter(t_MA, WP_OI_prev, WS_theta_top, WS_theta_bot, TABS_OI_prev, B_theta_air):
    # Set Point Hot Tank
    if t_MA >= param['CON']['WS_theta_e']:
        WS_theta_set = param['CON']['WS_theta_min']
    else:
        WS_theta_set =  param['CON']['WS_slope']*t_MA + (param['CON']['WS_theta_min']-param['CON']['WS_slope']*param['CON']['WS_theta_e'])
    if WS_theta_set > param['CON']['WS_theta_max']:
        WS_theta_set = param['CON']['WS_theta_max']
        
    # Hysterese der Wärmepumpe
    # if WP_OI_prev == 0 and WS_theta_top <= WS_theta_set - param['CON']['WP_Hyst']:
    if WP_OI_prev == 0 and WS_theta_top <= WS_theta_set:
        WP_OI = 1
    # elif WP_OI_prev == 1 and WS_theta_bot >= WS_theta_set: 
    elif WP_OI_prev == 1 and WS_theta_bot >= WS_theta_set + param['CON']['WP_Hyst']: 
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

# Startwerte
CON_WS_set, CON_WP_OI, CON_TABS_OI = [], [0], [0]

WP_COP, WP_Pel, WP_Qdot_prim,WP_Qdot_sec, WP_prim_theta_in, WP_prim_theta_out, WP_prim_Vpkt, WP_sek_theta_in, WP_sek_theta_out, WP_sek_Vpkt = [],[],[],[],[],[],[],[],[],[]

WS_theta = [np.linspace(20, 20, param['WS']['N']).tolist()]

PVT_Pel, PVT_theta_in, PVT_theta_out, PVT_Qpkt, PVT_theta_cell, PVT_G = [0],[df.loc[df.index[0], 't']],[df.loc[df.index[0], 't']],[0], [df.loc[df.index[0], 't']], [0]

B_theta_m, B_theta_s, B_theta_air, B_theta_operativ = [25], [25], [25], [25]

TABS_theta_in, TABS_theta_out, TABS_Vpkt, TABS_Qpkt = [], [], [], []

BAT_SOC = [1]
P_Grid = [0]

von = 0
bis = von + (24*4)*365
for i in df.index[von:bis]:
    b = B(B_phi_int = 90, 
          B_phi_HC = 0, 
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
   
a = pd.DataFrame({'G' : df['G'][von: bis],
                  'I_dif' : df['I_dif'][von: bis],
                  'I_dir' : df['I_dir'][von: bis],
                  't' : df['t'][von: bis],
                  't_MA' : df['t_MA'][von: bis],
                  'case' : df['case'][von: bis],
                  'sun_h' : df['sun_h'][von: bis],
                  'sun_az' : df['sun_az'][von: bis],
                  'B_theta_m': B_theta_m[1:],
                  'B_theta_s': B_theta_s[1:],
                  'B_theta_air': B_theta_air[1:],
                  'B_theta_operativ': B_theta_operativ[1:],

                  }, index=df.index[von: bis])




#%%
colors_map = {'Sommer': 'blue', 'Winter': 'red', np.nan: 'gray'}

plt.scatter(a['t_MA'], a['B_theta_air'], label = 'BAT_SOC', s = 1, c = [colors_map[val] for val in a['case']])
plt.plot([-10, 30],[26, 26], label = 'BAT_SOC')
plt.plot([-10, 30],[20, 20], label = 'BAT_SOC')
# plt.plot(a.index, [i/1000 for i in a['WP_Pel']], label = 'WP_Pel')
# plt.plot(a.index, [i/1000 for i in a['PVT_Pel']], label = 'PVT_Pel')
# plt.plot(a.index, [i/1000 for i in a['P_Grid']], label = 'P_Grid', alpha = 0.5)
plt.legend()
plt.xticks(rotation = 90)
plt.grid()










