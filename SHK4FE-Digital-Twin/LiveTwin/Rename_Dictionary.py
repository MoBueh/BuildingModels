# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 14:42:43 2023

@author: MoBueh
"""

rename = {
    'WS.DWD-----MEA.T-'                     : 'AMB_T',                          # Temperature from Weatherstation
    'STO.EL-INV---OUT-MEA.POW.EL-'          : 'Bat_Pel',                        # Electircal Power @ Inverter Output
    'GEN.HC-HP--R410a--CALC.POW.H-'         : 'WP_Qpkt',                        # Thermal Power of Heatpump
    'DIST.HC--MTR.HC-W.HC-RET-MEA.VF-'      : 'HM4_Vpkt',                       # Volume flow in l/h @ Heatmeter 4
    'STO.H---W.H-BOT-MEA.T-'                : 'WS_T',                           # Temperature @ Hot Tank
    'GEN.H-HP--W.H-RET-MEA.T-'              : 'WP_sek_T_RET',                   # Return Temperature at hot side of Heatpump
    'ODC----VER-MEA.GLIR-CALC'              : 'Pyr_ver',                        # Global Irradiation @ tilted Pyranometer
    'GEN.H--MTR.H-W.H-RET-MEA.VF-'          : 'HM10_Vpkt',                      # Volume flow in l/h @ Heatmeter 10
    'STO.EL-SCC---IN-MEA.CUR-'              : 'Bat_I',                          # Current from VarioTrack -> PVT but @ around 55V 
    'DIST.HC--MTR.HC-W.HC-RET-MEA.T-'       : 'HM4_T_RET',                      # Return Temperature @ Heatmeter 4
    'STO.EL-INV---OUT-MEA.VOLT-'            : 'INV_U_OUT',                      # Voltage @ Inverter Output
    'STO.EL-INV----MEA.VOLT-'               : 'BAT_U',                          # Voltage @ Inverter Input -> Battery
    'STO.EL-SCC----MEA.VOLT-'               : 'PVT_U_SCC',                      # Voltage @ SCC -> BAT
    'STO.EL-SCC----CALC.POW.EL-'            : 'PVT_Pel',                        # Electircall Power @ PVT-Cells
    'GEN.HC-SCOL-MTR.HC-BRINE-RET-CALC.T-'  : 'HM1_T_RET',                      # Return Temperature @ Heatmeter 1
    'ZONE---AIR.ID--MEA.T-'                 : 'Zone_T',                         # Inside Air Temperature
    'ODC---AIR.OD--CALC.T-'                 : 'AMB_T2',                         # Temperature measured from Sensor
    'GEN.H--MTR.H-W.H-SUP-CALC.T-'          : 'HM10_T_SUP',                     # Supply Temperature @ Heatmeter 10
    'STO.EL-INV---OUT-MEA.CUR-'             : 'INV_I_OUT',                      # Current @ Inverter Output
    # 'GEN.EL-PV----CALC.POW.EL-'             : 'PVT_Pel',                        # Electrical Power @ PVT
    'STO.EL-INV---IN-MEA.CUR-'              : 'INV_I_IN',                       # Current @ Inverter Input -> from Grid
    'DIST.HC--MTR.HC-W.HC-SUP-MEA.T-'       : 'HM4_T_SUP',                      # Volume flow in l/h @ Heatmeter 4,
    'GEN.C-HP--BRINE-SUP-CALC.T-'           : 'WP_prim_T_SUP',                  # Supply Temperature @ cold side of Heatpump
    'GEN.EL-PV----MEA.VOLT-'                : 'PVT_U',                          # Voltage @ PVT
    'STO.EL-INV---IN-MEA.VOLT-'             : 'INV_U_IN',                       # Voltage @ Inverter Input -> from Grid
    'GEN.HC-SCOL-MTR.HC-BRINE-SUP-CALC.T-'  : 'HM1_T_SUP',                      # Supply Temperature @ Heatmeter 1
    'STO.C---W.C-BOT-MEA.T-'                : 'KS_T',                           # Temperature @ Cold Tank 
    'GEN.HC-SCOL-MTR.HC-BRINE-RET-MEA.VF-'  : 'HM1_Vpkt',                       # Volume flow in l/h @ Heatmeter 1
    'GEN.HC-HP--R410a--CALC.COP.virt-'      : 'WP_COP',                         # COP of Heatpump
    'GEN.EL-PV----MEA.CUR-'                 : 'PVT_I',                          # Current in PVT
    'GEN.H-HP--W.H-SUP-MEA.T-'              : 'WP_sek_T_SUP',                   # Supply Temperature at hot side of Heatpump
    'GEN.C-HP--BRINE-RET-CALC.T-'           : 'WP_prim_T_RET',                  # Return Temperature @ cold side of Heatpump
    'STO.EL-INV---IN-MEA.POW.EL-'           : 'INV_P_IN',                       # Electirc Power @ Inverter Input -> from Grid
    'ODC----HOR-MEA.GLIR-CALC'              : 'Pyr_hor',                        # Global Irradiation @ horizontal Pyranometer
    'GEN.HC-HP.COMP--R410a--CALC.POW.EL-'   : 'WP_Pel',                         # Electircal ower consumption of Heatpump
    'GEN.H--MTR.H-W.H-RET-MEA.T-'           : 'HM10_T_RET',                     # Return Temperature @ Heatmeter 10
    'Warm_unten------'                      : 'Warm Unten',
    'Warm_oben------'                       : 'Warm Oben',
    'Kalt_unten------'                      : 'Kalt Unten',
    'Kalt_oben------'                       : 'Kalt Oben',
    }
