# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 15:13:39 2023

@author: 49157
"""
#%% Libraries
import math
# # #%%
# param = {}
# param['B1R1C'] = {
#     'UA' : 179.2,
#     'C_m' : 16500000,
#     }

#%% Modell

def B1R1C(theta_prev, theta_e, i_sol, phi_HC, phi_int, param):
    theta_i = theta_e + ((i_sol+phi_HC+phi_int)/param['UA']) + (theta_prev-theta_e-((i_sol+phi_HC+phi_int)/param['UA']))*math.exp(-60/(param['C_m']/param['UA']))
    # print(param['C_m'])
    return theta_i

#%% Testing


# GB_theta_i = [10]
# days = 10
# theta_e = [10]*24 + [20]*(24*days-24)

# for time in range(0,24*days,1):
#     # print(time)
#     GB = B1R1C(theta_e = theta_e[time],
#                 i_sol = 0,
#                 phi_HC = 0,
#                 phi_int = 0,
#                 theta_prev = GB_theta_i[-1],
#                 param = param['B1R1C']
#                 )
#     GB_theta_i.append(GB)

# import matplotlib.pyplot as plt
# plt.plot(range(0,len(GB_theta_i)),GB_theta_i)
