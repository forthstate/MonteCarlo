#!/usr/bin/python3
"""
visualisation code file ISING3D.PY
choose and uncomment plots necessary:
fig1 : output data from python3 file
fig2 : output data from fortran file 
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#Plot of E/M vs MCS step
"""
#---Load csv first Data file
data = pd.read_csv('ising3d_data.csv')

#---Figure initialisation E/M plot

fig1, axs1 = plt.subplots(1,2,figsize=(20, 12))

axs1[0].plot(data['Time'].values,data['E'].values,label='latticeEnergy')
axs1[0].plot(data['Time'].values,data['M'].values,label='latticeMagnetisation')
axs1[0].set_xlabel('MCS-step')
axs1[0].set_ylabel('E/M')
axs1[0].legend()

axs1[1].plot(data['Time'].values,data['E/N'].values,label='Energy/spin')
axs1[1].plot(data['Time'].values,data['M/N'].values,label='Magnetisation/spin')
axs1[1].set_xlabel('MCS-step')
axs1[1].set_ylabel('E/M-per-Spin')
axs1[1].legend()
"""
#Plot Cv vs Temperature plot
#---Load csv second Data file
data_cv = pd.read_csv('ising3d_Cv_data.csv')

#---Figure initialisation Cv plot
fig2, axs2 = plt.subplots(1,2,figsize=(20, 12))

axs2[0].plot(data_cv['Temp'].values,data_cv['E/N'].values,label='avgEnergy/spin')
axs2[0].set_xlabel('Temperature')
axs2[0].set_ylabel('averageEnergy')
axs2[0].legend()
axs2[0].grid()

axs2[1].plot(data_cv['Temp'].values,data_cv['Cv'].values,label='Cv')
axs2[1].set_xlabel('Temperature')
axs2[1].set_ylabel('SpecificHeat, Cv')
axs2[1].legend()
axs2[1].grid()

#set title as required
plt.title('ising3D, L=7')
plt.show()
