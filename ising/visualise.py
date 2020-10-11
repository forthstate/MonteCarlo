#!/usr/bin/python3
"""
visualisation code file ISING3D.PY
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#---Load csv Data file

data = pd.read_csv('ising3d_data.csv')

#---Figure initialisation

fig, axs = plt.subplots(1,2,figsize=(20, 6))

axs[0].plot(data['Time'].values,data['E'].values,label='latticeEnergy')
axs[0].plot(data['Time'].values,data['M'].values,label='latticeMagnetisation')
axs[0].set_xlabel('MCS-step')
axs[0].set_ylabel('E/M')
axs[0].legend()

axs[1].plot(data['Time'].values,data['E/N'].values,label='Energy/spin')
axs[1].plot(data['Time'].values,data['M/N'].values,label='Magnetisation/spin')
axs[1].set_xlabel('MCS-step')
axs[1].set_ylabel('E/M-per-Spin')
axs[1].legend()


plt.show()
