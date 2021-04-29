#!/usr/bin/python3
"""
3D-ising model, a computational physics exercise.
"""

import numpy as np
import random
import time
import csv

start = time.time()

#constants
J = 1 #J, exchange energy
KT = 4.0 #energy in units of J
MU = 0.33 #bohr magneton*gyromagnetic ratio
B = 0 #T, external magnetic field
   
#lattice
L = int(input('\nenter LATTICE side size: '))
N = L*L*L #number of spins in lattice

#initialization
lattice = np.ones((L,L,L)) #all spins UP

#initial energy
E=0.0 ; M=0
for i in range(L):
    for j in range(L):
        for k in range(L):
            x1=i-1; x2=i+1
            y1=j-1; y2=j+1
            z1=k-1; z2=k+1
            #periodic boundary conditions
            if i==0: x1=L-1
            if i==L-1: x2=0
            if j==0: y1=L-1
            if j==L-1: y2=0
            if k==0: z1=L-1
            if k==L-1: z2=0
            
            E+=-J*lattice[i,j,k]*(lattice[x1,j,k]+lattice[x2,j,k]+lattice[i,y1,k]\
               +lattice[i,y2,k]+lattice[i,j,z1]+lattice[i,j,z2])-MU*B*lattice[i,j,k]
            M+=lattice[i,j,k]
            
E *= 0.5
print('\nInitial Lattice Energy: ',E)
print('\nInitial Lattice Magnetisation: ',M)
#Metropolis---------------------------------------------------
iter = int(input('\nenter number of MCS iterations: '))

#write to file
filename = 'ising3d_data.csv'

with open(filename,'w') as file: 
    out = csv.writer(file)
    labels = ['Time','E','M','E/N','M/N']
    out.writerow(labels) #labels
    out.writerow([0,E,M,E/N,M/N])

    for t in range(1,iter+1):
        for n in range(L):
            for m in range(L):
                #random lattice nodes
                i = int(random.uniform(0,L))
                j = int(random.uniform(0,L))
                k = int(random.uniform(0,L))
            
                x1=i-1; x2=i+1
                y1=j-1; y2=j+1
                z1=k-1; z2=k+1
                #periodic boundary conditions
                if i==0: x1=L-1
                if i==L-1: x2=0
                if j==0: y1=L-1
                if j==L-1: y2=0
                if k==0: z1=L-1
                if k==L-1: z2=0
            
                #before spin flip
                Ei =-J*lattice[i,j,k]*(lattice[x1,j,k]+lattice[x2,j,k]+lattice[i,y1,k]\
                    +lattice[i,y2,k]+lattice[i,j,z1]+lattice[i,j,z2])-MU*B*lattice[i,j,k]
                #spin flip
                lattice[i,j,k] *= -1
                #after spin flip
                Ef =-J*lattice[i,j,k]*(lattice[x1,j,k]+lattice[x2,j,k]+lattice[i,y1,k]\
                    +lattice[i,y2,k]+lattice[i,j,z1]+lattice[i,j,z2])-MU*B*lattice[i,j,k]
                #energy difference
                dE = Ef-Ei
            
                if dE<=0:
                    E+=dE
                    M+=2*lattice[i,j,k]
                elif np.exp(-dE/KT)>=random.uniform(0,1):
                    E+=dE
                    M+=2*lattice[i,j,k]
                else:
                    lattice[i,j,k] *= -1
        #write data to file    
        out.writerow([t,E,M,E/N,M/N])

print('time taken: ',time.time()-start) 
             
