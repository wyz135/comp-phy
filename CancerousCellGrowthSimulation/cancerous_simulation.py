# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 15:36:24 2018

@author: WYZ
"""

import random as ran

# Rates of evolution
k0 = 0.0; k1 = 0.6; k2 = 0.30; k3 = 0.4; k4 = 0.4;
rates = [k0,k1,k2,k3,k4] # in terms of probabilities after a time period

def count_cancer(tissue):
  N = 0
  for i in range(0,len(tissue)):
    for j in range(0,len(tissue)):
      if tissue[i][j] == 'C':
        N += 1
  return N


def update_tissue(tissue,rates):
  # Update the prop values to determine which states will change
  profilerate = []
  for i in range(1,N-1):
    for j in range(1,N-1):
      if tissue[i][j] == 'N':
        if ran.random() < rates[0]:
          tissue[i][j] = 'C'
          
      elif tissue[i][j] == 'C':
        # Counting the number of cancerous neighbors
        q = 0
        for l in [-1,1]:
          if tissue[i+l][j] == 'C':
            q += 1
          if tissue[i][j+l] == 'C':
            q += 1
        # decide whether cancer cell will profilerate or be destroyed
        if ran.random() < rates[2]:
          tissue[i][j] = 'E'
        if ran.random() < rates[1]*(1-q/4):
          profilerate.append([i,j])
          
      elif tissue[i][j] =='E':
        if ran.random() < rates[3]:
          tissue[i][j] = 'D'
      elif tissue[i][j] == 'D':
        if ran.random() < rates[4]:
          tissue[i][j] = 'N'
  
  # Profileration of neighboring sites
  if len(profilerate) > 0:
    for site in profilerate:
      potential_sites = []
      for l in [[-1,0],[0,-1],[1,0],[0,1]]:
        if tissue[site[0]+l[0]][site[1]+l[1]] != 'C':
          potential_sites.append(l)
      
      if len(potential_sites) > 0:
        infected_site = ran.choice(potential_sites)
        tissue[site[0]+infected_site[0]][site[1]+infected_site[1]] = 'C'
    
  return tissue

xc= []
for ite in range(0,100):
  # Initialize tissue array of N by N with all cells normal
  N = 100
  tissue = []
  row = ['N']*N
  for n in range(0,N):
    tissue.append(list(row))

  # Randomly place N0 cancer cells
  N0 = 100
  coordinates = []
  for i in range(0,N):
    for j in range(0,N):
      coordinates.append([i,j])

  for i in range(0,N0):
    cancer_cell_replace = coordinates.pop(ran.randint(0,len(coordinates)-1))
    tissue[cancer_cell_replace[0]][cancer_cell_replace[1]] = 'C'

  # Running the simulation
  NC = [count_cancer(tissue)/(N**2)]
  days = [0]
  for i in range(1,501):
    tissue = update_tissue(tissue,rates)
    NC.append(count_cancer(tissue)/(N**2))
    days.append(i)

  # Obtaining values of xc
  for i in [200,250,300,350,400,450]:
    xc.append(NC[i])

# Plotting days vs number of cancer cells
import matplotlib.pyplot as plt
plt.plot(days,NC)
plt.xlabel('Days'); plt.ylabel('Number of cancer cells')
plt.show()

plt.xlabel('xc'); plt.ylabel('N')
plt.hist(xc)

import numpy as np
mean = np.mean(xc)
var = np.var(xc)
