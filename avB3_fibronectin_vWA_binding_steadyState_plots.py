# -*- coding: utf-8 -*-
"""
Created on Thu May  7 14:00:17 2020

@author: zeynep karagoz
"""

#Importing relevant packages
import tellurium as te # Python-based modeling environment for kinetic models
import roadrunner as rr # High-performance simulation and analysis library
import numpy as np # Scientific computing package
import random # Generate random numbers
import matplotlib.pylab as plt # Additional Python plotting utilities
import pandas as pd
import os
import math
import seaborn as sns
os.chdir(r'C:\Users\p70067970\Documents\Hudson_model_work\ligand_competition\figures\aVB3-FN-vWFA')
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 13}

matplotlib.rc('font', **font)
#%% 
# 2 ligand 1 integrin  
# ligands: fibronectin and von Willebrand Factor A (initial values adapted from kidney orgaoid iBAQ values)
# integrin avB3 (initial value from Hudson et al)
# here the rate for ligand binding-unbinding  is from Hudson et al (both ligands) 
# activation/inactivation is from  Cheng et al. 2020
# clustering forward / reverse is from Cheng et al 2020 (there is a problem here, the units in the paper are 1/s for both, in my equations, forward should be 1/Ms)
Ant_str = """
  model test # activation model 

  species i, I, $F, IF, $W, IW, C; 
  #inactive integrin, active integrin, fibronectin, integrin+fibronectin, vonWillebrand Factor A, integrin+vonWillebrand Factor A, clustered integrins respectively.
  
  #set initial values:
  i = 0.05; # integrin avB3
  I = 0; 
  F = 0.18   ; #fibronectin
  IF = 0;
  W = 0.33   ; #von Willebrand factor A
  IW = 0;  
  C = 0;

  J1: i -> I; k1*i - k2*I; # reaction; reaction rate law;   # activation step, k1 rate of activation, k2 rate of inactivation
  J2: I + $F -> IF; k3*I*F - k4*IF;                          # ligand binding step, k3 rate of fibronectin binding, k4 rate of dissociation
  J3: I + $W -> IW; k5*I*W - k6*IW;                          # alternative ligand binding step, k5 rate of vWFA binding, k6 rate of dissociation
  J4: IF + IF -> C; k7*IF^2 - k8*C;
  J5: IW + IW -> C; k7*IW^2 - k8*C;
  J6: IF + IW -> C; k7*IF*IW - k8*C;                         # clustering step, k7 rate of clustering, k8 rate of dissociation

  k1 = 23; k2 = 0.576 ; k3 = 1.6*10^8 ; k4 = 3.5*10^-1; k5 = 1.6*10^4 ; k6 = 2.3*10^-2; k7 = 1; k8 = 0.1; # assign constant values to global parameters
  end
  """
r2 = te.loada(Ant_str)
#%%

r2.conservedMoietyAnalysis = True
r2.steadyState()
#  3.2507367254238236e-15 --> this is a good number!
 
print(r2.getSteadyStateValuesNamedArray())
#      [IF],      [IW],         [I],        [i],      [W],         [F],        [C]
# [[ 0.0207107, 0.0207107, 6.48193e-08, 1.6233e-09, 0.459301, 0.000698937, 0.00428931]]

#these are the steady state values.
print(r2.getRatesOfChange())
# [ 2.11419424e-17  3.79470760e-19 -2.08708918e-17  0.00000000e+00  -5.42101086e-20] these are the rates of change of the 5 , all approaching 0 --> confirm steady state. 

#simulate for 3 mins!!! = 180s

#%%
#plot the integrin clusters in percent of the total 
r2 = te.loada(Ant_str)
result = pd.DataFrame(r2.simulate(0, 180 , 100 , ['time', 'i', 'I','F','W', 'IF', 'IW',  'C']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'vWA_bound', 'cluster'])

r2.reset()
r2.F = 0.46
r2.W = 0.50

result_old = pd.DataFrame(r2.simulate(0, 180 , 100 , ['time', 'i', 'I','F','W', 'IF', 'IW',  'C']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'vWA_bound', 'cluster'])
#create new column with total integrin amount at each time step
result = result.assign(Sum = result.inactive + result.active + result.F_bound + result.vWA_bound +2*result.cluster, Experiment="day18")
result = result.assign(percentCluster = 2*result.cluster / result.Sum,
                             percentActive = result.active / result.Sum,
                             percentInactive = result.inactive / result.Sum,
                             percentF_Bound = result.F_bound / result.Sum,
                             percentvWA_Bound = result.vWA_bound / result.Sum)

result_old = result_old.assign(Sum = result_old.inactive + result_old.active + result_old.F_bound + result_old.vWA_bound +2*result_old.cluster, Experiment="day25")
result_old = result_old.assign(percentCluster = 2*result_old.cluster / result_old.Sum,
                             percentActive = result_old.active / result_old.Sum,
                             percentInactive = result_old.inactive / result_old.Sum,
                             percentF_Bound = result_old.F_bound / result_old.Sum,
                             percentvWA_Bound = result_old.vWA_bound / result_old.Sum)

#make one df 
df2 = result.append(result_old)


#%% 
# inactive integrins absolute numbers
#plot with seaborn
g3=sns.relplot(x="time", y="inactive", kind="line", hue= "Experiment",data=df2)
#plt.ylim(0, 1)
g3.set_axis_labels("time (s)", "inactive integrins (nM)")
g3.fig.autofmt_xdate()
plt.savefig('inactiveIntegrins_day18_25.png')



#%%
# inactive integrins percent

#plot with seaborn
g4=sns.relplot(x="time", y="percentInactive", kind="line", hue= "Experiment",data=df2)
plt.ylim(0, 1)
g4.set_axis_labels("time (s)", "inactive / total integrins")
g4.fig.autofmt_xdate()
plt.savefig('Percent_inactiveIntegrins_day18_25.png')

#%%
# active integrins absolute numbers
#plot with seaborn
g5=sns.relplot(x="time", y="active", kind="line", hue= "Experiment",data=df2)
#plt.ylim(0, 1)
g5.set_axis_labels("time (s)", "active integrins (nM)")
g5.fig.autofmt_xdate()
plt.savefig('activeIntegrins_day18_25.png')



#%%
# active integrins percent

#plot with seaborn
g6=sns.relplot(x="time", y="percentActive", kind="line", hue= "Experiment",data=df2)
plt.ylim(0, 1)
g6.set_axis_labels("time (s)", "active / total integrins")
g6.fig.autofmt_xdate()
plt.savefig('Percent_activeIntegrins_day18_25.png')
#%% 
# F bound integrins absolute numbers
#plot with seaborn
g7=sns.relplot(x="time", y="F_bound", kind="line", hue= "Experiment",data=df2)
#plt.ylim(0, 1)
g7.set_axis_labels("time (s)", "fibronectin bound integrins (nM)")
g7.fig.autofmt_xdate()
plt.savefig('F_bound_Integrins_day18_25.png')



#%%
# F bound integrins percent

#plot with seaborn
g8=sns.relplot(x="time", y="percentF_Bound", kind="line", hue= "Experiment",data=df2)
plt.ylim(0, 1)
g8.set_axis_labels("time (s)", "Fibronectin bound / total integrins")
g8.fig.autofmt_xdate()
plt.savefig('Percent_F_bound_Integrins_day18_25.png')

#%% 
# vWFA bound integrins absolute numbers
#plot with seaborn
g9=sns.relplot(x="time", y="vWA_bound", kind="line", hue= "Experiment",data=df2)
#plt.ylim(0, 1)
g9.set_axis_labels("time (s)", "vWFA bound integrins (nM)")
g9.fig.autofmt_xdate()
plt.savefig('vWFA_bound_Integrins_day18_25.png')



#%%
# vWFA bound integrins percent

#plot with seaborn
g10=sns.relplot(x="time", y="percentvWA_Bound", kind="line", hue= "Experiment",data=df2)
plt.ylim(0, 1)
g10.set_axis_labels("time (s)", "vWFA bound / total integrins")
g10.fig.autofmt_xdate()
plt.savefig('Percent_vWFA_bound_Integrins_day18_25.png')

#%% 
# clustered integrins absolute
g11=sns.relplot(x="time", y="cluster", kind="line", hue= "Experiment",data=df2)

g11.set_axis_labels("time (s)", "clustered integrins")
g11.fig.autofmt_xdate()
plt.savefig('Clustered_Integrins_day18_25.png')

#%% 
# clustered integrins percent
g12=sns.relplot(x="time", y="percentCluster", kind="line", hue= "Experiment",data=df2)
plt.ylim(0, 1)
g12.set_axis_labels("time (s)", "clustered / total integrins")
g12.fig.autofmt_xdate()
plt.savefig('Percent_Clustered_Integrins_day18_25.png')

#%% 
# fibronectin absolute
g11=sns.relplot(x="time", y="F", kind="line", hue= "Experiment",data=df2)

g11.set_axis_labels("time (s)", "Fibronectin")
g11.fig.autofmt_xdate()
plt.savefig('Fibronectin_day18_25.png')
#%% 
# vWFA absolute
g11=sns.relplot(x="time", y="W", kind="line", hue= "Experiment",data=df2)

g11.set_axis_labels("time (s)", "vWFA")
g11.fig.autofmt_xdate()
plt.savefig('vWFA_day18_25.png')



