# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 16:19:51 2020
@author: zeynep karagoz
# -*- coding: utf-8 -*-
"""

#Importing relevant packages
import tellurium as te # Python-based modeling environment for kinetic models
import roadrunner as rr # High-performance simulation and analysis library
import numpy as np # Scientific computing package
import matplotlib.pylab as plt # Additional Python plotting utilities
import pandas as pd
import os
import seaborn as sns
os.chdir(r'C:\karagoz\01-RESEARCH\01-Projects\01-In_silico_modeling_of_Integrin_function\003-Ligand_competition_model\ligands_free\07_parameter_scan_L2')
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 13}

matplotlib.rc('font', **font)
#%% 
# 2 ligand 1 integrin  

# setting the ligands to EQUAL initial conditions
# differentiating between 3 types of clusters IF+IF, IW+IW and IF+IW 
# but the rates of clustering are the same between 3 clusters. only the identities are different


# ligands: fibronectin and von Willebrand Factor A (initial values adapted from kidney orgaoid iBAQ values)
# integrin avB3 (initial value from Hudson et al)
# here the rate for ligand binding-unbinding  is from Hudson et al (both ligands) 
# activation/inactivation is from  Yu et al. 2017
# clustering forward is numerically from Hudson, then reverse is estimated using the Kd from Yu et al .
Ant_str = """
  model test # activation model 

  species i, I, F, IF, W, IW, C1, C2, C3; 
  #inactive integrin, active integrin, fibronectin, integrin+fibronectin, vonWillebrand Factor A, integrin+vonWillebrand Factor A, clustered integrins respectively.
  
  #set initial values:
  i = 0.05; # integrin avB3 
  I = 0; 
  F = 0.84   ; #fibronectin (set the same as vWA)
  IF = 0;
  W = 0.50   ; #von Willebrand factor A 
  IW = 0;  
  C1 = 0;     # IF+IF cluster
  C2 = 0;     # IW + IW cluster
  C3 = 0;     # IF+IW cluster

  J1: i -> I; k1*i - k2*I; # reaction; reaction rate law;   # activation step, k1 rate of activation, k2 rate of inactivation
  J2: I + F -> IF; k3*I*F - k4*IF;                          # ligand binding step, k3 rate of fibronectin binding, k4 rate of dissociation
  J3: I + W -> IW; k5*I*W - k6*IW;                          # alternative ligand binding step, k5 rate of vWFA binding, k6 rate of dissociation
  J4: IF + IF -> C1; k7*IF^2 - k8*C1;
  J5: IW + IW -> C2; k7*IW^2 - k8*C2;
  J6: IF + IW -> C3; k7*IF*IW - k8*C3;                         # clustering step, k7 rate of clustering, k8 rate of dissociation

  k1 = 5*10^6; k2 = 10^8; k3 = 1.6*10^8 ; k4 = 3.5*10^-1; k5 = 1.6*10^4 ; k6 = 2.3*10^-2; k7 = 1.6*10^8; k8 = 0.5*10^7; # assign constant values to global parameters
  end
  """
#%% 
r2 = te.loada(Ant_str)

from tellurium import ParameterScan 

p_W= ParameterScan(r2) 

#  concentration of W-bound integrins with changing the initial amount of W 

p_W.endTime = 0.00001 
p_W.numberOfPoints = 100
p_W.startValue = 0.50
p_W.endValue = 2.97
p_W.value = 'W'
p_W.selection = ['IW']
p_W.ylabel="[L2-bound integrins]"
p_W.plotGraduatedArray()
#%%

# simulate for day25 initial conditions, gradually increase W
# store results in pandas DataFrame "result" : 
r2 = te.loada(Ant_str)
result_1 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time','F_bound', 'vWA_bound'])
result_1 = result_1.assign(grad="0.50")

#%%
# reset the model,
# change the initial ligand concentrations, simulate, store in DF

r2.reset()
r2.W = 0.77

result_2 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_2 = result_2.assign(grad="0.77")
#%%
r2.reset()
r2.W = 1.05

result_3 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_3 = result_3.assign(grad="1.05")

#%%
r2.reset()
r2.W = 1.32

result_4 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_4 = result_4.assign(grad="1.32")
#%%
r2.reset()
r2.W = 1.60

result_5 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_5 = result_5.assign(grad="1.60")

#%%
r2.reset()
r2.W = 1.87

result_6 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_6 = result_6.assign(grad="1.87")
#%%
r2.reset()
r2.W = 2.15

result_7 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_7 = result_7.assign(grad="2.15")

#%%
r2.reset()
r2.W = 2.42

result_8 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_8 = result_8.assign(grad="2.42")
#%%
r2.reset()
r2.W = 2.70

result_9 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_9 = result_9.assign(grad="2.70")
#%%
r2.reset()
r2.W = 2.97

result_10 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_10 = result_10.assign(grad="2.97")



#%%

#make one dataframe out of all results:
frames=[result_1, result_2,result_3,result_4,result_5,result_6,result_7,result_8,result_9,result_10 ]
df2 = pd.concat(frames)

df2.to_csv('L2_paramScan.csv', sep='\t', index=False)
