# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 15:02:47 2020

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
os.chdir(r'C:\karagoz\01-RESEARCH\01-Projects\01-In_silico_modeling_of_Integrin_function\003-Ligand_competition_model\Unique_clusters\figures_FN_vWF_equal_IC\parameterScan_k_values')
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

  species i, I, $F, IF, $W, IW, C1, C2, C3; 
  #inactive integrin, active integrin, fibronectin, integrin+fibronectin, vonWillebrand Factor A, integrin+vonWillebrand Factor A, clustered integrins respectively.
  
  #set initial values:
  i = 0.05; # integrin avB3 
  I = 0; 
  F = 0.33   ; #fibronectin (set the same as vWA)
  IF = 0;
  W = 0.33   ; #von Willebrand factor A 
  IW = 0;  
  C1 = 0;     # IF+IF cluster
  C2 = 0;     # IW + IW cluster
  C3 = 0;     # IF+IW cluster

  J1: i -> I; k1*i - k2*I; # reaction; reaction rate law;   # activation step, k1 rate of activation, k2 rate of inactivation
  J2: I + $F -> IF; k3*I*F - k4*IF;                          # ligand binding step, k3 rate of fibronectin binding, k4 rate of dissociation
  J3: I + $W -> IW; k5*I*W - k6*IW;                          # alternative ligand binding step, k5 rate of vWFA binding, k6 rate of dissociation
  J4: IF + IF -> C1; k7*IF^2 - k8*C1;
  J5: IW + IW -> C2; k7*IW^2 - k8*C2;
  J6: IF + IW -> C3; k7*IF*IW - k8*C3;                         # clustering step, k7 rate of clustering, k8 rate of dissociation

  k1 = 5*10^6; k2 = 10^8; k3 = 1.6*10^8 ; k4 = 3.5*10^-1; k5 = 1.6*10^4 ; k6 = 2.3*10^-2; k7 = 1.6*10^8; k8 = 0.5*10^7; # assign constant values to global parameters
  end
  """
#%% 
r2 = te.loada(Ant_str)

from tellurium import ParameterScan 

p_k5= ParameterScan(r2) 

#  concentration of W-bound integrins with changing the initial amount of W 

p_k5.endTime = 0.00001 
p_k5.numberOfPoints = 100
p_k5.startValue = 6400
p_k5.endValue = 25600
p_k5.value = 'k5'
p_k5.selection = ['IW']
p_k5.ylabel="[L2-bound integrins]"
p_k5.plotGraduatedArray()

#%%
r2 = te.loada(Ant_str)

r2.k5 = 6400

result_1 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_1 = result_1.assign(grad="6400")
#%%
r2.reset()
r2.k5 = 8533.33

result_2 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_2 = result_2.assign(grad="8533.33")
#%%
r2.reset()
r2.k5 = 10666.67

result_3 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_3 = result_3.assign(grad="10666.67")

#%%
r2.reset()
r2.k5 = 12800

result_4 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_4 = result_4.assign(grad="12800")
#%%
r2.reset()
r2.k5 = 14933.33

result_5 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_5 = result_5.assign(grad="14933.33")

#%%
r2.reset()
r2.k5 = 17066.67

result_6 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_6 = result_6.assign(grad="17066.67")
#%%
r2.reset()
r2.k5 = 19200

result_7 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_7 = result_7.assign(grad="19200")

#%%
r2.reset()
r2.k5 = 21333.33

result_8 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_8 = result_8.assign(grad="21333.33")
#%%
r2.reset()
r2.k5 = 23466.67

result_9 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_9 = result_9.assign(grad="23466.67")
#%%
r2.reset()
r2.k5 = 25600

result_10 = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'IF', 'IW']), columns=['time', 'F_bound', 'vWA_bound'])
result_10 = result_10.assign(grad="25600")



#%%

#make one dataframe out of day 25 results:
frames=[result_1, result_2,result_3,result_4,result_5,result_6,result_7,result_8,result_9,result_10 ]
df2 = pd.concat(frames)

df2.to_csv('k5_paramScan_IF_IW_results.csv', sep='\t', index=False)



