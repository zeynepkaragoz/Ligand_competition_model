# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 14:00:08 2020

@author: Zeynep Karagoz
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
import math
os.chdir(r'C:\karagoz\01-RESEARCH\01-Projects\01-In_silico_modeling_of_Integrin_function\003-Ligand_competition_model\ligands_free\08_parameter_sensitivity\01_integrins_in_excess')

#%% 
# 2 ligand 1 integrin  

#TEST1: what happens when integrin concentration is equal to total ligand amount
#TEST2: what happens when integrin concentration > total ligand amount


# differentiating between 3 types of clusters IF+IF, IW+IW and IF+IW 
# but the rates of clustering are the same between 3 clusters. only the identities are different


# ligands: fibronectin and von Willebrand Factor A (initial values adapted from kidney orgaoid iBAQ values)
# integrin avB3 (initial value from Hudson et al)
# here the rate for ligand binding-unbinding  is from Hudson et al (both ligands) 
# activation/inactivation is from  Yu et al. 2017
# clustering forward is numerically from Hudson, then reverse is estimated using the Kd from Yu et al .

# PARAMETER SENSITIVITY #

# for each k I run the model with k, k+0.2, k-0.2
# for each Ligand I run the model with L, L+0.2, L-0.2
# for i I run the model with i, i+0.2, i-0.2

# begin TEST1:
Ant_str = """
  model test # ligand competititon model 
  species i, a, F, IF, W, IW, C1, C2, C3; 
  #inactive integrin, active integrin, fibronectin, integrin+fibronectin, vonWillebrand Factor A, integrin+vonWillebrand Factor A, clustered integrins respectively.
  
  #set initial values:
  i = 0.51; # integrin avB3 
  a = 0; 
  F = 0.18   ; #fibronectin
  IF = 0;
  W = 0.33   ; #von Willebrand factor A 
  IW = 0;  
  C1 = 0;     # IF+IF cluster
  C2 = 0;     # IW + IW cluster
  C3 = 0;     # IF+IW cluster
  J1: i -> a; k1*i - k2*a; # reaction; reaction rate law;   # activation step, k1 rate of activation, k2 rate of inactivation
  J2: a + F -> IF; k3*a*F - k4*IF;                          # ligand binding step, k3 rate of fibronectin binding, k4 rate of dissociation
  J3: a + W -> IW; k5*a*W - k6*IW;                          # alternative ligand binding step, k5 rate of vWFA binding, k6 rate of dissociation
  J4: IF + IF -> C1; k7*IF^2 - k8*C1;
  J5: IW + IW -> C2; k7*IW^2 - k8*C2;
  J6: IF + IW -> C3; k7*IF*IW - k8*C3;                         # clustering step, k7 rate of clustering, k8 rate of dissociation
  k1 = 5*10^6; k2 = 10^8; k3 = 1.6*10^8 ; k4 = 3.5*10^-1; k5 = 1.6*10^4 ; k6 = 2.3*10^-2; k7 = 1.6*10^8; k8 = 0.5*10^7; # assign constant values to global parameters
  end
  """
r2 = te.loada(Ant_str)
#%%
r2.conservedMoietyAnalysis = True
r2.steadyState()
#  1.3972645263080818e-10 --> this is a good number!
 
print(r2.getSteadyStateValuesNamedArray())

#      [IF],         [IW],          [a],         [C3],          [i],     [W],     [C1],        [C2],      [F]
# [[ 0.0818642, -0.000124576, -5.41922e-10, -0.000326346, -1.08384e-08, 0.33045, 0.214456, 4.96614e-07, -0.33045]]

#these are the steady state values.
print(r2.getRatesOfChange())
# [-2.27366737e-13 -1.38054901e-10 -6.93804687e-18  2.27373675e-13 0.00000000e+00  8.47032947e-22] these are the rates of change of the 6 reactions, all approaching 0 --> confirm steady state. 


#%%
# simulate for day18 initial conditions,
# store results in pandas DataFrame "result" : 
r2 = te.loada(Ant_str)
result = pd.DataFrame(r2.simulate(0, 1 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered'])
result.to_csv('integrin_concentration_equal_to_ligand_simResults.csv', sep='\t', index=False)
print(result.columns)
for i, col in enumerate(result.columns[1:10]):
    plt.figure(i)
    sns.set_style("whitegrid")
    sns.relplot(x='time', y=col, kind='line' ,data=result)
    plt.ylabel(col+' uM')
    plt.xlabel("time (s)")
    plt.xticks(rotation=15)
# reset the model,
# change the initial ligand concentrations to day25, simulate,
# store results in pandas DataFrame "resul_old": 
r2.reset()
#%%
# begin TEST2:
Ant_str = """
  model test # ligand competititon model 
  species i, a, F, IF, W, IW, C1, C2, C3; 
  #inactive integrin, active integrin, fibronectin, integrin+fibronectin, vonWillebrand Factor A, integrin+vonWillebrand Factor A, clustered integrins respectively.
  
  #set initial values:
  i = 1; # integrin avB3 
  a = 0; 
  F = 0.18   ; #fibronectin
  IF = 0;
  W = 0.33   ; #von Willebrand factor A 
  IW = 0;  
  C1 = 0;     # IF+IF cluster
  C2 = 0;     # IW + IW cluster
  C3 = 0;     # IF+IW cluster
  J1: i -> a; k1*i - k2*a; # reaction; reaction rate law;   # activation step, k1 rate of activation, k2 rate of inactivation
  J2: a + F -> IF; k3*a*F - k4*IF;                          # ligand binding step, k3 rate of fibronectin binding, k4 rate of dissociation
  J3: a + W -> IW; k5*a*W - k6*IW;                          # alternative ligand binding step, k5 rate of vWFA binding, k6 rate of dissociation
  J4: IF + IF -> C1; k7*IF^2 - k8*C1;
  J5: IW + IW -> C2; k7*IW^2 - k8*C2;
  J6: IF + IW -> C3; k7*IF*IW - k8*C3;                         # clustering step, k7 rate of clustering, k8 rate of dissociation
  k1 = 5*10^6; k2 = 10^8; k3 = 1.6*10^8 ; k4 = 3.5*10^-1; k5 = 1.6*10^4 ; k6 = 2.3*10^-2; k7 = 1.6*10^8; k8 = 0.5*10^7; # assign constant values to global parameters
  end
  """
r2 = te.loada(Ant_str)
#%%
r2.conservedMoietyAnalysis = True
r2.steadyState()
#  5.112411681409751e-11 --> this is a good number!
 
print(r2.getSteadyStateValuesNamedArray())

#      [IF],         [IW],          [a],         [C3],          [i],      [W],     [C1],        [C2],       [F]
# [[ 0.11747, -7.19838e-05, -3.13241e-10, -0.000270589, -6.26482e-09, 0.330342, 0.441572, 1.65813e-07, -0.820342]]
#these are the steady state values.
print(r2.getRatesOfChange())
# [-4.54740412e-13  5.04372100e-11 -6.93889390e-18  4.54747351e-13 0.00000000e+00  0.00000000e+00] these are the rates of change of the 6 reactions, all approaching 0 --> confirm steady state. 


#%%
# simulate for day18 initial conditions,
# store results in pandas DataFrame "result" : 
r2 = te.loada(Ant_str)
result = pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered'])
result.to_csv('integrin_concentration_excess_to_ligand_simResults.csv', sep='\t', index=False)
print(result.columns)
for i, col in enumerate(result.columns[1:10]):
    plt.figure(i)
    sns.set_style("whitegrid")
    sns.relplot(x='time', y=col, kind='line' ,data=result)
    plt.ylabel(col+' uM')
    plt.xlabel("time (s)")
    plt.xticks(rotation=15)
# reset the model,
# change the initial ligand concentrations to day25, simulate,
# store results in pandas DataFrame "resul_old": 
r2.reset()
#%%
steady_states_all = pd.DataFrame(result.iloc[[99,]])
steady_states_all = steady_states_all.reset_index(drop=True)
## Index = 0 - SS with normal parameters, equal IC for ligands and 0.05 integrin IC
#%%
# reset the model,
# change the k1 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k1 = 5*10**6+(0.2*5*10**6)
## Index = 1 - SS with k1 increased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k1 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k1 = 5*10**6-(0.2*5*10**6)
## Index = 2 - SS with k1 decreased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k2 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k2 = 10**8+(0.2*10**8)
## Index = 3 - SS with k2 increased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k2 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k2 = 10**8-(0.2*10**8)
## Index = 4 - SS with k2 decreased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k3 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k3 = 1.6*10**8+(0.2*1.6*10**8)
## Index = 5 - SS with k3 increased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k3 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k3 = 1.6*10**8-(0.2*1.6*10**8)
## Index = 6 - SS with k3 decreased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k4 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k4 = 3.5*10**-1+(0.2*3.5*10**-1)
## Index = 7 - SS with k4 increased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k4 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k4 = 3.5*10**-1-(0.2*3.5*10**-1)
## Index = 8 - SS with k4 decreased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k5 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k5 = 1.6*10**4+(0.2*1.6*10**4)
## Index = 9 - SS with k5 increased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k5 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k5 = 1.6*10**4-(0.2*1.6*10**4)
## Index = 10 - SS with k5 decreased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k6 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k6 = 2.3*10**-2+(0.2*2.3*10**-2)
## Index = 11 - SS with k6 increased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k6 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k6 = 2.3*10**-2-(0.2*2.3*10**-2)
## Index = 12 - SS with k6 decreased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k7 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k7 = 1.6*10**8+(0.2*1.6*10**8)
## Index = 13 - SS with k7 increased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k7 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k7 = 1.6*10**8-(0.2*1.6*10**8)
## Index = 14 - SS with k7 decreased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k8 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k8 = 0.5*10**7+(0.2*0.5*10**7)
## Index = 15 - SS with k8 increased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the k8 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k8 = 0.5*10**7-(0.2*0.5*10**7)
## Index = 16 - SS with k8 decreased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the i to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.i = 1+(0.2*1)
## Index = 17 - SS with i increased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the i to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.i = 1-(0.2*1)
## Index = 18 - SS with i decreased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the F to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.F = 0.18+(0.2*0.18)
## Index = 19 - SS with F increased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the F to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.F = 0.18-(0.2*0.18)
## Index = 20 - SS with F decreased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the W to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.W = 0.33+(0.2*0.33)
## Index = 21 - SS with W increased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# reset the model,
# change the W to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.W = 0.33-(0.2*0.33)
## Index = 22 - SS with W decreased 20%
steady_states_all = steady_states_all.append(pd.DataFrame(r2.simulate(0, 0.03 , 100 , ['time', 'i', 'a','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered']).iloc[[99,]], ignore_index = True)
#%%
# add test condition explanation column
steady_states_all = steady_states_all.assign(test_condition = ['normal', 
                                                             'k1_up', 'k1_down', 
                                                             'k2_up', 'k2_down', 
                                                             'k3_up', 'k3_down', 
                                                             'k4_up', 'k4_down', 
                                                             'k5_up', 'k5_down', 
                                                             'k6_up', 'k6_down', 
                                                             'k7_up', 'k7_down', 
                                                             'k8_up', 'k8_down',
                                                             'i_up', 'i_down',
                                                             'F_up', 'F_down',
                                                             'W_up', 'W_down'])
# save dataframe to csv
steady_states_all.to_csv('parameterSensitivity_steadyStates_all_integrin_excess.csv', sep='\t', index=False)