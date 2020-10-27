# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 10:52:35 2020

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
os.chdir(r'C:\karagoz\01-RESEARCH\01-Projects\01-In_silico_modeling_of_Integrin_function\003-Ligand_competition_model\ligands_free\08_parameter_sensitivity')

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

# PARAMETER SENSITIVITY #

# for each k I run the model with k, k+0.2, k-0.2
# for each Ligand I run the model with L, L+0.2, L-0.2
# for i I run the model with i, i+0.2, i-0.2


Ant_str = """
  model test # ligand competititon model 
  species i, a, F, IF, W, IW, C1, C2, C3; 
  #inactive integrin, active integrin, fibronectin, integrin+fibronectin, vonWillebrand Factor A, integrin+vonWillebrand Factor A, clustered integrins respectively.
  
  #set initial values:
  i = 0.05; # integrin avB3 
  a = 0; 
  F = 0.18   ; #fibronectin (set the same as vWA)
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


## Index = 0 - SS with normal parameters, equal IC for ligands and 0.05 integrin IC
steadyStates_all = pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]'])

#%%
# reset the model,
# change the k1 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k1 = 5*10**6+(0.2*5*10**6)
r2.conservedMoietyAnalysis = True
## Index = 1 - SS with k1 increased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)

#%%
# reset the model,
# change the k1 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k1 = 5*10**6-(0.2*5*10**6)
r2.conservedMoietyAnalysis = True
## Index = 2 - SS with k1 decreased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)

#%%
# reset the model,
# change the k2 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k2 = 10**8+(0.2*10**8)
r2.conservedMoietyAnalysis = True
## Index = 3 - SS with k2 increased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the k2 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k2 = 10**8-(0.2*10**8)
r2.conservedMoietyAnalysis = True
## Index = 4 - SS with k2 decreased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the k3 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k3 = 1.6*10**8+(0.2*1.6*10**8)
r2.conservedMoietyAnalysis = True
## Index = 5 - SS with k3 increased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the k3 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k3 = 1.6*10**8-(0.2*1.6*10**8)
r2.conservedMoietyAnalysis = True
## Index = 6 - SS with k3 decreased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the k4 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k4 = 3.5*10**-1+(0.2*3.5*10**-1)
r2.conservedMoietyAnalysis = True
## Index = 7 - SS with k4 increased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the k4 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k4 = 3.5*10**-1-(0.2*3.5*10**-1)
r2.conservedMoietyAnalysis = True
## Index = 8 - SS with k4 decreased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the k5 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k5 = 1.6*10**4+(0.2*1.6*10**4)
r2.conservedMoietyAnalysis = True
## Index = 9 - SS with k5 increased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the k5 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k5 = 1.6*10**4-(0.2*1.6*10**4)
r2.conservedMoietyAnalysis = True
## Index = 10 - SS with k5 decreased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the k6 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k6 = 2.3*10**-2+(0.2*2.3*10**-2)
r2.conservedMoietyAnalysis = True
## Index = 11 - SS with k6 increased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the k6 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k6 = 2.3*10**-2-(0.2*2.3*10**-2)
r2.conservedMoietyAnalysis = True
## Index = 12 - SS with k6 decreased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the k7 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k7 = 1.6*10**8+(0.2*1.6*10**8)
r2.conservedMoietyAnalysis = True
## Index = 13 - SS with k7 increased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the k7 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k7 = 1.6*10**8-(0.2*1.6*10**8)
r2.conservedMoietyAnalysis = True
## Index = 14 - SS with k7 decreased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the k8 to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k8 = 0.5*10**7+(0.2*0.5*10**7)
r2.conservedMoietyAnalysis = True
## Index = 15 - SS with k8 increased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the k8 to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.k8 = 0.5*10**7-(0.2*0.5*10**7)
r2.conservedMoietyAnalysis = True
## Index = 16 - SS with k8 decreased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the i to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.i = 0.05+(0.2*0.05)
r2.conservedMoietyAnalysis = True
## Index = 17 - SS with i increased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the i to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.i = 0.05-(0.2*0.05)
r2.conservedMoietyAnalysis = True
## Index = 18 - SS with i decreased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the F to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.F = 0.18+(0.2*0.18)
r2.conservedMoietyAnalysis = True
## Index = 19 - SS with F increased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the F to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.F = 0.18-(0.2*0.18)
r2.conservedMoietyAnalysis = True
## Index = 20 - SS with F decreased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the W to be 0.20 more than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.W = 0.33+(0.2*0.33)
r2.conservedMoietyAnalysis = True
## Index = 21 - SS with W increased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# reset the model,
# change the W to be 0.20 less than normal value,
# append the DF with new steady state values: 
r2.reset()
r2.W = 0.33-(0.2*0.33)
r2.conservedMoietyAnalysis = True
## Index = 22 - SS with W decreased 20%
steadyStates_all = steadyStates_all.append(pd.DataFrame(r2.getSteadyStateValuesNamedArray(), columns=['[IF]','[IW]','[a]','[C3]','[i]', '[W]','[C1]','[C2]', '[F]']), ignore_index=True)
#%%
# add test condition explanation column
steadyStates_all = steadyStates_all.assign(test_condition = ['normal', 
                                                             'k1_up', 'k1_down', 'k2_up', 'k2_down', 'k3_up', 'k3_down', 'k4_up', 'k4_down', 'k5_up', 'k5_down', 'k6_up', 'k6_down', 'k7_up', 'k7_down', 'k8_up', 'k8_down',
                                                             'i_up', 'i_down',
                                                             'F_up', 'F_down',
                                                             'W_up', 'W_down'])
# save dataframe to csv
steadyStates_all.to_csv('parameterSensitivity_steadyStates_all.csv', sep='\t', index=False)
