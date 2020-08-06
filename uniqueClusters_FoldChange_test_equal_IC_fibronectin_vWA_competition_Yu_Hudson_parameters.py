# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 11:17:39 2020

@author: Zeynep Karagoz
"""
#Importing relevant packages
import tellurium as te # Python-based modeling environment for kinetic models
import roadrunner as rr # High-performance simulation and analysis library
import numpy as np # Scientific computing package
import matplotlib.pylab as plt # Additional Python plotting utilities
import pandas as pd
import os
import seaborn as sns
os.chdir(r'C:\karagoz\01-RESEARCH\01-Projects\01-In_silico_modeling_of_Integrin_function\003-Ligand_competition_model\Unique_clusters\figures_FN_vWF_foldChange_testing\FN_vWF_equal_foldChange')
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 13}

matplotlib.rc('font', **font)
#%% 
# 2 ligand 1 integrin  

# setting the ligands to EQUAL FOLD CHANGES between days 25 and 18
# differentiating between 3 types of clusters IF+IF, IT+IT and IF+IT 
# but the rates of clustering are the same between 3 clusters. only the identities are different


# ligands: fibronectin and von Willebrand Factor A (vWA) (initial values adapted from kidney orgaoid iBAQ values)
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
  W = 0.33   ; #vWA 
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
r2 = te.loada(Ant_str)

#k1 = 5*10^6; k2 = 10^8; k3 = 1.6*10^8; k4 = 3.5*10^-1; k5 = 1.6*10^8; k6 = 0.5*10^7; #rates of the ABC model adjusted for 60%clustering

#te.getODEsFromModel(r2)
#vJ1 = k1*i-k2*I
#vJ2 = k3*I*F-k4*IF
#vJ3 = k5*I*W-k6*IW
#vJ4 = k7*pow(IF,2)-k8*C1
#vJ5 = k7*pow(IW,2)-k8*C2
#vJ6 = k7*IF*IW-k8*C3

#di/dt = -vJ1
#dI/dt = vJ1 - vJ2 - vJ3
#dIF/dt = vJ2 - 2.0*vJ4 - vJ6
#dIW/dt = vJ3 - 2.0*vJ5 - vJ6
#dC1/dt = vJ4
#dC2/dt = vJ5
#dC3/dt = vJ6


#%%

r2.conservedMoietyAnalysis = True
r2.steadyState()
#  6.8141049931359555e-19 --> this is a good number!
 
print(r2.getSteadyStateValuesNamedArray())
#     [IF],        [IT],         [I],        [C3],         [i],      [C1],        [C2]
# [[ 0.0211892, 3.22444e-05, 1.40459e-10, 2.18634e-05, 2.80917e-09, 0.0143674, 3.32704e-08]]


#these are the steady state values.
print(r2.getRatesOfChange())
# [-2.84199747e-14 -7.53908050e-12 -1.70200683e-18  2.84217094e-14 0.00000000e+00  0.00000000e+00] these are the rates of change of the 6 reactions, all approaching 0 --> confirm steady state. 

#

#%%
# simulate for day18 initial conditions,
# store results in pandas DataFrame "result" : 
r2 = te.loada(Ant_str)
result = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'i', 'I','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered'])


# reset the model,
# change the initial ligand concentrations to day25, simulate,
# store results in pandas DataFrame "resul_old": 
r2.reset()
r2.F = 0.84
r2.W = 0.84

result_old = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'i', 'I','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'W_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered'])
#create new column with total integrin amount at each time step
result = result.assign(Sum = result.inactive + result.active + result.F_bound + result.W_bound +2*result.IF_IFclustered+2*result.IW_IWclustered+2*result.IF_IWclustered, Experiment="day18")
#result = result.assign(percentClustered = 2*result.clustered / result.Sum,
                             #percentActive = result.active / result.Sum,
                             #percentInactive = result.inactive / result.Sum,
                             #percentF_Bound = result.F_bound / result.Sum,
                             #percentvWA_Bound = result.vWA_bound / result.Sum)

result_old = result_old.assign(Sum = result_old.inactive + result_old.active + result_old.F_bound + result_old.W_bound +2*result_old.IF_IFclustered+2*result_old.IW_IWclustered+2*result_old.IF_IWclustered, Experiment="day25")
#result_old = result_old.assign(percentClustered = 2*result_old.clustered / result_old.Sum,
           
                  #percentActive = result_old.active / result_old.Sum,
                             #percentInactive = result_old.inactive / result_old.Sum,
                             #percentF_Bound = result_old.F_bound / result_old.Sum,
                             #percentvWA_Bound = result_old.vWA_bound / result_old.Sum)

#make one dataframe out of day18 and 25 results:
df2 = result.append(result_old)

df2.to_csv('F_vWA_equalFC_simResults.csv', sep='\t', index=False)