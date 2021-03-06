# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 11:00:59 2021

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
os.chdir(r'C:\karagoz\01-RESEARCH\01-Projects\01-In_silico_modeling_of_Integrin_function\003-Ligand_competition_model\ligands_free\09_no_competition')
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 13}

matplotlib.rc('font', **font)

#%% 
# 2 ligand 1 integrin  

# setting the ligands to different initial conditions
# differentiating between 3 types of clusters IF+IF, IW+IW and IF+IW 
# but the rates of clustering are the same between 3 clusters. only the identities are different

# set F IC to 0 to see the ligand bound integrins in a non-competitive case

Ant_str = """
  model test # activation model 

  species i, I, F, IF, W, IW, C1, C2, C3; 
  #inactive integrin, active integrin, fibronectin, integrin+fibronectin, vonWillebrand Factor A, integrin+vonWillebrand Factor A, clustered integrins respectively.
  
  #set initial values:
  i = 0.05; # integrin avB3 
  I = 0; 
  F = 0   ; #fibronectin 
  IF = 0;
  W = 0.33   ; #von Willebrand factor A 
  IW = 0;  
  C1 = 0;     # IF+IF cluster
  C2 = 0;     # IW + IW cluster
  C3 = 0;     # IF+IW cluster

  J1: i -> I; k1*i - k2*I; # reaction; reaction rate law;   # activation step, k1 rate of activation, k2 rate of inactivation
  J2: I + F -> IF; k3*I*F - k4*IF;                          # ligand binding step, k3 rate of fibronectin binding, k4 rate of dissociation
  J3: I + W -> IW; k5*I*W - k6*IW;                          # alternative ligand binding step, k5 rate of vWFA binding, k6 rate of dissociation
  J4: IF + IF -> C1; k7*IF^2 - k8*C1;
  J5: IW + IW -> C2; k7*IW^2 - k8*C2;
  #J6: IF + IW -> C3; k7*IF*IW - k8*C3;                         # clustering step, k7 rate of clustering, k8 rate of dissociation

  k1 = 5*10^6; k2 = 10^8; ; k3 = 1.6*10^8 ; k4 = 3.5*10^-1; k5 = 1.6*10^4 ; k6 = 2.3*10^-2; k7 = 1.6*10^8; k8 = 0.5*10^7; # assign constant values to global parameters
  end
  """
r2 = te.loada(Ant_str)


#%%

r2.conservedMoietyAnalysis = True
r2.steadyState()
#  2.7830757166182516e-07 --> this is a good number!
 
print(r2.getSteadyStateValuesNamedArray())

#      [IF],        [IW],         [I],         [C3],        [i],      [W],     [C1],      [C2],    [F]
# [[ 9.85657e-17, 0.021209, 1.08885e-07,  2.1777e-06,     0.280002,   0,         0,    0.0143943,    0]]

#these are the steady state values.
print(r2.getRatesOfChange())
# [-3.44980136e-17 -8.73114914e-11  3.44980105e-17  0.00000000e+00 0.00000000e+00] these are the rates of change of the 6 reactions, all approaching 0 --> confirm steady state. 

#%%
# simulate for day18 initial conditions,
# store results in pandas DataFrame "result" : 
r2 = te.loada(Ant_str)
result = pd.DataFrame(r2.simulate(0, 0.05 , 100 , ['time', 'i', 'I','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'vWA_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered'])
#%%

# reset the model,
# change the initial ligand concentrations to day25, simulate,
# store results in pandas DataFrame "resul_old": 
r2.reset()
r2.F = 0
r2.W = 0.50

result_old = pd.DataFrame(r2.simulate(0, 0.05 , 100 , ['time', 'i', 'I','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'vWA_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered'])
#create new column with total integrin amount at each time step
result = result.assign(Sum = result.inactive + result.active + result.F_bound + result.vWA_bound +2*result.IF_IFclustered+2*result.IW_IWclustered+2*result.IF_IWclustered, Experiment="day18")
#result = result.assign(percentClustered = 2*result.clustered / result.Sum,
                             #percentActive = result.active / result.Sum,
                             #percentInactive = result.inactive / result.Sum,
                             #percentF_Bound = result.F_bound / result.Sum,
                             #percentvWA_Bound = result.vWA_bound / result.Sum)

result_old = result_old.assign(Sum = result_old.inactive + result_old.active + result_old.F_bound + result_old.vWA_bound +2*result_old.IF_IFclustered+2*result_old.IW_IWclustered+2*result_old.IF_IWclustered, Experiment="day25")
#result_old = result_old.assign(percentClustered = 2*result_old.clustered / result_old.Sum,
           
                  #percentActive = result_old.active / result_old.Sum,
                             #percentInactive = result_old.inactive / result_old.Sum,
                             #percentF_Bound = result_old.F_bound / result_old.Sum,
                             #percentvWA_Bound = result_old.vWA_bound / result_old.Sum)

#make one dataframe out of day18 and 25 results:
df2 = result.append(result_old)
df2.to_csv('no_competition_for_L2_simResults.csv', sep='\t', index=False)

#%% loop over df2 to plot absolute numbers after 1 min of simulation 
# saves each plot in the current directory!!
print(df2.columns)
for i, col in enumerate(df2.columns[1:10]):
    plt.figure(i)
    sns.set_style("whitegrid")
    sns.relplot(x='time', y=col, kind='line', hue= 'Experiment' ,data=df2)
    plt.ylabel(col+' uM')
    plt.xlabel("time (s)")
    plt.xticks(rotation=15)
    #plt.savefig(str(i)+col+'_absolute_day18_25.png')


#%% 
# 2 ligand 1 integrin  

# setting the ligands to different initial conditions
# differentiating between 3 types of clusters IF+IF, IW+IW and IF+IW 
# but the rates of clustering are the same between 3 clusters. only the identities are different

# set W IC to 0 to see the ligand bound integrins in a non-competitive case

Ant_str = """
  model test # activation model 

  species i, I, F, IF, W, IW, C1, C2, C3; 
  #inactive integrin, active integrin, fibronectin, integrin+fibronectin, vonWillebrand Factor A, integrin+vonWillebrand Factor A, clustered integrins respectively.
  
  #set initial values:
  i = 0.05; # integrin avB3 
  I = 0; 
  F = 0.18   ; #fibronectin 
  IF = 0;
  W = 0   ; #von Willebrand factor A 
  IW = 0;  
  C1 = 0;     # IF+IF cluster
  C2 = 0;     # IW + IW cluster
  C3 = 0;     # IF+IW cluster

  J1: i -> I; k1*i - k2*I; # reaction; reaction rate law;   # activation step, k1 rate of activation, k2 rate of inactivation
  J2: I + F -> IF; k3*I*F - k4*IF;                          # ligand binding step, k3 rate of fibronectin binding, k4 rate of dissociation
  J3: I + W -> IW; k5*I*W - k6*IW;                          # alternative ligand binding step, k5 rate of vWFA binding, k6 rate of dissociation
  J4: IF + IF -> C1; k7*IF^2 - k8*C1;
  J5: IW + IW -> C2; k7*IW^2 - k8*C2;
  #J6: IF + IW -> C3; k7*IF*IW - k8*C3;                         # clustering step, k7 rate of clustering, k8 rate of dissociation

  k1 = 5*10^6; k2 = 10^8; ; k3 = 1.6*10^8 ; k4 = 3.5*10^-1; k5 = 1.6*10^4 ; k6 = 2.3*10^-2; k7 = 1.6*10^8; k8 = 0.5*10^7; # assign constant values to global parameters
  end
  """
r2 = te.loada(Ant_str)

#%%

r2.conservedMoietyAnalysis = True
r2.steadyState()
#  3.826949977119715e-06 --> this is a good number!
 
print(r2.getSteadyStateValuesNamedArray())

#      [IF],        [IW],         [I],         [C3],        [i],      [W],     [C1],      [C2],    [F]
# [[ 0.0212096, 3.12052e-38, 3.56893e-10, 7.13786e-09, -6.28998e-38, 0.13, 0.0143952, 1.58473e-38,    0]]

#these are the steady state values.
print(r2.getRatesOfChange())
# [-3.44980136e-17 -8.73114914e-11  3.44980105e-17  0.00000000e+00 0.00000000e+00] these are the rates of change of the 6 reactions, all approaching 0 --> confirm steady state. 

#%%
# simulate for day18 initial conditions,
# store results in pandas DataFrame "result" : 
r2 = te.loada(Ant_str)
result = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'i', 'I','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'vWA_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered'])
#%%

# reset the model,
# change the initial ligand concentrations to day25, simulate,
# store results in pandas DataFrame "resul_old": 
r2.reset()
r2.F = 0.46
r2.W = 0

result_old = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'i', 'I','F','W', 'IF', 'IW',  'C1', 'C2', 'C3']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'vWA_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered'])
#create new column with total integrin amount at each time step
result = result.assign(Sum = result.inactive + result.active + result.F_bound + result.vWA_bound +2*result.IF_IFclustered+2*result.IW_IWclustered+2*result.IF_IWclustered, Experiment="day18")
#result = result.assign(percentClustered = 2*result.clustered / result.Sum,
                             #percentActive = result.active / result.Sum,
                             #percentInactive = result.inactive / result.Sum,
                             #percentF_Bound = result.F_bound / result.Sum,
                             #percentvWA_Bound = result.vWA_bound / result.Sum)

result_old = result_old.assign(Sum = result_old.inactive + result_old.active + result_old.F_bound + result_old.vWA_bound +2*result_old.IF_IFclustered+2*result_old.IW_IWclustered+2*result_old.IF_IWclustered, Experiment="day25")
#result_old = result_old.assign(percentClustered = 2*result_old.clustered / result_old.Sum,
           
                  #percentActive = result_old.active / result_old.Sum,
                             #percentInactive = result_old.inactive / result_old.Sum,
                             #percentF_Bound = result_old.F_bound / result_old.Sum,
                             #percentvWA_Bound = result_old.vWA_bound / result_old.Sum)

#make one dataframe out of day18 and 25 results:
df2 = result.append(result_old)
df2.to_csv('no_competition_for_L1_simResults.csv', sep='\t', index=False)

#%% loop over df2 to plot absolute numbers after 1 min of simulation 
# saves each plot in the current directory!!
print(df2.columns)
for i, col in enumerate(df2.columns[1:10]):
    plt.figure(i)
    sns.set_style("whitegrid")
    sns.relplot(x='time', y=col, kind='line', hue= 'Experiment' ,data=df2)
    plt.ylabel(col+' uM')
    plt.xlabel("time (s)")
    plt.xticks(rotation=15)
    #plt.savefig(str(i)+col+'_absolute_day18_25.png')