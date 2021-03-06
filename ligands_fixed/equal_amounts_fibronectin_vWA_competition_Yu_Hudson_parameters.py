# -*- coding: utf-8 -*-
"""
Created on Mon May 25 09:08:39 2020

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
os.chdir(r'C:\karagoz\01-RESEARCH\01-Projects\01-In_silico_modeling_of_Integrin_function\003-Ligand_competition_model\figures\Yu_Hudson_params\aVB3-FN-vWFA\ligand_IC_equal')
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 13}

matplotlib.rc('font', **font)
#%% 
# 2 ligand 1 integrin  

# testing if the initial concentration of ligands have an effect on the steady state? 

# ligands: fibronectin and von Willebrand Factor A (initial values adapted from kidney orgaoid iBAQ values)
# integrin avB3 (initial value from Hudson et al)
# here the rate for ligand binding-unbinding  is from Hudson et al (both ligands) 
# activation/inactivation is from  Yu et al. 2017
# clustering forward is numerically from Hudson, then reverse is estimated using the Kd from Yu et al .
Ant_str = """
  model test # activation model 

  species i, I, $F, IF, $W, IW, C; 
  #inactive integrin, active integrin, fibronectin, integrin+fibronectin, vonWillebrand Factor A, integrin+vonWillebrand Factor A, clustered integrins respectively.
  
  #set initial values:
  i = 0.05; # integrin avB3 
  I = 0; 
  F = 0.33   ; #fibronectin (set the same as vWA)
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

  k1 = 5*10^6; k2 = 10^8; ; k3 = 1.6*10^8 ; k4 = 3.5*10^-1; k5 = 1.6*10^4 ; k6 = 2.3*10^-2; k7 = 1.6*10^8; k8 = 0.5*10^7; # assign constant values to global parameters
  end
  """
r2 = te.loada(Ant_str)

#k1 = 5*10^6; k2 = 10^8; k3 = 1.6*10^8; k4 = 3.5*10^-1; k5 = 1.6*10^8; k6 = 0.5*10^7; #rates of the ABC model adjusted for 60%clustering

#te.getODEsFromModel(r2)
#vJ1 = k1*i-k2*I 
#vJ2 = k3*I*F-k4*IF 
#vJ3 = k5*I*W-k6*IW
#vJ4 = k7*pow(IF,2)-k8*C 
#vJ5 = k7*pow(IW,2)-k8*C 
#vJ6 = k7*IF*IW-k8*C 

#di/dt = -vJ1 
#dI/dt = vJ1 - vJ2 - vJ3 
#dIF/dt = vJ2 - 2.0*vJ4 - vJ6 
#dIW/dt = vJ3 - 2.0*vJ5 - vJ6 
#dC/dt = vJ4 + vJ5 + vJ6

#%%

r2.conservedMoietyAnalysis = True
r2.steadyState()
#  6.8141049931359555e-19 --> this is a good number!
 
print(r2.getSteadyStateValuesNamedArray())
#     [IF],     [IW],         [I],         [i],        [C]
#[[ 0.0163967, 0.0163967, 2.12321e-10, 4.24643e-09, 0.00860328]]
#these are the steady state values.
print(r2.getRatesOfChange())
# [ 2.77319834e-13 -2.77318750e-13 -4.55364912e-18  3.46944695e-18] these are the rates of change of the 5 , all approaching 0 --> confirm steady state. 

#simulate for 3 mins!!! = 180s

#%%
# simulate for day18 initial conditions,
# store results in pandas DataFrame "result" : 
r2 = te.loada(Ant_str)
result = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'i', 'I','F','W', 'IF', 'IW',  'C']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'vWA_bound', 'clustered'])


# reset the model,
# change the initial ligand concentrations to day25, simulate,
# store results in pandas DataFrame "resul_old": 
r2.reset()
r2.F = 0.46
r2.W = 0.50

result_old = pd.DataFrame(r2.simulate(0, 0.00001 , 100 , ['time', 'i', 'I','F','W', 'IF', 'IW',  'C']), columns=['time', 'inactive', 'active','F','W', 'F_bound', 'vWA_bound', 'clustered'])
#create new column with total integrin amount at each time step
result = result.assign(Sum = result.inactive + result.active + result.F_bound + result.vWA_bound +2*result.clustered, Experiment="day18")
result = result.assign(percentClustered = 2*result.clustered / result.Sum,
                             percentActive = result.active / result.Sum,
                             percentInactive = result.inactive / result.Sum,
                             percentF_Bound = result.F_bound / result.Sum,
                             percentvWA_Bound = result.vWA_bound / result.Sum)

result_old = result_old.assign(Sum = result_old.inactive + result_old.active + result_old.F_bound + result_old.vWA_bound +2*result_old.clustered, Experiment="day25")
result_old = result_old.assign(percentClustered = 2*result_old.clustered / result_old.Sum,
                             percentActive = result_old.active / result_old.Sum,
                             percentInactive = result_old.inactive / result_old.Sum,
                             percentF_Bound = result_old.F_bound / result_old.Sum,
                             percentvWA_Bound = result_old.vWA_bound / result_old.Sum)

#make one dataframe out of day18 and 25 results:
df2 = result.append(result_old)

#%% loop over df2 to plot absolute numbers after 1 min of simulation 
# saves each plot in the current directory!!

for i, col in enumerate(df2.columns[1:8]):
    plt.figure(i)
    sns.relplot(x='time', y=col, kind='line', hue= 'Experiment' ,data=df2)
    plt.ylabel(col+' uM')
    plt.xticks(rotation=45)
    plt.savefig(col+'_absolute_day18_25.png')

#%% loop over df2 to plot percentage of each integrin-related species after 1 min of simulation
# saves each plot in the current directory!!

for i, col in enumerate(df2.columns[10:16]):
    plt.figure(i)
    sns.relplot(x='time', y=col, kind='line', hue= 'Experiment' ,data=df2)
    plt.ylim(0, 1)
    plt.ylabel(col+' integrins')
    plt.xticks(rotation=45)
    plt.savefig(col+'_day18_25.png')
