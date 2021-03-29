# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 22:09:01 2021

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
os.chdir(r'C:\karagoz\01-RESEARCH\01-Projects\01-In_silico_modeling_of_Integrin_function\003-Ligand_competition_model\ligands_free\10_3LigandsCompetition')
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 13}

matplotlib.rc('font', **font)


#%% 
# 3 ligands 1 integrin  

# setting the ligands to different initial conditions
# differentiating between 3 types of clusters IF+IF, IW+IW and IF+IW 
# but the rates of clustering are the same between 3 clusters. only the identities are different


# ligands: fibronectin and von Willebrand Factor A (initial values adapted from kidney orgaoid iBAQ values) and a third ligand with lowest affinity
# integrin avB3 (initial value from Hudson et al)
# here the rate for ligand binding-unbinding  is from Hudson et al (both ligands) 
# activation/inactivation is from  Yu et al. 2017
# clustering forward is numerically from Hudson, then reverse is estimated using the Kd from Yu et al .
Ant_str = """
  model test # activation model 

  species i, I, F, IF, W, IW, G, IG, C1, C2, C3, C4, C5, C6; 
  #inactive integrin, active integrin, fibronectin, integrin+fibronectin, vonWillebrand Factor A, integrin+vonWillebrand Factor A, clustered integrins respectively.
  
  #set initial values:
  i = 0.05; # integrin avB3 
  I = 0; 
  F = 0.18   ; #fibronectin 
  IF = 0;
  W = 0.33   ; #von Willebrand factor A 
  G = 0.90;
  IW = 0;  
  C1 = 0;     # IF+IF cluster
  C2 = 0;     # IW + IW cluster
  C3 = 0;     # IF+IW cluster
  C4 = 0;     # IG + IG cluster
  C5 = 0;     # IW + IG cluster 
  C6 = 0;     # IF + IG cluster


  J1: i -> I; k1*i - k2*I; # reaction; reaction rate law;   # activation step, k1 rate of activation, k2 rate of inactivation
  J2: I + F -> IF; k3*I*F - k4*IF;                          # ligand binding step, k3 rate of fibronectin binding, k4 rate of dissociation
  J3: I + W -> IW; k5*I*W - k6*IW;                          # alternative ligand binding step, k5 rate of vWFA binding, k6 rate of dissociation
  
  J4: IF + IF -> C1; k7*IF^2 - k8*C1;
  J5: IW + IW -> C2; k7*IW^2 - k8*C2;
  J6: IF + IW -> C3; k7*IF*IW - k8*C3;                         # clustering step, k7 rate of clustering, k8 rate of dissociation

  J7: I + G -> IG; k9*I*G - k10*IG;
  J8: IG + IG -> C4; k7*IG*IG - k8*C4;
  J9 : IW + IG -> C5; k7*IW*IG - k8*C5;
  J10 : IF + IG -> C6; k7*IF*IG - k8*C6;

  k1 = 5*10^6; k2 = 10^8; ; k3 = 1.6*10^8 ; k4 = 3.5*10^-1; k5 = 1.6*10^4 ; k6 = 2.3*10^-2; k7 = 1.6*10^8; k8 = 0.5*10^7; k9 = 1.6*10; k10 = 2.3*10^-2 # assign constant values to global parameters
  end
  """
r2 = te.loada(Ant_str)
#%%
r2.conservedMoietyAnalysis = True
r2.steadyState()
#  3.5906061591688018e-06 --> this is a good number!
 
print(r2.getSteadyStateValuesNamedArray())

#     [IF],        [IW],        [IG],         [I],        [C3],         [i],        [C6],        [C5],       [W],      [F],      [C1],       [C4],        [C2],     [G]
# [[ 0.0211576, 8.16089e-05, 2.22662e-07, 3.55641e-10, 5.52528e-05, 7.11282e-09, 1.50752e-07, 5.81479e-10, 0.329863, 0.130138, 0.0143246, 1.58652e-12, 2.1312e-07, 0.9]]
#these are the steady state values.

print(r2.getRatesOfChange())
#[ 2.77555756e-17 -3.75459663e-11  1.56067455e-10  3.26736342e-23  0.00000000e+00  0.00000000e+00 -2.77555756e-17  0.00000000e+00  0.00000000e+00  0.00000000e+00]


#%%
# simulate for day18 initial conditions,
# store results in pandas DataFrame "result" : 
r2 = te.loada(Ant_str)
result = pd.DataFrame(r2.simulate(0, 0.000016 , 100 , ['time', 'i', 'I','F', 'W', 'G', 'IF', 'IW', 'IG', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6' ]), 
                      columns=['time', 'inactive', 'active','F','W', 'G','F_bound', 'vWA_bound','G_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered', 'IG_IGclustered', 'IW_IGclustered', 'IF_IGclustered' ])


# reset the model,
# change the initial ligand concentrations to day25, simulate,
# store results in pandas DataFrame "resul_old": 
r2.reset()
r2.F = 0.46
r2.W = 0.50
r2.G = 2.7

result_old = pd.DataFrame(r2.simulate(0, 0.000016 , 100 , ['time', 'i', 'I','F', 'W', 'G', 'IF', 'IW', 'IG', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6']),
                          columns=['time', 'inactive', 'active','F','W', 'G','F_bound', 'vWA_bound','G_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered', 'IG_IGclustered', 'IW_IGclustered', 'IF_IGclustered'])
#create new column with total integrin amount at each time step
result = result.assign(Sum = result.inactive + result.active + result.F_bound + result.vWA_bound + result.G_bound +2*result.IF_IFclustered+2*result.IW_IWclustered+2*result.IF_IWclustered+2*result.IG_IGclustered+2*result.IW_IGclustered+2*result.IF_IGclustered, Experiment="day18")
#result = result.assign(percentClustered = 2*result.clustered / result.Sum,
                             #percentActive = result.active / result.Sum,
                             #percentInactive = result.inactive / result.Sum,
                             #percentF_Bound = result.F_bound / result.Sum,
                             #percentvWA_Bound = result.vWA_bound / result.Sum)

result_old = result_old.assign(Sum = result_old.inactive + result_old.active + result_old.F_bound + result_old.vWA_bound + result_old.G_bound + 2*result_old.IF_IFclustered+2*result_old.IW_IWclustered+2*result_old.IF_IWclustered+2*result_old.IG_IGclustered+2*result_old.IW_IGclustered+2*result_old.IF_IGclustered, Experiment="day25")
#result_old = result_old.assign(percentClustered = 2*result_old.clustered / result_old.Sum,
           
                  #percentActive = result_old.active / result_old.Sum,
                             #percentInactive = result_old.inactive / result_old.Sum,
                             #percentF_Bound = result_old.F_bound / result_old.Sum,
                             #percentvWA_Bound = result_old.vWA_bound / result_old.Sum)

#make one dataframe out of day18 and 25 results:
df2 = result.append(result_old)
df2.to_csv('3LigandsModel_lowAff_simResults.csv', sep='\t', index=False)

#%% 
# 3 ligands 1 integrin  

# setting the ligands to different initial conditions
# differentiating between 3 types of clusters IF+IF, IW+IW and IF+IW 
# but the rates of clustering are the same between 3 clusters. only the identities are different


# ligands: fibronectin and von Willebrand Factor A (initial values adapted from kidney orgaoid iBAQ values) and a third ligand with lowest affinity
# integrin avB3 (initial value from Hudson et al)
# here the rate for ligand binding-unbinding  is from Hudson et al (both ligands) 
# activation/inactivation is from  Yu et al. 2017
# clustering forward is numerically from Hudson, then reverse is estimated using the Kd from Yu et al .
Ant_str = """
  model test # activation model 

  species i, I, F, IF, W, IW, G, IG, C1, C2, C3, C4, C5, C6; 
  #inactive integrin, active integrin, fibronectin, integrin+fibronectin, vonWillebrand Factor A, integrin+vonWillebrand Factor A, clustered integrins respectively.
  
  #set initial values:
  i = 0.05; # integrin avB3 
  I = 0; 
  F = 0.33   ; #fibronectin 
  IF = 0;
  W = 0.33   ; #von Willebrand factor A 
  G = 0.33;
  IW = 0;  
  C1 = 0;     # IF+IF cluster
  C2 = 0;     # IW + IW cluster
  C3 = 0;     # IF+IW cluster
  C4 = 0;     # IG + IG cluster
  C5 = 0;     # IW + IG cluster 
  C6 = 0;     # IF + IG cluster


  J1: i -> I; k1*i - k2*I; # reaction; reaction rate law;   # activation step, k1 rate of activation, k2 rate of inactivation
  J2: I + F -> IF; k3*I*F - k4*IF;                          # ligand binding step, k3 rate of fibronectin binding, k4 rate of dissociation
  J3: I + W -> IW; k5*I*W - k6*IW;                          # alternative ligand binding step, k5 rate of vWFA binding, k6 rate of dissociation
  
  J4: IF + IF -> C1; k7*IF^2 - k8*C1;
  J5: IW + IW -> C2; k7*IW^2 - k8*C2;
  J6: IF + IW -> C3; k7*IF*IW - k8*C3;                         # clustering step, k7 rate of clustering, k8 rate of dissociation

  J7: I + G -> IG; k9*I*G - k10*IG;
  J8: IG + IG -> C4; k7*IG*IG - k8*C4;
  J9 : IW + IG -> C5; k7*IW*IG - k8*C5;
  J10 : IF + IG -> C6; k7*IF*IG - k8*C6;

  k1 = 5*10^6; k2 = 10^8; ; k3 = 1.6*10^8 ; k4 = 3.5*10^-1; k5 = 1.6*10^4 ; k6 = 2.3*10^-2; k7 = 1.6*10^8; k8 = 0.5*10^7; k9 = 1.6*10^10; k10 = 2.3*10^-1 # assign constant values to global parameters
  end
  """
r2 = te.loada(Ant_str)
#%%
r2.conservedMoietyAnalysis = True
r2.steadyState()
#  3.5906061591688018e-06 --> this is a good number!
 
print(r2.getSteadyStateValuesNamedArray())

#      [IF],      [IW],      [IG],         [I],        [C3],         [i],        [C6],        [C5],       [W],      [F],      [C1],      [C4],        [C2],      [G]
# [[ 0.000163164, 2.485e-07, 0.0211056, 1.08248e-12, 1.29748e-09, 2.16496e-11, 0.000110197, 1.67832e-07, 0.33, 0.329725, 8.51919e-07, 0.0142543, 1.97606e-12, 0.280275]]
#these are the steady state values.

print(r2.getRatesOfChange())
#[-7.88702428e-11 -9.24282909e-11  1.13575815e-13 -4.06575815e-20 -8.67361738e-19  2.71050543e-20 -1.13686838e-13  1.11022302e-16  0.00000000e+00 -1.35525272e-20]
#%%
# simulate for day18 initial conditions,
# store results in pandas DataFrame "result" : 
r2 = te.loada(Ant_str)
result = pd.DataFrame(r2.simulate(0, 0.000016 , 100 , ['time', 'i', 'I','F', 'W', 'G', 'IF', 'IW', 'IG', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6' ]), 
                      columns=['time', 'inactive', 'active','F','W', 'G','F_bound', 'vWA_bound','G_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered', 'IG_IGclustered', 'IW_IGclustered', 'IF_IGclustered' ])


# reset the model,
# change the initial ligand concentrations to day25, simulate,
# store results in pandas DataFrame "resul_old": 
r2.reset()
r2.F = 0.50
r2.W = 0.50
r2.G = 0.50

result_old = pd.DataFrame(r2.simulate(0, 0.000016 , 100 , ['time', 'i', 'I','F', 'W', 'G', 'IF', 'IW', 'IG', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6']),
                          columns=['time', 'inactive', 'active','F','W', 'G','F_bound', 'vWA_bound','G_bound', 'IF_IFclustered', 'IW_IWclustered', 'IF_IWclustered', 'IG_IGclustered', 'IW_IGclustered', 'IF_IGclustered'])
#create new column with total integrin amount at each time step
result = result.assign(Sum = result.inactive + result.active + result.F_bound + result.vWA_bound + result.G_bound +2*result.IF_IFclustered+2*result.IW_IWclustered+2*result.IF_IWclustered+2*result.IG_IGclustered+2*result.IW_IGclustered+2*result.IF_IGclustered, Experiment="day18")
#result = result.assign(percentClustered = 2*result.clustered / result.Sum,
                             #percentActive = result.active / result.Sum,
                             #percentInactive = result.inactive / result.Sum,
                             #percentF_Bound = result.F_bound / result.Sum,
                             #percentvWA_Bound = result.vWA_bound / result.Sum)

result_old = result_old.assign(Sum = result_old.inactive + result_old.active + result_old.F_bound + result_old.vWA_bound + result_old.G_bound + 2*result_old.IF_IFclustered+2*result_old.IW_IWclustered+2*result_old.IF_IWclustered+2*result_old.IG_IGclustered+2*result_old.IW_IGclustered+2*result_old.IF_IGclustered, Experiment="day25")
#result_old = result_old.assign(percentClustered = 2*result_old.clustered / result_old.Sum,
           
                  #percentActive = result_old.active / result_old.Sum,
                             #percentInactive = result_old.inactive / result_old.Sum,
                             #percentF_Bound = result_old.F_bound / result_old.Sum,
                             #percentvWA_Bound = result_old.vWA_bound / result_old.Sum)

#make one dataframe out of day18 and 25 results:
df2 = result.append(result_old)
df2.to_csv('3LigandsModel_highAff_simResults.csv', sep='\t', index=False)