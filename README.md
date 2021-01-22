# Ligand_competition_model

Use scripts in Ligands_free folder to simulate the model, create data and plot the results used in the manuscript. 

1) simulations for different test conditions: \
use \
1- 01_different_IC_different_FC_ligands_free.py for ligand competition different initial conditions \
2- 02_equal_IC_different_FC_ligands_free.py for ligand competition equal initial conditions \
3- 03_equal_IC_equal_FC_ligands_free.py for ligand competition equal initial conditions, equal fold change \
4- 04_equal_IC_high_FC_ligands_free.py for ligand competition equal intial condition, high fold change for L2 \
5- 05_different_IC_equal_BR_ligands_free.py for ligand competition equal initial conditions, equal binding rates \
all these scripts export the simulation results in csv files.  

read these files in: \
06_plots_simulations.R to plot the results.  

2) parameter scan of L2 initial condition (manual) \
use \
07_paramScan_manual_L2.py to do the parameter scan for L2=0.5 to L2=2.97 \
exports to csv \
then use \
06_plots_simulations.R to plot the results.  

3) Parameter Sensitivity analysis of the original ligand competition model with excess ligands compared to integrins: \
use \
08_parameterSensitivity_calculations.py \
then use \
09_parameterSensitivity_formula_plot.R to plot.  

4) Parameter Sensitivity of the ligand competition model with excess integrins compared to ligands: \
use \
10_parameterSensitivity_calculations_integrin_excess.py  
then the script below plots the supplementary figure 2 \
11_parameterSensitivity_formula_plot_integrins_excess.R \

5) extras: \
12_no_competition.py : simulates the model without competition    
equations.ipynb : prints out the differential equations of the model  

