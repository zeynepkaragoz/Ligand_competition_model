# Ligand_competition_model

1) simulations for different test conditions: \
use \
1- uniqueClusters_fibronectin_vWA_competition_Yu_Hudson_parameters.py for ligand competition different IC \
2- uniqueClusters_equal_IC_fibronectin_vWA_competition_Yu_Hudson_parameters.py for ligand competition equal IC \
3- uniqueClusters_FoldChange_test_equal_IC_fibronectin_vWA_competition_Yu_Hudson_parameters.py for ligand competition equal IC, equal FC \
4- uniqueClusters_equal_IC_fibronectin_THBS_competition_Yu_Hudson_parameters.py for ligand competition equal IC, high FC \
5- uniqueClusters_equal_bindingRates_equal_IC_fibronectin_vWA_competition_Yu_Hudson_parameters.py for ligand competition equal IC, equal BR \
all these scripts export the simulation results in csv files. \
\
read these files in: \
plots_for_publication_v2.R to plot the results. \
\
2) parameter scan of L2 initial condition (manual) \
use \
paramScan_manual_vWA.py to do the parameter scan for L2=0.5 to L2=2.97 \
exports to csv \
then use \
plots_for_publication_v2.R to plot the results. \
\ 
3) uncertainty using util from tellurium: \
use \
paramScan_utility_saveCSV.py for uncertainty analysis, saves as csv. \
then use \
paramScan_utility_plots_for_publication.R to plot. \




