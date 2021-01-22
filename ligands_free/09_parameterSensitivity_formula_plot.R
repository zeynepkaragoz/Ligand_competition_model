###### calculate parameter sensitivity using Sun et al. formula ######

# general formula: 
#
# PS = ( abs(concentration_k_up - concentration_normal)/concentration_normal ) / (delta_k/k)
# (the same applies for k_down)
#
# in our case delta_k/k is 20% = 0.2
# 
# PS to up and down should be calculated for each molecular species (i, a, IF, IW, C1, C2, C3) and for each parameter (k1-8, i, F, W)

###### read the csv with calculated steady states ######

setwd("C:/karagoz/01-RESEARCH/01-Projects/01-In_silico_modeling_of_Integrin_function/003-Ligand_competition_model/ligands_free/08_parameter_sensitivity/")

#install + load packages

pkg<-c("magrittr", "ggplot2", "dplyr", "cowplot", 'tidyr')
lapply(pkg, require, character.only = TRUE)


# read file
steadyState_all <- read.csv("parameterSensitivity_steadyStates_all.csv", header = TRUE, sep="\t")

###### calculate parameter sensitivity ######

# for ups #
PS_values_ups <- data.frame()
for(i in 1:11){
PS_values_ups[i,1:9] <- (abs(steadyState_all[grep("up", steadyState_all$test_condition),][i,1:9]-steadyState_all[1,1:9])/ steadyState_all[1,1:9])/0.2
}
PS_values_ups$param_name <- c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "i", "L1", "L2")
PS_values_ups$up_down <- "up"

PS_values_ups_long <- gather(PS_values_ups, "mol_species","PS_value", 1:9)

# for downs # 
PS_values_downs <- data.frame()
for(i in 1:11){
  PS_values_downs[i,1:9] <- (abs(steadyState_all[grep("down", steadyState_all$test_condition),][i,1:9]-steadyState_all[1,1:9])/ steadyState_all[1,1:9])/0.2
}
PS_values_downs$param_name <- c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "i", "L1", "L2")
PS_values_downs$up_down <- "down"

PS_values_downs_long <- gather(PS_values_downs, "mol_species","PS_value", 1:9)

# append the two dfs 

PS_values_all <- rbind(PS_values_ups_long, PS_values_downs_long)
PS_values_all$param_name <- factor(PS_values_all$param_name, levels = c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "i", "L1", "L2"))
PS_values_all$up_down <- factor(PS_values_all$up_down, levels=c("up", "down"))
PS_values_all$mol_species <- substring(PS_values_all$mol_species,3)
PS_values_all$mol_species <- factor(PS_values_all$mol_species, levels = c("i.", "a.", "IF.", "IW.", "C1.","C2.","C3.","F.","W."))

#ready to plot

###### plot the parameter sensitivities ######
# x = param name
# y = PS value 
# grid cols = ups/downs 
# grid rows = molecular species

molecule_names <- c(i.="Inactive", 
                    a.="Active", 
                    IF.="L1-bound", 
                    IW.="L2-bound", 
                    C1.="L1-L1 cluster",
                    C2.="L2-L2 cluster",
                    C3.="L1-L2 cluster",
                    F.= "L1",
                    W.="L2")
up_down_label = c(up="Parameter increased 20%", 
                  down="Parameter decreased 20%")

ggplot(PS_values_all, aes(x=param_name, y=PS_value, fill=param_name, alpha=up_down)) + 
  geom_col()+
  geom_hline(yintercept = 1, color= "red")+
  facet_grid(rows = vars(mol_species), cols = vars(up_down), labeller = labeller(mol_species = molecule_names, up_down=up_down_label) )+
  #scale_y_continuous(limits = c(0,2))+
  scale_alpha_manual("up_dow", values = c(0.9, 0.5))+
  labs(x="", y="Parameter sensitivity")+
  theme_bw(base_size = 15)+
  theme(legend.position = "none", 
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18)) 


#exclude L1 - L2 sensitivities
ggplot(PS_values_all[!(PS_values_all$mol_species %in% names(molecule_names[8:9])),], aes(x=param_name, y=PS_value, fill=param_name, alpha=up_down)) + 
  geom_col()+
  geom_hline(yintercept = 1, color= "red")+
  facet_grid(rows = vars(mol_species), cols = vars(up_down), labeller = labeller(mol_species = molecule_names, up_down=up_down_label) )+
  scale_alpha_manual("up_dow", values = c(0.9, 0.5))+
  labs(x="", y="Parameter sensitivity")+
  theme_bw(base_size = 15)+
  theme(legend.position = "none", 
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=12)) 

#save pdf, landscape, 16*13
ggsave("05_figure5_parameter_sensitivity.pdf", 
       width = 18,
       height = 25,
       units = "cm")



#plot for only L1 & L2-bound
ggplot(PS_values_all[PS_values_all$mol_species %in% c("IF.", "IW.") & PS_values_all$up_down == "down",], aes(x=param_name, y=PS_value, fill=param_name, alpha=up_down)) + 
  geom_col()+
  geom_hline(yintercept = 1, color= "red")+
  facet_grid(rows = vars(mol_species), cols = vars(up_down), labeller = labeller(mol_species = molecule_names, up_down=up_down_label) )+
  #scale_y_continuous(limits = c(0,2))+
  scale_alpha_manual("up_dow", values = c(0.9, 0.5))+
  theme_bw(base_size = 15)+
  theme(legend.position = "none", 
        axis.text = element_text(size =15, face = "bold"),
        axis.text.x = element_text(angle = 15),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size=15)) +
  labs(x="", y="Parameter sensitivity")
