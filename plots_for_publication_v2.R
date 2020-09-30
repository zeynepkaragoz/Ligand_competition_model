setwd("C:/karagoz/01-RESEARCH/01-Projects/01-In_silico_modeling_of_Integrin_function/003-Ligand_competition_model/Unique_clusters/pretty_plots")

#install + load packages

pkg<-c("magrittr", "ggplot2", "dplyr", "cowplot")
#install.packages(pkg)
lapply(pkg, require, character.only = TRUE)


# read files
equalIC_res <- read.csv("../figures_FN_vWF_equal_IC/equalIC_simResults.csv", header = TRUE, sep="\t")
diffIC_res <- read.csv("../figures_FN_vWF/diffIC_simResults.csv", header = TRUE, sep = "\t")

THBS_equalIC_res <- read.csv("../figures_FN_THBS_equal_IC/F_THBS_equalIC_simResults.csv", header = TRUE, sep="\t")
colnames(THBS_equalIC_res) <- colnames(diffIC_res)
equalFC_res <- read.csv("../figures_FN_vWF_foldChange_testing/FN_vWF_equal_foldChange/F_vWA_equalFC_simResults.csv", header = TRUE, sep="\t")
colnames(equalFC_res) <- colnames(diffIC_res)
#equal_bindingRate_res <- read.csv("../figures_FN_vWF_foldChange_testing/FN_vWF_equal_FC_equal_bindingRates/F_vWA_equalFC_equal_bindingRates_simResults.csv", header = TRUE, sep="\t")
#colnames(equal_bindingRate_res) <- colnames(diffIC_res)

equal_bindingRate_res <- read.csv("../figures_FN_vWF_equalBR_differentIC/F_vWA_equal_bindingRates_differentIC_differentFC_simResults.csv", header = TRUE, sep="\t")
colnames(equal_bindingRate_res) <- colnames(diffIC_res)

# concat. all dfs in one:
equalIC_res$test <- "Equal IC"
diffIC_res$test <- "Different IC"
THBS_equalIC_res$test <- "Equal IC, High FC"
equalFC_res$test <- "Equal IC, Equal FC"
equal_bindingRate_res$test <- "Different IC, Equal BR"

full_DF <- rbind(equalIC_res, diffIC_res, THBS_equalIC_res, equalFC_res)

factors <- c("Different IC","Equal IC","Equal IC, Equal FC",  "Equal IC, High FC" )

full_DF$test <- factor(full_DF$test, levels=factors)


# get legend 
legend <- get_legend(ggplot(full_DF,aes(x=time, y=F_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-bound Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "top", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) +
  scale_y_continuous( limits=c(0,0.022)) + 
  facet_grid(rows = vars(test)))

##### Figure 1: L1 & L2 bound integrins, all tests ###################

L1_bound_grid <- ggplot(full_DF,aes(x=time, y=F_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-bound Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) +
  scale_y_continuous( limits=c(0,0.022)) + 
  facet_grid(rows = vars(test))

L2_bound_grid <- ggplot(full_DF,aes(x=time, y=vWA_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0,1.5e-5))+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) +
  facet_grid(rows = vars(test))

eq_BR_L1_bound <- ggplot(equal_bindingRate_res,aes(x=time, y=F_bound, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0, 0.015))+
  labs(x = "Time (s)",
       y = expression(paste("L1-bound Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size= 14, face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size =13, face = "bold"))+
  facet_grid(rows = vars(test))

eq_BR_L2_bound <- ggplot(equal_bindingRate_res,aes(x=time, y=vWA_bound, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0, 0.015))+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size= 14, face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size =13, face = "bold")) + 
  facet_grid(rows = vars(test))

top_row <- plot_grid(L1_bound_grid, L2_bound_grid, ncol = 2, labels = c("A", "B"))
bot_row <- plot_grid(eq_BR_L1_bound, eq_BR_L2_bound, ncol=2, labels = c("C", "D"))

plot_grid(legend,top_row, bot_row, nrow=3, ncol = 1,rel_heights = c(0.4,4,1.4))
#save pdf 15x9




##### Figure 2: L2 bound integrin Initial Concentration test ########
#install.packages("colorspace")
library(colorspace)
paramScan_L2 <- read.csv("../figures_FN_vWF_equal_IC/vWA_paramScan.csv", header = TRUE, sep="\t")
paramScan_L2$grad <- as.character(paramScan_L2$grad)
ggplot(paramScan_L2,aes(x=time, y=vWA_bound, group=grad)) + 
  geom_line(aes(color=grad),size=1.5) +
  scale_color_discrete_diverging(palette="Green-Orange")+
  scale_y_continuous(limits = c(0,1.2e-5))+
  geom_hline(yintercept = 2.977955e-06, color="black", size=1.5, linetype="dotted" )+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Integrin ( ", mu,"M)",sep = "")),
       color= "L2 initial concentration:")+
  theme_bw(base_size = 13)+
  theme(legend.position = "right", 
        legend.title = element_text(size=15), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) 
# save pdf 6x9




##### Figure 3: Integrin clusters, all tests
L1_L1_cluster_grid <- ggplot(full_DF,aes(x=time, y=IF_IFclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-L1 Bound Clusters ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(angle = 15)) +
  scale_y_continuous( limits=c(0,0.022)) + 
  facet_grid(rows = vars(test))

L2_L2_cluster_grid <- ggplot(full_DF,aes(x=time, y=IW_IWclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L2-L2 Bound Clusters ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(angle = 15)) +
  scale_y_continuous( limits=c(0,4e-9)) + 
  facet_grid(rows = vars(test))

L1_L2_cluster_grid <- ggplot(full_DF,aes(x=time, y=IF_IWclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-L2 Bound Clusters ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(angle = 15)) +
  scale_y_continuous( limits=c(0, 8e-6)) + 
  facet_grid(rows = vars(test))

eq_BR_L1_L1cluster <- ggplot(equal_bindingRate_res,aes(x=time, y=IF_IFclustered, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0, 0.008))+
  labs(x = "Time (s)",
       y = expression(paste("L1-L1 Bound Cluster ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size= 14, face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size =13, face = "bold"),
        axis.text.x = element_text(angle = 15))+
  facet_grid(rows = vars(test))

eq_BR_L2_L2cluster <- ggplot(equal_bindingRate_res,aes(x=time, y=IW_IWclustered, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0, 0.008))+
  labs(x = "Time (s)",
       y = expression(paste("L2-L2 Bound Cluster ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size= 14, face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size =13, face = "bold"),
        axis.text.x = element_text(angle = 15))+
  facet_grid(rows = vars(test))

eq_BR_L1_L2cluster <- ggplot(equal_bindingRate_res,aes(x=time, y=IF_IWclustered, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0, 0.008))+
  labs(x = "Time (s)",
       y = expression(paste("L1-L2 Bound Cluster ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size= 14, face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size =13, face = "bold"),
        axis.text.x = element_text(angle = 15))+
  facet_grid(rows = vars(test))

cluster_top_row <- plot_grid(L1_L1_cluster_grid, L2_L2_cluster_grid, L1_L2_cluster_grid, labels = c("A", "B", "C"), ncol=3)
cluster_bottom_row <- plot_grid(eq_BR_L1_L1cluster, eq_BR_L2_L2cluster, eq_BR_L1_L2cluster, labels = c("D", "E", "F"), ncol = 3)
plot_grid(legend,cluster_top_row, cluster_bottom_row, nrow = 3, rel_heights = c(0.2, 4, 1.4))
#save pdf 16x13