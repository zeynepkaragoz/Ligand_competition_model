setwd("C:/karagoz/01-RESEARCH/01-Projects/01-In_silico_modeling_of_Integrin_function/003-Ligand_competition_model/ligands_free/06_plots_simulations_out")

#install + load packages

pkg<-c("magrittr", "ggplot2", "dplyr", "cowplot")
#install.packages(pkg)
lapply(pkg, require, character.only = TRUE)


# read files
diffIC_diffFC_res <- read.csv("../01_different_IC_different_FC/diffIC_diffFC_simResults.csv", header = TRUE, sep="\t")

equalIC_diffFC_res <- read.csv("../02_equal_IC_different_FC/equalIC_diffFC_simResults.csv", header = TRUE, sep = "\t")

equalIC_equalFC_res <- read.csv("../03_equal_IC_equal_FC/equalIC_equalFC_simResults.csv", header = TRUE, sep="\t")
colnames(equalIC_equalFC_res) <- colnames(diffIC_diffFC_res)

equalIC_highFC_res <- read.csv("../04_equal_IC_high_FC/equalIC_highFC_simResults.csv", header = TRUE, sep="\t")
colnames(equalIC_highFC_res) <- colnames(diffIC_diffFC_res)

diffIC_equalBR_res <- read.csv("../05_different_IC_equal_BR/differentIC_equalBR_simResults.csv", header = TRUE, sep="\t")
colnames(diffIC_equalBR_res) <- colnames(diffIC_diffFC_res)

# concat. all dfs in one:

diffIC_diffFC_res$test <- "Different IC, Different FC"
equalIC_diffFC_res$test <- "Equal IC, Different FC"
equalIC_equalFC_res$test <- "Equal IC, Equal FC"
equalIC_highFC_res$test <- "Equal IC, High FC"
diffIC_equalBR_res$test <- "Different IC, Equal BR"

full_DF <- rbind(diffIC_diffFC_res, equalIC_diffFC_res, equalIC_equalFC_res, equalIC_highFC_res)

factors <- c("Different IC, Different FC","Equal IC, Different FC","Equal IC, Equal FC",  "Equal IC, High FC" )

full_DF$test <- factor(full_DF$test, levels=factors)

# define labels for tests to plot

test_labels <- c("Different IC, Different FC" = "1", 
                 "Equal IC, Different FC"="2", 
                 "Equal IC, Equal FC"="3", 
                 "Equal IC, High FC"="4", 
                 "Different IC, Equal BR"="5")
                   
# get legend 
legend <- get_legend(ggplot(full_DF,aes(x=time, y=F_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-bound Integrin ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "top", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(3, units = "cm"),
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
       y = expression(paste("L1-bound Integrin ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
  facet_grid(rows = vars(test),labeller = labeller(test = test_labels))+
  theme(legend.position = "none", 
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18)) +
  scale_y_continuous( limits=c(0,0.022)) 

L2_bound_grid <- ggplot(full_DF,aes(x=time, y=vWA_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0,1.5e-5))+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Integrin ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
  facet_grid(rows = vars(test), labeller = labeller(test = test_labels))+
  theme(legend.position = "none", 
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18)) #change here!!!

eq_BR_L1_bound <- ggplot(diffIC_equalBR_res,aes(x=time, y=F_bound, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0, 0.015))+
  labs(x = "Time (s)",
       y = expression(paste("L1-bound Integrin ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
  facet_grid(rows = vars(test), labeller = labeller(test = test_labels))+
  theme(legend.position = "none", 
        axis.text = element_text(size =16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18))

eq_BR_L2_bound <- ggplot(diffIC_equalBR_res,aes(x=time, y=vWA_bound, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0, 0.015))+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Integrin ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
  facet_grid(rows = vars(test), labeller = labeller(test = test_labels))+
  theme(legend.position = "none", 
        axis.text = element_text(size =16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18)) 

top_row <- plot_grid(L1_bound_grid,NULL, 
                     L2_bound_grid, 
                     ncol = 3, 
                     rel_widths = c(1,0.1,1), 
                     labels = c("A","", "B"), 
                     label_size = 20)
bot_row <- plot_grid(eq_BR_L1_bound,NULL, 
                     eq_BR_L2_bound, ncol=3, 
                     rel_widths = c(1,0.1,1),  
                     labels = c("C","", "D"), 
                     label_size = 20)

plot_grid(legend,top_row,NULL, bot_row, nrow=4, ncol = 1,rel_heights = c(0.4,4,0.2,1.4))

ggsave("01_figure2_ligand_bound_integrins.pdf",
       height = 19,
       width = 12,
       units = "in")





##### Figure 2: L2 bound integrin Initial Concentration test ########
#install.packages("colorspace")
library(colorspace)
paramScan_L2 <- read.csv("../07_parameter_scan_L2/L2_paramScan.csv", header = TRUE, sep="\t")
paramScan_L2$grad <- as.character(paramScan_L2$grad)
ggplot(paramScan_L2,aes(x=time, y=vWA_bound, group=grad)) + 
  geom_line(aes(color=grad),size=1.5) +
  scale_color_discrete_diverging(palette="Green-Orange", guide = guide_legend(reverse = TRUE))+
  scale_y_continuous(limits = c(0,1.2e-5))+
  geom_hline(yintercept = 2.977955e-06, color="black", size=1.5, linetype="dotted" )+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Integrin ( nM )",sep = "")),
       color= "L2 initial concentration (nM):")+
  theme_bw(base_size = 13)+
  theme(legend.position = "right", 
        legend.title = element_text(size=15, face="bold"), 
        legend.text = element_text(size = 16,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size =16, face = "bold"),
        axis.title = element_text(size = 20, face = "bold")) 

ggsave("02_figure3_L2_paramScan.pdf", 
       width = 10,
       height = 8,
       units = "in",)


##### Figure 3: Integrin clusters, all tests ######
L1_L1_cluster_grid <- ggplot(full_DF,aes(x=time, y=IF_IFclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-L1 Bound Clusters ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
  scale_y_continuous( limits=c(0,0.022)) + 
  facet_grid(rows = vars(test),labeller = labeller(test = test_labels))+
  theme(legend.position = "none", 
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18)) 

L2_L2_cluster_grid <- ggplot(full_DF,aes(x=time, y=IW_IWclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L2-L2 Bound Clusters ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
  scale_y_continuous( limits=c(0,4e-9)) + 
  facet_grid(rows = vars(test),labeller = labeller(test = test_labels))+
  theme(legend.position = "none", 
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18)) 

L1_L2_cluster_grid <- ggplot(full_DF,aes(x=time, y=IF_IWclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-L2 Bound Clusters ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
  facet_grid(rows = vars(test),labeller = labeller(test = test_labels))+
  theme(legend.position = "none", 
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18)) +
  scale_y_continuous( limits=c(0, 8e-6)) 

eq_BR_L1_L1cluster <- ggplot(diffIC_equalBR_res,aes(x=time, y=IF_IFclustered, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0, 0.008))+
  labs(x = "Time (s)",
       y = expression(paste("L1-L1 Bound Cluster ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
  facet_grid(rows = vars(test),labeller = labeller(test = test_labels))+
  theme(legend.position = "none", 
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18))

eq_BR_L2_L2cluster <- ggplot(diffIC_equalBR_res,aes(x=time, y=IW_IWclustered, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0, 0.008))+
  labs(x = "Time (s)",
       y = expression(paste("L2-L2 Bound Cluster ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
  facet_grid(rows = vars(test),labeller = labeller(test = test_labels))+
  theme(legend.position = "none", 
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18))

eq_BR_L1_L2cluster <- ggplot(diffIC_equalBR_res,aes(x=time, y=IF_IWclustered, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0, 0.008))+
  labs(x = "Time (s)",
       y = expression(paste("L1-L2 Bound Cluster ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
  facet_grid(rows = vars(test),labeller = labeller(test = test_labels))+
  theme(legend.position = "none", 
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18))

cluster_top_row <- plot_grid(L1_L1_cluster_grid,NULL, 
                             L2_L2_cluster_grid,NULL, 
                             L1_L2_cluster_grid, 
                             rel_widths = c(1,0.1,1,0.1,1) ,
                             labels = c("A","", "B","", "C"), 
                             label_size = 20,
                             nrow=1)

cluster_bottom_row <- plot_grid(eq_BR_L1_L1cluster,NULL, 
                                eq_BR_L2_L2cluster,NULL, 
                                eq_BR_L1_L2cluster, 
                                rel_widths = c(1,0.1,1,0.1,1),
                                labels = c("D","", "E","", "F"),
                                label_size = 20,
                                nrow = 1)

plot_grid(legend,cluster_top_row, NULL, cluster_bottom_row, nrow = 4, rel_heights = c(0.2, 4,0.1, 1.4))

ggsave("03_figure4_integrin_clusters1.pdf", 
       width = 17,
       height = 19,
       units = "in")


######### Fig S active-inactive integrins ########

inactive_grid <- ggplot(full_DF,aes(x=time, y=inactive, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("Inactive Integrin ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
 # scale_y_continuous( limits=c(0,0.022)) + 
  facet_grid(rows = vars(test),labeller = labeller(test = test_labels))+
  theme(legend.position = "none", 
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18)) 

active_grid <- ggplot(full_DF,aes(x=time, y=active, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  #scale_y_continuous(limits = c(0,1.5e-5))+
  labs(x = "Time (s)",
       y = expression(paste("Active Integrin ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
  facet_grid(rows = vars(test),labeller = labeller(test = test_labels))+
  theme(legend.position = "none", 
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18)) 

eq_BR_inactive <- ggplot(diffIC_equalBR_res,aes(x=time, y=inactive, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  #scale_y_continuous(limits = c(0, 0.015))+
  labs(x = "Time (s)",
       y = expression(paste("Inactive Integrin ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
  facet_grid(rows = vars(test),labeller = labeller(test = test_labels))+
  theme(legend.position = "none", 
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18))

eq_BR_active <- ggplot(diffIC_equalBR_res,aes(x=time, y=active, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  #scale_y_continuous(limits = c(0, 0.015))+
  labs(x = "Time (s)",
       y = expression(paste("Active Integrin ( nM )",sep = "")))+
  theme_bw(base_size = 13)+
  facet_grid(rows = vars(test),labeller = labeller(test = test_labels))+
  theme(legend.position = "none", 
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size=18)) 

top_row <- plot_grid(inactive_grid, NULL,active_grid, nrow = 1, labels = c("A","", "B"), rel_widths = c(1,0.1,1))
bot_row <- plot_grid(eq_BR_inactive, NULL, eq_BR_active, nrow=1, labels = c("C","", "D"),  rel_widths = c(1,0.1,1))

plot_grid(legend,top_row, NULL, bot_row, nrow=4, ncol = 1,rel_heights = c(0.4,4,0.2,1.4))

ggsave("04_figureS1_active_inactive_integrins1.pdf",
       height = 19,
       width = 12,
       units = "in")

#### Figures for presentation ####

# no ligand competition #

# L1 initial concentration = 0, L2 IC = 0.33 #

# read files

no_competition_L2_res <- read.csv("../09_no_competition/no_competition_for_L2_simResults.csv", header = TRUE, sep="\t")

no_competition_L1_res <- read.csv("../09_no_competition/no_competition_for_L1_simResults.csv", header = TRUE, sep="\t")

diffIC_diffFC_res <- read.csv("../01_different_IC_different_FC/diffIC_diffFC_simResults.csv", header = TRUE, sep="\t")

diffIC_equalBR_res <- read.csv("../05_different_IC_equal_BR/differentIC_equalBR_simResults.csv", header = TRUE, sep="\t")
colnames(diffIC_equalBR_res) <- colnames(diffIC_diffFC_res)

legend2 <- get_legend(ggplot(no_competition_L2_res,aes(x=time, y=vWA_bound, group=Experiment)) + 
             geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
             scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
             scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
             labs(x = "Time (s)",
                  y = expression(paste("L2-bound Integrin ( nM )",sep = "")))+
             theme_bw(base_size = 13)+
             theme(legend.position = "top", 
                   legend.title = element_blank(), 
                   legend.text = element_text(size = 12), 
                   legend.key.size = unit(3, units = "cm"),
                   axis.text.x = element_text(angle = 30),
                   axis.text = element_text(size = 13),
                   axis.title = element_text(size = 13)) +
             scale_y_continuous( limits=c(0,0.022))) 

no_competition_L2_bound_day18 <- ggplot(no_competition_L2_res[no_competition_L2_res$Experiment == "day18",],aes(x=time, y=vWA_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Integrin ( nM )",sep = "")), title = "No ligand competition", subtitle = "L2 initial concentration: 0.33 nM")+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.key.size = unit(3, units = "cm"),
        axis.text.x = element_text(angle = 30),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) +
  scale_y_continuous( limits=c(0,0.022))

no_competition_L2_bound_day25 <- ggplot(no_competition_L2_res,aes(x=time, y=vWA_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Integrin ( nM )",sep = "")), title = "No ligand competition", subtitle = "L2 initial concentration: 0.50 nM")+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.key.size = unit(3, units = "cm"),
        axis.text.x = element_text(angle = 30),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) +
  scale_y_continuous( limits=c(0,0.022))

no_competition_L1_bound_day18 <- ggplot(no_competition_L1_res[no_competition_L1_res$Experiment == "day18",],aes(x=time, y=F_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-bound Integrin ( nM )",sep = "")), title = "No ligand competition", subtitle = "L1 initial concentration: 0.18 nM")+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.key.size = unit(3, units = "cm"),
        axis.text.x = element_text(angle = 30),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) +
  scale_y_continuous( limits=c(0,0.022))

no_competition_L1_bound_day25 <- ggplot(no_competition_L1_res,aes(x=time, y=F_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-bound Integrin ( nM )",sep = "")), title = "No ligand competition", subtitle = "L1 initial concentration: 0.46 nM")+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.key.size = unit(3, units = "cm"),
        axis.text.x = element_text(angle = 30),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) +
  scale_y_continuous( limits=c(0,0.022))

diff_IC_diffFC_L1_bound_day18 <- ggplot(diffIC_diffFC_res[diffIC_diffFC_res$Experiment == "day18",],aes(x=time, y=F_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-bound Integrin ( nM )",sep = "")), title = "Ligand competition", subtitle = "L1 initial concentration: 0.18 nM")+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.key.size = unit(3, units = "cm"),
        axis.text.x = element_text(angle = 30),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) +
  scale_y_continuous( limits=c(0,0.022))

diff_IC_diffFC_L1_bound_day25 <- ggplot(diffIC_diffFC_res,aes(x=time, y=F_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-bound Integrin ( nM )",sep = "")), title = "Ligand competition", subtitle = "L1 initial concentration: 0.46 nM")+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.key.size = unit(3, units = "cm"),
        axis.text.x = element_text(angle = 30),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) +
  scale_y_continuous( limits=c(0,0.022))

diff_IC_diffFC_L2_bound_day18 <- ggplot(diffIC_diffFC_res[diffIC_diffFC_res$Experiment == "day18",],aes(x=time, y=vWA_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Integrin ( nM )",sep = "")), title = "Ligand competition", subtitle = "L2 initial concentration: 0.33 nM")+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.key.size = unit(3, units = "cm"),
        axis.text.x = element_text(angle = 30),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) #+
  #scale_y_continuous( limits=c(0,0.022))

diff_IC_diffFC_L2_bound_day25 <- ggplot(diffIC_diffFC_res,aes(x=time, y=vWA_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Integrin ( nM )",sep = "")), title = "Ligand competition", subtitle = "L2 initial concentration: 0.50 nM")+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.key.size = unit(3, units = "cm"),
        axis.text.x = element_text(angle = 30),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) #+
  #scale_y_continuous( limits=c(0,0.022))

equalBR_L1_bound_day18 <- ggplot(diffIC_equalBR_res[diffIC_equalBR_res$Experiment == "day18",],aes(x=time, y=F_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-bound Integrin ( nM )",sep = "")), title = "Ligand competition with equal binding rates", subtitle = "L1 initial concentration: 0.18 nM")+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.key.size = unit(3, units = "cm"),
        axis.text.x = element_text(angle = 30),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) +
  scale_y_continuous( limits=c(0,0.022))

equalBR_L1_bound_day25 <- ggplot(diffIC_equalBR_res,aes(x=time, y=F_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-bound Integrin ( nM )",sep = "")),title = "Ligand competition with equal binding rates",  subtitle = "L1 initial concentration: 0.46 nM")+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.key.size = unit(3, units = "cm"),
        axis.text.x = element_text(angle = 30),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) +
  scale_y_continuous( limits=c(0,0.022))

equalBR_L2_bound_day18 <- ggplot(diffIC_equalBR_res[diffIC_equalBR_res$Experiment == "day18",],aes(x=time, y=vWA_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Integrin ( nM )",sep = "")), title = "Ligand competition with equal binding rates", subtitle = "L2 initial concentration: 0.33 nM")+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.key.size = unit(3, units = "cm"),
        axis.text.x = element_text(angle = 30),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) +
  scale_y_continuous( limits=c(0,0.022))

equalBR_L2_bound_day25 <- ggplot(diffIC_equalBR_res,aes(x=time, y=vWA_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Integrin ( nM )",sep = "")), title = "Ligand competition with equal binding rates",  subtitle = "L2 initial concentration: 0.50 nM")+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.key.size = unit(3, units = "cm"),
        axis.text.x = element_text(angle = 30),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) +
  scale_y_continuous( limits=c(0,0.022))
?ggsave

ggsave("no_comp_L1_bound_day18.pdf",
       plot = no_competition_L1_bound_day18,
       width = 10,
       height = 10,
       units = "cm")

ggsave("no_comp_L1_bound_day25.pdf",
       plot = no_competition_L1_bound_day25,
       width = 10,
       height = 10,
       units = "cm")

ggsave("no_comp_L2_bound_day18.pdf",
       plot = no_competition_L2_bound_day18,
       width = 10,
       height = 10,
       units = "cm")

ggsave("no_comp_L2_bound_day25.pdf",
       plot = no_competition_L2_bound_day25,
       width = 10,
       height = 10,
       units = "cm")

ggsave("diffIC_diff_FC_L1_bound_day18.pdf",
       plot = diff_IC_diffFC_L1_bound_day18,
       width = 10,
       height = 10,
       units = "cm")

ggsave("diffIC_diff_FC_L1_bound_day25.pdf",
       plot = diff_IC_diffFC_L1_bound_day25,
       width = 10,
       height = 10,
       units = "cm")

ggsave("diffIC_diff_FC_L2_bound_day18.pdf",
       plot = diff_IC_diffFC_L2_bound_day18,
       width = 10,
       height = 10,
       units = "cm")

ggsave("diffIC_diff_FC_L2_bound_day25.pdf",
       plot = diff_IC_diffFC_L2_bound_day25,
       width = 10,
       height = 10,
       units = "cm")

ggsave("equal_BR_L1_bound_day18.pdf",
       plot = equalBR_L1_bound_day18,
       width = 15,
       height = 15,
       units = "cm")

ggsave("equal_BR_L1_bound_day25.pdf",
       plot = equalBR_L1_bound_day25,
       width = 15,
       height = 15,
       units = "cm")

ggsave("equal_BR_L2_bound_day18.pdf",
       plot = equalBR_L2_bound_day18,
       width = 15,
       height = 15,
       units = "cm")

ggsave("equal_BR_L2_bound_day25.pdf",
       plot = equalBR_L2_bound_day25,
       width = 15,
       height = 15,
       units = "cm")
