# plotting simulation results 
# ligand competition model with unique clusters
# 2 set of results: 1) different initial condition for ligands 
#                   2) equal initial conditions for ligands. 

setwd("C:/karagoz/01-RESEARCH/01-Projects/01-In_silico_modeling_of_Integrin_function/003-Ligand_competition_model/Unique_clusters/pretty_plots")

#install + load packages
  
pkg<-c("magrittr", "ggplot2", "dplyr", "cowplot")
#install.packages(pkg)
lapply(pkg, require, character.only = TRUE)


# read files
equalIC_res <- read.csv("../figures_FN_vWF_equal_IC/equalIC_simResults.csv", header = TRUE, sep="\t")
diffIC_res <- read.csv("../figures_FN_vWF/diffIC_simResults.csv", header = TRUE, sep = "\t")

THBS_equalIC_res <- read.csv("../figures_FN_THBS_equal_IC/F_THBS_equalIC_simResults.csv", header = TRUE, sep="\t")
# plot


p <- ggplot(diffIC_res,aes(x=time, y=inactive, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("Inactive Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 15)+
  theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.text = element_text(size= 14, face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size =13, face = "bold"))
legend <- get_legend(p)

diff_inactive <- ggplot(diffIC_res,aes(x=time, y=inactive, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("Inactive Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size= 14, face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size =13, face = "bold"))

        
diff_active <- ggplot(diffIC_res,aes(x=time, y=active, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0,0.002))+
  labs(x = "Time (s)",
       y = expression(paste("Active Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size= 14, face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size= 13, face = "bold"))

diff_Fbound <- ggplot(diffIC_res,aes(x=time, y=F_bound, group=Experiment)) + 
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
  scale_y_continuous( limits=c(0,0.022))

diff_Wbound <- ggplot(diffIC_res,aes(x=time, y=vWA_bound, group=Experiment)) + 
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
        axis.title = element_text(size = 13, face = "bold")) 

diff_F_F_cluster <- ggplot(diffIC_res,aes(x=time, y=IF_IFclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-bound Cluster ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) 

diff_W_W_cluster <- ggplot(diffIC_res,aes(x=time, y=IW_IWclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0,4e-9))+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Cluster ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) 


diff_F_W_cluster <- ggplot(diffIC_res,aes(x=time, y=IF_IWclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0,8e-6))+
  labs(x = "Time (s)",
       y = expression(paste("L1-L2-bound Cluster ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) 


plot_grid(diff_inactive, diff_active, legend, diff_Fbound, diff_Wbound, NULL, diff_F_F_cluster, diff_W_W_cluster, diff_F_W_cluster,
  ncol=3, labels = c("A", "B", "","C", "D", "", "E", "F", "G"))
# export as pdf, 14x18


# equal initial condition for ligands 

eq_inactive <- ggplot(equalIC_res,aes(x=time, y=inactive, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("Inactive Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size= 14, face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size =13, face = "bold"))


eq_active <- ggplot(equalIC_res,aes(x=time, y=active, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0,0.002))+
  labs(x = "Time (s)",
       y = expression(paste("Active Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size= 14, face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size= 13, face = "bold"))

eq_Fbound <- ggplot(equalIC_res,aes(x=time, y=F_bound, group=Experiment)) + 
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
  scale_y_continuous( limits=c(0,0.022))

eq_Wbound <- ggplot(equalIC_res,aes(x=time, y=vWA_bound, group=Experiment)) + 
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
        axis.title = element_text(size = 13, face = "bold")) 

eq_F_F_cluster <- ggplot(equalIC_res,aes(x=time, y=IF_IFclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-bound Cluster ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) 

eq_W_W_cluster <- ggplot(equalIC_res,aes(x=time, y=IW_IWclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0,4e-9))+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Cluster ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) 


eq_F_W_cluster <- ggplot(equalIC_res,aes(x=time, y=IF_IWclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0,8e-6))+
  labs(x = "Time (s)",
       y = expression(paste("L1-L2-bound Cluster ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) 


plot_grid(eq_inactive, eq_active, legend, eq_Fbound, eq_Wbound, NULL, eq_F_F_cluster, eq_W_W_cluster, eq_F_W_cluster,
          ncol=3, labels = c("A", "B", "","C", "D", "", "E", "F", "G"))
#export as pdf 14x18



# plot equal amount of fibronectin and THBS 
# but THBS increases 9fold on day 25 while fibronectin increases only 2.5fold

THBSeq_inactive <- ggplot(THBS_equalIC_res,aes(x=time, y=inactive, group=Experiment)) + 
  geom_line(aes(color= Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("Inactive Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size= 14, face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size =13, face = "bold"))


THBSeq_active <- ggplot(THBS_equalIC_res,aes(x=time, y=active, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0,0.002))+
  labs(x = "Time (s)",
       y = expression(paste("Active Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size= 14, face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size= 13, face = "bold"))

THBSeq_Fbound <- ggplot(THBS_equalIC_res,aes(x=time, y=F_bound, group=Experiment)) + 
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
  scale_y_continuous( limits=c(0,0.022))

THBSeq_Tbound <- ggplot(THBS_equalIC_res,aes(x=time, y=T_bound, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0,1.5e-5))+
  labs(x = "Time (s)",
       y = expression(paste("L3-bound Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) 

THBSeq_F_F_cluster <- ggplot(THBS_equalIC_res,aes(x=time, y=IF_IFclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("L1-bound Cluster ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) 

THBSeq_T_T_cluster <- ggplot(THBS_equalIC_res,aes(x=time, y=IT_ITclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0,4e-9))+
  labs(x = "Time (s)",
       y = expression(paste("L3-bound Cluster ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) 


THBSeq_F_T_cluster <- ggplot(THBS_equalIC_res,aes(x=time, y=IF_ITclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  scale_y_continuous(limits = c(0,8e-6))+
  labs(x = "Time (s)",
       y = expression(paste("L1-L3-bound Cluster ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) 


plot_grid(THBSeq_inactive, THBSeq_active, legend, THBSeq_Fbound, THBSeq_Tbound, NULL, THBSeq_F_F_cluster, THBSeq_T_T_cluster, THBSeq_F_T_cluster,
          ncol=3, labels = c("A", "B", "","C", "D", "", "E", "F", "G"))

#change layout! 
# plot species independently, but for all conditions 
# A = FN & vWA different initial conditions 
# B = FN & vWA equal initial conditions 
# C = FN & vWA with 9 fold increase in day 25, equal initial conditions

# inactive integrin concentration 
plot_grid(diff_inactive, eq_inactive, THBSeq_inactive, legend, ncol=4, labels=c("A", "B", "C", ""), rel_widths = c(1,1,1,0.3))
# save pdf 4.8 x 20

# active integrin concentration

plot_grid(diff_active, eq_active, THBSeq_active, legend, ncol=4, labels=c("A", "B", "C", ""), rel_widths = c(1,1,1,0.3))

# Fibronectin-bound integrin

plot_grid(diff_Fbound, eq_Fbound, THBSeq_Fbound,legend, ncol=4, labels=c("A", "B", "C", ""), rel_widths = c(1,1,1,0.3) )

#vWA-bound integrin

plot_grid(diff_Wbound, eq_Wbound, THBSeq_Tbound, legend, ncol=4, labels=c("A", "B", "C", ""), rel_widths = c(1,1,1,0.3) )


# F clusters
plot_grid(diff_F_F_cluster, eq_F_F_cluster, THBSeq_F_F_cluster, legend, ncol=4, labels=c("A", "B", "C", ""), rel_widths = c(1,1,1,0.3))

#vWA clusters
plot_grid(diff_W_W_cluster, eq_W_W_cluster, THBSeq_T_T_cluster, legend, ncol=4, labels=c("A", "B", "C", ""), rel_widths = c(1,1,1,0.3))

# mixed cluster
plot_grid(diff_F_W_cluster, eq_F_W_cluster, THBSeq_F_T_cluster, legend, ncol=4, labels=c("A", "B", "C", ""), rel_widths = c(1,1,1,0.3))


#### paper v1 layout

# fig1: ligand bound integrins with experimental initial conditions
# A: L1-bound integrin
# B: L2-bound integrin

plot_grid(diff_Fbound, diff_Wbound,legend, ncol=3, labels = c("A", "B", ""), rel_widths = c(1,1,0.3))
# save pdf 4.8 x 14

# fig2: ligand bound integrins with equal initial conditions 
# A: L1-bound I (2.5-fold change)
# B: L2-bound I (1.5-fold change)
# C: L1-bound I
# D: L3-bound I (9-fold change) 

plot_grid(eq_Fbound, eq_Wbound, legend, THBSeq_Fbound, THBSeq_Tbound, ncol=3, nrow = 2, labels = c("A", "B", "", "C", "D"), rel_widths = c(1,1,0.3,1,1) )
# save pdf 9.6 x 14

# fig3: integrin cluster composition reflects the ligand competition
# equal IC

# A: L1-bound clusters
# B: L2-bound clusters
# C: mixed L1-L2 clusters
# D: L1-bound clusters
# E: L3-bound clusters
# F: mixed L1-L3 bound clusters

plot_grid(eq_F_F_cluster, eq_W_W_cluster, eq_F_W_cluster, legend, THBSeq_F_F_cluster, THBSeq_T_T_cluster, THBSeq_F_T_cluster, ncol=4, nrow=2, labels = c("A", "B", "C", "", "D", "E", "F"), rel_widths = c(1,1,1,0.3,1,1,1))
#save pdf 9.6 x 20


#figure 4: param scan of L2
# read files
install.packages("colorspace")
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
