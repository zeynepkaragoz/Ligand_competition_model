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
       y = expression(paste("Fibronectin-bound Integrin ( ", mu,"M)",sep = "")))+
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
  labs(x = "Time (s)",
       y = expression(paste("vWF-bound Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) +
  scale_y_continuous( limits=c(0,6e-6))

diff_F_F_cluster <- ggplot(diffIC_res,aes(x=time, y=IF_IFclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("F-bound Cluster ( ", mu,"M)",sep = "")))+
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
  labs(x = "Time (s)",
       y = expression(paste("vWF-bound Cluster ( ", mu,"M)",sep = "")))+
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
  labs(x = "Time (s)",
       y = expression(paste("F-vWF-bound Cluster ( ", mu,"M)",sep = "")))+
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
       y = expression(paste("Fibronectin-bound Integrin ( ", mu,"M)",sep = "")))+
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
  labs(x = "Time (s)",
       y = expression(paste("vWF-bound Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) +
  scale_y_continuous( limits=c(0,6e-6))

eq_F_F_cluster <- ggplot(equalIC_res,aes(x=time, y=IF_IFclustered, group=Experiment)) + 
  geom_line(aes(color=Experiment, linetype=Experiment),size=1.5) +
  scale_linetype_manual(labels=c("Day 18", "Day 25"), values=c("solid", "dotted"))+
  scale_color_manual(labels=c("Day 18", "Day 25"), values = c("gray" , "#FC4E07"))+
  labs(x = "Time (s)",
       y = expression(paste("F-bound Cluster ( ", mu,"M)",sep = "")))+
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
  labs(x = "Time (s)",
       y = expression(paste("vWF-bound Cluster ( ", mu,"M)",sep = "")))+
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
  labs(x = "Time (s)",
       y = expression(paste("F-vWF-bound Cluster ( ", mu,"M)",sep = "")))+
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
       y = expression(paste("Fibronectin-bound Integrin ( ", mu,"M)",sep = "")))+
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
  labs(x = "Time (s)",
       y = expression(paste("THBS-bound Integrin ( ", mu,"M)",sep = "")))+
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
       y = expression(paste("F-bound Cluster ( ", mu,"M)",sep = "")))+
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
  labs(x = "Time (s)",
       y = expression(paste("THBS-bound Cluster ( ", mu,"M)",sep = "")))+
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
  labs(x = "Time (s)",
       y = expression(paste("F-THBS-bound Cluster ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)+
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) 


plot_grid(THBSeq_inactive, THBSeq_active, legend, THBSeq_Fbound, THBSeq_Tbound, NULL, THBSeq_F_F_cluster, THBSeq_T_T_cluster, THBSeq_F_T_cluster,
          ncol=3, labels = c("A", "B", "","C", "D", "", "E", "F", "G"))
