setwd("C:/karagoz/01-RESEARCH/01-Projects/01-In_silico_modeling_of_Integrin_function/003-Ligand_competition_model/Unique_clusters/figures_FN_vWF_equal_IC/parameterScan_k_values/")

#install + load packages

pkg<-c("magrittr", "ggplot2", "dplyr", "cowplot", "tidyverse")
#install.packages(pkg)
lapply(pkg, require, character.only = TRUE)

# read files
equalIC_k5_scan <- read.csv("k5_paramScan_IF_IW_results.csv", header = TRUE, sep="\t")



#install.packages("colorspace")
library(colorspace)
equalIC_k5_scan$grad = as.character(equalIC_k5_scan$grad)
ggplot(equalIC_k5_scan,aes(x=time, y=vWA_bound, group=grad)) + 
  geom_line(aes(color=grad),size=1.5) +
  scale_color_discrete_diverging(palette="Green-Orange")+
  labs(x = "Time (s)",
       y = expression(paste("L2-bound Integrin ( ", mu,"M)",sep = "")),
       color= "k5 =")+
  theme_bw(base_size = 13)+
  theme(legend.position = "right", 
        legend.title = element_text(size=15), 
        legend.text = element_text(size = 14,face="bold"), 
        legend.key.size = unit(1.2, units = "cm"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")) 
# save pdf 6x9

# plot k5 vs L2-bound integrin

DF <- dplyr::group_by(equalIC_k5_scan, grad) %>% 
  dplyr::summarise(steady_state=max(vWA_bound)) 

ggplot(DF,aes(x=as.numeric(grad), y=steady_state)) + 
  geom_line(size=1.5) +
  #scale_color_discrete_diverging(palette="Green-Orange")+
  labs(x = "k5",y = expression(paste("L2-bound Integrin ( ", mu,"M)",sep = "")))+
  theme_bw(base_size = 13)#+
#  theme(legend.position = "right", 
#        legend.title = element_text(size=15), 
 #       legend.text = element_text(size = 14,face="bold"), 
  #      legend.key.size = unit(1.2, units = "cm"),
   #     axis.text = element_text(size = 13, face = "bold"),
    #    axis.title = element_text(size = 13, face = "bold")) 