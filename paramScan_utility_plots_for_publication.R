setwd("C:/karagoz/01-RESEARCH/01-Projects/01-In_silico_modeling_of_Integrin_function/003-Ligand_competition_model/Unique_clusters/figures_FN_vWF_equal_IC/uncertainty/")

#install + load packages

pkg<-c("magrittr", "ggplot2", "dplyr", "cowplot", "tidyverse")
install.packages(pkg)
lapply(pkg, require, character.only = TRUE)

# read uncertainty results all in a list:

myfiles = list.files(pattern="*.csv", full.names=TRUE)
data_list <- lapply(myfiles, read.csv)
names(data_list) <- myfiles #name the individual dataframes in the list with appropriate filenames


ordered_file_names <- c("./k1_i.csv","./k1_a.csv", "./k1_IF.csv", "./k1_IW.csv","./k1_C1.csv", "./k1_C2.csv", "./k1_C3.csv",
                        "./k2_i.csv","./k2_a.csv", "./k2_IF.csv", "./k2_IW.csv","./k2_C1.csv","./k2_C2.csv", "./k2_C3.csv",
                        "./k3_i.csv","./k3_a.csv", "./k3_IF.csv", "./k3_IW.csv","./k3_C1.csv", "./k3_C2.csv", "./k3_C3.csv",
                        "./k4_i.csv","./k4_a.csv", "./k4_IF.csv","./k4_IW.csv", "./k4_C1.csv", "./k4_C2.csv", "./k4_C3.csv",
                        "./k5_i.csv","./k5_a.csv",  "./k5_IF.csv", "./k5_IW.csv", "./k5_C1.csv", "./k5_C2.csv", "./k5_C3.csv",
                        "./k6_i.csv","./k6_a.csv",  "./k6_IF.csv", "./k6_IW.csv","./k6_C1.csv", "./k6_C2.csv", "./k6_C3.csv",
                        "./k7_i.csv", "./k7_a.csv",  "./k7_IF.csv", "./k7_IW.csv", "./k7_C1.csv", "./k7_C2.csv", "./k7_C3.csv",
                        "./k8_i.csv", "./k8_a.csv", "./k8_IF.csv", "./k8_IW.csv", "./k8_C1.csv", "./k8_C2.csv", "./k8_C3.csv")
data_list <- data_list[ordered_file_names]


# plot 
# using ggplot ribbon plot

plots_list_k1 <- list()
for (n in ordered_file_names[1:7]){
  plots_list_k1[[n]] <- ggplot(data_list[[n]], aes(x=time)) +  
  geom_ribbon(aes(ymin=lower.limit, ymax=upper.limit), fill="pink",alpha=0.5)+
  geom_line(aes(y=mean), size=1) +
  labs(x = "", y = "")+
  theme_bw(base_size = 15)+
  theme(axis.text = element_text(size =13, face = "bold"))
}
plots_list_k1[["./k1_C1.csv"]] #works

plots_list_k2 <- list()
for (n in ordered_file_names[8:14]){
  plots_list_k2[[n]] <- ggplot(data_list[[n]], aes(x=time)) +  
    geom_ribbon(aes(ymin=lower.limit, ymax=upper.limit), fill="pink1",alpha=0.5)+
    geom_line(aes(y=mean), size=1) +
    labs(x = "", y = "")+
    theme_bw(base_size = 15)+
    theme(axis.text = element_text(size =13, face = "bold"))
}
plots_list_k2[["./k2_C1.csv"]] #works

plots_list_k3 <- list()
for (n in ordered_file_names[15:21]){
  plots_list_k3[[n]] <- ggplot(data_list[[n]], aes(x=time)) +  
    geom_ribbon(aes(ymin=lower.limit, ymax=upper.limit), fill="pink2",alpha=0.5)+
    geom_line(aes(y=mean), size=1) +
    labs(x = "", y = "")+
    theme_bw(base_size = 15)+
    theme(axis.text = element_text(size =13, face = "bold"))
}
plots_list_k3[["./k3_IW.csv"]]

plots_list_k4 <- list()
for (n in ordered_file_names[22:28]){
  plots_list_k4[[n]] <- ggplot(data_list[[n]], aes(x=time)) +  
    geom_ribbon(aes(ymin=lower.limit, ymax=upper.limit), fill="pink3",alpha=0.5)+
    geom_line(aes(y=mean), size=1) +
    labs(x = "", y = "")+
    theme_bw(base_size = 15)+
    theme(axis.text = element_text(size =13, face = "bold"))
}
plots_list_k4[["./k4_a.csv"]]

plots_list_k5 <- list()
for (n in ordered_file_names[29:35]){
  plots_list_k5[[n]] <- ggplot(data_list[[n]], aes(x=time)) +  
    geom_ribbon(aes(ymin=lower.limit, ymax=upper.limit), fill="pink3",alpha=0.5)+
    geom_line(aes(y=mean), size=1) +
    labs(x = "", y = "")+
    theme_bw(base_size = 15)+
    theme(axis.text = element_text(size =13, face = "bold"))
}
plots_list_k5[["./k5_IW.csv"]]

plots_list_k6 <- list()
for (n in ordered_file_names[36:42]){
  plots_list_k6[[n]] <- ggplot(data_list[[n]], aes(x=time)) +  
    geom_ribbon(aes(ymin=lower.limit, ymax=upper.limit), fill="plum",alpha=0.5)+
    geom_line(aes(y=mean), size=1) +
    labs(x = "", y = "")+
    theme_bw(base_size = 15)+
    theme(axis.text = element_text(size =13, face = "bold"))
}
plots_list_k6[["./k6_IW.csv"]]


plots_list_k7 <- list()
for (n in ordered_file_names[43:49]){
  plots_list_k7[[n]] <- ggplot(data_list[[n]], aes(x=time)) +  
    geom_ribbon(aes(ymin=lower.limit, ymax=upper.limit), fill="plum",alpha=0.5)+
    geom_line(aes(y=mean), size=1) +
    labs(x = "", y = "")+
    theme_bw(base_size = 15)+
    theme(axis.text = element_text(size =13, face = "bold"))
}
plots_list_k7[["./k7_C2.csv"]]

plots_list_k8 <- list()
for (n in ordered_file_names[50:56]){
  plots_list_k8[[n]] <- ggplot(data_list[[n]], aes(x=time)) +  
    geom_ribbon(aes(ymin=lower.limit, ymax=upper.limit), fill="plum1", alpha=0.5)+
    geom_line(aes(y=mean), size=1) +
    labs(x = "", y = "")+
    theme_bw(base_size = 15)+
    theme(axis.text = element_text(size =13, face = "bold"))
}
plots_list_k8[["./k8_C3.csv"]]


# cowplot all 
# rows: different parameters k1-k8
# cols: each model variable 

row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="k1", angle = 90, size=10) + theme_void() 

row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="k2", angle = 90, size=10) + theme_void() 

row3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="k3", angle = 90, size=10) + theme_void() 

row4 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="k4", angle = 90, size=10) + theme_void() 

row5 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="k5", angle = 90, size=10) + theme_void() 

row6 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="k6", angle = 90, size=10) + theme_void() 

row7 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="k7", angle = 90, size=10) + theme_void() 

row8 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="k8", angle = 90, size=10) + theme_void() 

row_list1 <- list(row1,row2, row3,row4)
row_list2 <- list(row5,row6,row7,row8)

col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="inactive integrin", size=10) + theme_void() 

col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="active integrin", size=10) + theme_void()

col3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="L1-bound integrin", size=10) + theme_void() 

col4 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="L2-bound integrin", size=10) + theme_void()

col5 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="L1-integrin cluster", size=10) + theme_void() 

col6 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="L2-integrin cluster", size=10) + theme_void()

col7 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="mixed integrin cluster", size=10) + theme_void()

col_list <- list(col1, col2, col3, col4, col5, col6, col7)


plot_grid(plot_grid(plotlist = row_list1, nrow=4), 
          plot_grid(plot_grid(plotlist = col_list, nrow = 1),
          plot_grid(plotlist = plots_list_k1, nrow = 1),
          plot_grid(plotlist = plots_list_k2, nrow = 1),
          plot_grid(plotlist = plots_list_k3, nrow = 1),
          plot_grid(plotlist = plots_list_k4, nrow = 1), ncol = 1, nrow=4, rel_heights = c(0.2,1,1,1,1)), ncol = 2, rel_widths = c(0.01,1))


plot_grid(plot_grid(plotlist = row_list2, nrow=4),
          plot_grid(plot_grid(plotlist = plots_list_k5, nrow= 1),
          plot_grid(plotlist = plots_list_k6, nrow=1),
          plot_grid(plotlist = plots_list_k7,nrow = 1),
          plot_grid(plotlist = plots_list_k8, nrow = 1), ncol = 1, nrow = 4), ncol = 2, rel_widths = c(0.01,1))
# save both as pdf: 20x38