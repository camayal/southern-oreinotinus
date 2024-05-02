### This script is derived from Nora's script for env analyses in Oreinotinus paper
### adapted for Bolivian history

#clean previous run
rm(list = ls())


# # install packages
# install.packages("tidyverse")
# install.packages("ggplot2")
# install.packages("factoextra")
# install.packages("FactoMineR")
# install.packages("ggsci")
# install.packages("ggpubr")
# install.packages("corrplot")
# install.packages("caret")
# install.packages("ggforce")
# install.packages("RColorBrewer")
# install.packages("ggrepel)" 
# install.packages("viridis")
# install.packages("ggridges")
# install.packages("extrafont")
# install.packages("rstatix")
# install.packages("svglite")


# load packages
library(tidyverse)
library(ggplot2)
library(factoextra)
library(FactoMineR)
library(ggsci)
library(ggpubr)
library(corrplot)
library(caret)
library(ggforce)
library(RColorBrewer)
library(ggrepel) 
library(viridis)
library(ggridges)
library(extrafont)
library(rstatix)
library(cowplot)
library(svglite)
loadfonts(device = "win")


#### set-up ####
# set workdir
setwd("C:\\Users\\camay\\Dropbox\\Postdocs_research\\Edwards-lab\\bolivia\\env_analyses\\")

# bring in data
data <- read.csv("Different_classifications\\fulldist_categoriesByCAML_withChelsey.csv", header = T, sep = ",", na.strings = "")



## SET CLASSIFICATION COLUMN TO TEST
# CLASSIFICATION = "clas1_cladeBased"
CLASSIFICATION = "clas2_cladeBasedExtended"
# CLASSIFICATION = "clas3_spFromTree"
# CLASSIFICATION = "clas4_structureK3"
# CLASSIFICATION = "clas5_structureK4"

  # select only columns of interest (environmental variables, elevation, and groups)
  data <- data[, c(CLASSIFICATION,paste0(colnames(data)[18:36]),"elevation")]


  # put all variables in correct format
  data[,1] <- as.factor(data[,1])

  number_columns = ncol(data)
  for (i in 2:number_columns) {
    data[,i] <- as.numeric(data[,i])
  }

  # remove rows with NA values
  # notice that here a lot of data is being removed because no all points are selected in the classification
  data <- na.omit(data)
  any(is.na(data))

  # standardize variable names for easier plotting
  vec <- c("group", "bio01", "bio02", "bio03", "bio04", "bio05", "bio06", 
          "bio07", "bio08", "bio09", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", 
          "bio16", "bio17", "bio18", "bio19", "eleva")
  colnames(data) <- vec


  # not sure if this is doing something in my matrices
  data[,1] <- droplevels(data[,1])


  # make a theme
  theme_clim <- function(){
    theme_bw() +
      theme(text = element_text(family = "Helvetica"),
            axis.text = element_text(size = 20), 
        #    axis.text.x = element_text(face = "italic", vjust = 0.9, hjust = 1, angle = 30),
            axis.title = element_text(size = 25),
    #        axis.title.y = element_text(vjust = 0.4),
            #           axis.line.x = element_line(color = "black"), 
            #           axis.line.y = element_line(color = "black"),
            #           panel.border = element_line(color = "black"),
    #        axis.title.x = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.major.x = element_blank(),                                          
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),  
            panel.border = element_blank(),
            plot.margin = unit(c(1, 1, 1, 1), units = , "cm"))
            #           plot.title = element_text(size = 18, vjust = 1, hjust = 0),
            #           legend.text = element_text(size = 12),          
            #           legend.title = element_blank(),                              
            #           legend.position = c(0.95, 0.15), 
            #           legend.key = element_blank(),
            #           legend.background = element_rect(color = "black", 
            #                                            fill = "transparent", 
            #                                            size = 2, linetype = "blank"),
            # strip.text = element_text(size = 10, color = "black", face = "bold.italic"),
            # strip.background = element_rect(color = "white", fill = "white", size = 1))
  }

  # load Flat Violin geom
  source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

  # run a quick loop to determine which variables best differentiat hart, juc, and laut
  # bio1, bio2, bio4, bio7, bio12, bio20

  # make a vector of titles for the graph
  col_names <- c("Group",
                "01 - Mean Annual Temp. (°C)",
                "02 - Mean Diurnal Air Range (°C)",
                "03 - Isothermality (°C)",
                "04 - Temp. Seasonality\n(standard deviation °C/100)",
                "05 - Max Temp. of Warmest Month (°C)",
                "06 - Min Temp. of Coldest Month (°C)",
                "07 - Temp. Annual Range (°C)",
                "08 - Mean Temp. of Wettest Quarter (°C",
                "09 - Mean Temp. of Driest Quarter (°C)",
                "10 - Mean Temp. of Warmest Quarter (°C)",
                "11 - Mean Temp. of Coldest Quarter (°C)",
                "12 - Annual Precip. (mm/year)",
                "13 - Precip. of Wettest Month\n(mm/month)",
                "14 - Precip. of Driest Month\n(mm/month)",
                "15 - Precip. Seasonality\n(coefficient of variation)",
                "16 - Precip. of Wettest Quarter\n(mm/quarter)",
                "17 - Precip. of Driest Quarter\n(mm/quarter)",
                "18 - Precip. of Warmest Quarter\n(mm/quarter)",
                "19 - Precip. of Coldest Quarter\n(mm/quarter)",
                "Elevation (m)"
                )
  
  
  col_names <- c("Group",
                 "Chelsa01\nMean Annual Temp.(°C)",
                 "Chelsa02\nMean Diurnal Air (°C)",
                 "Chelsa03\nIsothermality (°C)",
                 "Chelsa04\nTemp. Seasonality (°C/100)",
                 "Chelsa05\nMax Temp. of Warmest Mth. (°C)",
                 "Chelsa06\nMin Temp. of Coldest Mth. (°C)",
                 "Chelsa07\nTemp. Annual Range (°C)",
                 "Chelsa08\nMean Temp. of Wettest Qtr. (°C",
                 "Chelsa09\nMean Temp. of Driest Qtr. (°C)",
                 "Chelsa10\nMean Temp. of Warmest Qtr. (°C)",
                 "Chelsa11\nMean Temp. of Coldest Qtr. (°C)",
                 "Chelsa12\nAnnual Precip. (mm/year)",
                 "Chelsa13\nPrecip. Wettest Mth.(mm/mth)",
                 "Chelsa14\nPrecip. Driest Mth. (mm/mth)",
                 "Chelsa15\nPrecip. Seasonality (coeff. var.)",
                 "Chelsa16\nPrecip. Wettest Qtr. (mm/qtr)",
                 "Chelsa17\nPrecip. Driest Qtr.(mm/qtr)",
                 "Chelsa18\nPrecip. Warmest Qtr.(mm/qtr)",
                 "Chelsa19\nPrecip. Coldest Qtr.(mm/qtr)",
                 "Elevation (m)"
  )


  # make x-axis labels
  group_labels <- unique(data[["group"]])



# create a plot for each variable including elevation (ONE FIGURE WITH ALL VARIABLES)

#define a function to plot with i as only parameter
plotting <- function(i){
     ggplot(data = data, aes(x = group, y = data[,i], fill = group)) +
     geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.8) +
     geom_point(aes(y = data[,i], color = group), 
                position = position_jitter(width = 0.15), size = 3, alpha = 0.3) +
     geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
     labs(y = paste0(col_names[i], "\n"), x = "\nGroups") +
     guides(fill = "none", color = "none") +
     scale_fill_manual(values = c("#66C2A4","#FB8D61","#8D9FCA","#E789C3","#A6D753","#FFD92E","#E4C493")) +
     scale_colour_manual(values = c("#66C2A4","#FB8D61","#8D9FCA","#E789C3","#A6D753","#FFD92E","#E4C493")) +
    #  scale_x_discrete(labels = group_labels) +
    # stat_compare_means(comparison=my_comparisons, label="p.format", method="wilcox.test", p.adjust.method = "bonferroni") +
    #  stat_compare_means(comparisons = list(c("hartwegii", "jucundum"), 
    #                                        c("hartwegii", "lautum"), 
    #                                        c("jucundum", "lautum")), size = 7, label = "p.signif") +
     theme_clim()



}



#use map to pass i from 2 to ncolumns into the function above created
plot_list <- map(c(2:number_columns),plotting)

#merge all plots in one page
unified_plot = ggarrange(plotlist = plot_list)


# create folders if do not exist
dir.create(file.path(getwd(), "Results"))
dir.create(file.path(getwd(), "Results", CLASSIFICATION))

#save as png or svg
output_file = paste0("Results\\",CLASSIFICATION,"\\",CLASSIFICATION,".svg")
ggsave(plot = unified_plot, 
      filename = output_file,
      width = 90,
      height = 70,
      units = "cm",
      )

output_file = paste0("Results\\",CLASSIFICATION,"\\",CLASSIFICATION,".png")
ggsave(plot = unified_plot, 
       filename = output_file,
       width = 90,
       height = 70,
       units = "cm",
)


#### plot in one single line

#use map to pass i from 2 to ncolumns into the function above created
plot_list <- map(c(2:number_columns),plotting)


#merge all plots in one page
unified_plot = ggarrange(plotlist = plot_list, nrow = 1)


# create folders if do not exist
dir.create(file.path(getwd(), "Results"))
dir.create(file.path(getwd(), "Results", CLASSIFICATION))

#save as png or svg
output_file = paste0("Results\\",CLASSIFICATION,"\\",CLASSIFICATION,"_line.png")
ggsave(plot = unified_plot, 
      filename = output_file,
      width = 50000,
      height = 2000,
      units = "px",
      limitsize = FALSE
      )




############## TEMP CODE FOR PRESENTATION AND SLIDES PLOTS


#define a function to plot with i as only parameter
plotting <- function(i){
  ggplot(data = data, aes(x=factor(group, level=c('E7', 'E6', 'E5', 'E4', 'E3', 'E2', 'E1')), y = data[,i], fill = group)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.8) +
    geom_point(aes(y = data[,i], color = group), 
               position = position_jitter(width = 0.15), size = 3, alpha = 0.3) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    labs(y = paste0(col_names[i], "\n"), x = "\nClades") +
    guides(fill = "none", color = "none") +
    scale_fill_manual(values = c("#66C2A4","#FB8D61","#8D9FCA","#E789C3","#A6D753","#FFD92E","#E4C493")) +
    scale_colour_manual(values = c("#66C2A4","#FB8D61","#8D9FCA","#E789C3","#A6D753","#FFD92E","#E4C493")) +
    #  scale_x_discrete(labels = group_labels) +
    # stat_compare_means(comparison=my_comparisons, label="p.format", method="wilcox.test", p.adjust.method = "bonferroni") +
     stat_compare_means(comparisons = list(c("E7", "E6"),
                                           c("E7", "E5"),
                                           c("E7", "E4"),
                                           c("E7", "E3"),
                                           c("E7", "E2"),
                                           c("E7", "E1")), size = 7, label = "p.signif") +
    # stat_compare_means(comparisons=comparisons, size = 4, label="p.signif", method="wilcox.test", p.adjust.method = "bonferroni", hide.ns = TRUE) +
    
    # stat_compare_means(ref.group = "E7", size = 4, label="p.signif", method="wilcox.test", p.adjust.method = "bonferroni", hide.ns = TRUE) +
    # stat_compare_means(ref.group = ".all.", size = 4, label="p.signif", method="wilcox.test", p.adjust.method = "bonferroni") +  ## this is not useful because it compare if it is significant repect with the whole mean
    
  
    
    theme_clim()+
    coord_flip()
  
  
  
}

comparisons = list(c("E7", "E6"),
                  c("E7", "E5"),
                  c("E7", "E4"),
                  c("E7", "E3"),
                  c("E7", "E2"),
                  c("E7", "E1")
                   )



#I found that non correlated vars (incuding latitude) are:
# features = ['CHELSA_01', 
#             'CHELSA_02', 
#             # 'CHELSA_03', #removed for high corr. with lat
#             # 'CHELSA_04', #removed for high corr. with lat
#             'CHELSA_05',
#             'CHELSA_06', 
#             # 'CHELSA_07', #removed for high corr. with lat 
#             # 'CHELSA_08', 
#             'CHELSA_09', 
#             'CHELSA_10',
#             # 'CHELSA_11', 
#             'CHELSA_12', 
#             'CHELSA_13', 
#             'CHELSA_14', 
#             # 'CHELSA_15', #removed for high corr. with lat
#             # 'CHELSA_16', 
#             # 'CHELSA_17', 
#             'CHELSA_18', 



#testing how plot looks
bio = 1
plotting(bio+1)
bio = 5
plotting(bio+1)
bio = 9
plotting(bio+1)
bio = 14
plotting(bio+1)


# create a plot for each variable including elevation (ONE FIGURE PER VARIABLE)
for (i in 2:number_columns) {

 plot <- plotting(i)

   #create folders if do not exist
   dir.create(file.path(getwd(), "Results"))
   dir.create(file.path(getwd(), "Results", CLASSIFICATION))

   #save each figure in a folder
   ggsave(plot,
         filename = paste0("Results\\",CLASSIFICATION,"\\",CLASSIFICATION,"_",colnames(data)[i], ".svg"),
         width = 25,
         height = 25,
         units = "cm",
         dpi = 300)


}


