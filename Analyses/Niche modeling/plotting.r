library("dismo")
library("raster")
library("sp")
library("maptools")
library("scales")
library("maps")
library("ENMTools")


setwd("C:/CAML/Viburnum_bolivia/maxent_withoutR")




#get shapefile of some countries
data(wrld_simpl)
shape <- subset(wrld_simpl, NAME=="Peru" | NAME=="Bolivia" | NAME=="Chile" | NAME=="Argentina" )



classifications = c("clas1_cladeBased",
                    "clas2_cladeBasedExtended",
                    "clas3_spFromTree",
                    "clas4_structureK3",
                    "clas5_structureK4")


classifications = c("clas1")


for (classification in classifications){
  cat("Processing ", classification, "...\n")
  
  files <- list.files(paste0("./rectangular_studyarea/", classification, "/"), pattern="_rectangular_studyarea.asc$", full.names = T) 
  
  #set up colors
  i_col = 0
  main_colors <- hue_pal()(length(files))                       
  
  
  svg(paste0("./", classification, "_merged.svg"))
  #Plot initial map
  plot(shape, xlim=c(-80,-60), ylim=c(-27,-8), main=classification)
  # map("world", xlim=c(-80,-60), ylim=c(-27,-8), main=classification)
  
  for (file in files){
      i_col = i_col + 1
      # load(file)
  
      ## crop and mask continent
      # cropped_sA <- crop(studyArea_rectangle, extent(shape))
      # masked_sA <- mask(cropped_sA, shape)
  
      # run prediction with maxent
      # ped2 <- predict(mod, masked_sA)
      
      #load raster prediction produced by maxent
      pred <- read.table(file, skip = 6, header = FALSE, sep = " ")
      #replace missing data from maxent notation into NA
      pred[pred == -9999] <- NA
      #convert to R raster
      pred <- as.raster(as.matrix(pred))
      
      
      #create palette and plot map (piled)
      crp = colorRampPalette(c("white", main_colors[i_col]))( 100 )
      plot(pred, add=T, alpha= 0.9, zlim=0.3, col=crp, legend=F)
  
  }
  
  
  legend("topright",
         legend = mapply(basename, files),
         fill = main_colors,   
         bty="n",
         border="white",
         cex=0.9
         )
  
  dev.off() 
}
