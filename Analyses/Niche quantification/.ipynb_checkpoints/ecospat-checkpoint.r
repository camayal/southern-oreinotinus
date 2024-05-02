## Following: https://plantarum.ca/2021/07/29/ecospat/

# install.packages('ecospat')


# 
#   #####                                     
#  #     # #       ####  #####    ##   #      
#  #       #      #    # #    #  #  #  #      
#  #  #### #      #    # #####  #    # #      
#  #     # #      #    # #    # ###### #      
#  #     # #      #    # #    # #    # #      
#   #####  ######  ####  #####  #    # ###### 
#                                             
# 
#General execution (needed always)


library(ecospat)
library(raster)
library(maptools)
library(ade4)
library(hash)



setwd("C:/Users/camay/Dropbox/Postdocs_research/Edwards-lab/bolivia/ecospat/")


# Load env layers
clim_list <- list.files("C:/CAML/Viburnum_bolivia/bioclim/", pattern=".tif$", full.names = T) # '..' leads to the path above the folder where the .rmd file is located

# stacking the bioclim variables to process them at one go 
wclim <- raster::stack(clim_list)

# Subset 
occ_full = read.csv("../nichmodeling/nichmodeling_1/data/incaseedum_categoriesByCAML_withChelsey.csv")
coordinates(occ_full) <- ~ decimalLongitude + decimalLatitude


data(wrld_simpl) # load the maptools worldmap

# par(mar = c(0,0, 0, 0))
# plot(wrld_simpl, border = "gray80")
# points(occ_full, pch = 16, col = 2, cex = 0.3)
# 
# 
# 
# par(mar = c(0,0, 3, 1))
# plot(wclim[[1]], main = "bio1")

#extract occ_full values for the climate rasters
#this produce a spatialpoint object
lsOccs <- cbind(occ_full, extract(wclim, occ_full))


# test = subset(occ_full, occ_full$clas2_cladeBasedExtended == "E1")
# 
# 
# plot(test)

## END General execution



# 
#  ###                                                                                   
#   #  #    # #####  # #    # # #####  #    #   ##   #         ##### ######  ####  ##### 
#   #  ##   # #    # # #    # # #    # #    #  #  #  #           #   #      #        #   
#   #  # #  # #    # # #    # # #    # #    # #    # #           #   #####   ####    #   
#   #  #  # # #    # # #    # # #    # #    # ###### #           #   #           #   #   
#   #  #   ## #    # #  #  #  # #    # #    # #    # #           #   #      #    #   #   
#  ### #    # #####  #   ##   # #####   ####  #    # ######      #   ######  ####    #   
#                                                                                        
# 
### Begining of single comparison test


#splitting data
# create categories for comparison
c1 = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == "E1")
c2 = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == "E2")
c3 = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == "E3")
c4 = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == "E4")
c5 = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == "E5")
c6 = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == "E6")
c7 = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == "E7")



par(mar = c(1, 0, 0, 0))
plot(wrld_simpl, axes = FALSE)
points(c1, pch = 16, col = 'red', cex = 0.5)
points(c7, pch = 16, col = 'darkgreen', cex = 0.5)


##### test with c1 c2 only, plan is do this with all clades
## Crop Climate Layers:
c1EnvR <- crop(wclim, c1)
c2EnvR <- crop(wclim, c2)

## Extract values to matrix:
c1EnvM <- getValues(c1EnvR)
c2EnvM <- getValues(c2EnvR)

## Clean out missing values:
c1EnvM <- c1EnvM[complete.cases(c1EnvM), ]
c2EnvM <- c2EnvM[complete.cases(c2EnvM), ]

## Combined global environment:
globalEnvM <- rbind(c1EnvM, c2EnvM)



#niche quantification
pca.clim <- dudi.pca(globalEnvM, center = TRUE,
                     scale = TRUE, scannf = FALSE, nf = 2)
global.scores <- pca.clim$li

c1LS.scores <-
  suprow(pca.clim,
         data.frame(c1)[, colnames(globalEnvM)])$li   
c2LS.scores <-
  suprow(pca.clim,
         data.frame(c1)[, colnames(globalEnvM)])$li

c1Env.scores <- suprow(pca.clim, c1EnvM)$li
c2Env.scores <- suprow(pca.clim, c2EnvM)$li



c1Grid <- ecospat.grid.clim.dyn(global.scores,
                                    c1Env.scores,
                                    c1LS.scores)

c2Grid <- ecospat.grid.clim.dyn(global.scores,
                                      c2Env.scores, 
                                      c2LS.scores)

ecospat.plot.niche.dyn(c1Grid, 
                       c2Grid, 
                       quant = 0.05, 
                       colZ1="#66C2A4", 
                       colZ2="#FB8D61", 
                       col.unf ="#66C2A4",
                       col.exp = "#FB8D61",
                       col.stab = "#eeeeee",
                       interest=1) 

# Explaination of the plot: The resulting plot shows us the environmental conditions present in Eurasia (inside the green line) 
# and North America (inside the red line). The green area represents environments occupied by Lythrum salicaria in Eurasia, 
# but not in North America, the red area shows environments occupied in North America and not Eurasia, and the blue area shows 
# environments occupied in both ranges. We can also see that there are a few areas in Eurasia with environments not present in 
# North America, and vice versa. However, for the most part, Lythrum salicara doesn’t occur in this environments (except for a 
# tiny bit of green in the center of the plot).

# More notes about this plot:
# stability  ## both present
# expansion  ## only sp2 present
# unfilling  ## only sp1 present


# Plot on map
geoGrid <- expand.grid(longitude =
                         seq(-77, -60, length.out = 500),
                       latitude =
                         seq(-7, -20, length.out = 500))


GeoGrid_c1 <- ecospat.grid.clim.dyn(geoGrid, geoGrid, coordinates(c1))
GeoGrid_c2 <- ecospat.grid.clim.dyn(geoGrid, geoGrid, coordinates(c2))



ecospat.plot.niche.dyn(GeoGrid_c1, GeoGrid_c2, quant = 0)
plot(wrld_simpl, add = TRUE)


### End of single comparison test



# 
#     #                                                                                              
#    # #   #      #          ####   ####  #    # #####  # #    #   ##   ##### #  ####  #    #  ####  
#   #   #  #      #         #    # #    # ##  ## #    # # ##   #  #  #    #   # #    # ##   # #      
#  #     # #      #         #      #    # # ## # #####  # # #  # #    #   #   # #    # # #  #  ####  
#  ####### #      #         #      #    # #    # #    # # #  # # ######   #   # #    # #  # #      # 
#  #     # #      #         #    # #    # #    # #    # # #   ## #    #   #   # #    # #   ## #    # 
#  #     # ###### ######     ####   ####  #    # #####  # #    # #    #   #   #  ####  #    #  ####  
#                                                                                                    
# 
### create all possible combination in a programatic way for CLADES

# Define groups
groups = c("E1", "E2", "E3", "E4", "E5", "E6", "E7")

#set hash to put colors according to the other materials
colors <- hash()
colors[["E6"]] = "#FFD92E"
colors[["E4"]] = "#E789C3"
colors[["E3"]] = "#8D9FCA"
colors[["E5"]] = "#A6D753"
colors[["E2"]] = "#FB8D61"
colors[["E1"]] = "#66C2A4"
colors[["E7"]] = "#E4C493"



if(!file.exists("./1-niche_quantification/byclades")) dir.create("./1-niche_quantification/byclades")

# Get all possible pairs
all_comparisons = combn(groups, 2, simplify=F)

# Iterate over each pair and do the niche quantification
for (comparison in all_comparisons){
  
  name = paste(comparison, collapse = "_")
  
  gA = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == comparison[1])
  gB = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == comparison[2])
  
  
  ## Crop Climate Layers:
  gAEnvR <- crop(wclim, gA)
  gBEnvR <- crop(wclim, gB)
  
  ## Extract values to matrix:
  gAEnvM <- getValues(gAEnvR)
  gBEnvM <- getValues(gBEnvR)
  
  ## Clean out missing values:
  gAEnvM <- gAEnvM[complete.cases(gAEnvM), ]
  gBEnvM <- gBEnvM[complete.cases(gBEnvM), ]
  
  ## Combined global environment:
  globalEnvM <- rbind(gAEnvM, gBEnvM)
  
  
  
  #niche quantification
  pca.clim <- dudi.pca(globalEnvM, center = TRUE,
                       scale = TRUE, scannf = FALSE, nf = 2)
  
  global.scores <- pca.clim$li
  
  gALS.scores <-
    suprow(pca.clim,
           data.frame(gA)[, colnames(globalEnvM)])$li   
  gBLS.scores <-
    suprow(pca.clim,
           data.frame(gB)[, colnames(globalEnvM)])$li
  
  gAEnv.scores <- suprow(pca.clim, gAEnvM)$li
  gBEnv.scores <- suprow(pca.clim, gBEnvM)$li
  
  
  
  gAGrid <- ecospat.grid.clim.dyn(global.scores,
                                  gAEnv.scores,
                                  gALS.scores)
  
  gBGrid <- ecospat.grid.clim.dyn(global.scores,
                                  gBEnv.scores, 
                                  gBLS.scores)
  


  #save plto as svg
  svg(paste0("./1-niche_quantification/byclades/niche_quantification_", name, ".svg"))

  ecospat.plot.niche.dyn(gAGrid, 
                         gBGrid, 
                         quant = 0.05, 
                         colZ1=colors[[comparison[1]]],
                         colZ2=colors[[comparison[2]]],
                         col.unf =colors[[comparison[1]]],
                         col.exp = colors[[comparison[2]]],
                         col.stab = "#eeeeee",
                         interest=1,
                         title=name,
                         transparency = 10) 
  
  # More notes about this plot:
  # stability  ## both present
  # expansion  ## only sp2 present
  # unfilling  ## only sp1 present
  

  dev.off()
  
}





### create all possible combination in a programatic way for STRUCTURE GROUPS

# Define groups
groups = c("s1k4", "s2k4", "s3k4", "s4k4", "s5k4")


#set hash to put colors according to the other materials



colors <- hash()
colors[["s1k4"]] = hue_pal()(4)[1]
colors[["s2k4"]] = hue_pal()(4)[2]
colors[["s3k4"]] = hue_pal()(4)[3]
colors[["s4k4"]] = hue_pal()(4)[4]
colors[["s5k4"]] = hue_pal()(4)[5]



if(!file.exists("./1-niche_quantification/byGeneticGroups")) dir.create("./1-niche_quantification/byGeneticGroups")

# Get all possible pairs
all_comparisons = combn(groups, 2, simplify=F)

# Iterate over each pair and do the niche quantification
for (comparison in all_comparisons){
  
  name = paste(comparison, collapse = "_")
  
  gA = subset(lsOccs, lsOccs$clas5_structureK4 == comparison[1])
  gB = subset(lsOccs, lsOccs$clas5_structureK4 == comparison[2])
  
  
  ## Crop Climate Layers:
  gAEnvR <- crop(wclim, gA)
  gBEnvR <- crop(wclim, gB)
  
  ## Extract values to matrix:
  gAEnvM <- getValues(gAEnvR)
  gBEnvM <- getValues(gBEnvR)
  
  ## Clean out missing values:
  gAEnvM <- gAEnvM[complete.cases(gAEnvM), ]
  gBEnvM <- gBEnvM[complete.cases(gBEnvM), ]
  
  ## Combined global environment:
  globalEnvM <- rbind(gAEnvM, gBEnvM)
  
  
  
  #niche quantification
  pca.clim <- dudi.pca(globalEnvM, center = TRUE,
                       scale = TRUE, scannf = FALSE, nf = 2)
  
  global.scores <- pca.clim$li
  
  gALS.scores <-
    suprow(pca.clim,
           data.frame(gA)[, colnames(globalEnvM)])$li   
  gBLS.scores <-
    suprow(pca.clim,
           data.frame(gB)[, colnames(globalEnvM)])$li
  
  gAEnv.scores <- suprow(pca.clim, gAEnvM)$li
  gBEnv.scores <- suprow(pca.clim, gBEnvM)$li
  
  
  
  gAGrid <- ecospat.grid.clim.dyn(global.scores,
                                  gAEnv.scores,
                                  gALS.scores)
  
  gBGrid <- ecospat.grid.clim.dyn(global.scores,
                                  gBEnv.scores, 
                                  gBLS.scores)
  
  
  
  #save plto as svg
  svg(paste0("./1-niche_quantification/byGeneticGroups/niche_quantification_", name, ".svg"))
  
  ecospat.plot.niche.dyn(gAGrid, 
                         gBGrid, 
                         quant = 0.05, 
                         colZ1=colors[[comparison[1]]],
                         colZ2=colors[[comparison[2]]],
                         col.unf =colors[[comparison[1]]],
                         col.exp = colors[[comparison[2]]],
                         col.stab = "#eeeeee",
                         interest=1,
                         title=name,
                         transparency = 10) 
  
  # More notes about this plot:
  # stability  ## both present
  # expansion  ## only sp2 present
  # unfilling  ## only sp1 present
  
  
  dev.off()
  
}


### create all possible combination in a programatic way for taxonomic groups

# Define groups
groups = c("incarum", "seemenii", "new_name_2")


#set hash to put colors according to the other materials



colors <- hash()
colors[["incarum"]] = hue_pal()(3)[1]
colors[["seemenii"]] = hue_pal()(3)[2]
colors[["new_name_2"]] = hue_pal()(3)[3]




if(!file.exists("./1-niche_quantification/byTaxonomic")) dir.create("./1-niche_quantification/byTaxonomic")

# Get all possible pairs
all_comparisons = combn(groups, 2, simplify=F)

# Iterate over each pair and do the niche quantification
for (comparison in all_comparisons){
  
  name = paste(comparison, collapse = "_")
  
  gA = subset(lsOccs, lsOccs$clas3_spFromTree == comparison[1])
  gB = subset(lsOccs, lsOccs$clas3_spFromTree == comparison[2])
  
  
  ## Crop Climate Layers:
  gAEnvR <- crop(wclim, gA)
  gBEnvR <- crop(wclim, gB)
  
  ## Extract values to matrix:
  gAEnvM <- getValues(gAEnvR)
  gBEnvM <- getValues(gBEnvR)
  
  ## Clean out missing values:
  gAEnvM <- gAEnvM[complete.cases(gAEnvM), ]
  gBEnvM <- gBEnvM[complete.cases(gBEnvM), ]
  
  ## Combined global environment:
  globalEnvM <- rbind(gAEnvM, gBEnvM)
  
  
  
  #niche quantification
  pca.clim <- dudi.pca(globalEnvM, center = TRUE,
                       scale = TRUE, scannf = FALSE, nf = 2)
  
  global.scores <- pca.clim$li
  
  gALS.scores <-
    suprow(pca.clim,
           data.frame(gA)[, colnames(globalEnvM)])$li   
  gBLS.scores <-
    suprow(pca.clim,
           data.frame(gB)[, colnames(globalEnvM)])$li
  
  gAEnv.scores <- suprow(pca.clim, gAEnvM)$li
  gBEnv.scores <- suprow(pca.clim, gBEnvM)$li
  
  
  
  gAGrid <- ecospat.grid.clim.dyn(global.scores,
                                  gAEnv.scores,
                                  gALS.scores)
  
  gBGrid <- ecospat.grid.clim.dyn(global.scores,
                                  gBEnv.scores, 
                                  gBLS.scores)
  
  
  
  #save plto as svg
  svg(paste0("./1-niche_quantification/byTaxonomic/niche_quantification_", name, ".svg"))
  
  ecospat.plot.niche.dyn(gAGrid, 
                         gBGrid, 
                         quant = 0.05, 
                         colZ1=colors[[comparison[1]]],
                         colZ2=colors[[comparison[2]]],
                         col.unf =colors[[comparison[1]]],
                         col.exp = colors[[comparison[2]]],
                         col.stab = "#eeeeee",
                         interest=1,
                         title=name,
                         transparency = 10) 
  
  # More notes about this plot:
  # stability  ## both present
  # expansion  ## only sp2 present
  # unfilling  ## only sp1 present
  
  
  dev.off()
  
}





###### Trying to plot all niches in one single plot
# the idea is getting the grid with ecospat and after that just plotin the w elements in the grid object
#done by hand -_- quick solution


colors <- hash()
colors[["S1k3"]] = "#F88F8A40"
colors[["S2k3"]] = "#33B8FF40"
colors[["S3k3"]] = "#33CA8540"
#weird bug if I put alpha in plot(), so I am adding alpha here (xxxxxx40)


gA = subset(lsOccs, lsOccs$clas4_structureK3 == groups[1])
gB = subset(lsOccs, lsOccs$clas4_structureK3 == groups[2])
gC = subset(lsOccs, lsOccs$clas4_structureK3 == groups[3])


## Crop Climate Layers:
gAEnvR <- crop(wclim, gA)
gBEnvR <- crop(wclim, gB)
gCEnvR <- crop(wclim, gC)

## Extract values to matrix:
gAEnvM <- getValues(gAEnvR)
gBEnvM <- getValues(gBEnvR)
gCEnvM <- getValues(gCEnvR)


## Clean out missing values:
gAEnvM <- gAEnvM[complete.cases(gAEnvM), ]
gBEnvM <- gBEnvM[complete.cases(gBEnvM), ]
gCEnvM <- gCEnvM[complete.cases(gCEnvM), ]

## Combined global environment:
globalEnvM <- rbind(gAEnvM,gBEnvM,gCEnvM)



#niche quantification
pca.clim <- dudi.pca(globalEnvM, center = TRUE,
                     scale = TRUE, scannf = FALSE, nf = 2)

global.scores <- pca.clim$li

gALS.scores <-
  suprow(pca.clim,
         data.frame(gA)[, colnames(globalEnvM)])$li  
gBLS.scores <-
  suprow(pca.clim,
         data.frame(gB)[, colnames(globalEnvM)])$li
gCLS.scores <-
  suprow(pca.clim,
         data.frame(gC)[, colnames(globalEnvM)])$li


gAEnv.scores <- suprow(pca.clim, gAEnvM)$li
gBEnv.scores <- suprow(pca.clim, gBEnvM)$li
gCEnv.scores <- suprow(pca.clim, gCEnvM)$li



gAGrid <- ecospat.grid.clim.dyn(global.scores,
                                gAEnv.scores,
                                gALS.scores)
gBGrid <- ecospat.grid.clim.dyn(global.scores,
                                gBEnv.scores,
                                gBLS.scores)
gCGrid <- ecospat.grid.clim.dyn(global.scores,
                                gCEnv.scores,
                                gCLS.scores)


svg(paste0("./1-niche_quantification/byGeneticGroups/niche_quantification_ALLinONE.svg"))
plot(gAGrid$w, col=c(colors[[groups[1]]]), breaks=c(0.1,1.0))
plot(gBGrid$w, col=c(colors[[groups[2]]]), breaks=c(0.1,1.0), add=T )
plot(gCGrid$w, col=c(colors[[groups[3]]]), breaks=c(0.1,1.0), add=T )
dev.off()


###now for clades, by hand again _-_ 

groups = c("E1", "E2", "E3", "E4", "E5", "E6", "E7")

#set hash to put colors according to the other materials
colors <- hash()
colors[["E6"]] = "#FFD92E40"
colors[["E4"]] = "#E789C340"
colors[["E3"]] = "#8D9FCA40"
colors[["E5"]] = "#A6D75340"
colors[["E2"]] = "#FB8D6140"
colors[["E1"]] = "#66C2A440"
colors[["E7"]] = "#E4C49340"
#weird bug if I put alpha in plot(), so I am adding alpha here (xxxxxx40)


gA = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == groups[1])
gB = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == groups[2])
gC = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == groups[3])
gD = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == groups[4])
gE = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == groups[5])
gF = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == groups[6])
gG = subset(lsOccs, lsOccs$clas2_cladeBasedExtended == groups[7])

## Crop Climate Layers:
gAEnvR <- crop(wclim, gA)
gBEnvR <- crop(wclim, gB)
gCEnvR <- crop(wclim, gC)
gDEnvR <- crop(wclim, gD)
gEEnvR <- crop(wclim, gE)
gFEnvR <- crop(wclim, gF)
gGEnvR <- crop(wclim, gG)

## Extract values to matrix:
gAEnvM <- getValues(gAEnvR)
gBEnvM <- getValues(gBEnvR)
gCEnvM <- getValues(gCEnvR)
gDEnvM <- getValues(gDEnvR)
gEEnvM <- getValues(gEEnvR)
gFEnvM <- getValues(gFEnvR)
gGEnvM <- getValues(gGEnvR)

## Clean out missing values:
gAEnvM <- gAEnvM[complete.cases(gAEnvM), ]
gBEnvM <- gBEnvM[complete.cases(gBEnvM), ]
gCEnvM <- gCEnvM[complete.cases(gCEnvM), ]
gDEnvM <- gDEnvM[complete.cases(gDEnvM), ]
gEEnvM <- gEEnvM[complete.cases(gEEnvM), ]
gFEnvM <- gFEnvM[complete.cases(gFEnvM), ]
gGEnvM <- gGEnvM[complete.cases(gGEnvM), ]

## Combined global environment:
globalEnvM <- rbind(gAEnvM,gBEnvM,gCEnvM,gDEnvM,gEEnvM,gFEnvM,gGEnvM)



#niche quantification
pca.clim <- dudi.pca(globalEnvM, center = TRUE,
                     scale = TRUE, scannf = FALSE, nf = 2)

global.scores <- pca.clim$li

gALS.scores <-
  suprow(pca.clim,
         data.frame(gA)[, colnames(globalEnvM)])$li  
gBLS.scores <-
  suprow(pca.clim,
         data.frame(gB)[, colnames(globalEnvM)])$li
gCLS.scores <-
  suprow(pca.clim,
         data.frame(gC)[, colnames(globalEnvM)])$li
gDLS.scores <-
  suprow(pca.clim,
         data.frame(gD)[, colnames(globalEnvM)])$li
gELS.scores <-
  suprow(pca.clim,
         data.frame(gE)[, colnames(globalEnvM)])$li
gFLS.scores <-
  suprow(pca.clim,
         data.frame(gF)[, colnames(globalEnvM)])$li
gGLS.scores <-
  suprow(pca.clim,
         data.frame(gG)[, colnames(globalEnvM)])$li


gAEnv.scores <- suprow(pca.clim, gAEnvM)$li
gBEnv.scores <- suprow(pca.clim, gBEnvM)$li
gCEnv.scores <- suprow(pca.clim, gCEnvM)$li
gDEnv.scores <- suprow(pca.clim, gDEnvM)$li
gEEnv.scores <- suprow(pca.clim, gEEnvM)$li
gFEnv.scores <- suprow(pca.clim, gFEnvM)$li
gGEnv.scores <- suprow(pca.clim, gGEnvM)$li



gAGrid <- ecospat.grid.clim.dyn(global.scores,
                                gAEnv.scores,
                                gALS.scores)
gBGrid <- ecospat.grid.clim.dyn(global.scores,
                                gBEnv.scores,
                                gBLS.scores)
gCGrid <- ecospat.grid.clim.dyn(global.scores,
                                gCEnv.scores,
                                gCLS.scores)
gDGrid <- ecospat.grid.clim.dyn(global.scores,
                                gDEnv.scores,
                                gDLS.scores)
gEGrid <- ecospat.grid.clim.dyn(global.scores,
                                gEEnv.scores,
                                gELS.scores)
gFGrid <- ecospat.grid.clim.dyn(global.scores,
                                gFEnv.scores,
                                gFLS.scores)
gGGrid <- ecospat.grid.clim.dyn(global.scores,
                                gGEnv.scores,
                                gGLS.scores)


svg(paste0("./1-niche_quantification/byclades/niche_quantification_ALLinONE.svg"))
plot(gAGrid$w, col=c(colors[[groups[1]]]), breaks=c(0.1,1.0))
plot(gBGrid$w, col=c(colors[[groups[2]]]), breaks=c(0.1,1.0), add=T )
plot(gCGrid$w, col=c(colors[[groups[3]]]), breaks=c(0.1,1.0), add=T )
plot(gDGrid$w, col=c(colors[[groups[4]]]), breaks=c(0.1,1.0), add=T )
plot(gEGrid$w, col=c(colors[[groups[5]]]), breaks=c(0.1,1.0), add=T )
plot(gFGrid$w, col=c(colors[[groups[6]]]), breaks=c(0.1,1.0), add=T )
plot(gGGrid$w, col=c(colors[[groups[7]]]), breaks=c(0.1,1.0), add=T )
dev.off()






