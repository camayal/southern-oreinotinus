### Testing if we have isolation by distance effect.




setwd("C:/Users/camay/Dropbox/Postdocs_research/Edwards-lab/bolivia/isolation_by_distance/")


library(ggplot2)
library(poppr)
library(adegenet)


#trying to fit the ugeno format produced by ipyrad, I converted it into a csv and load it as table
total_dat = read.table("./data/bolivia_history.ugeno.csv", 
                       sep = ",", 
                       # nrows=1000,
                       na.strings = "9")

#adegenet require individuals in rows and each column a marker. I need to transpose my table
total_dat = t(total_dat)


#convert into genind object (adegenet native format)
total_dat <- df2genind(total_dat, ploidy=1, sep="", pop=c(1:46))
total_dat


total_dat_genpop <- genind2genpop(total_dat)
total_dat_genpop


## Load lat/long info
coords <- read.csv("./data/coords.csv", row.names=1)

## Calculate distances
geo_dist_total <- dist(coords, method = "euclidean") #the "dist"function then supposedly estimates relative distances (Euclidean).
# geo_dist_total <- dist(coords, method = "maximum")
gen_dist_total <- dist.genpop(total_dat_genpop, method = 2 ) ## method = 2 is Nei's
# gen_dist_total
# geo_dist_total

## do a mantel test and plot the thing
ibd_total <- mantel.randtest(gen_dist_total, geo_dist_total, nrepet = 1e6)
sink("./mantel_result.txt")
ibd_total
sink()

svg(paste0("./ibd_bars.svg"))
plot(ibd_total)
title("Isolation by distance plot bars")
dev.off()



#This clearly indicates that there is no IBD in our data

svg(paste0("./ibd_plot.svg"))
plot(geo_dist_total, gen_dist_total)
dist_lm <- lm(as.vector(gen_dist_total) ~ as.vector(geo_dist_total))
abline(dist_lm, col="red", lty=2)
title("Isolation by distance plot")
dev.off()

#another plot this time coloring local densities, showing which location in the plot
#has the most points
library(MASS)
dens <- kde2d(as.vector(geo_dist_total), as.vector(gen_dist_total), n=300)
myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))

svg(paste0("./ibd_density.svg"))
plot(geo_dist_total, gen_dist_total, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(dist_lm)
title("Isolation by distance plot with density")
dev.off()










################
################
################


### trying putting clades as pops
total_dat = read.table("./data/bolivia_history.ugeno.csv", 
                       sep = ",", 
                       # nrows=1000,
                       na.strings = "9")

#adegenet require individuals in rows and each column a marker. I need to transpose my table
total_dat = t(total_dat)

#load clade info
groups = read.csv("./data/groups.csv", row.names=1)
total_dat <- df2genind(total_dat, ploidy=1, sep="", pop=groups$clade)
total_dat


total_dat_genpop <- genind2genpop(total_dat)
total_dat_genpop


## Load lat/long info

coords <- read.csv("./data/groups_coords.csv", row.names=1)

## Calculate distances
geo_dist_total <- dist(coords, method = "euclidean") #the "dist"function then supposedly estimates relative distances (Euclidean).
gen_dist_total <- dist.genpop(total_dat_genpop, method = 2 ) ## method = 2 is Nei's
# gen_dist_total
# geo_dist_total

## do a mantel test and plot the thing
ibd_total <- mantel.randtest(gen_dist_total, geo_dist_total)
sink("./mantel_result_byclades.txt")
ibd_total
sink()


svg(paste0("./ibd_bars_byclades.svg"))
plot(ibd_total)
dev.off()

ibd_total

svg(paste0("./ibd_plot_byclades.svg"))
plot(geo_dist_total, gen_dist_total)
dist_lm <- lm(as.vector(gen_dist_total) ~ as.vector(geo_dist_total))
abline(dist_lm, col="red", lty=2)
dev.off()




######################
#######################
##########################3

#Testing with less missing data



#trying to fit the ugeno format produced by ipyrad, I converted it into a csv and load it as table
total_dat = read.table("./minimum_missing_data/bolivian_history_incaseedum_mincov0.75.csv", 
                       sep = ",", 
                       # nrows=1000,
                       na.strings = "9", header = T, row.names = 1)

#adegenet require individuals in rows and each column a marker. I need to transpose my table
# total_dat = t(total_dat)


#convert into genind object (adegenet native format)
total_dat <- df2genind(total_dat, ploidy=1, sep="", pop=c(1:46))
total_dat


total_dat_genpop <- genind2genpop(total_dat)
total_dat_genpop


## Load lat/long info
coords <- read.csv("./data/coords.csv", row.names=1)

## Calculate distances
geo_dist_total <- dist(coords, method = "euclidean") #the "dist"function then supposedly estimates relative distances (Euclidean).
# geo_dist_total <- dist(coords, method = "maximum")
gen_dist_total <- dist.genpop(total_dat_genpop, method = 2 ) ## method = 2 is Nei's
# gen_dist_total
# geo_dist_total

## do a mantel test and plot the thing
ibd_total <- mantel.randtest(gen_dist_total, geo_dist_total, nrepet = 1e6)
sink("./mantel_result_mincov0.75.txt")
ibd_total
sink()

svg(paste0("./ibd_bars_mincov0.75.svg"))
plot(ibd_total, main="mincov = 0.75")
dev.off()


svg(paste0("./ibd_plot_mincov0.75.svg"))
plot(geo_dist_total, gen_dist_total)
dist_lm <- lm(as.vector(gen_dist_total) ~ as.vector(geo_dist_total))
abline(dist_lm, col="red", lty=2)
title("Isolation by distance plot bars mincov = 0.75")
dev.off()



################################
################################
################################
################################
################################



#trying to fit the ugeno format produced by ipyrad, I converted it into a csv and load it as table
total_dat = read.table("./minimum_missing_data/bolivian_history_incaseedum_mincov0.9.csv", 
                       sep = ",", 
                       # nrows=1000,
                       na.strings = "9", header = T, row.names = 1)

#adegenet require individuals in rows and each column a marker. I need to transpose my table
# total_dat = t(total_dat)


#convert into genind object (adegenet native format)
#added pop here, indicating that each individual is one independent pop
#maybe I need to improve that, maybe each pop is a clade or a genogroup
total_dat <- df2genind(total_dat, ploidy=1, sep="", pop=c(1:46))
total_dat


total_dat_genpop <- genind2genpop(total_dat)
total_dat_genpop


## Load lat/long info
coords <- read.csv("./data/coords.csv", row.names=1)

## Calculate distances
geo_dist_total <- dist(coords, method = "euclidean") #the "dist"function then supposedly estimates relative distances (Euclidean).
# geo_dist_total <- dist(coords, method = "maximum")
gen_dist_total <- dist.genpop(total_dat_genpop, method = 2 ) ## method = 2 is Nei's
# gen_dist_total
# geo_dist_total

## do a mantel test and plot the thing
ibd_total <- mantel.randtest(gen_dist_total, geo_dist_total, nrepet = 1e6)
sink("./mantel_result_mincov0.90.txt")
ibd_total
sink()

svg(paste0("./ibd_bars_mincov0.90.svg"))
plot(ibd_total, main="Mantel test for isolation by distance (p=0.018)", 
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 2)
dev.off()


svg(paste0("./ibd_plot_mincov0.90.svg"))
plot(geo_dist_total, gen_dist_total)
dist_lm <- lm(as.vector(gen_dist_total) ~ as.vector(geo_dist_total))
abline(dist_lm, col="red", lty=2)
title("Isolation by distance")
dev.off()



#another plot this time coloring local densities, showing which location in the plot
#has the most points
library(MASS)
dens <- kde2d(as.vector(geo_dist_total), as.vector(gen_dist_total), n=300, lims = c(-1, 20, 0, 1))
myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))

svg(paste0("./ibd_density0.90.svg"))
plot(geo_dist_total, gen_dist_total, pch=20, cex=.5,
     main="Kernel density estimate - Isolation by distance", 
     xlab="Geographic dist. (Euclidian)", 
     ylab="Genetic dist. (Nei's distance)",
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 2)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(dist_lm, lty=2)
dev.off()






#trying to fit the ugeno format produced by ipyrad, I converted it into a csv and load it as table
total_dat = read.table("./minimum_missing_data/bolivian_history_incaseedum_mincov1.0.csv", 
                       sep = ",", 
                       # nrows=1000,
                       na.strings = "9", header = T, row.names = 1)

#adegenet require individuals in rows and each column a marker. I need to transpose my table
# total_dat = t(total_dat)


#convert into genind object (adegenet native format)
total_dat <- df2genind(total_dat, ploidy=1, sep="", pop=c(1:46))
total_dat


total_dat_genpop <- genind2genpop(total_dat)
total_dat_genpop


## Load lat/long info
coords <- read.csv("./data/coords.csv", row.names=1)

## Calculate distances
geo_dist_total <- dist(coords, method = "euclidean") #the "dist"function then supposedly estimates relative distances (Euclidean).
# geo_dist_total <- dist(coords, method = "maximum")
gen_dist_total <- dist.genpop(total_dat_genpop, method = 2 ) ## method = 2 is Nei's
# gen_dist_total
# geo_dist_total

## do a mantel test and plot the thing
ibd_total <- mantel.randtest(gen_dist_total, geo_dist_total, nrepet = 1e6)
sink("./mantel_result_mincov1.0.txt")
ibd_total
sink()

svg(paste0("./ibd_bars_mincov1.0.svg"))
plot(ibd_total, main="mincov = 1.0")
dev.off()


svg(paste0("./ibd_plot_mincov1.0.svg"))
plot(geo_dist_total, gen_dist_total)
dist_lm <- lm(as.vector(gen_dist_total) ~ as.vector(geo_dist_total))
abline(dist_lm, col="red", lty=2)
title("Isolation by distance plot bars")
dev.off()


####################
######################
###################3333
#Now mincov = 0.5 (almost the original missing data) bt with imputation in ipa
#imputation was done using kmeans with k=3
#k=3 shows a similar structure than using non imputation.
#increasing k (7,10, etc), makes groups more distant, but I think this is an artifact.

#trying to fit the ugeno format produced by ipyrad, I converted it into a csv and load it as table
total_dat = read.table("./minimum_missing_data/bolivian_history_incaseedum_mincov0.5_IMPUTED.csv", 
                       sep = ",", 
                       # nrows=1000,
                       na.strings = "9", header = T, row.names = 1)

#adegenet require individuals in rows and each column a marker. I need to transpose my table
# total_dat = t(total_dat)


#convert into genind object (adegenet native format)
#added pop here, indicating that each individual is one independent pop
#maybe I need to improve that, maybe each pop is a clade or a genogroup
total_dat <- df2genind(total_dat, ploidy=1, sep="", pop=c(1:46))
total_dat


total_dat_genpop <- genind2genpop(total_dat)
total_dat_genpop


## Load lat/long info
coords <- read.csv("./data/coords.csv", row.names=1)

## Calculate distances
geo_dist_total <- dist(coords, method = "euclidean") #the "dist"function then supposedly estimates relative distances (Euclidean).
# geo_dist_total <- dist(coords, method = "maximum")
gen_dist_total <- dist.genpop(total_dat_genpop, method = 2 ) ## method = 2 is Nei's
# gen_dist_total
# geo_dist_total

## do a mantel test and plot the thing
ibd_total <- mantel.randtest(gen_dist_total, geo_dist_total, nrepet = 1e6)
sink("./mantel_result_mincov0.5_IMP.txt")
ibd_total
sink()

svg(paste0("./ibd_bars_mincov0.5_IMP.svg"))
plot(ibd_total, main="mincov = 0.5 IMPUTED")
dev.off()


svg(paste0("./ibd_plot_mincov0.5_IMP.svg"))
plot(geo_dist_total, gen_dist_total)
dist_lm <- lm(as.vector(gen_dist_total) ~ as.vector(geo_dist_total))
abline(dist_lm, col="red", lty=2)
title("Isolation by distance plot bars mincov = 0.5 imputed")
dev.off()






#trying to fit the ugeno format produced by ipyrad, I converted it into a csv and load it as table
total_dat = read.table("./minimum_missing_data/bolivian_history_incaseedum_mincov0.75_IMPUTED.csv", 
                       sep = ",", 
                       # nrows=1000,
                       na.strings = "9", header = T, row.names = 1)

#adegenet require individuals in rows and each column a marker. I need to transpose my table
# total_dat = t(total_dat)


#convert into genind object (adegenet native format)
#added pop here, indicating that each individual is one independent pop
#maybe I need to improve that, maybe each pop is a clade or a genogroup
total_dat <- df2genind(total_dat, ploidy=1, sep="", pop=c(1:46))
total_dat


total_dat_genpop <- genind2genpop(total_dat)
total_dat_genpop


## Load lat/long info
coords <- read.csv("./data/coords.csv", row.names=1)

## Calculate distances
geo_dist_total <- dist(coords, method = "euclidean") #the "dist"function then supposedly estimates relative distances (Euclidean).
# geo_dist_total <- dist(coords, method = "maximum")
gen_dist_total <- dist.genpop(total_dat_genpop, method = 2 ) ## method = 2 is Nei's
# gen_dist_total
# geo_dist_total

## do a mantel test and plot the thing
ibd_total <- mantel.randtest(gen_dist_total, geo_dist_total, nrepet = 1e6)
sink("./mantel_result_mincov0.75_IMP.txt")
ibd_total
sink()

svg(paste0("./ibd_bars_mincov0.75_IMP.svg"))
plot(ibd_total, main="mincov = 0.75 IMPUTED")
dev.off()


svg(paste0("./ibd_plot_mincov0.75_IMP.svg"))
plot(geo_dist_total, gen_dist_total)
dist_lm <- lm(as.vector(gen_dist_total) ~ as.vector(geo_dist_total))
abline(dist_lm, col="red", lty=2)
title("Isolation by distance plot bars mincov = 0.75 imputed")
dev.off()







