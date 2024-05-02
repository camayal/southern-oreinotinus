## IMPORTING LIBS

library("raster")
library("dismo")
library("rgeos")

# prepare folders for data input and output
if(!file.exists("C:/CAML/Viburnum_bolivia/")) dir.create("C:/CAML/Viburnum_bolivia/")
if(!file.exists("C:/CAML/Viburnum_bolivia/bioclim")) dir.create("C:/CAML/Viburnum_bolivia/bioclim")
if(!file.exists("C:/CAML/Viburnum_bolivia/studyarea")) dir.create("C:/CAML/Viburnum_bolivia/studyarea")
if(!file.exists("./output")) dir.create("./output")


require(utils)
# download climate data from worldclim.org
# utils::download.file(url="https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_bio.zip",
#                      destfile=paste0("./data/bioclim/wc2.1_10m_bio.zip")) 
# utils::unzip("./data/bioclim/wc2.1_10m_bio.zip",exdir="./data/bioclim/")

# this part was replaced by chelsea layers, bioclim is down.

urls = c("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio10_1981-2010_V.2.1.tif", 
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio11_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio12_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio13_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio14_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio15_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio16_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio17_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio18_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio19_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio1_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio2_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio3_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio4_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio5_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio6_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio7_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio8_1981-2010_V.2.1.tif",
         "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio9_1981-2010_V.2.1.tif"
)

# download climate data 
for (url in urls){
  
  dest = strsplit(url, split = "/", fixed=T)[[1]][11]
  print(dest)
  
  destfile = paste0("C:/CAML/Viburnum_bolivia/bioclim/", dest)
  
  if(!file.exists(destfile)){
    
    utils::download.file(url=url, destfile=destfile)
  }
}


## END SETUP


## START MODELING
## config ONLY RUN ONCE PER SESSION

buffer = 2 # in degrees to determine seudoabsences


setwd("C:/Users/camay/Dropbox/Postdocs_research/Erika-lab/bolivia/nichmodeling/nichmodeling_1/")


# Load env layers
clim_list <- list.files("C:/CAML/Viburnum_bolivia/bioclim/", pattern=".tif$", full.names = T) # '..' leads to the path above the folder where the .rmd file is located

# stacking the bioclim variables to process them at one go 
clim <- raster::stack(clim_list)

# Subset 
occ_full = read.csv("./data/incaseedum_categoriesByCAML_withChelsey.csv")


# create a common study area (same that model for all points as one species)
## Run ONLY ONCE
# modeling the niche only having a very small region (e.g. locality with one pop
# and for that one clade is not good for comparison).
# Instead I will maintain the same range for all to see if there is some
# overlap or not
occ_temp = occ_full
#convert into spatial
coordinates(occ_temp) <- ~ decimalLongitude + decimalLatitude
# define crs for this spatial dataframe
crs(occ_temp) <- "+proj=utm +zone=1 +datum=WGS84"
# this creates a n-decimal-degree buffer around the occurrence data 
occ_temp_buff <- buffer(occ_temp,buffer) 
# crop study area to a manageable extent (rectangle shaped)
studyArea_rectangle <- crop(clim, extent(occ_temp_buff))  
# the 'study area' created by extracting the buffer area from the raster stack
studyArea <- mask(studyArea_rectangle, occ_temp_buff)

# save the new study area rasters as ascii but first create folder
path_studyarea = paste0("C:/CAML/Viburnum_bolivia/studyarea/")
if(!file.exists(path_studyarea)){ 
  dir.create(path_studyarea)
  writeRaster(studyArea,
              # a series of names for output files
              filename=paste0(path_studyarea, names(studyArea),".asc"), 
              format="ascii",  ## the output format
              bylayer=TRUE, ## this will save a series of layers
              overwrite=T)
}

path_studyarea = paste0("C:/CAML/Viburnum_bolivia/rectangular_studyarea/")
if(!file.exists(path_studyarea)){ 
  dir.create(path_studyarea)
  writeRaster(studyArea_rectangle,
              # a series of names for output files
              filename=paste0(path_studyarea, names(studyArea_rectangle),".asc"), 
              format="ascii",  ## the output format
              bylayer=TRUE, ## this will save a series of layers
              overwrite=T)
}


