#########################
## SNOW COVER DURATION ##
#########################

# Info ----
# Daily snow cover at 20 m from Sentinel-2 and MODIS 
# using random forest algorithm
# Contact person: Chiara Richiardi - richiardi@iia.cnr.it

# Build dataset with the following variables:
# Topographic features (Elevation, slope, aspect)
# MODIS Snow
# DOY
# Lat, Lon
# Year

# APPROACH: 
# Regression RF of MODIS FSC + regression RF on Sentinel-2 NDSI 
# to get daily NDSI & binary snow/snow free maps

# Libraries ----
# Installing packages
#install.packages("caTools")    # For sampling the dataset
#install.packages("ranger")    # For implementing random forest algorithm
#install.packages('lubridate')
#install.packages("writexl")
#install.packages("sp")
#install.packages("sf")
#install.packages("fields")
#install.packages("MASS")
#install.packages("parsnip")
#install.packages("tuneRanger")

# Loading package
library(caTools)
library(ranger)
library(devtools)
library(raster)
library(gdalUtils)
library(rgdal)
library(sf)
library(data.table)
library(lubridate)
library(sp)
library(parallel)
library(doParallel)
library(fields)
library(MASS)
library(dplyr)

numCores <- detectCores() 
registerDoParallel(numCores)

# Set folders ----

yy <- as.list(2015:2021)   #set the desired year (from 2015 to 2021)

root <- getwd() # SET THE 
L2folder <- file.path(paste0(root,"/Sentinel2/",tl,"/Level2"))
L2_yy <- file.path(paste0(L2folder,"/",yy))

rfFolder <- file.path(paste0(root,"/SCD/Approccio2"))
dir.create(rfFolder,showWarnings = FALSE)

modisFolder <- file.path(root,"/MODIS/MOD10A1")
modisGF <- file.path(paste0(modisFolder,"/MOD10A1_GapFilled"))

### STEP 1: Regression RF of MODIS FSC #########################################
# Parameters ----

#specify the number of points to be randomly sampled from the good indices each day. 
sub.size <- 25 #% of pixels per altitudinal band selected for training the regression RF

# Load and organize training data ## ----

# Topographic features and AOI
setwd(root)
aoi <- shapefile(paste0("TAGLIO_",sito,"_",tl,"_S2_32632.shp"))
dem <- setMinMax(raster::raster(paste0("DEM_",sito,"_modis_buffer.tif")))
#slope <- raster::raster(paste0("Slope_",sito,"_modis_buffer.tif"))
#aspect <- raster::raster(paste0("Aspect_",sito,"_modis_buffer.tif"))
slope <- terrain(dem, opt="slope", unit="degrees", neighbors=8)
aspect <- terrain(dem, opt="aspect", unit="degrees", neighbors=8)


# MODIS
modis <- raster("D:/ABRESO/PNGP/MODIS/MOD10A1/173130042/MOD10A1_A2000090_h18v04_061_2020040195629_MOD_Grid_Snow_500m_NDSI_Snow_Cover_de15620c.tif")
aoi_modis <- buffer(aoi, width=1000) %>% 
  spTransform(aoi, CRSobj = crs(modis))

# Lat, Lon
locs <- xyFromCell(dem, 1:ncell(dem), spatial=FALSE)
locs.x <- locs[,1]
locs.y <- locs[,2]

# Build the training dataset)----
# get list of DOYs from 2015 to 2021
all.days <- as.Date(as.Date("2015-07-30"):as.Date("2021-12-31"), origin="1970-01-01")
n.days <- length(all.days)

# from every MODIS build a dataset for the gap filling
columns <- c("MODIS_FSC", "Elevation", "Slope", "Aspect","X","Y","DOY","Year")
df <- matrix(ncol = length(columns)) 
colnames(df) = columns
data.stack <- stack()

setwd(modisFolder)
n <- 1

for (i in all.days[1:n.days]) {
  
  d <- all.days[n]
  doy <- strftime(d, format = "%j")    #convert date in doy
  y <- strftime(d, format = "%Y")
  setwd(modisFolder)
  if (nchar(as.character(doy)) == 3){
    modis.name <- list.files(pattern=glob2rx(paste0("MOD10A1_A",y,doy,"_*_NDSI_Snow_Cover*.tif")),recursive = TRUE)
  } else if (nchar(as.character(doy)) == 2){
    modis.name <- list.files(pattern=glob2rx(paste0("MOD10A1_A",y,"0",doy,"_*_NDSI_Snow_Cover*.tif")),recursive = TRUE)
  } else if (nchar(as.character(doy)) == 1){
    modis.name <- list.files(pattern=glob2rx(paste0("MOD10A1_A",y,"00",doy,"_*_NDSI_Snow_Cover*.tif")),recursive = TRUE)
  }
  
  if (length(modis.name)>=1){
    modis <- raster(modis.name[1])                          # Get MODIS FSC raster
    modis[modis > 100] <- NA                                # 0-100: fractional snow cover, set flags as no data
    modis <- crop(modis, aoi_modis)                         # Crop on aoi
    modis.stack <- stack(modis,dem,slope,aspect)            # create stack for each snow cover scene
    modis.data <- as.matrix(modis.stack)                    # convert in matrix
    modis.data <- cbind(modis.data,locs.x,locs.y,doy,y)     # add predictive variables
    colnames(modis.data) = columns                          # assign column names
    modis.data <- na.omit(modis.data)                       # omit pixels with no data
    
    # select random pixel 
    train.size <- round((nrow(modis.data)*sub.size/100),0)
    if (train.size < nrow(modis.data) & train.size != 1){
      set.seed(489)
      train.data <- as.matrix(modis.data[sample(nrow(modis.data), train.size), ])
      df <- rbind(df,train.data)
      print(paste0("doy ", doy, " of year ", y, " ok - loop progress ", n,"/",n.days, ", df rows: ",nrow(df)))
      n <- n+1  
    } else {
      print(paste0("doy ", doy, " of year ", y, " all cloudy - loop progress ", n,"/",n.days))
      n <- n+1
      next
    }
  } else {
    print(paste0("doy ", doy, " of year ", y, " missing - loop progress ", n,"/",n.days))
    n <- n+1
    next
  }
}
temp <- df
df <-    #Delete first empty row
  #df <- df[-1,]
  df <- df %>% 
  df[-1,] %>% 
  as.data.frame(df) %>% 
  mutate_if(is.character,as.numeric)


setwd(rfFolder)
#df <- na.omit(df)
#df <- apply(df, 2, as.integer) %>% # convert al values in integer to decrease size
#as.data.frame
write.matrix(df, file = paste0("Training_dataset_MODIS_FSC_regression_subset",sub.size,".csv"), sep = ",")
#df <- read.csv(paste0(rfFolder,"/Training_dataset_MODIS_NDSI_regression_subset",sub.size,".csv"), header = TRUE, sep = ",") 

# Free memory
rm(sc.stack,sc.data,train.data,dem_band)
rm(aoi,dem,slope,aspect,locs,locs.x,locs.y,modis) #lulc
gc()


# Growing and saving the regression forest ----

print(Sys.time())
ranger.modis.regression <- ranger(MODIS_FSC ~ ., 
                                  data = modisData,
                                  num.trees = 500,   
                                  mtry = 7, #round(((length(columns)-1)/3),digits = 0)
                                  min.node.size = 3,
                                  #importance = 'impurity_corrected',
                                  splitrule = "variance",
                                  save.memory = TRUE,
                                  verbose = TRUE,
                                  num.threads = numCores)

print(Sys.time())

print(object.size(ranger.modis.regression),units="GB")
ranger.modis.regression$prediction.error

print(Sys.time())
save(ranger.modis.regression ,file=paste0(rfFolder,"/ranger.modis.FSC.regression_sub.size",as.character(sub.size),".Rda"))
print(Sys.time())
#rm(ranger.regression)

#variable importance
feat.imp <- importance(ranger.classifier)

names.imp <- names(feat.imp)

feat.imp <- sort(feat.imp)

png(paste0(rfFolder,"/MODIS_importance.png"),width=2000,height=800)
par(mfrow=c(1,2),mai=c(1,4,1,1))
barplot(feat.imp/1e6,xlab='Importance (scaled)',las=1,cex.lab=2.5, cex.axis=2, cex.names=2.5,horiz=TRUE)
dev.off()


# Build the gap filled MODIS dataset ----
load(paste0(rfFolder,"/ranger.modis.FSC.regression_sub.size",as.character(sub.size),".Rda"))

# get list of DOYs from 2015 to 2021
all.days <- as.Date(as.Date("2015-07-30"):as.Date("2021-12-31"), origin="1970-01-01")
n.days <- length(all.days)

# from every MODIS build a dataset for the gap filling
columns <- c("MODIS", "Elevation", "Slope", "Aspect","X","Y","DOY","Year")
df.SCD <- matrix(ncol = length(columns)) 
colnames(df.SCD) = columns
SCD.stack <- stack()

rr <- raster("D:/ABRESO/PNGP/MODIS/MOD10A1/173130042/MOD10A1_A2000090_h18v04_061_2020040195629_MOD_Grid_Snow_500m_NDSI_Snow_Cover_de15620c.tif")
aoi_modis <- buffer(aoi, width=1000) %>% 
  spTransform(aoi, CRSobj = crs(modis))

setwd(modisFolder)

modis_output <- file.path(paste0(root,"/MODIS/MOD10A1/MOD10A1_GapFilled"))

n <- 1

for (i in all.days[1:n.days]) {
  
  d <- all.days[n]
  doy <- strftime(d, format = "%j")    #convert date in doy
  y <- strftime(d, format = "%Y")
  setwd(modisFolder)
  if (nchar(as.character(doy)) == 3){
    modis.name <- list.files(pattern=glob2rx(paste0("MOD10A1_A",y,doy,"_*_NDSI_Snow_Cover*.tif")),recursive = TRUE)
  } else if (nchar(as.character(doy)) == 2){
    modis.name <- list.files(pattern=glob2rx(paste0("MOD10A1_A",y,"0",doy,"_*_NDSI_Snow_Cover*.tif")),recursive = TRUE)
  } else if (nchar(as.character(doy)) == 1){
    modis.name <- list.files(pattern=glob2rx(paste0("MOD10A1_A",y,"00",doy,"_*_NDSI_Snow_Cover*.tif")),recursive = TRUE)
  }
  
  #Find if there is a MODIS scene and if not create an empty raster
  if (length(modis.name)>=1){
    modis <- raster(modis.name[1])                          # Get MODIS FSC raster
    modis[modis > 100] <- NA                                # 0-100: fractional snow cover, set flags as no data
    modis <- crop(modis, aoi_modis)                         # Crop on aoi
    modis.stack <- stack(modis,dem,slope,aspect)            # create stack for each snow cover scene
    modis.data <- as.matrix(modis.stack)                    # convert in matrix
    modis.data <- cbind(modis.data,locs.x,locs.y,doy,y)     # add predictive variables
    colnames(modis.data) = columns                          # assign column names
    modis.data <- na.omit(modis.data)                       # omit pixels with no data
  } else {
    modis <- raster(ncol=ncol(rr), nrow=nrow(rr), ext=extent(rr), crs = crs(rr))  #create an empty raster
    modis <- crop(modis, aoi_modis)
    values(modis) <- NA    # assign all value as NA
  }
  
  SCD.stack <- stack(modis,dem,slope,aspect)    # create stack for each snow cover scene
  SCD.data <- as.matrix(SCD.stack) #%>%     # convert in matrix
  SCD.data <- cbind(SCD.data,locs.x,locs.y,doy,y)  
  df.SCD <- rbind(df.SCD,SCD.data)
  df.SCD <- df.SCD[-1,]    #Delete first empty row
  df.SCD <- as.data.frame(df.SCD)
  colnames(df.SCD) = columns
  df.SCD <- mutate_if(df.SCD, is.character,as.numeric) # Convert as numeric
  df.SCD <- subset(df.SCD, !is.na(df.SCD$Elevation))
  regressionvalues <- df.SCD[is.na(df.SCD$MODIS), ]
  
  ## 
  ## Predict snow data
  ##
  
  if (nrow(regressionvalues) >=1){
    ranger.temp <- predict(ranger.modis.regression,
                           data = regressionvalues,
                           predict.all = FALSE,
                           num.trees = ranger.modis.regression$num.trees,
                           type = "response",   #terminalNodes
                           se.method = "infjack",
                           seed = 567,
                           num.threads = numCores,
                           verbose = TRUE,
                           inbag.counts = NULL)
    
    # Add retrieved values to the s2 raster
    df.SCD$MODIS[is.na(df.SCD$MODIS)] <- ranger.temp$predictions
    df.SCD <- data.frame(lapply(df.SCD, function(x) ifelse(!is.na(as.numeric(x)), as.numeric(x), x)))
    coordinates(df.SCD) <- ~ X + Y
    # coerce to SpatialPixelsDataFrame
    gridded(df.SCD) <- TRUE
    # coerce to raster
    r <- raster(df.SCD)
  } else {
    r <- modis
  }
  crs(r) <- crs(modis)
  
  # Save the gap filled NDSI raster
  setwd(modis_output)
  writeRaster(r, paste0("MOD10A1_GF_A",y,doy,"_h18v04_061_MOD_Grid_Snow_500m_NDSI_SnowCover.tif"), datatype='INT2S', format = "GTiff", overwrite=TRUE, option="COMPRESS=LZW", prj=TRUE)
  
  # Empty dataframe for the next scene
  df.SCD <- matrix(ncol = length(columns)) 
  colnames(df.SCD) = columns
  SCD.stack <- stack()
  
  print(paste0("doy ", doy, " of year ", y, " ok - loop progress ", n,"/",n.days))
  n <- n+1  
}




### STEP 2: Regression RF of Sentinel-2 NDSI ###################################
# Parameters ----

#specify the number of points to be randomly sampled from the good indices each day. 
sub.size <- 3 #% of pixels per altitudinal band selected for training the regression RF

#specify the cloud threshold for the snow cover ground truth
ct <- 10

#specify altitudinal band threshold (in meters)
bh <- 100

# Load and organize training data ----

# Topographic features and AOI 
setwd(root)
aoi <- shapefile(paste0("TAGLIO_",sito,"_",tl,"_S2_32632.shp"))
dem <- setMinMax(raster::raster(paste0("DTM_",sito,"_32632_20m.tif")))
#slope <- raster::raster(paste0("Slope_",sito,"_32632_20m.tif"))
#aspect <- raster::raster(paste0("Aspect_",sito,"_32632_20m.tif"))
slope <- terrain(dem, opt="slope", unit="degrees", neighbors=8)
aspect <- terrain(dem, opt="aspect", unit="degrees", neighbors=8)

# Convert 
dem <- round(dem, digits = 0)
slope <- round(slope, digits = 0)
aspect <- round(aspect, digits = 0)
TRI <- (round(TRI, digits = 1))*10

# Lat, Lon 
locs <- xyFromCell(dem, 1:ncell(dem), spatial=FALSE)
locs.x <- locs[,1]
locs.y <- locs[,2]

# MODIS
modisEX <- raster("D:/ABRESO/PNGP/MODIS/MOD10A1/173130042/MOD10A1_A2000090_h18v04_061_2020040195629_MOD_Grid_Snow_500m_NDSI_Snow_Cover_de15620c.tif")
aoi_modis <- spTransform(aoi, crs(modisEX))

# Sentinel-2 Snow Cover (build the training dataset)----
sc.names <- list()
columns <- c("SNOW", "Elevation", "Slope", "Aspect","MODIS","X","Y","DOY","Year") #"LULC"
df <- matrix(ncol = length(columns)) 
colnames(df) = columns
sc.stack <- stack()

for (y in yy){
  
  L2_yy_out <- file.path(paste0(L2folder,"/",y,"/output_",tl))
  setwd(L2_yy_out)
  
  sc <- list.files(pattern=glob2rx(paste0("Snow&Cloud_mask_20m_",y,"*TP*.tif")))
  num_sc <- length(sc)
  n <- 1
  
  for (s in sc[1:num_sc]){
    setwd(L2_yy_out)
    sr <- raster(s)
    total_cloud_cells <- sum(sr [!is.na(sr)]  == 205)
    total_cells <- ncell(sr[sr != 254])
    total_cloud_fraction <- round(total_cloud_cells*100/total_cells, 0)
    
    if (total_cloud_fraction < ct){
      sr[sr > 100] <- NA
      sc.names <- append(sc.names,s)
      sc.date <- substr(s,start=21,stop=28)
      sc.doy <- strptime(x=as.character(sc.date),format="%Y%m%d")
      doy <- yday(as.character.Date(sc.doy))    # get DOY
      m <- substr(s,start=25,stop=26)
      g <- substr(s,start=27,stop=28)
      ndsi <- raster::raster(list.files(pattern=glob2rx(paste0("NDSI_20m_",y,m,g,"*.tif"))))
      
      setwd(modisGF)
      if (nchar(as.character(doy)) == 3){
        modis.name <- list.files(pattern=glob2rx(paste0("MOD10A1_GF_A",y,doy,"_*_NDSI_SnowCover.tif")),recursive = TRUE)
      } else if (nchar(as.character(doy)) == 2){
        modis.name <- list.files(pattern=glob2rx(paste0("MOD10A1_GF_A",y,"0",doy,"_*_NDSI_SnowCover.tif")),recursive = TRUE)
      } else if (nchar(as.character(doy)) == 1){
        modis.name <- list.files(pattern=glob2rx(paste0("MOD10A1_GF_A",y,"00",doy,"_*_NDSI_SnowCover.tif")),recursive = TRUE)
      }
      
      modis <- raster::raster(modis.name[1])
      crs(modis) <- crs(modisEX)
      # Modis pre-processing
      modis <- crop(modis, aoi_modis)
      # check CRS
      if (compareCRS(sr, modis) == FALSE){
        modis <- projectRaster(modis, crs = crs(sr), res=xres(sr))
      }
      # check resolution
      rres <- xres(modis)/xres(sr)
      if (rres != 1){
        modis <- raster::resample(modis, sr, method="ngb")
      }
      # check extent
      if (extent(modis) != extent(sr)){
        modis <- crop(modis, sr)
      }
      if (extent(modis) != extent(sr)){
        modis <- setExtent(modis, extent(sr), keepres=FALSE, snap=FALSE)
      }
      
      sc.stack <- stack(ndsi,dem,slope,aspect,modis)  # create stack for each snow cover scene
      sc.data <- as.data.frame(sc.stack)    # convert in dataframe
      sc.data <- cbind(sc.data,locs.x,locs.y,doy,y)
      colnames(sc.data) = columns
      sc.data <- na.omit(sc.data)    # omit pixels with no data
      
      # select random pixel for each DEM altitudinal bands
      total_dem_cells <- ncell(dem)
      DEMmin <- round(minValue(dem),-2) - bh
      zmax <- DEMmax <- round(maxValue(dem),-2) + bh
      zmin <- zmax - bh
      size <- 0   # counter
      
      while (zmin >= DEMmin){
        
        dem_band <- overlay(dem, sr, fun = function(x,y) {return(x > zmin & x <= zmax & !is.na(y))})
        dem_band_pixels <- sum(dem_band[!is.na(dem_band)])
        train.size <- round((dem_band_pixels*sub.size/100),0)   # take % from every altitudinal band
        train.data <- sc.data[sc.data$Elevation <= zmax & sc.data$Elevation > zmin, ]
        set.seed(489)
        train.data <- train.data[sample(nrow(train.data), train.size), ]
        train.data <- as.matrix(train.data)
        df <- rbind(df,train.data)
        
        zmax <- zmax - bh
        zmin <- zmin - bh
        size <- size + train.size
      }
      tot.size <- round((size*100/total_cells),0)
      print(paste0(tot.size,"% of the scene added to training data"))
      #df <- rbind(df,train.data)
      setwd(L2_yy_out)
    }
    print(paste0(n,"/", num_sc, " done"))
    n <- n+1
  }
  msg1 <- paste0("Year ",y," completed")
  print(msg1)
}

# Free memory
#rm(sc.stack,sc.data,train.data,dem_band)
#rm(aoi,dem,slope,aspect,locs,locs.x,locs.y,sr,modis) #lulc
gc()

setwd(rfFolder)

df <-
  na.omit() %>%
  as.data.frame() %>% 
  mutate_if(is.character,as.numeric)

setwd(rfFolder)
write.matrix(df, file = paste0("Training_dataset_S2.NDSI_modis.FSC_regression_subset",sub.size,".csv"), sep = ",")

#snowData <- read.csv(paste0(rfFolder,"/Training_dataset_NDSI_regression_subset",sub.size,".csv"), header = TRUE, sep = ",")
snowData <- as.data.frame(df)
snowData <- na.omit(snowData)
rm(df)
gc()

# Growing and saving the regression random forest ----
setwd(rfFolder)
#snowData <- read.csv(paste0("Training_dataset_S2.NDSI_modis.FSC_regression_subset",sub.size,".csv"), sep = ",")
snowData <- na.omit(snowData)

n.variables <- length(columns)-1 # s2

print(Sys.time())
ranger.s2.regression <- ranger(SNOW ~ ., 
                            data = snowData,
                            num.trees = 100,   # 100
                            mtry = round((n.variables/3),digits = 0),
                            #min.node.size = 5,   #1
                            #importance = 'impurity_corrected',
                            splitrule = "variance",
                            save.memory = TRUE,
                            verbose = TRUE,
                            num.threads = numCores)

print(Sys.time())

print(object.size(ranger.s2.regression),units="GB")
ranger.s2.regression$prediction.error

print(Sys.time())
save(ranger.s2.regression,file=paste0(rfFolder,"/ranger.s2.ndsi.MODIS.FSC.regression_sub.size",as.character(sub.size),".Rda"))
print(Sys.time())
#rm(ranger.regression)
#load(paste0(rfFolder,"/ranger.regression_sub.size",as.character(sub.size),".Rda"))


#variable importance
feat.imp <- importance(ranger.s2.regression)
names.imp <- names(feat.imp)
feat.imp <- sort(feat.imp)

png(paste0(rfFolder,"/importance_S2.ndsi_modis_FSC_regression.png"),width=2000,height=800)
par(mfrow=c(1,2),mai=c(1,4,1,1))
barplot(feat.imp/1e6,xlab='Importance (scaled)',las=1,cex.lab=2.5, cex.axis=2, cex.names=2.5,horiz=TRUE)
dev.off()


# Build the dataset for the regression ----

sr <- raster("D:/ABRESO/PNGP/Sentinel2/T32TLR/Level2/2015/output_T32TLR/Snow&Cloud_mask_20m_20150730_S2A_L2A_T32TLR_R108_TA103016_TP103016.tif")
setwd(modisFolder)
m <- 1
t <- FALSE

# Find the first cloud free day to extract the first water mask ----
for (y in yy) {
  L2_yy_out <- file.path(paste0(L2folder,"/",y,"/output_",tl))
  setwd(L2_yy_out)
  s2.names <- list.files(pattern=glob2rx("Snow&Cloud_mask_20m_*.tif"))
  n.scenes <- length(s2.names)
  
  for (s in s2.names) {
    s2.mask <- raster(s)
    date <- substr(s,21,28)
    total_cloud_cells <- sum(s2.mask[!is.na(s2.mask)] == 205)
    total_cells <- ncell(s2.mask[s2.mask != 254])
    total_cloud_fraction <- (total_cloud_cells*100)/total_cells
    total_NA_cells <- sum(s2.mask[!is.na(s2.mask)] == 254)
    total_NA_fraction <- (total_NA_cells*100)/ncell(s2.mask)
    
    if (total_cloud_fraction < ct & total_NA_fraction == 0) {
      t <- TRUE
      B03 <- raster::raster(list.files(pattern=glob2rx(paste0("B3_10m_",as.character(date),"*.tif"))))
      B08 <- raster::raster(list.files(pattern=glob2rx(paste0("B8_10m_",as.character(date),"*.tif"))))
      ## Water bodies extraction ##
      NDWI <- (B03 - B08) / (B03 + B08)
      w1 <- 0.5 # water threshold
      waterbodies <- calc(NDWI, fun = function(x) {return(x >= w1)})
      # check resolution
      rres <- xres(waterbodies)/xres(sr)
      if (rres != 1){
        waterbodies <- raster::resample(waterbodies, sr, method="bilinear")
      }
      # check extent
      if (extent(waterbodies) != extent(sr)){
        waterbodies <- crop(waterbodies, sr)
      }
      if (extent(waterbodies) != extent(sr)){
        waterbodies <- setExtent(waterbodies, extent(sr), keepres=FALSE, snap=FALSE)
      }
      ## end ##
      print(paste0("date ", as.character(date), " used for water mask"))
      break
      
    } else {
      print(paste0("date ", as.character(date), " cloudy or NA - loop progress ", m,"/",n.scenes))
      m <- m + 1
    }
  }
  
  if (t == TRUE){
    break
  }
}

#  Main loop ----

load(paste0(rfFolder,"/ranger.s2.ndsi.MODIS.FSC.regression_sub.size",as.character(sub.size),".Rda"))

# get list of DOYs from 2015 to 2021
all.days <- as.Date(as.Date("2015-07-30"):as.Date("2021-12-31"), origin="1970-01-01")
n.days <- length(all.days)

# from every MODIS build a dataset for the gap filling
columns <- c("SNOW", "Elevation", "Slope", "Aspect","MODIS","X","Y","DOY","Year")
df.SCD <- matrix(ncol = length(columns)) 
colnames(df.SCD) = columns
SCD.stack <- stack()

# Begin the main loop
n <- 1

for (i in all.days[n:n.days]) {
  
  d <- all.days[n]
  doy <- strftime(d, format = "%j")    #convert date in doy
  y <- strftime(d, format = "%Y")
  setwd(modisGF)
  if (nchar(as.character(doy)) == 3){
    modis.name <- list.files(pattern=glob2rx(paste0("MOD10A1_GF_A",y,doy,"_*_NDSI_SnowCover.tif")),recursive = TRUE)
  } else if (nchar(as.character(doy)) == 2){
    modis.name <- list.files(pattern=glob2rx(paste0("MOD10A1_GF_A",y,"0",doy,"_*_NDSI_SnowCover.tif")),recursive = TRUE)
  } else if (nchar(as.character(doy)) == 1){
    modis.name <- list.files(pattern=glob2rx(paste0("MOD10A1_GF_A",y,"00",doy,"_*_NDSI_SnowCover.tif")),recursive = TRUE)
  }
  
  modis <- raster(modis.name[1])    #get the MODIS scene
  # check CRS
  if (compareCRS(sr, modis) == FALSE){
    modis <- projectRaster(modis, crs = crs(sr), res=xres(sr))
  }
  # check resolution
  rres <- xres(modis)/xres(sr)
  if (rres != 1){
    modis <- raster::resample(modis, sr, method="ngb")
  }
  # check extent
  if (extent(modis) != extent(sr)){
    modis <- crop(modis, sr)
  }
  if (extent(modis) != extent(sr)){
    modis <- setExtent(modis, extent(sr), keepres=FALSE, snap=FALSE)
  }
  
  #Find if there is a sentinel2 scene and if not create an empty raster
  y <- format(d, format = "%Y")
  m <- format(d, format = "%m")
  g <- format(d, format = "%d")
  L2_yy_out <- file.path(paste0(L2folder,"/",y,"/output_",tl))
  setwd(L2_yy_out)
  s2.date <- paste0(y,m,g)
  s2.name <- list.files(pattern=glob2rx(paste0("Snow&Cloud_mask_20m_",s2.date,"*.tif")))
  
  if (length(s2.name)>0){
    s2.mask <- raster(s2.name)
    total_cloud_cells <- sum(s2.mask[!is.na(s2.mask)] == 205)
    total_cells <- ncell(s2.mask[s2.mask != 254])
    total_cloud_fraction <- (total_cloud_cells*100)/total_cells
    #total_NA_cells <- sum(s2.mask[!is.na(s2.mask)] == 254)
    #total_NA_fraction <- (total_NA_cells*100)/ncell(s2.mask)
    
    s2.mask[s2.mask == 205] <- NA #Set cloud as NA
    s2.mask[s2.mask == 254] <- NA #Set No data as NA
    total_NA_cells <- ncell(s2.mask[is.na(s2.mask)])
    
    if (total_cloud_fraction < ct & (total_NA_cells > 0)) {
      s2 <- raster::raster(list.files(pattern=glob2rx(paste0("NDSI_20m_",y,m,g,"*.tif"))))
        
    } else if (total_cloud_fraction >= ct) {
      s2 <- raster(ncol=ncol(sr), nrow=nrow(sr), ext=extent(sr), crs = crs(sr))  #create an empty raster
      values(s2) <- NA    # assign all value as NA
      
    } else if (total_cloud_fraction == 0 & total_NA_fraction == 0) {
      r <- s2
      writeRaster(r, paste0("A2_NDSI_20m_",doy,"_",y,".tif"), datatype='INT2S', format = "GTiff", overwrite=TRUE, option="COMPRESS=LZW")
      snow2 <- s2.mask
      snow2[snow2 == 100] <- 1
      writeRaster(snow2, paste0("A2_Snow_mask_regression_20m_",doy,"_",y,".tif"), datatype='INT2S', format = "GTiff", overwrite=TRUE, option="COMPRESS=LZW")
      
      print(paste0("doy ", doy, " of year ", y, " ok - loop progress ", n,"/",n.days))
      n <- n+1  
      next
    }
    
  } else {
    s2 <- raster(ncol=ncol(sr), nrow=nrow(sr), ext=extent(sr), crs = crs(sr))  #create an empty raster
    values(s2) <- NA    # assign all value as NA
    
  }
  
  SCD.stack <- stack(s2,dem,slope,aspect,modis)    # lulc   create stack for each snow cover scene
  SCD.data <- as.matrix(SCD.stack) %>%     # convert in matrix
    cbind(locs.x,locs.y,doy,y)  
  
  df.SCD <- rbind(df.SCD,SCD.data)
  df.SCD <- df.SCD[-1,]    #Delete first empty row
  df.SCD <- apply(df.SCD, 2, as.integer) %>% # convert all values in integer to decrease size
    as.data.frame
  colnames(df.SCD) = columns
  regressionvalues <- df.SCD[is.na(df.SCD$SNOW), ]
  
  ## 
  ## Predict snow data
  ##
  
  ranger.temp <- predict(ranger.s2.regression ,
                         data = regressionvalues,
                         predict.all = FALSE,
                         num.trees = ranger.s2.regression$num.trees,
                         type = "response",   #terminalNodes
                         se.method = "infjack",
                         seed = 567,
                         num.threads = numCores,
                         verbose = TRUE,
                         inbag.counts = NULL)
  
  
  # Add retrieved values to the s2 raster
  df.SCD$SNOW[is.na(df.SCD$SNOW)] <- ranger.temp$predictions
  
  # from dataframe to raster desired raster
  r <- raster(nrow=nrow(s2), ncol=ncol(s2), ext=extent(s2), crs=crs(s2))
  values(r) <- as.numeric(df.SCD$SNOW)
  
  # Save the gap filled NDSI raster
  writeRaster(r, paste0("A2_NDSI_20m_",doy,"_",y,".tif"), datatype='INT2S', format = "GTiff", overwrite=TRUE, option="COMPRESS=LZW")
  
  ### SNOWLINE ELEVATION ###
  
  bh <- 100    # elevation band
  DEMmin <- round(minValue(dem),-2) - bh
  zmax <- DEMmax <- round(maxValue(dem),-2) + bh
  zmin <- zmax - bh
  p1 <- 4000   # NDSI threshold 1
  p2 <- 2000   # NDSI threshold 2  # 1500
  p3 <- 0.1    # snow pixels limit (0.1%)
  
  #calcola estensione banda altitudinale libera da nuvole
  snow1 <- r
  snow1 <- calc(snow1, fun = function(x) {return(x >= p1)})
  
  dem_band <- calc(dem, fun = function(x) {return(x > zmin & x <= zmax)})
  dem_band_pixel <- sum(dem_band[!is.na(dem_band)])
  
  snow_band <- overlay(snow1,
                       dem_band,
                       fun=function(x,y) {return(x == 1 & y ==1)})
  
  snowy_band_pixel <- sum(snow_band[!is.na(snow_band)]==1)
  snow_pixels <- (snowy_band_pixel*100)/dem_band_pixel
  
  while (snow_pixels > p3 & zmin >= DEMmin) {
    
    dem_band <- calc(dem, fun = function(x) {return(x > zmin & x <= zmax)})
    dem_band_pixel <- sum(dem_band[!is.na(dem_band)])
    
    snow_band <- overlay(snow1,
                         dem_band,
                         fun=function(x,y) {return(x == 1 & y ==1)})
    
    snowy_band_pixel <- sum(snow_band[!is.na(snow_band)]==1)
    snow_pixels <- (snowy_band_pixel*100)/dem_band_pixel
    
    zmax <- zmax - bh
    zmin <- zmin - bh
    
  }
  
  #definisci limite inferiore neve come due bande inferiori a quella precedentemente trovata
  zs <- zmin - (2*bh)
  
  snow2 <- overlay(r,
                   dem,
                   waterbodies,
                   fun=function(x,y,z) {return(x >= p2 & y >= zs & z != 1)})
  
  # Save the gap filled raster
  writeRaster(snow2, paste0("A2_Snow_mask_regression_20m_",doy,"_",y,".tif"), datatype='INT2S', format = "GTiff", overwrite=TRUE, option="COMPRESS=LZW")
  
  #### END Snowline elevation ###
  
  # Empty dataframe for the next scene
  df.SCD <- matrix(ncol = length(columns)) 
  colnames(df.SCD) = columns
  SCD.stack <- stack()
  
  print(paste0("doy ", doy, " of year ", y, " ok - loop progress ", n,"/",n.days))
  n <- n+1  
}


