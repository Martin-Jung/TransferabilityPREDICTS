library(gdalUtils)
if( is.null( getOption("gdalUtils_gdalPath") ) ){
  gdal_setInstallation("C:/Program Files/GDAL/",rescan = T)
}
library(rgdal)
library(raster)
library(zoo)
library(stringr)
library(RStoolbox)
library(foreach);library(doParallel)
source("000_HelperFunction.R")
source("../../Scripts/Icarus/R/miscellaneous/latlong2UTMzone.R")
rasterOptions(progress = "text")
extractPath = "C:/LocalStorage/LocalFolder/LandsatPREDICTS_AnnualStack" # Where are the extracts per study stored ?
extractPath.Predictors = "C:/LocalStorage/LocalFolder/LandsatPREDICTSclassification_Predictors" # Where are the extracted Predictors per study saved ?
outPath = "Z:/Private/mj291/predicts_stacks/";if(!dir.exists(outPath)) dir.create(outPath)
outPath.predictors = "Z:/Private/mj291/predicts_stacksPredictors/";if(!dir.exists(outPath.predictors)) dir.create(outPath.predictors)

tempdir = "C:/temp"

buffersize = 1000 # THe buffer size used
overw = FALSE # Should new files be created if already existing?
smoothGapF = TRUE # Should gaps as well as outliers be smoothed

# ---------------------------- #
# Load 2006 dataset and retain only those with a latitude coordinate
predicts <- readRDS("../../Data/PREDICTS_v1/sites.rds") %>% dplyr::select(SS,SSBS,Longitude,Latitude,
                                                                          Sample_midpoint,Sample_start_earliest,Sample_end_latest,
                                                                          Predominant_land_use,Habitat_as_described,Use_intensity)
predicts$midyear <- lubridate::year(lubridate::ymd(predicts$Sample_midpoint))
predicts$startyear <- lubridate::year(lubridate::ymd(predicts$Sample_start_earliest))
predicts <- predicts[which(!is.na(predicts$Latitude)),]

#### LS extraction -- Load and subset a studies raster layer per site #### 
# Load each study raster layer
# Crop to buffered site per study
bands = c('blue','green','red','nir','swir1','swir2','lst')
years = seq(1982,2015)
nn <- paste0(
  bands,"_",as.vector(sapply(years,simplify = T, function(x) rep(x,7 )))
)

# Load all files
ll <- list.files(extractPath,pattern = "*.tif",full.names = T,ignore.case = T)
ll <- ll[grep("*aux.xml",ll,invert = T)] # Remove aux.xml jic
ll <- ll[grep("SiteComposite",ll,invert = T)] # Remove aux.xml jic
# Get the study names from the files
ll.SS <- str_split(str_split(basename(ll),"-",simplify = T)[,1],"Composite_",simplify = T)[,2]
ll.SS <- str_remove_all(ll.SS,"\\.tif")

registerDoParallel(cores = parallel::detectCores()-2)

# Now loop through each study id and build the raster
for(fname in unique(ll.SS) ){
  myLog("Processing study: ",fname)
  # Get Study name from file name
  sub.study <- subset(predicts,SS ==  fname)
  if(nrow(sub.study)==0) {print("Study not in public PREDICTS dataset");next() }
  
  # Get all files with that study name
  ll.ex <- ll[which(ll.SS == fname)]
  
  if(length(ll.ex)==1){
    # Build the stack
    ss <- stack(ll.ex)
    ss.path = ll.ex
  } else {
    myLog("--- Multiple files. Mosaic stack (VRT)...")
    o = gdalbuildvrt(ll.ex,output.vrt = paste0(tempdir,"/",fname,".vrt"),verbose = F,overwrite = T)
    if(!is.null(o)) stop("Mosaicing did not suceed...")
    ss <- stack(paste0(tempdir,"/",fname,".vrt"))
    ss.path = paste0(tempdir,"/",fname,".vrt")
  }
  # Correct file names
  stopifnot(nlayers(ss) == length(nn))
  names(ss) <- nn
  ss <- setZ(ss,z = str_split(nn,"_",simplify = T)[,2] ) 
  
  # Now for each study site
  for( site in unique(sub.study$SSBS) ){
    # Correct special symbols in SSBS name
    site <- sub("/","_",site)
    if( file.exists(paste0(outPath,site,".nc")) & !overw) {next()}
    
    myLog("--> ",site)
    sub.site <- subset(sub.study,SSBS == site)
    # Get UTM zone projection
    zone = CRS(paste0("+proj=utm +zone=",latlong2UTMzone(lon = sub.site$Longitude,lat = sub.site$Latitude )," +datum=WGS84 +units=m +no_defs"))
    # Make spatial file
    coordinates(sub.site) <- ~Longitude+Latitude
    proj4string(sub.site) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    # Transform to a metric meter based projection
    sub.site <- spTransform(sub.site,CRSobj = zone )
    buf <- rgeos::gBuffer(sub.site,quadsegs = 50,width = buffersize)
    buf <- spTransform(buf, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") )# Transform back
    suppressWarnings(
      writeOGR(SpatialPolygonsDataFrame(Sr = buf,data = sub.site@data,match.ID = F),
               dsn = tempdir,layer = "site",driver = "ESRI Shapefile",overwrite_layer = TRUE)
    )
    # Now crop the layer stack
    o <- gdalwarp(
      srcfile = ss.path,dstfile = paste0(outPath,site,".tif"),
      co = "COMPRESS = LZW",output_Raster=TRUE,multi=TRUE,
      q= T, cutline = paste0(tempdir,"/site.shp"),crop_to_cutline = TRUE,overwrite = TRUE
    )
    if(is.null(o)) {stop("Cropped to cutline of site did not work!") }
    
    o[o == 0] <- NA # Overwrite NA
    # Now rename and format and save to spatial-temporal frame
    names(o) <- nn
    o <- setZ(o,z = str_split(nn,"_",simplify = T)[,2] ) 
    
    if(smoothGapF){
      require(forecast)
      newbrick <- list() # Ugly hack to get the order back
      # Temporally gapfilling 
      for(b in bands){
        print(b)
        sub <- o[[which(str_detect(names(o),b))]]
        if(nrow(sub) == 0) stop("No bands found")
        
        a1 = Sys.time()
        oo <- foreach(i = 1:nrow(sub),
                      .combine='rbind', 
                      .multicombine =TRUE,
                      .errorhandling = 'stop',
                      .packages=c('raster','zoo','lubridate','forecast','stringr'),
                      .verbose = FALSE) %dopar% {
                        # Get all values for a specific row
                        val <- raster::getValues(sub,i)
                        # Now apply per row and time series
                        return( 
                          base::t.default( apply(val, 1, function(y) preparePixel(y,years)  ) )
                        )
                      }
        print( Sys.time()-a1 )
        sub[] <-oo
        
        # Interpolate single missing years
        sub <- approxNA(sub,NArule = NA)
        sub <- setZ(sub,z = years ) 
        
        newbrick[[b]] <- sub
      }
      rm(o) # Delete the original
      # Reconstruct the stack per year
      for( y in sort(years)){
        if(!exists('o')) {
          o <- brick(
            newbrick[[bands[1]]][[paste0(bands[1],"_",y)]],
            newbrick[[bands[2]]][[paste0(bands[2],"_",y)]],
            newbrick[[bands[3]]][[paste0(bands[3],"_",y)]],
            newbrick[[bands[4]]][[paste0(bands[4],"_",y)]],
            newbrick[[bands[5]]][[paste0(bands[5],"_",y)]],
            newbrick[[bands[6]]][[paste0(bands[6],"_",y)]],
            newbrick[[bands[7]]][[paste0(bands[7],"_",y)]]
          )
        } else {
          o <- addLayer(o,newbrick[[bands[1]]][[paste0(bands[1],"_",y)]])
          o <- addLayer(o,newbrick[[bands[2]]][[paste0(bands[2],"_",y)]])
          o <- addLayer(o,newbrick[[bands[3]]][[paste0(bands[3],"_",y)]])
          o <- addLayer(o,newbrick[[bands[4]]][[paste0(bands[4],"_",y)]])
          o <- addLayer(o,newbrick[[bands[5]]][[paste0(bands[5],"_",y)]])
          o <- addLayer(o,newbrick[[bands[6]]][[paste0(bands[6],"_",y)]])
          o <- addLayer(o,newbrick[[bands[7]]][[paste0(bands[7],"_",y)]])
        }
      }
      names(o) <- nn
      o <- setZ(o,z = str_split(nn,"_",simplify = T)[,2] ) 
    }
    # Use only the year of sampling
    o[[grep(sub.site$midyear,names(o))]]    
    writeRaster(x = o[[grep(sub.site$midyear,names(o))]],
                paste0(outPath,site,".nc"),
                format = "CDF",NAflag = -9999,
                varname = "Surface reflectance",varunit = "bandwith",
                xname = "Longitude", yname = "Latitude", zname = "Time", zunit = "year",
                force_v4=TRUE, compression=7,
                overwrite= TRUE
    )
    # Clean up
    file.remove(paste0(outPath,site,".tif"))
    try(lapply(list.files(tempdir,pattern = "site",full.names = T), file.remove),silent=TRUE)
  }
}

# ----------Site composites---------- #
# Repeat for the SiteComposites only
ll <- list.files(extractPath,pattern = "*.tif",full.names = T,ignore.case = T)
ll <- ll[grep("*aux.xml",ll,invert = T)] # Remove aux.xml 
ll <- ll[grep("SiteComposite",ll,invert = F)] # Site comps
# Get the study names from the files
ll.SSBS <- str_split(str_split(basename(ll),"-",simplify = T)[,1],"SiteComposite_",simplify = T)[,2]
ll.SSBS <- str_remove_all(ll.SSBS,"\\.tif")

# Now loop through each study id and build the raster
for(site in unique(ll.SSBS) ){
  if( file.exists(paste0(outPath,site,".nc")) & !overw) {print("Already existing...");next()}
  myLog("--> ",site)
  sub.site <- subset(predicts,SSBS == site)
  if(nrow(sub.site)==0) {print("Site not in public PREDICTS dataset");next() }
  
  # Get all files with that study name
  ll.ex <- ll[which(ll.SSBS == site)]
  ss.path = ll.ex
  # Get UTM zone projection
  zone = CRS(paste0("+proj=utm +zone=",latlong2UTMzone(lon = sub.site$Longitude,lat = sub.site$Latitude )," +datum=WGS84 +units=m +no_defs"))
  # Make spatial file
  coordinates(sub.site) <- ~Longitude+Latitude
  proj4string(sub.site) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  # Transform to a metric meter based projection
  sub.site <- spTransform(sub.site,CRSobj = zone )
  buf <- rgeos::gBuffer(sub.site,quadsegs = 50,width = buffersize)
  buf <- spTransform(buf, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") )# Transform back
  suppressWarnings(
    writeOGR(SpatialPolygonsDataFrame(Sr = buf,data = sub.site@data,match.ID = F),
             dsn = tempdir,layer = "site",driver = "ESRI Shapefile",overwrite_layer = TRUE)
  )
  # Now crop the layer stack
  o <- gdalwarp(
    srcfile = ss.path,dstfile = paste0(outPath,site,".tif"),
    co = "COMPRESS = LZW",output_Raster=TRUE,multi=TRUE,
    q= T, cutline = paste0(tempdir,"/site.shp"),crop_to_cutline = TRUE,overwrite = TRUE
  )
  if(is.null(o)) {stop("Cropped to cutline of site did not work!") }
  
  o[o == 0] <- NA # Overwrite NA
  # Now rename and format and save to spatial-temporal frame
  names(o) <- nn
  o <- setZ(o,z = str_split(nn,"_",simplify = T)[,2] ) 
  
  if(smoothGapF){
    require(forecast)
    newbrick <- list() # Ugly hack to get the order back
    # Temporally gapfilling 
    for(b in bands){
      print(b)
      sub <- o[[which(str_detect(names(o),b))]]
      
      a1 = Sys.time()
      oo <- foreach(i = 1:nrow(sub),
                    .combine='rbind', 
                    .multicombine =TRUE,
                    .errorhandling = 'stop',
                    .packages=c('raster','zoo','lubridate','forecast','stringr'),
                    .verbose = FALSE) %dopar% {
                      # Get all values for a specific row
                      val <- raster::getValues(sub,i)
                      # Now apply per row and time series
                      return( 
                        base::t.default( apply(val, 1, function(y) preparePixel(y,years)  ) )
                      )
                    }
      print( Sys.time()-a1 )
      sub[] <-oo
      
      # Interpolate single missing years
      sub <- approxNA(sub,NArule = NA)
      sub <- setZ(sub,z = years ) 
      
      newbrick[[b]] <- sub
    }
    rm(o) # Delete the original
    # Reconstruct the stack per year
    for( y in sort(years)){
      if(!exists('o')) {
        o <- brick(
          newbrick[[bands[1]]][[paste0(bands[1],"_",y)]],
          newbrick[[bands[2]]][[paste0(bands[2],"_",y)]],
          newbrick[[bands[3]]][[paste0(bands[3],"_",y)]],
          newbrick[[bands[4]]][[paste0(bands[4],"_",y)]],
          newbrick[[bands[5]]][[paste0(bands[5],"_",y)]],
          newbrick[[bands[6]]][[paste0(bands[6],"_",y)]],
          newbrick[[bands[7]]][[paste0(bands[7],"_",y)]]
        )
      } else {
        o <- addLayer(o,newbrick[[bands[1]]][[paste0(bands[1],"_",y)]])
        o <- addLayer(o,newbrick[[bands[2]]][[paste0(bands[2],"_",y)]])
        o <- addLayer(o,newbrick[[bands[3]]][[paste0(bands[3],"_",y)]])
        o <- addLayer(o,newbrick[[bands[4]]][[paste0(bands[4],"_",y)]])
        o <- addLayer(o,newbrick[[bands[5]]][[paste0(bands[5],"_",y)]])
        o <- addLayer(o,newbrick[[bands[6]]][[paste0(bands[6],"_",y)]])
        o <- addLayer(o,newbrick[[bands[7]]][[paste0(bands[7],"_",y)]])
      }
    }
    names(o) <- nn
    o <- setZ(o,z = str_split(nn,"_",simplify = T)[,2] ) 
  }
  writeRaster(x = o,
              paste0(outPath,site,".nc"),
              format = "CDF",NAflag = -9999,
              varname = "Surface reflectance",varunit = "bandwith",
              xname = "Longitude", yname = "Latitude", zname = "Time", zunit = "year",
              force_v4=TRUE, compression=7,
              overwrite= TRUE
  )
  # Clean up
  file.remove(paste0(outPath,site,".tif"))
  try(lapply(list.files(tempdir,pattern = "site",full.names = T), file.remove),silent=TRUE)
  
}

doParallel::stopImplicitCluster()
stop("DONE!")

# ------------------- #
#### Clip and prepare the predictors ####
# Do the same for the additional predictors precip, 

pred.names <- c('precip_median','precip_range','elevation','slope')

nn = c(paste0(pred.names[1:2],rep(years,each=2)), c(pred.names[3],pred.names[4]) )

# Load all files
ll <- list.files(extractPath.Predictors,pattern = "*.tif",full.names = T,ignore.case = T)
ll <- ll[grep("*aux.xml",ll,invert = T)] # Remove aux.xml jic
# Get the study names from the files
ll.SS <- str_split(str_split(basename(ll),"-",simplify = T)[,1],"CompositePred_",simplify = T)[,2]
ll.SS <- str_remove_all(ll.SS,"\\.tif")

registerDoParallel(cores = parallel::detectCores()-2)
# Now loop through each study id and build the raster
for(fname in unique(ll.SS) ){
  myLog("Processing study: ",fname)
  # Get Study name from file name
  sub.study <- subset(predicts,SS ==  fname)
  if(nrow(sub.study)==0) {print("Study not in public PREDICTS dataset");next() }
  
  # Get all files with that study name
  ll.ex <- ll[which(ll.SS == fname)]
  
  if(length(ll.ex)==1){
    # Build the stack
    ss <- stack(ll.ex)
    ss.path = ll.ex
  } else {
    myLog("--- Multiple files. Mosaic stack (VRT)...")
    o = gdalbuildvrt(ll.ex,output.vrt = paste0(tempdir,"/",fname,".vrt"),verbose = F,overwrite = T)
    if(!is.null(o)) stop("Mosaicing did not suceed...")
    ss <- stack(paste0(tempdir,"/",fname,".vrt"))
    ss.path = paste0(tempdir,"/",fname,".vrt")
  }
  # Correct file names
  stopifnot(nlayers(ss) == length(nn))
  names(ss) <- nn
  
  # Now for each study
  for( site in unique(sub.study$SSBS) ){
    if( file.exists(paste0(outPath.predictors,site,".nc")) & overw) {next()}
    myLog("--> ",site)
    sub.site <- subset(sub.study,SSBS == site)
    # Get UTM zone projection
    zone = CRS(paste0("+proj=utm +zone=",latlong2UTMzone(lon = sub.site$Longitude,lat = sub.site$Latitude )," +datum=WGS84 +units=m +no_defs"))
    # Make spatial file
    coordinates(sub.site) <- ~Longitude+Latitude
    proj4string(sub.site) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    # Transform to a metric meter based projection
    sub.site <- spTransform(sub.site,CRSobj = zone )
    buf <- rgeos::gBuffer(sub.site,quadsegs = 50,width = buffersize)
    buf <- spTransform(buf, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") )# Transform back
    suppressWarnings(
      writeOGR(SpatialPolygonsDataFrame(Sr = buf,data = sub.site@data,match.ID = F),
               dsn = tempdir,layer = "site",driver = "ESRI Shapefile",overwrite_layer = TRUE)
    )
    # Now crop the layer stack
    o <- gdalwarp(
      srcfile = ss.path,dstfile = paste0(outPath.predictors,site,".tif"),
      co = "COMPRESS = LZW",output_Raster=TRUE,multi=TRUE,
      q= T, cutline = paste0(tempdir,"/site.shp"),crop_to_cutline = TRUE,overwrite = TRUE
    )
    if(is.null(o)) {stop("Cropped to cutline of site did not work!") }
    # Now rename and format and save to spatial-temporal frame
    names(o) <- nn
    
    # For all precip layers assign the most common value
    # Neccessary because of the coarse resolution
    o.mask <- o$elevation
    o.mask[o.mask == 0] <- NA
    for(l in grep('precip',nn) ){
      o[[l]][] <- cellStats(o[[l]],Mode)
      o[[l]] <- raster::mask(o[[l]],o.mask,updatevalue=NA)
    }
    
    writeRaster(x = o,
                paste0(outPath.predictors,site,".nc"),
                format = "CDF",NAflag = -9999,
                varname = "Predictors for stack",varunit = "mm and m",
                xname = "Longitude", yname = "Latitude", zname = "Time", zunit = "year",
                force_v4=TRUE, compression=7,
                overwrite= TRUE
    )
    # Clean up
    file.remove(paste0(outPath.predictors,site,".tif"))
    try(lapply(list.files(tempdir,pattern = "site",full.names = T), file.remove),silent=TRUE)
  }
}
doParallel::stopImplicitCluster()
stop("DONE!")







#### Testing ####
sub <- stack("resSaves/AD1_2006__Diekoetter 1  90.tif")
# Correct file names
stopifnot(nlayers(ss) == length(nn))
names(ss) <- nn
ss <- setZ(ss,z = str_split(nn,"_",simplify = T)[,2] ) 


# Write the result
writeRaster(ss,"out",overwrite=TRUE)
