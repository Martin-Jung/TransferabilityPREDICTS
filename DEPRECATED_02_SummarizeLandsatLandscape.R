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
outPath = "Z:/Private/mj291/predicts_landscapesummary/";if(!dir.exists(outPath)) dir.create(outPath)

tempdir = "C:/temp"

buffersize = 1000 # THe buffer size used
overw = FALSE # Should new files be created if already existing?
smoothGapF = FALSE # Should gaps as well as outliers be smoothed

# ---------------------------- #
# Load 2006 dataset and retain only those with a latitude coordinate
predicts <- readRDS("../../Data/PREDICTS_v1/sites.rds") %>% dplyr::select(SS,SSBS,Longitude,Latitude,
                                                                          Sample_midpoint,Sample_start_earliest,Sample_end_latest,
                                                                          Predominant_land_use,Habitat_as_described,Use_intensity)
predicts$midyear <- lubridate::year(lubridate::ymd(predicts$Sample_midpoint))
predicts$startyear <- lubridate::year(lubridate::ymd(predicts$Sample_start_earliest))
predicts <- predicts[which(!is.na(predicts$Latitude)),]


# ---------------------------- #
# Coefficient of variation
cv <- function(x) { ( sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)  )  }

# Requires a raster layer
# Focal metrics
# Proportion of missing data | Mean and cv of NDVI |
# Mean spectral values per site
# Shannon diversity index and RAO's Q
# Multidimensional 
summarize.landscape <- function(o){
  require(glcm);require(vegan)
  # Source the spectral rao's Q metric
  #source("https://raw.githubusercontent.com/mattmar/spectralrao/master/spectralrao.r")
  
  # Calculate Vegetation indices
  sv <- RStoolbox::spectralIndices(o,
                                   blue = grep("blue",names(o)),green = grep("green",names(o)),red = grep("red",names(o)),
                                   nir = grep("nir",names(o)),swir2 = grep("swir1",names(o)),swir3 = grep("swir2",names(o)),
                                   indices = c("NDVI","EVI2","NDWI2")
                                   )
  # raster PCA over all bands
  try( pc <- RStoolbox::rasterPCA(o), silent = T )
  
  # Calculate graylevel texture matrix texture
  tex <- glcm(sv$NDVI,window = c(7,7),  shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),na_opt = "any")
  # non-linear exponential relationship between dissimi and entropy. To catch large values use those
  
  # --- #
  # Calculate stats
  band.means <- data.frame( t(cellStats(o,"mean",na.rm = T)))
  names(band.means) <- paste0(str_to_upper(str_split(names(band.means),"_",simplify = T)[,1]), "_mean")
  band.cv <- data.frame( t( apply(values(o), 2, cv) ))
  names(band.cv) <- paste0(str_to_upper(str_split(names(band.cv),"_",simplify = T)[,1]), "_cv")
  
  # Add vegetation indices
  out <- data.frame(
    band.means,band.cv,
    missing = cellStats(sv$NDVI,"countNA"),
    NDVI_mean = cellStats(sv$NDVI,stat = "mean",na.rm = T),
    NDVI_cv = cv(values(sv$NDVI)),
    EVI2_mean = cellStats(sv$EVI2,stat = "mean",na.rm = T),
    EVI2_cv = cv(values(sv$EVI2)),
    NDWI_mean = cellStats(sv$NDWI2,stat = "mean",na.rm = T),
    NDWI_cv = cv(values(sv$NDWI2))
    
  )
  rm(band.means,band.cv)
  
  # Add PC results
  if(exists("pc")){
    pc.means <- data.frame( t( cellStats( pc$map[[c("PC1","PC2")]],"mean") ) )
    names(pc.means) <-  paste0(str_to_upper( names(pc.means)[grep("PC",names(pc.means))] ), "_mean")
    pc.cv <- data.frame( t( apply(values(pc$map[[c("PC1","PC2")]]), 2, cv) ))
    names(pc.cv) <- paste0(str_to_upper(str_split(names(pc.cv),"_",simplify = T)[,1]), "_cv")
    
    out <- data.frame(out, pc.means,pc.cv)
    rm(pc.means,pc.cv)
  } else {
    out$PC1_mean <- NA; out$PC2_mean <- NA; out$PC1_cv <- NA; out$PC2_cv <- NA; 
  }

  # Add texture results
  tex.mean <- data.frame( t( cellStats( tex[[c("glcm_entropy","glcm_dissimilarity")]] ,"mean") )  )
  names(tex.mean) <-  paste0( "NDVI_",names(tex.mean), "_mean")
  tex.cv <- data.frame( t( apply(values( tex[[c("glcm_entropy","glcm_dissimilarity")]]  ), 2, cv) ))
  names(tex.cv) <- paste0( "NDVI_",names(tex.cv), "_cv")
  
  out <- data.frame(out, tex.mean,tex.cv)
  rm(tex.mean,tex.cv)
  
  return(out)  
}

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
  if(file.exists( paste0(outPath,str_replace_all(fname," ","_"),".rds") )){
    print("Study already computed. Skip");next()
  } else {
    myLog("Processing study: ",fname)
  }
  # Get Study name from file name
  sub.study <- subset(predicts,SS ==  fname)
  if(nrow(sub.study)==0) {print("Study not in public PREDICTS dataset");next() }
  
  # Get all files with that study name
  ll.ex <- ll[which(ll.SS == fname)]
  if(length(ll.ex)==0) { print("No files for study found"); next() }
    
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
  
  study.output <- data.frame() # The output dataframe
  
  # Now for each study site
  for( site in as.character( unique(sub.study$SSBS) ) ){
    # Correct special symbols in SSBS name
    site <- sub("/","_",site)
    if( file.exists(paste0(outPath,site,".rds")) & !overw) {next()}
    
    myLog("--> ",site)
    sub.site <- subset(sub.study,SSBS == site)
    if(nrow(sub.site)==0) { next() }
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

    # Temporal grid cell forecasting. Takes considerably longer
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
    
    # Now subset to target year and summarize
    # Use only the year of sampling
    o <- o[[grep(sub.site$midyear,names(o))]]   
    o <- o * 0.0001# Correct unit
    o[o > 1] <- NA # Reflectance should be below and above this value. Thus set to NA
    
    # Calculate landscape statistics at 1000m
    # --- #
    # Check if bands have no data at all somewhere
    if( any(is.na(cellStats(o,stat = mean,na.rm=T)) ) ){
      # Dummy dataset'
      out1 <- data.frame(
      "BLUE_mean" = NA, "GREEN_mean" = NA, "RED_mean" = NA, "NIR_mean" = NA,                   
      "SWIR1_mean"  = NA, "SWIR2_mean" = NA, "LST_mean" = NA, "BLUE_cv" = NA,                    
      "GREEN_cv" = NA, "RED_cv" = NA, "NIR_cv" = NA,"SWIR1_cv" = NA,                 
      "SWIR2_cv" = NA, "LST_cv" = NA,"missing"  = 100000, "NDVI_mean" = NA,            
      "NDVI_cv"   = NA,"EVI2_mean" = NA,  "EVI2_cv" = NA,"NDWI_mean" = NA,                
      "NDWI_cv"  = NA, "PC1_mean" = NA, "PC2_mean" = NA, "PC1_cv"= NA,                   
      "PC2_cv"= NA,  "NDVI_glcm_entropy_mean"  = NA,"NDVI_glcm_dissimilarity_mean" = NA, "NDVI_glcm_entropy_cv" = NA,       
      "NDVI_glcm_dissimilarity_cv" = NA
    ) 
    
    } else {
    out1 <- summarize.landscape(o)
    }
    out1$scale <- 1000
    # --- #
    
    sub.site <- subset(sub.study,SSBS == site)
    # Get UTM zone projection
    zone = CRS(paste0("+proj=utm +zone=",latlong2UTMzone(lon = sub.site$Longitude,lat = sub.site$Latitude )," +datum=WGS84 +units=m +no_defs"))
    # Make spatial file
    coordinates(sub.site) <- ~Longitude+Latitude
    proj4string(sub.site) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    # Transform to a metric meter based projection
    sub.site <- spTransform(sub.site,CRSobj = zone )
    buf <- rgeos::gBuffer(sub.site,quadsegs = 50,width = 500)
    buf <- spTransform(buf, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") )# Transform back
    
    # At 500 m
    o <- raster::crop(o,buf);o <- raster::mask(o, buf )
    # --- #
    try(out2 <- summarize.landscape(o),silent = T)
    if(!exists("out2") ){
      out2 <- out1
      out2[names(out2)] <- NA
    }
    out2$scale <- 500
    # --- #
    out <- rbind(out1,out2)
    out$SS <- as.character( sub.site$SS )
    out$SSBS <- as.character( sub.site$SSBS )
    out$startyear <- sub.site$startyear
    out$midyear <- sub.site$midyear
    
    # Combine with previous site estimates
    study.output <- rbind( study.output, out )
    rm(out,out1,out2) # Clean up
    
    file.remove(paste0(outPath,site,".tif"))
    suppressMessages( try(lapply(list.files(tempdir,pattern = "site",full.names = T), file.remove),silent=TRUE) )
  }
  # Save study output
  saveRDS(study.output,file = paste0(outPath,str_replace_all(fname," ","_"),".rds") )
  rm(study.output,ss)
}

stopImplicitCluster();stop("Finished with real composites!")

#### ----------Site composites---------- ####
# Repeat for the SiteComposites only
ll <- list.files(extractPath,pattern = "*.tif",full.names = T,ignore.case = T)
ll <- ll[grep("*aux.xml",ll,invert = T)] # Remove aux.xml 
ll <- ll[grep("SiteComposite",ll,invert = F)] # Site comps
# Get the study names from the files
ll.SSBS <- str_split(str_split(basename(ll),"-",simplify = T)[,1],"SiteComposite_",simplify = T)[,2]
ll.SSBS <- str_remove_all(ll.SSBS,"\\.tif")

# Get study ids
studies <- predicts %>% dplyr::select(SS,SSBS) %>% distinct() %>% dplyr::filter(SSBS %in% unique(ll.SSBS))

registerDoParallel(cores = parallel::detectCores()-2)
print(paste0("Processing ",n_distinct(studies$SS)," studies per site"))

# Now loop through each study id and build the raster
for(fname in unique(studies$SS) ){
  if(file.exists( paste0(outPath,str_replace_all(fname," ","_"),".rds") )){
    print("Study already computed. Skip");next()
  } else {
    myLog("Processing study: ",fname)
  }
  sub.study <- subset(predicts,SS == fname)
  # Save study output  
  study.output <- data.frame() # The output dataframe

  for(  site in as.character( unique(sub.study$SSBS) )  ){
    # Get all files with that study name
    ll.ex <- ll[which(ll.SSBS == site)]
    if(length(ll.ex)==0) { next()  }
    ss.path = ll.ex
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
    
    # Now subset to target year and summarize
    # Use only the year of sampling
    o <- o[[grep(sub.site$midyear,names(o))]]   
    o <- o * 0.0001# Correct unit
    o[o > 1] <- NA # Reflectance should be below and above this value. Thus set to NA
    
    # Calculate landscape statistics at 1000m
    # --- #
    # Check if bands have no data at all somewhere
    if( any(is.na(cellStats(o,stat = mean,na.rm=T)) ) ){
      # Dummy dataset'
      out1 <- data.frame(
        "BLUE_mean" = NA, "GREEN_mean" = NA, "RED_mean" = NA, "NIR_mean" = NA,                   
        "SWIR1_mean"  = NA, "SWIR2_mean" = NA, "LST_mean" = NA, "BLUE_cv" = NA,                    
        "GREEN_cv" = NA, "RED_cv" = NA, "NIR_cv" = NA,"SWIR1_cv" = NA,                 
        "SWIR2_cv" = NA, "LST_cv" = NA,"missing"  = 100000, "NDVI_mean" = NA,            
        "NDVI_cv"   = NA,"EVI2_mean" = NA,  "EVI2_cv" = NA,"NDWI_mean" = NA,                
        "NDWI_cv"  = NA, "PC1_mean" = NA, "PC2_mean" = NA, "PC1_cv"= NA,                   
        "PC2_cv"= NA,  "NDVI_glcm_entropy_mean"  = NA,"NDVI_glcm_dissimilarity_mean" = NA, "NDVI_glcm_entropy_cv" = NA,       
        "NDVI_glcm_dissimilarity_cv" = NA
      ) 
      
    } else {
      out1 <- summarize.landscape(o)
    }
    out1$scale <- 1000
    # --- #
    
    sub.site <- subset(sub.study,SSBS == site)
    # Get UTM zone projection
    zone = CRS(paste0("+proj=utm +zone=",latlong2UTMzone(lon = sub.site$Longitude,lat = sub.site$Latitude )," +datum=WGS84 +units=m +no_defs"))
    # Make spatial file
    coordinates(sub.site) <- ~Longitude+Latitude
    proj4string(sub.site) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    # Transform to a metric meter based projection
    sub.site <- spTransform(sub.site,CRSobj = zone )
    buf <- rgeos::gBuffer(sub.site,quadsegs = 50,width = 500)
    buf <- spTransform(buf, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") )# Transform back
    
    # At 500 m
    o <- raster::crop(o,buf);o <- raster::mask(o, buf )
    # --- #
    try(out2 <- summarize.landscape(o),silent = T)
    if(!exists("out2") ){
      out2 <- out1
      out2[names(out2)] <- NA
    }
    out2$scale <- 500
    # --- #
    out <- rbind(out1,out2)
    out$SS <- as.character( sub.site$SS )
    out$SSBS <- as.character( sub.site$SSBS )
    out$startyear <- sub.site$startyear
    out$midyear <- sub.site$midyear
    
    # Combine with previous site estimates
    study.output <- rbind( study.output, out )
    rm(out,out1,out2) # Clean up
    
    file.remove(paste0(outPath,site,".tif"))
    suppressMessages( try(lapply(list.files(tempdir,pattern = "site",full.names = T), file.remove),silent=TRUE) )
  }
  # Save output
  saveRDS(study.output,file = paste0(outPath,str_replace_all(fname," ","_"),".rds") )
  
  # Clean up
  file.remove(paste0(outPath,site,".tif"))
  try(lapply(list.files(tempdir,pattern = "site",full.names = T), file.remove),silent=TRUE)
  
}

doParallel::stopImplicitCluster()
stop("DONE!")
