## Lake Okeechobee Algae Dynamics
## Created by: Paul Julian (pjulian@evergladesfoundation.org)
## Created on: 2024-06-26

## Data retrieval procedures

## Libraries
library(AnalystHelper)
library(EVERSpatDat)
library(reshape2)
library(plyr)
library(sf)
library(raster)
library(terra)

wd="C:/Julian_LaCie/_GitHub/LOK_AlgaeEval_manu"
paths=paste0(wd,c("/Plots/","/Data/","/src/","/_documents/"))

#Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
data.path=paths[2]


GIS.path.gen="C:/Julian_LaCie/_GISData"

# CRS
nad83.pro=st_crs("EPSG:4269")
NAD83_HARN <- st_crs("EPSG:2881")
utm17=st_crs("EPSG:26917")
wgs84=st_crs("EPSG:4326")

# Functions
source(paste0(wd,"/src/func.R"))

# Data --------------------------------------------------------------------
dates <- date.fun(c("1974-10-01","2023-09-30"))

dates2 <- dates # for the sonde data
dates2[1] <- date.fun("2016-10-01")

## Reg Schedules/Periods/Seasons -------------------------------------------
schs.all <- data.frame(
  name = c(
    "Alternate Regulation Schedule",
    "Interim Regulation Schedule",
    "Interim Regulation Schedule",
    "Run25",
    "WSE",
    "LORS08",
    "LOSOM"
  ),
  yr = c(1965, 1972, 1978, 1994, 1999, 2008, 2024)
)

schs.all.ts <- do.call(rbind, Map(function(nm, i) {
  data.frame(CY = schs.all$yr[i]:(schs.all$yr[i + 1] - 1), name = nm)
}, schs.all$name[-nrow(schs.all)], seq_len(nrow(schs.all) - 1)))

rownames(schs.all.ts) <- NULL

## Orignal LOSOM periods
CY.periods <- data.frame(
  CY = c(1973:1979, 1980:1999, 2000:2008, 2009:2023),
  Period = rep(c("73-79", "80-99", "00-08", "09-23"), 
               times = c(length(1973:1979), length(1980:1999), length(2000:2008), length(2009:2023)))
)
# Last period originally 2009 - 2019 "09-19"

## Hurricane years idenitified by during initial LOSOM modeling effort
CY.periods$walker.hurr <- as.integer(CY.periods$CY %in% c(2006:2008, 2018:2019, 2022))

# Define seasons in one step using vector recycling
month.season <- data.frame(
  month = 1:12,
  Season = rep(c("Winter", "Spring", "Summer", "Fall"), times = c(3, 1, 4, 4))
)

# Define LOSOM hydro seasons
hydro.sea <- data.frame(
  month = 1:12,
  LOSOM.hydrosea = rep(NA_character_, 12),
  FL.hydrosea = ifelse(1:12 %in% 5:10, "Wet", "Dry"),
  hydrosea2 = rep(NA_character_, 12)
)

# Assign LOSOM.hydrosea using direct logical indexing
hydro.sea$LOSOM.hydrosea[hydro.sea$month %in% c(11, 12, 1, 2)] <- "EarlyDry"
hydro.sea$LOSOM.hydrosea[hydro.sea$month %in% 3:5] <- "LateDry"
hydro.sea$LOSOM.hydrosea[hydro.sea$month %in% 6:10] <- "Wet"

# Assign hybrid hydro seasons
hydro.sea$hydrosea2[hydro.sea$month %in% 5:7] <- "EarlyWet"
hydro.sea$hydrosea2[hydro.sea$month %in% 8:10] <- "LateWet"
hydro.sea$hydrosea2[hydro.sea$month %in% c(11, 12, 1)] <- "EarlyDry"
hydro.sea$hydrosea2[hydro.sea$month %in% 2:4] <- "LateDry"

## GIS Data ----------------------------------------------------------------
data("LOK")

## Load all files in geopackage
# Path to your GeoPackage
gpkg_file <- paste0(data.path,"LOK.gpkg")

# List all layers
layers <- st_layers(gpkg_file)$name

# Read all layers into a named list
for (layer_name in layers) {
  assign(layer_name, st_read(gpkg_file, layer = layer_name))
}

# Sanity Check
plot(lok_zones[1])
plot(LOK_sites_all[1])
plot(AL.hurdat.LOK[1])
Walker.sites
site.zones
LOK.hur

## Hurricane Data Lookup ---------------------------------------------------
CY.periods$Hurr <- with(CY.periods,ifelse(CY%in%unique(LOK.hur$CY),1,0))
CY.periods$maj.Hurr <- with(CY.periods,ifelse(CY%in%unique(subset(LOK.hur,storm.cat!="TS")$CY),1,0))

WY.periods <- data.frame(WY=seq(WY(dates[1],"Fed"),WY(dates[2],"Fed"),1))
WY.periods$Hurr <- with(WY.periods,ifelse(WY%in%unique(LOK.hur$WY),1,0))
WY.periods$maj.Hurr <- with(WY.periods,ifelse(WY%in%unique(subset(LOK.hur,storm.cat!="TS")$WY),1,0))

## Flow --------------------------------------------------------------------
# Download each flow DBKEY
if(!any(list.files(data.path) == "LOK_Q.csv")) {
  
  pb <- txtProgressBar(min = 0, max = length(flow.dbkeys$DBKEY), style = 3)
  
  LOK_q_list <- lapply(seq_along(flow.dbkeys$DBKEY), function(i) {
    dbkey <- flow.dbkeys$DBKEY[i]
    
    # Attempt to retrieve data and handle errors
    result <- tryCatch({
      tmp <- DBHYDRO_daily(dates[1], dates[2], dbkey)
      tmp$DBKEY <- as.character(dbkey)
      tmp
    }, error = function(e) {
      message(sprintf("Skipping DBKEY %s due to error: %s", dbkey, e$message))
      NULL  # Return NULL if there's an error
    })
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
    
    return(result)
  })
  
  # Close the progress bar
  close(pb)
  
  # Remove NULL values from the list
  LOK_q_list <- Filter(Negate(is.null), LOK_q_list)
  
  # Merge the combined data with the dbkeys data frame
  dbvars <- c("DBKEY", "Suwatershed", "GenBasin", "Basin", "Qsite", "WQSite", "flow_dir","pref")
  LOK.q.da <- do.call(rbind, LOK_q_list)|>
    merge(flow.dbkeys[,dbvars], by = "DBKEY")|> 
    mutate(Date.EST = date.fun(Date),
           WY = WY(Date,"Fed"))
  
  write.csv(LOK.q.da, paste0(data.path, "LOK_Q.csv"), row.names = FALSE)
}

## Stage -------------------------------------------------------------------
comp.dbkey <- data.frame(DBKEY=c("06832","00268"),Priority=c("P2","P1"))

if(!any(list.files(data.path) == "LOK_STG.csv")) {
  pb <- txtProgressBar(min = 0, max = length(comp.dbkey$DBKEY), style = 3)
  
  lok.stg_list <- lapply(seq_along(comp.dbkey$DBKEY), function(i) {
    dbkey <- comp.dbkey$DBKEY[i]
    
    # Attempt to retrieve data and handle errors
    result <- tryCatch({
      tmp <- DBHYDRO_daily(dates[1], dates[2], dbkey)
      tmp$DBKEY <- as.character(dbkey)
      tmp
    }, error = function(e) {
      message(sprintf("Skipping DBKEY %s due to error: %s", dbkey, e$message))
      NULL  # Return NULL if there's an error
    })
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
    
    return(result)
  })
  
  # Close the progress bar
  close(pb)
  
  # Remove NULL values from the list
  lok.stg_list <- Filter(Negate(is.null), lok.stg_list)
  
  # Merge the combined data with the dbkeys data frame
  lok.stg <- do.call(rbind, lok.stg_list)
  
  write.csv(lok.stg, paste0(data.path, "LOK_STG.csv"), row.names = FALSE)
}

## WQ ----------------------------------------------------------------------
params <- data.frame(test_num = c(21,20,18,80,
                                  25,23,
                                  61,179,112,178,
                                  7,11,
                                  16),
                     param    = c("TKN","NH4","NOx","TN",
                                  "TP","SRP",
                                  "Chla","Chla","Chla","Chla",
                                  "TempC","SD",
                                  "TSS"))

if(!any(list.files(data.path) == "LOK_WQ.csv")) {
  pb <- txtProgressBar(min = 0, max = length(LOK_sites_all$STATION), style = 3)
  
  lok.wq_list <- lapply(seq_along(LOK_sites_all$STATION), function(i) {
    site <- LOK_sites_all$STATION[i]
    
    # Attempt to retrieve data and handle errors
    result <- tryCatch({
      tmp <- DBHYDRO_WQ(dates[1], dates[2], site, params$test_num)
      tmp
    }, error = function(e) {
      message(sprintf("Skipping STATION_ID %s due to error: %s", site, e$message))
      NULL  # Return NULL if there's an error
    })
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
    
    return(result)
  })
  
  # Close the progress bar
  close(pb)
  
  # Remove NULL values from the list
  lok.wq_list <- Filter(Negate(is.null), lok.wq_list)
  
  # Merge the combined data with the dbkeys data frame
  lok.wq <- do.call(rbind, lok.wq_list)
  
  
  write.csv(lok.wq, paste0(data.path, "LOK_WQ.csv"), row.names = FALSE)
}

# Structures WQ (for Load analysis) ------------------------------------------
wq.str.sites <- unique(flow.dbkeys$WQSite) 
params.str <- subset(params,!(param%in%c("TempC","Chla","SD")))

if(!any(list.files(data.path) == "LOK_STR_WQ.csv")) {
  pb <- txtProgressBar(min = 0, max = length(wq.str.sites), style = 3)
  
  str.wq.dat.list <- lapply(seq_along(wq.str.sites), function(i) {
    site <- wq.str.sites[i]
    
    # Attempt to retrieve data and handle errors
    result <- tryCatch({
      tmp <- DBHYDRO_WQ(dates[1], dates[2], site, params.str$test_num)
      tmp
    }, error = function(e) {
      message(sprintf("Skipping STATION_ID %s due to error: %s", site, e$message))
      NULL  # Return NULL if there's an error
    })
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
    
    return(result)
  })
  
  # Close the progress bar
  close(pb)
  
  # Remove NULL values from the list
  str.wq.dat.list <- Filter(Negate(is.null), str.wq.dat.list)
  
  # Merge the combined data with the dbkeys data frame
  str.wq.dat <- do.call(rbind, str.wq.dat.list)
  
  
  write.csv(str.wq.dat, paste0(data.path, "LOK_STR_WQ.csv"), row.names = FALSE)
}


# RSMBN -------------------------------------------------------------------
## Preps DSSRip
# options(dss_override_location="C:\\projects\\dssrip\\monolith")
# options(dss_config_filename="C:\\projects\\dssrip\\dssrip2.config")
# options(dss_default_config="monolith-win-x86_64")
# options(dss_allowed_states="untested")
# options(dssrip_debug=T)
# library(dssrip)


## LOSOM -------------------------------------------------------------------
LOSOM.RSM.datpath <- "./rsmbn/" # local path, downloaded from SMMS 
LOSOM.RSM.datpath <-"D:/RSM_Data/LOSOM/Iteration_3/Model_Output/"

LOSOM.alts  <-  c("NA25f","PA25"); # FWO and Preferred alternative


### Stage -------------------------------------------------------------------
if(!any(list.files(data.path) == "RSMBN_LOK.csv")) {
  lok.stg.loc <- "LOK"
  LOK.stg.list <- vector("list", length(LOSOM.alts) * length(lok.stg.loc))
  index <- 1  # Initialize an index for the list
  
  for (j in seq_along(LOSOM.alts)) {
    dss_out <- opendss(file.path(LOSOM.RSM.datpath, LOSOM.alts[j], "RSMBN_output.dss"))
    cat.dat <- pathsToDataFrame(getCatalogedPathnames(dss_out), simplify = TRUE)
    
    paths <- subset(cat.dat, LOCATION == lok.stg.loc & PARAMETER == "STAGE")$PATH
    
    tmp <- getFullTSC(dss_out, paths) |> data.frame()
    tmp$Date <- date.fun(rownames(tmp))
    rownames(tmp) <- NULL
    tmp$SITE <- lok.stg.loc
    tmp$Alt <- LOSOM.alts[j]
    tmp$proj <- "LOSOM"
    
    # Store the result in the preallocated list
    LOK.stg.list[[index]] <- tmp
    index <- index + 1
    
    print(j)
  }
  
  # Combine all results into a single data frame at the end
  LOK.LOSOM.stg <- do.call(rbind, LOK.stg.list)
}

### Discharge ---------------------------------------------------------------
Q.inflow.sites <- c('S65E','FEC','TOTAL_ISTOK','S77BF','S4BP','S3','S2',
                    'C12ABP','C12BP',"C10BP",'C4ABP','S236','P5WPS','S308BF',
                    'TCNSQ','S154','S135'); # missing from water budget file is LKTFPL, MDS, S271BK
Q.outflow.sites <-c('NELKSH_WS_QWS','NLKSH_WS_QWS','S77','S4_WS','S271','C12A','C12','C10',
                    'S352',"S351",'C4A','C3','S354','BRIGHTON_WS','S308')

LOK.wb.sites <- rbind(
  data.frame(SITE = Q.inflow.sites,dir = "inflow"),
  data.frame(SITE = Q.outflow.sites,dir = "outflow")
)

if(!any(list.files(data.path) == "RSMBN_LOK_Q.csv")) {
  
  # Preallocate a list to store results
  LOK.Q.list <- vector("list", length(LOSOM.alts) * nrow(LOK.wb.sites))
  index <- 1  # Initialize an index for the list
  
  for (j in seq_along(LOSOM.alts)) {
    dss_out <- opendss(file.path(LOSOM.RSM.datpath, LOSOM.alts[j], "RSMBN_output.dss"))
    cat.dat <- pathsToDataFrame(getCatalogedPathnames(dss_out), simplify = TRUE)
    
    for (i in seq_len(nrow(LOK.wb.sites))) {
      paths <- subset(cat.dat, LOCATION == LOK.wb.sites$SITE[i] & PARAMETER == "FLOW")$PATH
      
      # Skip if no paths are found
      if (length(paths) == 0) next
      
      tmp <- getFullTSC(dss_out, paths) |> data.frame()
      tmp$Date <- date.fun(rownames(tmp))
      rownames(tmp) <- NULL
      tmp$SITE <- LOK.wb.sites$SITE[i]
      tmp$Alt <- LOSOM.alts[j]
      tmp$proj <- "LOSOM"
      
      # Store the result in the preallocated list
      LOK.Q.list[[index]] <- tmp
      index <- index + 1
    }
    print(j)
  }
  
  LOK.LOSOM.Q <- do.call(rbind, LOK.Q.list)|>
    merge(LOK.wb.sites,"SITE")
}

## LOCAR -------------------------------------------------------------------
## Gathering file paths
LOCAR.RSM.datpath <-  "C:/Julian_LaCie/_GitHub/LOCAR_Eval/Data/RSMBN/"

alts.all <-  sapply(strsplit(list.files(LOCAR.RSM.datpath),"LOCAR_|_RSMBN"),"[",2)
out.files <-  list.files(paste0(LOCAR.RSM.datpath,list.files(LOCAR.RSM.datpath),"/",alts.all),full.names = T)
out.files <-  out.files[grep("output",out.files)]

LOCAR.alts = c("PA_FWOLL","LCR1","ECB23L")


### Stage -------------------------------------------------------------------

if(!any(list.files(data.path) == "RSMBN_LOK.csv")) {
  LOK.stg.list <- vector("list", length(LOCAR.alts) * length(lok.stg.loc))
  index <- 1  # Initialize an index for the list
  
  for (j in seq_along(LOCAR.alts)) {
    dss_out <- opendss(file.path(out.files[grep(LOCAR.alts[j],out.files)], "RSMBN_output.dss"))
    cat.dat <- pathsToDataFrame(getCatalogedPathnames(dss_out), simplify = TRUE)
    
    paths <- subset(cat.dat, LOCATION == lok.stg.loc & PARAMETER == "STAGE")$PATH
    
    tmp <- getFullTSC(dss_out, paths) |> data.frame()
    tmp$Date <- date.fun(rownames(tmp))
    rownames(tmp) <- NULL
    tmp$SITE <- lok.stg.loc
    tmp$Alt <- LOCAR.alts[j]
    tmp$proj <- "LOCAR"
    
    # Store the result in the preallocated list
    LOK.stg.list[[index]] <- tmp
    index <- index + 1
    
    print(j)
  }
  # Combine all results into a single data frame at the end
  LOK.LOCAR.stg <- do.call(rbind, LOK.stg.list)
  
  LOK.stg <- rbind(LOK.LOSOM.stg,subset(LOK.LOCAR.stg,Alt!="ECB23L"))
  write.csv(LOK.stg,paste0(data.path,"RSMBN_LOK.csv"),row.names = F)
}


### Discharge ---------------------------------------------------------------

## adjusted for changes in Indian Prairie Basin Refinement
Q.inflow.sites <- c(Q.inflow.sites,"S84","C40C41LSCOUT")
Q.outflow.sites <-c(Q.outflow.sites,'G207G208')

LOK.wb.sites <- rbind(
  data.frame(SITE = Q.inflow.sites,dir = "inflow"),
  data.frame(SITE = Q.outflow.sites,dir = "outflow")
)

if(!any(list.files(data.path) == "RSMBN_LOK_Q.csv")) {
  
  # Preallocate a list to store results
  LOK.Q.list <- vector("list", length(LOCAR.alts) * nrow(LOK.wb.sites))
  index <- 1  # Initialize an index for the list
  
  for (j in seq_along(LOCAR.alts)) {
    dss_out <- opendss(file.path(out.files[grep(LOCAR.alts[j],out.files)], "RSMBN_output.dss"))
    cat.dat <- pathsToDataFrame(getCatalogedPathnames(dss_out), simplify = TRUE)
    
    for (i in seq_len(nrow(LOK.wb.sites))) {
      paths <- subset(cat.dat, LOCATION == LOK.wb.sites$SITE[i] & PARAMETER == "FLOW")$PATH
      
      # Skip if no paths are found
      if (length(paths) == 0) next
      
      tmp <- getFullTSC(dss_out, paths) |> data.frame()
      tmp$Date <- date.fun(rownames(tmp))
      rownames(tmp) <- NULL
      tmp$SITE <- LOK.wb.sites$SITE[i]
      tmp$Alt <- LOCAR.alts[j]
      tmp$proj <- "LOCAR"
      
      # Store the result in the preallocated list
      LOK.Q.list[[index]] <- tmp
      index <- index + 1
    }
    print(j)
  }
  
  # Combine all results into a single data frame at the end
  LOK.LOCAR.Q <- do.call(rbind, LOK.Q.list)|>
    merge(LOK.wb.sites,"SITE")
  
  LOK.RSMBN.Q <- rbind(LOK.LOSOM.Q,subset(LOK.LOCAR.Q,Alt!="ECB23L"))
  write.csv(LOK.RSMBN.Q,paste0(export.path,"RSMBN_LOK_Q.csv"),row.names = F)
}
