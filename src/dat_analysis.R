## Lake Okeechobee Algae Dynamics - Data Analysis
## Created by: Paul Julian (pjulian@evergladesfoundation.org)
## Created on: 2024-06-26

library(AnalystHelper)
library(EVERSpatDat)
library(reshape2)
library(plyr)
library(sf)
library(raster)
library(terra)

## GAM
library(mgcv)
library(gratia)
library(ggplot2)
library(viridis)
library(DHARMa)

theme_set(theme_minimal(base_size = 16)+
            theme_bw()+
            theme(panel.border = element_rect("black",fill=NA,linewidth=1)))


# Other stats
library(vegan)
library(REdaS)

## Data viz/reporting
library(flextable)

# Paths
wd="C:/Julian_LaCie/_GitHub/LOK_AlgaeEval_manu"
paths=paste0(wd,c("/Plots/","/Data/","/src/","/_documents/"))

#Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
data.path=paths[2]

# Stored on local drive
GIS.path.gen="C:/Julian_LaCie/_GISData"

# CRS
nad83.pro=st_crs("EPSG:4269")
NAD83_HARN <- st_crs("EPSG:2881")
utm17=st_crs("EPSG:26917")
wgs84=st_crs("EPSG:4326")


# Functions
source(paste0(wd,"/src/func.R"))


# load(paste0(export.path,"LOK_AlgaeEval_analysis.RData"))

# Data --------------------------------------------------------------------
dates <- date.fun(c("1974-10-01","2023-09-30"))

dates2 <- dates # for the sonde data
dates2[1] <- date.fun("2016-10-01")

## GIS ---------------------------------------------------------------------
data("LOK")
bath <- terra::rast(paste0(data.path,"LOK_bathy.tif"))|>
  terra::project(utm17$wkt)

gpkg_file <- paste0(data.path,"LOK.gpkg")

# List all layers
st_layers(gpkg_file)$name

lok_zones <- st_read(gpkg_file,"lok_zones")
LOK_sites_all <- st_read(gpkg_file,"LOK_sites_all")
AL.hurdat.LOK <- st_read(gpkg_file,"AL.hurdat.LOK")


LOK.hur <- st_drop_geometry(AL.hurdat.LOK)|>as.data.frame()
storm.cat.val <- c("TS",paste0("Cat_",1:5))
LOK.hur$storm.cat <- storm.cat.val[findInterval(LOK.hur$maxWind,c(0,74,96,111,130,157))]
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

WY.periods <- data.frame(WY=seq(WY(dates[1],"Fed"),WY(dates[2],"Fed"),1))
WY.periods$Hurr <- with(WY.periods,ifelse(WY%in%unique(LOK.hur$WY),1,0))
WY.periods$maj.Hurr <- with(WY.periods,ifelse(WY%in%unique(subset(LOK.hur,storm.cat!="TS")$WY),1,0))

## Stage -------------------------------------------------------------------
comp.dbkey=data.frame(DBKEY=c("06832","00268"),Priority=c("P2","P1"))

lok.stg <- read.csv(paste0(data.path,"LOK_stg.csv"),colClasses = c(rep("character",3),"numeric",rep("character",3)))
lok.stg <- merge(lok.stg,comp.dbkey,"DBKEY")
lok.stg$Date <- date.fun(lok.stg$Date,tz="America/New_York")

LakeO.xtab <- dcast(lok.stg,Date~Priority,value.var="Data.Value",mean)
LakeO.xtab$STG29 <- with(LakeO.xtab,ifelse(is.na(P1)==T,P2,P1))

# Fill for continuous timeseries (just dates)
fill <- data.frame(Date=date.fun(seq(min(LakeO.xtab$Date),max(LakeO.xtab$Date),"1 days"),tz="America/New_York"))

min.stg=11.5
LakeO.xtab <- merge(LakeO.xtab,fill,"Date",all.y=T)|>
  mutate(
    Date.EST = date.fun(Date),
    DOWY = hydro.day(Date.EST,"Fed"),
    CY = as.numeric(format(Date.EST,"%Y")),
    month = as.numeric(format(Date.EST,"%m")),
    recess_30day = c(rep(NA,30),diff(STG29,lag=30)),
    delta.min.stg = pmax(STG29-min.stg,0,na.rm=T))

lok.stg <- merge(LakeO.xtab,month.season,"month")
lok.stg <- merge(lok.stg,CY.periods,"CY")
lok.stg <- lok.stg[order(lok.stg$Date.EST),]
lok.stg$WY <- WY(lok.stg$Date.EST,"Fed")

## Lake Stage - Volume - Area ----------------------------------------------
# Pre-calculate constants
cell_area_m2 <- res(bath)[1] * res(bath)[2]  # Cell area in m²
m2.to.acres # m2_to_acres <- function(x) x / 4046.86
m3.to.acft <- function(x) x / 1233.48
acft.to.m3 <- function(x) x * 1233.48

# Function to compute inundation, area, and volume

compute_metrics <- function(z) {
  if(class(bath)=="SpatRaster"){
    inundate <- bath < z
    depth <- (z - bath)*inundate
    m.depth <- ifel(depth>0,depth,NA)
    
    volume <- ft.to.m(depth) * cell_area_m2  # m³
    area <- inundate * cell_area_m2  # m²
    
    Tarea <- global(area, fun = "sum", na.rm = TRUE)|>as.numeric();  # m2
    Tvol <- global(volume, fun = "sum", na.rm = TRUE)|>as.numeric(); # m3
    # mean.depth <- global(m.depth,fun="mean",na.rm=T)|>as.numeric(); # ft
  }
  
  data.frame(
    STG29.ft = z,
    depth.ft = m3.to.acft(Tvol)/m2.to.acres(Tarea),
    area.acres = m2.to.acres(Tarea),
    volume.acft = m3.to.acft(Tvol),
    STG29.m = ft.to.m(z),
    depth.m = Tvol/Tarea,
    area.m2 = Tarea,
    volume.m3 = Tvol
  )
}


# Apply the function across all stages
if(!any(list.files(data.path) == "LOK_StgAreaVol.csv")) {
  stg.seq <- seq(8, 19, 0.1)
  LOK.stg.area.vol <- do.call(rbind, lapply(stg.seq, compute_metrics))
  write.csv(LOK.stg.area.vol,paste0(data.path,"LOK_StgAreaVol.csv"),row.names = F)
}else{
  LOK.stg.area.vol <- read.csv(paste0(data.path,"LOK_StgAreaVol.csv")) 
}

# Interpolation functions
get_LOK_area <- approxfun(LOK.stg.area.vol$STG29.ft, LOK.stg.area.vol$area.acres, rule = 2)    # Volume -> Area
get_LOK_stage <- approxfun(LOK.stg.area.vol$volume.acft, LOK.stg.area.vol$STG29.ft, rule = 2)  # Volume -> Stage
get_LOK_volume <- approxfun(LOK.stg.area.vol$STG29.ft, LOK.stg.area.vol$volume.acft, rule = 2) # Stage -> Volume
get_LOK_depth <- approxfun(LOK.stg.area.vol$STG29.ft, LOK.stg.area.vol$depth.ft, rule = 2) # Stage -> mean Depth

lok.stg <- lok.stg|>
  mutate(
    volume.acft = get_LOK_volume(STG29),
    area.acres = get_LOK_area(STG29),
    depth = get_LOK_depth(STG29)
  )

mon.vol <- ddply(lok.stg,c("CY","month","WY"),summarise,
              mean.vol.km3=mean(acft.to.m3(volume.acft)/1e9,na.rm=T),
              mean.area.km2=mean(acres.to.m2(area.acres)/1e6,na.rm=T),
              mean.depth.m = mean(ft.to.m(depth),na.rm=T),
              mean.STG = mean(STG29,na.rm=T))

WY.vol <- ddply(lok.stg,c("WY"),summarise,
             mean.vol.km3=mean(acft.to.m3(volume.acft)/1e9,na.rm=T),
             mean.area.km2=mean(acres.to.m2(area.acres)/1e6,na.rm=T),
             mean.depth.m = mean(ft.to.m(depth),na.rm=T),
             mean.STG = mean(STG29,na.rm=T))

## Discharge ---------------------------------------------------------------
flow.dbkeys <- st_read( paste0(data.path,"LOK.gpkg"), layer = "flow.dbkeys")

flow.dat <- read.csv(paste0(data.path,"LOK_Q.csv"))
flow.dat$Date <- date.fun(flow.dat$Date,tz="America/New_York")
flow.dat$Date.EST <- date.fun(flow.dat$Date.EST)

flow.xtab <- dcast(flow.dat,
                   Date+WY+Suwatershed+GenBasin+Qsite+flow_dir+Basin+WQSite~pref,
                   value.var = "Data.Value",sum,na.rm=T)|>
  mutate(
    Date.EST = date.fun(Date),
    fflow.cfs = ifelse(is.na(P1),P2,P1)*flow_dir,
    direct = ifelse(fflow.cfs<0,"Outflow","Inflow"),
    month = as.numeric(format(Date,"%m")),
    CY = as.numeric(format(Date,"%Y"))
  )

flow.xtab2 <- dcast(flow.xtab,WY+CY+month+Date+direct~Qsite,
                    value.var = "fflow.cfs",
                    fun.aggregate=function(x) sum(abs(x),na.rm=T)
)

idvars <- c("WY", "CY", "month","Date", "direct")
struct.vars <- c("C10", "C12", "C12A",  "C41H78", "C4A", 
                 "CU5A", "CV5", "FISHCR", "G207", "G208","G33", "G34", "G74", "G75", 
                 "G76", "INDUST", "L8.441", "S127", "S129", "S131", "S133", "S135", 
                 "S154","S154C", "S191", "S2", "S236", "S3", "S308", "S351", "S352", "S354", 
                 "S4", "S65", "S65E", "S65EX1", "S68", "S71", "S72", "S77", "S84", 
                 "S84X", "S65E.tot", "BasinOut", "S84.tot", "S72.adj", "S71.adj")
struct.vars <- struct.vars[!(struct.vars%in%c("S65","S68", # upstream basins
                                              "S65E.tot", # no need for aggregated 
                                              "S84.tot", # no need for aggregated
                                              "S71.adj","S72.adj","BasinOut") # not calculating adjusted values
)]

QSite.meta <- ddply(flow.dbkeys,c("GenBasin","Suwatershed","Basin","WQSite","Qsite"),summarise,N.val = N.obs(DBKEY))
QSite.meta <- QSite.meta[-ncol(QSite.meta)]
QSite.meta <- QSite.meta|>
  rbind(subset(QSite.meta,Qsite=="S84")|>
          mutate(Qsite="S84.tot"))|>
  rbind(subset(QSite.meta,Qsite=="S65E")|>
          mutate(Qsite="S65E.tot"))

flow.xtab.melt <- melt(flow.xtab2[,c(idvars,struct.vars)],id.vars=idvars,
                       variable.name="Qsite",value.name="flow")|>
  merge(QSite.meta,by="Qsite",all.x=T)

flow.mon.sum <- dcast(flow.xtab.melt,Suwatershed+GenBasin+Qsite+WQSite+WY+CY+month~direct,
                      value.var="flow",
                      fun.aggregate=function(x) sum(abs(cfs.to.acftd(x)),na.rm=T))
flow.WY.sum <- ddply(flow.mon.sum,c("WY"),summarise,
                     inflow.Q = sum(Inflow,na.rm=T),
                     outflow.Q = sum(Outflow,na.rm=T))

flow.mon.CY.sum <- ddply(flow.mon.sum,c("CY","month","WY"),summarise,
                         inflow.Q = sum(Inflow,na.rm=T),
                         outflow.Q = sum(Outflow,na.rm=T))

# In-lake Water Quality ---------------------------------------------------
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

LOK_sites_all <- st_read(paste0(data.path,"LOK.gpkg"),"LOK_sites_all")

WQ.dat <- read.csv(paste0(data.path,"LOK_WQ.csv"))
WQ.dat$Date <- date.fun(WQ.dat$Date,tz="America/New_York")

WQ.dat <- WQ.dat|>
  mutate(
    Date.EST = date.fun(Date),
    CY = as.numeric(format(Date.EST,"%Y")),
    month = as.numeric(format(Date.EST,"%m")),
    WY = WY(Date.EST,"Fed")
  )

WQ.dat <- merge(WQ.dat,st_drop_geometry(LOK_sites_all),by.x="Station.ID",by.y="STATION",all.x=T)|>
  merge(params,by.x="Test.Number",by.y="test_num")

LOK.sites.all.coords <- cbind(data.frame(Station.ID=LOK_sites_all$STATION),
                           UTMX=st_coordinates(LOK_sites_all)[,1],
                           UTMY=st_coordinates(LOK_sites_all)[,2])

vars.val <- c("Station.ID", "walker.Region", "EcoZone3", "bath.NGVD29", "Date.EST", 
              "Chla", "NH4", "NOx", "SD", "SRP", "TempC", "TKN", "TN", "TP")
WQ.dat.xtab <- dcast(WQ.dat,
                  Station.ID+walker.Region+EcoZone3+bath.NGVD29+Date.EST~param,
                  value.var="HalfMDL",mean)[,vars.val]|>
  subset(is.na(Chla)==F)|>
  mutate(TN=TN_Combine(NOx,TKN,TN),
         DIN=pmax(NH4,psum_fun(NOx,NH4,na.rm=T)),
         CY=as.numeric(format(Date.EST,"%Y")),
         month=as.numeric(format(Date.EST,"%m")),
         DOY=as.numeric(format(Date.EST,"%j")),
         WY=WY(Date.EST,"Fed"),
         FLWY=WY(Date.EST),
         hydro.season=FL.Hydroseason(Date.EST))|>
  merge(CY.periods,'CY')|>
  merge(month.season,"month")|>
  merge(lok.stg[,c("Date.EST","STG29","volume.acft")],"Date.EST",all.x=T)|>
  mutate(z.ft=pmax(STG29-bath.NGVD29,0,na.rm=T))|>
  merge(LOK.sites.all.coords,"Station.ID")



## sample screening
structure.sites <- unique(subset(LOK_sites_all,zone=="struct")$Station.ID)
chla.sea.screen <- dcast(subset(WQ.dat.xtab,!(Station.ID%in%structure.sites)),
                         Station.ID+WY~hydro.season,value.var="Chla",
                         fun.aggregate = function(x) N.obs(x))
chla.sea.screen$Tsamp <- rowSums(chla.sea.screen[,c("A_Wet","B_Dry")],na.rm=T)
chla.sea.screen$scrn.val <- with(chla.sea.screen,ifelse(A_Wet>0&B_Dry>0&Tsamp>=4,1,0))
chla.WY.screen <- ddply(subset(chla.sea.screen,scrn.val==1),"Station.ID",
                        summarise,N.WY=N.obs(WY),minWY=min(WY),maxWY=max(WY))
# check how many samples screend samples 
# nrow(subset(chla.sea.screen,scrn.val==0))
# (nrow(subset(chla.sea.screen,scrn.val==0))/nrow(chla.sea.screen))*100

WQ.dat.xtab <- merge(WQ.dat.xtab,chla.sea.screen[,c("Station.ID","WY","scrn.val")],c("Station.ID","WY"))


# Structure Water Quality -------------------------------------------------
str.wq.dat <-  read.csv(paste0(data.path,"LOK_STR_WQ.csv"))
str.wq.dat$Date <- date.fun(str.wq.dat$Date,tz="America/New_York")
str.wq.dat$Date.EST <- date.fun(str.wq.dat$Date)

str.wq.dat <- merge(str.wq.dat,params,
                    by.x = "Test.Number",by.y = "test_num")

# Daily mean WQ removing LAB project code.
str.wq.dat.xtab <- subset(str.wq.dat,Project.Code!="LAB"&
                            Collection.Method == "G")|>
  dcast(Station.ID+Date.EST~param,value.var="HalfMDL",mean,na.rm=T)
## Check if TKN field is present, if not add 
if(sum(names(str.wq.dat.xtab)%in%c("TKN"))==0){str.wq.dat.xtab$TKN <- NA}
str.wq.dat.xtab$TN <- with(str.wq.dat.xtab,TN_Combine(NOx,TKN,TN))
str.wq.dat.xtab$DIN <- with(str.wq.dat.xtab,NOx+NH4)
head(str.wq.dat.xtab)

str.wq.dat.xtab$CY <- as.numeric(format(str.wq.dat.xtab$Date.EST,"%Y"))


## Flow and WQ pairing - Load Data -----------------------------------------
ALIAS.vals <- ddply(flow.dbkeys,"Qsite",summarise,N.val=N.obs(DBKEY))
date.fill.frame <- expand.grid(Date.EST=date.fun(seq(dates[1],dates[2],"1 days")),
                               Qsite=ALIAS.vals$Qsite)
flow.xtab.melt$Date.EST <- date.fun(flow.xtab.melt$Date)

# check flow and wq site pairing
ddply(flow.dbkeys,c("Qsite","WQSite"),summarise,N.val=N.obs(Qsite))|>
  ddply("WQSite",summarise,N.val=N.obs(Qsite))

flow.vars <- c("Qsite","WQSite","GenBasin","Basin","Suwatershed","Date.EST","WY","direct","flow")
wq.vars <- c("Date.EST","Station.ID","TP","TN",'TSS')

flow.wq <- merge(flow.xtab.melt[,flow.vars],
                 str.wq.dat.xtab[,wq.vars],
                 by.x=c("Date.EST","WQSite"),
                 by.y=c("Date.EST","Station.ID"),all.x=T)
flow.wq <- merge(date.fill.frame,flow.wq,c("Date.EST","Qsite"),all.y=T)
flow.wq <- flow.wq[order(flow.wq$Qsite,flow.wq$Date.EST),]

flow.wq$TP.int <- with(flow.wq,ave(TP,Qsite,FUN = function(x)dat.interp(x)))
flow.wq$TPLoad.kg <- with(flow.wq,Load.Calc.kg(abs(flow),TP.int))
flow.wq$TN.int <- with(flow.wq,ave(TN,Qsite,FUN = function(x)dat.interp(x)))
flow.wq$TNLoad.kg <- with(flow.wq,Load.Calc.kg(abs(flow),TN.int))
flow.wq$TSS.int <- with(flow.wq,ave(TSS,Qsite,FUN = function(x)dat.interp(x)))
flow.wq$TSSLoad.kg <- with(flow.wq,Load.Calc.kg(abs(flow),TSS.int))
flow.wq$month <- as.numeric(format(flow.wq$Date.EST,"%m"))
flow.wq$CY <- as.numeric(format(flow.wq$Date.EST,"%Y"))

## Monthly
Qdat <- dcast(flow.wq,CY+month+WY~direct,value.var = "flow",fun.aggregate = function(x) sum(cfs.to.m3d(x),na.rm=T))
Qdat.names <- c(names(Qdat)[!(names(Qdat)%in%c("Inflow","Outflow"))],paste(c("Inflow","Outflow"),"Q.m3",sep="."))
colnames(Qdat) <- Qdat.names

TPLdat <- dcast(flow.wq,CY+month+WY~direct,value.var = "TPLoad.kg",sum,na.rm=T)
TPLdat.names <- c(names(TPLdat)[!(names(TPLdat)%in%c("Inflow","Outflow"))],paste(c("Inflow","Outflow"),"TPL.kg",sep="."))
colnames(TPLdat) <- TPLdat.names

TNLdat <- dcast(flow.wq,CY+month+WY~direct,value.var = "TNLoad.kg",sum,na.rm=T)
TNLdat.names <- c(names(TNLdat)[!(names(TNLdat)%in%c("Inflow","Outflow"))],paste(c("Inflow","Outflow"),"TNL.kg",sep="."))
colnames(TNLdat) <- TNLdat.names

flow.TP.mon.CY  <- merge(Qdat,TPLdat,c("CY", "month","WY"))|>
  merge(mon.vol,c("CY", "month","WY"),all.x=T)|>
  mutate(Date.monCY = date.fun(paste(CY,month,"01",sep="-")))

fill.val <- data.frame(Date.monCY = seq(min(flow.TP.mon.CY$Date.monCY),max(flow.TP.mon.CY$Date.monCY),"1 month"))|>
  mutate(
    CY = as.numeric(format(Date.monCY,"%Y")),
    month = as.numeric(format(Date.monCY,"%m")),
    WY = WY(Date.monCY,"Fed")
  )
flow.TP.mon.CY <- merge(flow.TP.mon.CY,fill.val,c("CY", "month","WY","Date.monCY"),all.y=T)
flow.TP.mon.CY$WRT.yr  <- with(flow.TP.mon.CY,ifelse(Outflow.Q.m3==0,(mean.vol.km3/0.0001)/12,
                                                     (mean.vol.km3/(Outflow.Q.m3/1e9))/12))
flow.TP.mon.CY <- flow.TP.mon.CY[order(flow.TP.mon.CY$CY,flow.TP.mon.CY$month),]

# Quick plot
ggplot(flow.TP.mon.CY,aes(x=Date.monCY,y=WRT.yr))+
  geom_line()+
  geom_point(shape=21,fill="grey")+
  scale_y_log10()+
  labs(x = "Date (CY)",
       y = "Residence Time (years)")

## WY
Qdat.WY <- dcast(flow.wq,WY~direct,value.var = "flow",fun.aggregate = function(x) sum(cfs.to.m3d(x),na.rm=T))
Qdat.names <- c(names(Qdat.WY)[!(names(Qdat.WY)%in%c("Inflow","Outflow"))],paste(c("Inflow","Outflow"),"Q.m3",sep="."))
colnames(Qdat.WY) <- Qdat.names

TPLdat.WY <- dcast(flow.wq,WY~direct,value.var = "TPLoad.kg",sum,na.rm=T)
TPLdat.names <- c(names(TPLdat.WY)[!(names(TPLdat.WY)%in%c("Inflow","Outflow"))],paste(c("Inflow","Outflow"),"TPL.kg",sep="."))
colnames(TPLdat.WY) <- TPLdat.names

TNLdat.WY <- dcast(flow.wq,WY~direct,value.var = "TNLoad.kg",sum,na.rm=T)
TNLdat.names <- c(names(TNLdat.WY)[!(names(TNLdat.WY)%in%c("Inflow","Outflow"))],paste(c("Inflow","Outflow"),"TNL.kg",sep="."))
colnames(TNLdat.WY) <- TNLdat.names

flow.TP.WY  <- merge(Qdat.WY,TPLdat.WY,c("WY"))|>
  merge(TNLdat.WY,c("WY"))|>
  merge(WY.vol,c("WY"),all.x=T)

flow.TP.WY$WRT.yr  <- with(flow.TP.WY,(mean.vol.km3/(Outflow.Q.m3/1e9)))

ggplot(flow.TP.WY,aes(x=WY,y=WRT.yr))+
  geom_line()+
  geom_point(shape=21,fill="grey")+
  scale_y_log10()+
  labs(x = "Water Year",
       y = "Residence Time (years)")


# Sonde Data ---------------------------------------------
PC.dat <- read.csv(paste0(data.path,"LOK_sonde.csv"))
PC.dat$Date <- date.fun(PC.dat$Date,tz="America/New_York")

PC.dat.xtab <- dcast(PC.dat,SITE+Region+Date~param,value.var = "Data.Value",mean)|>
  mutate(Date.EST=date.fun(Date),
         month=as.numeric(format(Date.EST,"%m")),
         CY=as.numeric(format(Date.EST,"%Y")),
         DOY=as.numeric(format(Date.EST,"%j")),
         DOWY=hydro.day(Date.EST,"Fed"),
         WY=WY(Date.EST,'Fed'),
         decMonth=dec.month(Date.EST),
         decWY=decimal.WY(Date.EST,"Fed"),
         decCY=lubridate::decimal_date(Date.EST)
  )
head(PC.dat.xtab)

PC.dat.xtab=merge(PC.dat.xtab,lok.stg[,c("Date.EST","STG29")],"Date.EST")
PC.dat.xtab$STATION=with(PC.dat.xtab,ifelse(substr(SITE,1,1)=="L",paste(SITE,"S",sep="-"),paste0(SITE,"S")))

## some data have to be screened out
plot(FTURB~PU,PC.dat.xtab);abline(0,1,col="red")
points(FTURB~PU, subset(PC.dat.xtab,SITE=="L001"&WY==2021&DOY<20),pch=19,col="red")
subset(PC.dat.xtab,PU>800)

# Analyses ----------------------------------------------------------------
## https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/mgcv-parallel.html
library(parallel)
detectCores()
nc <- 6
cl <- makeCluster(nc)


## NMDS Data --------------------------------------------------------------------
lok.stg2  <-  lok.stg|>
  mutate(FLWY = WY(Date.EST),
         FedWY = WY(Date.EST,"Fed"))

## Flow, TP and TN load, stage, volume, area, depth, WRT
flow.TP.WY

## In lake WQ
WQ.dat.xtab <-  WQ.dat.xtab|>
  mutate(FLWY = WY(Date.EST),
         FedWY = WY(Date.EST,"Fed"),
         z.m = ifelse(z.ft<0.328,0.1,ft.to.m(z.ft)),
         TP.mgm2=(TP*1000)*z.m, # convert mg/L to mg/m3 then multiply by depth m = mg/m2
         TN.mgm2=(TN*1000)*z.m # convert mg/L to mg/m3 then multiply by depth m = mg/m2
  )

range(WQ.dat.xtab$z.m)
vars <- c("Station.ID","EcoZone3","month","CY",'Date.EST',"TN","TP","DIN","SRP","Chla","TempC","z.m")
WQ.melt <-  melt(subset(WQ.dat.xtab,scrn.val==1)[,vars],id.vars=vars[1:5])|>
  mutate(FLWY = WY(Date.EST),
         FedWY = WY(Date.EST,"Fed"))
head(WQ.melt)

WQ.melt.WY  <-  ddply(WQ.melt,c("Station.ID","EcoZone3","FedWY","variable"),
                   summarise,
                   mean.val = mean(value,na.rm=T),
                   GM.val = exp(mean(log(value),na.rm=T)))|>
  dcast(EcoZone3+FedWY~variable,value.var = "GM.val",mean,na.rm=T)

## In lake Load
M_in.site.WY <- ddply(WQ.dat.xtab,c("Station.ID","EcoZone3","FedWY"),
                     summarise,
                     M_in.TP.mgm2 = mean(TP.mgm2,na.rm=T),
                     M_in.TN.mgm2 = mean(TN.mgm2,na.rm=T))


lok_zones$area  <-  st_area(lok_zones)
LOK.EcoZone3.area = ddply(lok_zones,"EcoZone3",summarise,Tarea = sum(as.numeric(area)))

M_in.WY <- ddply(M_in.site.WY,c("EcoZone3","FedWY"),
                summarise,
                M_in.TP.mgm2 = mean(M_in.TP.mgm2,na.rm=T),
                M_in.TN.mgm2 = mean(M_in.TN.mgm2,na.rm=T))|>
  merge(LOK.EcoZone3.area,'EcoZone3')|>
  mutate(M_in.TP = (M_in.TP.mgm2*Tarea)*1e-6,
         M_in.TN = (M_in.TN.mgm2*Tarea)*1e-6)

ggplot(M_in.WY,aes(x=FedWY,y=M_in.TP,color=EcoZone3))+
  geom_point(shape=19)+
  geom_line(linewidth = 0.75)+
  scale_y_continuous(label = function(x)x/1e3)+
  labs(y="TP In-Lake Load (x10\u00B3 kg)",
       x="Federal Water Year",
       color = "Ecological\nRegions")

WQ.hydro.WY <-  WQ.melt.WY|>
  merge(flow.TP.WY,by.x="FedWY",by.y="WY")|>
  merge(M_in.WY[,c("FedWY","EcoZone3","M_in.TP","M_in.TN")],c("FedWY","EcoZone3"))

### NMDS --------------------------------------------------------------------
# https://r.qcbs.ca/workshop09/book-en/nonmetric-multidimensional-scaling.html
# https://ourcodingclub.github.io/tutorials/ordination/
## PCA is more focused on the dimensions themselves,
## and seek to maximize explained variance, whereas MDS is 
## more focused on relations among the scaled objects.
## from http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/122-multidimensional-scaling-essentials-algorithms-and-r-code/#google_vignette

vars <- c("FedWY", "EcoZone3",
       "TN", "TP", "DIN", "SRP", "Chla", "TempC", 
       "mean.STG", "mean.vol.km3","WRT.yr", 
       "Inflow.Q.m3","Inflow.TNL.kg", "Inflow.TPL.kg", 
       "M_in.TP", "M_in.TN")

ord.dat <-  na.omit(WQ.hydro.WY[,vars])
ord.dat$EcoZone3 <- as.factor(ord.dat$EcoZone3)
KMOS(ord.dat[,vars[3:length(vars)]])

# Run the NMDS
## checking other distance methods
ord1 <- metaMDS(ord.dat[,vars[3:length(vars)]],distance = "gower", k = 2)
ord2 <- metaMDS(ord.dat[,vars[3:length(vars)]],distance = "kulczynski", k = 2)
ord3 <- metaMDS(ord.dat[,vars[3:length(vars)]],distance = "jaccard", k = 2)
ord4 <- metaMDS(ord.dat[,vars[3:length(vars)]],distance = "bray", k = 2)

data.frame(stress = c(ord1$stress,ord2$stress,ord3$stress,ord4$stress),
           method = c(ord1$distmethod,ord2$distmethod,ord3$distmethod,ord4$distmethod))

ord <- metaMDS(ord.dat[,vars[3:length(vars)]],distance = "kulczynski", k = 2)
ord$stress;

envfit(ord, ord.dat[,vars[3:length(vars)]])

stress.diag <- goeveg::dimcheckMDS(ord.dat[,vars[3:length(vars)]])
shep.plot <-  stressplot(ord)|>as.data.frame()

shep <- MASS::Shepard(metaMDSredist(ord),ord$points)
shep.vals  <-  stress.fit.fun(ord)

library(patchwork)
dim_vals <- 1:6
stress_df <- data.frame(Dimension = dim_vals, Stress = stress.diag)

# Stress and Shepard Plots
p1 <- ggplot(stress_df, aes(x = Dimension, y = Stress)) +
  geom_hline(yintercept = seq(0, 0.4, 0.1), color = "grey80", linetype = "dotted") +
  geom_vline(xintercept = dim_vals, color = "grey80", linetype = "dotted") +
  geom_line(color = "black", size = 1.5) +
  geom_point(shape = 21, fill = "indianred1", color = "black", size = 4) +
  geom_text(data = stress_df[1:3, ],
            aes(label = round(Stress, 3)),
            hjust = -0.2, vjust = 0, size = 3) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(breaks = dim_vals, limits = c(1, 6)) +
  scale_y_continuous(breaks = seq(0, 0.4, 0.1), limits = c(0, 0.4)) +
  labs(x = "Dimension", y = "Stress") +
  # theme_minimal(base_family = "serif") +
  theme(plot.title = element_text(hjust = 1)) +
  ggtitle("A")

# Shepard plot
p2 <- ggplot(as.data.frame(shep), aes(x = x, y = y)) +
  # geom_hline(yintercept = seq(0, 0.5, 0.1), color = "grey80", linetype = "dotted") +
  # geom_vline(xintercept = seq(0, 0.35, 0.1), color = "grey80", linetype = "dotted") +
  geom_point(alpha = 0.5, color = "dodgerblue1", size = 0.8) +
  geom_line(aes(y = yf),color = "red", size = 1.2) +
  scale_x_continuous(breaks = seq(0, 0.35, 0.1), limits = c(0, 0.35)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) +
  annotate("text", x = 0, y = 0.4875,
           label = paste0(shep.vals$Fit, "\nFit R²: ", shep.vals$R2),
           hjust = 0, size = 3) +
  labs(x = "Ordination Distance", y = "Observed Dissimilarity") +
  # theme_minimal(base_family = "serif") +
  theme(plot.title = element_text(hjust = 1)) +
  ggtitle("B")

# Combine using patchwork
p1 + p2

# NMDS plot
cols <- MetBrewer::met.brewer("Hiroshige", 5)
# Site scores
sites_df <- as.data.frame(scores(ord, display = "sites"))
sites_df$EcoZone3 <- ord.dat$EcoZone3

# Species scores
species_df <- as.data.frame(scores(ord, display = "species"))
species_df$label <- c("TN", "TP", "DIN", "SRP", "Chl-a", "Temp", "Mean Stage", 
                      "Mean Vol.", "WRT", "Inflow Q", "Inflow TN Load", 
                      "Inflow TP Load", "In Lake\nTP Load", "In Lake\nTN Load")

ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(data = sites_df, aes(x = NMDS1, y = NMDS2, fill = EcoZone3),
             shape = 21, color = "grey30", size = 2, stroke = 0.2) +
  scale_fill_manual(values = cols) +
  stat_ellipse(type = "sd", level = 0.95)+
  geom_segment(data = species_df, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "indianred1", linewidth = 1) +
  geom_text(data = species_df, aes(x = NMDS1, y = NMDS2, label = label),
            color = "black", fontface = "bold", size = 3) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2") +
  coord_cartesian(xlim = c(-0.25, 0.35), ylim = c(-0.3, 0.3))+
  theme(legend.position = "right",
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 8)))

### ANOSIM ------------------------------------------------------------------
vars2 <- c("TN", "TP", "DIN", "SRP", "Chla", "TempC", "mean.STG",
           "M_in.TP", "M_in.TN","WRT.yr","mean.vol.km3","Inflow.TPL.kg","Inflow.TNL.kg")
ano = anosim(ord.dat[,vars2],
             grouping = ord.dat$EcoZone3,
             distance = "kulczynski")
plot(ano)
summary(ano)

R.values <- with(ano, data.frame(R = c(statistic, perm) ) )
R.values$Type <- c("actual", rep("perm", length(R.values$R) - 1))
plot(density(R.values$R));abline(v=R.values[R.values$Type == "actual" , "R"],col="red")

## Pairwise 
pwrslt <-  anosim.pw(ord.dat[,vars2],
                   grouping = ord.dat$EcoZone3,
                   sim.method = "kulczynski")|>
  mutate(R1 = factor(sapply(strsplit(pairs,".vs."),"[",1),levels=unique(as.character(ord.dat$EcoZone3))),
         R2 = factor(sapply(strsplit(pairs,".vs."),"[",2),levels=unique(as.character(ord.dat$EcoZone3))),
         R_pval = ifelse(round(anosimR,3)==0,
                         paste0(format(round(anosimR,4),scientific = F),"\n(",ifelse(p.adj<=(0.01),"<0.01",format(p.adj,digits=3)),")"),
                         paste0(round(anosimR,3),"\n(",ifelse(p.adj<=(0.01),"<0.01",format(p.adj,digits=3)),")"))
  )
pwrslt
pwrslt.long <- rbind(data.frame(R1 ="Global",
                             R2="Global",
                             anosimR=ano$statistic,
                             p.adj=ano$signif),
                  pwrslt[,c("R1","R2","anosimR","p.adj")])

pwrslt.long|>
  flextable()|>
  merge_h()|>
  colformat_double(j="anosimR",i=~round(anosimR,3)==0,digits=4)|>
  colformat_double(j="anosimR",i=~round(anosimR,3)!=0,digits=3)|>
  # colformat_double(j=3,digits=4)|>
  compose(j="p.adj",i=~p.adj<(0.05),value=as_paragraph("< 0.05"))|>
  compose(j="p.adj",i=~p.adj<=(0.01),value=as_paragraph("< 0.01"))|>
  compose(j="p.adj",i=~p.adj>0.05,value=as_paragraph(format(round(p.adj,2),digits=2)))|>
  italic(j="p.adj",i=~p.adj<0.05)|>bold(j="p.adj",i=~p.adj<0.05)|>
  flextable::width(width=c(1,1,1,0.75))|>
  padding(padding=1.25,part="all")|>
  flextable::align(j=1:4,align="center",part="all")|>
  set_header_labels(
    "R1"="Value 1",
    "R2"="Value 2",
    "anosimR"="ANOSIM R",
    "p.adj"="\u03C1-value"
  )|>
  font(fontname="Times New Roman",part="all")|>
  fontsize(size=12,part="all")|>
  bold(part="header") # |>print("docx")

pwrslt.xtab  <-  pwrslt|>
  dcast(R1~R2,value.var = "R_pval",fun.aggregate = function(x) x[1])
colmax <- ifelse(is.na(pwrslt.xtab),"lightgrey",NA)

pwrslt.xtab|>
  flextable()|>
  bg(bg=colmax)|>
  flextable::align(j=2:5,align="center",part="all")|>
  compose(i=1,j=1,as_paragraph("Nearshore"))|>
  compose(i=2,j=1,as_paragraph("Pelagic"))|>
  compose(i=3,j=1,as_paragraph("Littoral North"))|>
  compose(i=4,j=1,as_paragraph("Littoral South"))|>
  flextable::width(width=rep(1,5))|>
  padding(padding=1.25,part="all")|>
  hline(border=officer::fp_border(color="white"))|>
  vline(border=officer::fp_border(color="white"))|>
  set_header_labels(
    "R1"="",
    "pelagic"="Pelagic",
    "Littoral_North"="Littoral North",
    "Littoral_South"="Littoral South",
    "Littoral_West"="Littoral West"
  )

## Chla spatio-temporal GAM ------------------------------------------------
ctrl = gam.control(trace = T)

WQ.dat.xtab2 <- subset(WQ.dat.xtab,zone!="struct"&CY%in%seq(2008,2023,1)&scrn.val==1)
WQ.dat.xtab2$decMonth <- dec.month(WQ.dat.xtab2$Date.EST)
WQ.dat.xtab2$EcoZone3 <- with(WQ.dat.xtab2,as.factor(ifelse(EcoZone2 == "Pelagic_North","pelagic",as.character(EcoZone2))))
WQ.dat.xtab2$decCY <- lubridate::decimal_date(WQ.dat.xtab2$Date.EST)

yr.k  <-  16
mon.k <- 50
utm.k <- 60
stg.k <- 15

chla.m3<-bam(Chla ~
               s(STG29,k=15) +
               s(decMonth,bs="cc",k=mon.k) +
               s(CY,k=yr.k) +
               s(UTMX,UTMY,bs="ds",m=c(1,0.5),k=utm.k) + 
               ti(CY, decMonth,bs=c("tp","cc"),k=c(yr.k,mon.k+20)) +
               ti(UTMX,UTMY,CY,d=c(2,1),bs=c("ds","tp"),m = list(c(1, 0.5), NA),k=c(utm.k,yr.k)) +
               ti(UTMX,UTMY,decMonth,d=c(2,1),bs=c("ds","cc"), m = list(c(1, 0.5), NA),k=c(utm.k,mon.k+30)),
             data=WQ.dat.xtab2,discrete = T,nthreads=c(6,1),
             family=gaussian(link="log"),
             control = ctrl)

chla.m3.sum <- summary(chla.m3);
chla.m3.sum

layout(matrix(1:4,2,2,byrow=T));gam.check(chla.m3,pch=21,col="lightblue",bg="grey")
dev.off();plot(chla.m3,pages=1)

res  <- recalculateResiduals(simulateResiduals(chla.m3), group = na.omit(WQ.dat.xtab2[,c(names(chla.m3$model),"Date.EST")])$Date.EST)
testTemporalAutocorrelation(res,time=unique(na.omit(WQ.dat.xtab2[,c(names(chla.m3$model),"Date.EST")])$Date.EST))

range(WQ.dat.xtab2$Chla,na.rm=T); # range of observed values
range(chla.m3$fitted.values,na.rm=T);# range of predicted (fitted) values


# individual term residual check
pred.org <- (predict(chla.m3,type="terms"))
chla.m3.partial.resids<-pred.org+residuals(chla.m3)
ncol(chla.m3.partial.resids)
layout(matrix(1:8,2:4))
for(i in 1:ncol(chla.m3.partial.resids)){hist(chla.m3.partial.resids[,i],main=smooths(chla.m3)[i])}  

acf(residuals(chla.m3))
pacf(residuals(chla.m3))

draw(chla.m3)

### Chla GAM devriative -----------------------------------------------------
reg.ext <- extent(LOK)

smooths(chla.m3)
var.vals <- c("stg","decMon","CY","UTM","CYdecMon","CYUTM","decMonUTM")
chla.m3.stg.d <- dev.data.val(chla.m3,1,n=400,var.names=var.vals)
chla.m3.CY.d <- dev.data.val(chla.m3,3,n=400,var.names=var.vals)
chla.m3.decMon.d <- dev.data.val(chla.m3,2,n=400,var.names=var.vals)

range(subset(chla.m3.stg.d, is.na(dsig.incr)==F)$STG29,na.rm=T)
range(subset(chla.m3.stg.d, is.na(dsig.decr)==F)$STG29,na.rm=T)

## Phyco GAM Data -----------------------------------------------------------
LOK.PC.sites.coords <- data.frame(
  SITE = c("POLESOUT1", "POLESOUT3", "L001", "L006", "LZ40", "L005"), 
  UTMX = c(510009.358377801, 513237.031031199, 520473.178839964, 
           521579.963180089, 520952.59610862, 502741.140514964),
  UTMY = c(2989284.20365265, 2985712.5026207, 3001835.22052724, 
                       2966800.10696748, 2975577.34428215, 2981642.49996558)
)
PC.dat.xtab=merge(PC.dat.xtab,LOK.PC.sites.coords,"SITE")

##
PC.dat.xtab2=subset(PC.dat.xtab,is.na(PU)==F&is.na(Region)==F)
PC.dat.xtab2$Region.f=as.factor(PC.dat.xtab2$Region)
PC.dat.xtab2

ggplot(PC.dat.xtab2,aes(x = DOY,y = PU,color = as.factor(WY)))+
  geom_point()+
  facet_wrap(~SITE)

dev.off();plot(PU~DOY,PC.dat.xtab2,type="n")
points(PU~DOY,subset(PC.dat.xtab2,SITE=="L001"&WY==2021),
       pch=21,bg=ifelse(DOY<20,"red","grey"),cex=1.25)
points(PU~DOY,subset(PC.dat.xtab2,SITE=="L001"&WY==2021),
       pch=19,col=ifelse(PU>800,NA,"grey"),cex=1.25)

PC.dat.xtab2[PC.dat.xtab2$SITE=="L001"&PC.dat.xtab2$WY==2021&PC.dat.xtab2$DOY<20,"PU"] <- NA

### PC GAM -----------------------------------------------------------------
set.seed(123)
PC.m0 <- bam(PU~
            s(STG29,bs="tp",k=20)+
            s(DOY,bs="cc",k=50)+
            s(UTMY,UTMX,bs="ds",m=c(1,0.5),k=20)+
            ti(DOY,decCY,k=c(20,10))+
            s(decCY,k=20)+
            ti(UTMX,UTMY,decCY,d=c(2,1),bs=c("ds","tp"),m = list(c(1, 0.5), NA),k=c(15,40)),
          data=PC.dat.xtab2,
          discrete = T,nthreads=c(10,1),family= tw(),
          control = ctrl)

layout(matrix(1:4,2,2,byrow=T));gam.check(PC.m0,pch=21,col="lightblue",bg="grey");
PC.m0.sum <-  summary(PC.m0); 
PC.m0.sum

dev.off();plot(PC.m0,pages=1)

res  <- recalculateResiduals(simulateResiduals(PC.m0), group = na.omit(WQ.dat.xtab2[,c(names(PC.m0$model),"Date.EST")])$Date.EST)
testTemporalAutocorrelation(res,time=unique(na.omit(WQ.dat.xtab2[,c(names(PC.m0$model),"Date.EST")])$Date.EST))

range(WQ.dat.xtab2$Chla,na.rm=T); # range of observed values
range(PC.m0$fitted.values,na.rm=T);# range of predicted (fitted) values

# individual term residual check
pred.org <- (predict(PC.m0,type="terms"))
PC.m0.partial.resids<-pred.org+residuals(PC.m0)
ncol(PC.m0.partial.resids)
layout(matrix(1:8,2:4))
for(i in 1:ncol(PC.m0.partial.resids)){hist(PC.m0.partial.resids[,i],main=smooths(PC.m0)[i])}  

acf(residuals(PC.m0))
pacf(residuals(PC.m0))

draw(PC.m0)

### PC Derivative analysis -----------------------------------------------------

dev.vars <- c("stg","DOY","UTM","decCY_DOY","decCY","UTM_decCY")
PC.m.stg.d <- dev.data.val(PC.m0,1,n=400,var.names=dev.vars)
PC.m.DOY.d <- dev.data.val(PC.m0,2,n=400,var.names=dev.vars)
PC.m.decCY.d <- dev.data.val(PC.m0,5,n=400,var.names=dev.vars)

## Interogating the results
pracma::findpeaks(PC.m.DOY.d$fit.DOY)
as.Date(212,origin=as.Date("2023-01-01"))
as.Date(mean(c(335,367)),origin=as.Date("2023-01-01"))

range(PC.m.decCY.d$dsig.incr,na.rm=T)
range(subset(PC.m.decCY.d,is.na(dsig.incr)==F)$decCY,na.rm=T)|>
  lubridate::date_decimal()

dec.events <-  with(PC.m.stg.d,rle(ifelse(is.na(dsig.decr)==T,0,1)))
events <-  rep(dec.events$values==1,dec.events$lengths)
PC.m.stg.d[events,]

PC.m.stg.d$event_id <-  rep(seq_along(dec.events$values),dec.events$lengths)
tmp  <-  ddply(PC.m.stg.d,"event_id",summarise,
            N.decline = N.obs(dsig.decr),
            min.val = min(STG29,na.rm=T),
            max.val = max(STG29,na.rm=T),
            min.fit = min(dsig.decr,na.rm=T),
            max.fit = max(dsig.decr,na.rm=T))

plot(min.val~event_id,tmp)
points(max.val~event_id,tmp)

### Chla & PC Response curves ---------------------------------------------------------
## Average response curves for chlorophyll and Phycocyanin models
## curve that limits the time frame consistent with Chla 
smooths(chla.m3)
crit.t <- qt(0.025, df.residual(chla.m3), lower.tail = FALSE)
chla.stg.pdat <- expand.grid(
  STG29=seq(chla.m3$var.summary$STG29[1],chla.m3$var.summary$STG29[3],length.out = 100),
  decMonth=chla.m3$var.summary$decMonth[2],
  CY=chla.m3$var.summary$CY[2],
  UTMX=chla.m3$var.summary$UTMX[2],
  UTMY=chla.m3$var.summary$UTMY[2]
)
Chl.stg.resp.dat <- predict(chla.m3,chla.stg.pdat,type="response",se=T)|>
  mutate(
    resp.UCI=(fit+(crit.t*se.fit)),
    resp.LCI=(fit-(crit.t*se.fit)),
    fit=(fit)
  )|>
  as.data.frame()
chla.stg.pdat <- cbind(chla.stg.pdat,Chl.stg.resp.dat)

crit.t <- qt(0.025, df.residual(PC.m0), lower.tail = FALSE)
PC.stg.pdat <- expand.grid(
  STG29=seq(PC.m0$var.summary$STG29[1],PC.m0$var.summary$STG29[3],length.out = 100),
  DOY=PC.m0$var.summary$DOY[2],
  UTMX=PC.m0$var.summary$UTMX[2],
  UTMY=PC.m0$var.summary$UTMY[2],
  decCY = PC.m0$var.summary$decCY[2]
)
PC.stg.resp.dat <- predict(PC.m0,PC.stg.pdat,type="response",se=T)|>
  mutate(
    resp.UCI=(fit+(crit.t*se.fit)),
    resp.LCI=(fit-(crit.t*se.fit)),
    fit=(fit)
  )|>
  as.data.frame()
PC.stg.pdat=cbind(PC.stg.pdat,PC.stg.resp.dat)


## Frequency Analysis - Mixed Models Approach (LMM) --------------------------
## see Walker (2020) for more information - https://swampthingecology.org/LOK_AlgaeEval/LOK_AlgaeEvalTool.html

lok.stg.sea.mean2 <- ddply(lok.stg,c("WY","Season"),summarise,
                        mean.STG=mean(STG29,na.rm=T),
                        mean.recess=mean(recess_30day,na.rm=T),
                        mean.delta.min=mean(delta.min.stg,na.rm=T),
                        N.val=N.obs(STG29))

monCY.flow <- ddply(flow.mon.sum,c("CY","month","WY"),summarise,
                    inflow.Q=sum(cfs.to.km3d(Inflow),na.rm=T),
                    outflow.Q=sum(cfs.to.km3d(Outflow),na.rm=T))
monCY.flow$Qin.3mon.sum <- zoo::rollapply(monCY.flow$inflow.Q,
                                          3,sum,na.rm=T,align="right",fill=NA)
monCY.flow$Qin.6mon.sum <- zoo::rollapply(monCY.flow$inflow.Q,
                                          6,sum,na.rm=T,align="right",fill=NA)
monCY.flow$Qin.12mon.sum <- zoo::rollapply(monCY.flow$inflow.Q,
                                           12,sum,na.rm=T,align="right",fill=NA)
monCY.flow <- merge(monCY.flow,month.season,"month",all.x=T,sort=FALSE)
monCY.flow <- monCY.flow[order(monCY.flow$CY,monCY.flow$month),]

lok.q.sea.mean2 <- ddply(monCY.flow,c("WY","Season"),summarise,
                      mean.Qin=mean(inflow.Q,na.rm=T),
                      total.Qin=sum(inflow.Q,na.rm=T),
                      mean.3m.sum=mean(Qin.3mon.sum,na.rm=T),
                      mean.6m.sum=mean(Qin.6mon.sum,na.rm=T),
                      mean.12m.sum=mean(Qin.12mon.sum,na.rm=T))

site.sea.mean <- ddply(subset(WQ.dat.xtab,!(Station.ID%in%structure.sites)),
                    c("Station.ID","Season","EcoZone3","WY","UTMX","UTMY"),
                    summarise,
                    # N.sites = N.obs(unique(Station.ID)),
                    N.Chla.vals = N.obs(Chla),
                    mean.val = mean(Chla,na.rm=T),
                    DIN = mean(DIN,na.rm=T),
                    SRP = mean(SRP,na.rm=T),
                    Temp = mean(TempC,na.rm=T),
                    f20 = sum(Chla>20,na.rm=T),
                    NB20 = N.Chla.vals-f20,
                    f40 = sum(Chla>40,na.rm=T),
                    NB40 = N.Chla.vals-f40)
site.sea.mean.sum <- ddply(site.sea.mean,c("Season","EcoZone3","WY"),summarise,
                          N.vals = N.obs(mean.val),
                          mean.val = mean(mean.val,na.rm=T),
                          mean.f20 = mean(f20/N.Chla.vals,na.rm=T),
                          mean.f40 = mean(f40/N.Chla.vals,na.rm=T),
                          sum.f20 = sum(f20,na.rm=T),
                          sum.NB20 = sum(NB20,na.rm=T),
                          sum.f40 = sum(f40,na.rm=T),
                          sum.NB40 = sum(NB40,na.rm=T),
                          Tcount.Chla = sum(N.Chla.vals),
                          mean.DIN = mean(DIN,na.rm=T),
                          mean.SRP = mean(SRP, na.rm=T),
                          mean.Temp = mean(Temp,na.rm=T))

melt(site.sea.mean.sum[,c("WY","Season","EcoZone3","mean.val","mean.f20","mean.f40")],
     id.vars = c("WY","Season","EcoZone3"))|>
  merge(data.frame(variable =c("mean.val","mean.f20","mean.f40"),
                    param = c("Chl-a (ug/L)","f20 (prop)","f40 (prop)")))|>
  subset(Season=="Summer")|>
  ggplot(aes(x=WY,y=value,color=EcoZone3))+
  geom_point(shape=19)+
  geom_line(linewidth = 0.75)+
  facet_wrap(param~EcoZone3,nrow=3,scales="free_y")


### Seasonal mods (site season stats)-----------------------------------------
WY.vals.scn <- seq(2000,2023,1)[!(seq(2000,2023,1)%in%subset(WY.periods,maj.Hurr==1)$WY)]

site.sea.mean.sum <- site.sea.mean.sum|>
  merge(lok.stg.sea.mean2,c("WY","Season"))|>
  merge(lok.q.sea.mean2,c("WY","Season"))

site.sea.mean.sum.su <-  subset(site.sea.mean.sum,Season=="Summer"&WY%in%WY.vals.scn)|>
  mutate(EcoZone3.f=as.factor(EcoZone3))
### Mean Chlorophyll --------------------------------------------------------
mod_site.seamean1.stg <- gam(mean.val ~ mean.delta.min  +
                               s(EcoZone3.f,bs="re")+
                               s(EcoZone3.f,mean.delta.min,bs="re"),
                             data = site.sea.mean.sum.su,method = 'REML',family=tw())
stg.chl.me <-  summary(mod_site.seamean1.stg)
stg.chl.me

gratia::variance_comp(mod_site.seamean1.stg); 
gammit::extract_fixed(mod_site.seamean1.stg)

coef(mod_site.seamean1.stg)
# SE sqrt(diag(vcov(mod_site.seamean1.stg)))
layout(matrix(1:4,2,2));gam.check(mod_site.seamean1.stg,pch=21,col="lightblue",bg="grey")
dev.off();testResiduals(simulateResiduals(mod_site.seamean1.stg))

stg.mod.fe  <-  gam(mean.val ~ mean.delta.min,
                 data = site.sea.mean.sum.su,method = 'REML',family=tw())
summary(stg.mod.fe)

#### Cross Validation --------------------------------------------------------
# k-fold approach

set.seed(83)
nfolds = 5
case.folds = rep(1:nfolds, length.out = nrow(site.sea.mean.sum.su)) # divide the cases as evenly as possible
case.folds <- sample(case.folds) # randomly permute the order
# bandwidths <- (1:5)/10 # Evenly space bandwidths from 0.1 to 0.5 (only for npreg in example)
kfold.result=data.frame()
for(fold in 1:nfolds){
  train  <-  site.sea.mean.sum.su[case.folds!=fold,]
  test <- site.sea.mean.sum.su[case.folds==fold,]
  # Fit training
  mod <- gam(mean.val ~ mean.delta.min  +
              s(EcoZone3.f,bs="re")+
              s(EcoZone3.f,mean.delta.min,bs="re"),
            data = train,method = 'REML',family=tw())
  sum.mod <- summary(mod)
  devexpl <- sum.mod$dev.expl
  disp <- sum.mod$dispersion
  
  pred <- predict(mod, newdata = test)
  eval <- model_fit_params2(test$mean.val,exp(pred))
  eval <- cbind(eval,data.frame(test.N = N.obs(test$mean.val),
                              train.N=N.obs(train$mean.val),
                              mod.devexpl = devexpl, mod.disp = disp,fold=fold))
  kfold.result <- rbind(kfold.result,eval)
}
kfold.result

apply(kfold.result,2,mean)

### f20  --------------------------------------------------------------------
site.sea.mean.sum.su$tot.f20 <- with(site.sea.mean.sum.su,sum.f20/Tcount.Chla)

f20.mod.fe <- gam(cbind(sum.f20,sum.NB20) ~ mean.delta.min,
                 data = site.sea.mean.sum.su,method = 'REML',family=binomial(link = "logit"))
summary(f20.mod.fe)

mod_site.sea.f20 <- gam(cbind(sum.f20,sum.NB20) ~ mean.delta.min +
                          s(EcoZone3.f,bs="re")+
                          s(EcoZone3.f,mean.delta.min,bs="re"),
                        data = site.sea.mean.sum.su,method = 'REML',family= binomial(link = "logit"))
f20.chl.me <- summary(mod_site.sea.f20);
f20.chl.me
gratia::variance_comp(mod_site.sea.f20); 

layout(matrix(1:4,2,2));gam.check(mod_site.sea.f20,pch=21,col="lightblue",bg="grey")
dev.off();testResiduals(simulateResiduals(mod_site.sea.f20))

set.seed(83)
nfolds = 5
case.folds <- rep(1:nfolds, length.out = nrow(site.sea.mean.sum.su)) # divide the cases as evenly as possible
case.folds <- sample(case.folds) # randomly permute the order
kfold.f20.result=data.frame()
for(fold in 1:nfolds){
  train <- site.sea.mean.sum.su[case.folds!=fold,]
  test <- site.sea.mean.sum.su[case.folds==fold,]
  # Fit training
  mod <- gam(cbind(sum.f20,sum.NB20) ~ mean.delta.min +
              s(EcoZone3.f,bs="re")+
              s(EcoZone3.f,mean.delta.min,bs="re"),
            data = train,method = 'REML',family= binomial(link = "logit"))
  sum.mod <- summary(mod)
  devexpl <- sum.mod$dev.expl
  disp <- sum.mod$dispersion
  
  pred <- predict(mod, newdata = test ,type="response")
  eval <- model_fit_params2(test$tot.f20,pred)
  eval <- cbind(eval,data.frame(test.N = N.obs(test$sum.f20),
                              train.N=N.obs(train$sum.f20),
                              mod.devexpl = devexpl, mod.disp = disp,fold=fold))
  kfold.f20.result <- rbind(kfold.f20.result,eval)
}
kfold.f20.result
apply(kfold.f20.result,2,mean)

#### f40  --------------------------------------------------------------------
site.sea.mean.sum.su$tot.f40 <-  with(site.sea.mean.sum.su,sum.f40/Tcount.Chla)

f40.mod.fe <-  gam(cbind(sum.f40,sum.NB40) ~ mean.delta.min ,
                 data = site.sea.mean.sum.su,method = 'REML',family=binomial(link = "logit"))
summary(f40.mod.fe)

mod_site.sea.f40 <- gam(cbind(sum.f40,sum.NB40) ~ mean.delta.min +
                          s(EcoZone3.f,bs="re")+
                          s(EcoZone3.f,mean.delta.min,bs="re"),
                        data = site.sea.mean.sum.su,method = 'REML',family= binomial(link = "logit"))
f40.chl.me <-  summary(mod_site.sea.f40);f40.chl.me
gratia::variance_comp(mod_site.sea.f40); 

dev.off();testResiduals(simulateResiduals(mod_site.sea.f40,refit=T));# dispersion test good
dev.off();testResiduals(simulateResiduals(mod_site.sea.f40,refit=F));# dispersion test p<0.05 (higher dispersion rate?)

set.seed(83)
nfolds  <-  5
case.folds  <-  rep(1:nfolds, length.out = nrow(site.sea.mean.sum.su)) # divide the cases as evenly as possible
case.folds <- sample(case.folds) # randomly permute the order
kfold.f40.result=data.frame()
for(fold in 1:nfolds){
  train <-  site.sea.mean.sum.su[case.folds!=fold,]
  test <-  site.sea.mean.sum.su[case.folds==fold,]
  # Fit training
  mod <-  gam(cbind(sum.f40,sum.NB40) ~ mean.delta.min +
              s(EcoZone3.f,bs="re")+
              s(EcoZone3.f,mean.delta.min,bs="re"),
            data = train,method = 'REML',family= binomial(link = "logit"))
  sum.mod  <-  summary(mod)
  devexpl <-  sum.mod$dev.expl
  disp <-  sum.mod$dispersion
  
  pred  <-  predict(mod, newdata = test,type="response")
  eval  <-  model_fit_params2(test$tot.f40,pred)
  eval <- cbind(eval,data.frame(test.N = N.obs(test$sum.f40),
                              train.N=N.obs(train$sum.f40),
                              mod.devexpl = devexpl, mod.disp = disp,fold=fold))
  kfold.f40.result  <- rbind(kfold.f40.result,eval)
}
kfold.f40.result
apply(kfold.f40.result,2,mean)

## WRT - Chla analysis -----------------------------------------------------
idvars <- c("Station.ID", "WY", "Date.EST", "EcoZone3", 
            "CY", "month", "hydro.season", "scrn.val")
val.vars <- c("TKN", "NH4", "NOx", "TN", "TP", "SRP", "Chla", "TempC", "SD", "DIN") 
idvars[!(idvars%in%names(WQ.dat.xtab))]; #double check varables are in data.frame name
val.vars[!(val.vars%in%names(WQ.dat.xtab))]; #double check varables are in data.frame name
lake.wq.melt <- melt(WQ.dat.xtab[,c(idvars,val.vars)],id.vars=idvars)
lake.wq.melt <- merge(lake.wq.melt,
                      month.season,
                      "month")

lake.wq.melt.site <- subset(lake.wq.melt,scrn.val==1)|>
  ddply(c("Station.ID", "WY","EcoZone3", "CY", "month","Season","variable"),
        summarise,mean.val = mean(value,na.rm=T))

flow.TP.sea.mean <- melt(flow.TP.mon.CY,id.vars = c("CY","month","WY","Date.monCY"))|>
  merge(month.season,"month")|>
  dcast(WY+Season~variable,value.var="value",mean,na.rm=T)

lake.wq.mon.xtab2 <- dcast(lake.wq.melt.site,EcoZone3+WY+Season~variable,
                           value.var = "mean.val",mean,na.rm=T)|>
  merge(flow.TP.sea.mean,c("WY","Season"))
lake.wq.mon.xtab2$EcoZone3.f <- as.factor(lake.wq.mon.xtab2$EcoZone3)

ggplot(lake.wq.mon.xtab2,aes(x=WY,y=Chla,color=EcoZone3))+
  geom_point()+
  geom_smooth(method='loess')+
  scale_y_continuous(trans = scales::log_trans(),breaks=c(1,10,100,500),labels=scales::comma) +
  facet_wrap(EcoZone3~Season,nrow=5)

## LOK HABAM ---------------------------------------------------------------
# Summer Only Model
WYs.val <- seq(2000,2023,1)
WYs.val <-WYs.val[!(WYs.val%in%subset(WY.periods,maj.Hurr==1)$WY)]
test.WYs.val <- seq(1977,1999,1)
test.WYs.val <-test.WYs.val[!(test.WYs.val%in%subset(WY.periods,maj.Hurr==1)$WY)]

wq.hydro.sum.mean <- subset(lake.wq.mon.xtab2,WY%in%WYs.val&Season=="Summer")
wq.hydro.sum.mean.test <- subset(lake.wq.mon.xtab2,WY%in%test.WYs.val&Season=="Summer")

vars <- c("Chla","TP","DIN", "EcoZone3.f","mean.depth.m","Inflow.Q.m3","mean.vol.km3","WRT.yr")

### Chla --------------------------------------------------------------------
Chla_HGAM_mod_sum.null <- gam(Chla~1, data = wq.hydro.sum.mean,family=tw(),method="REML")

Chla_HGAM_mod_sum3 <- gam(Chla ~ WRT.yr + TP + DIN + 
                            ti(TP,DIN,by=EcoZone3.f,k=6)+
                            s(mean.depth.m,k=19,by=EcoZone3.f)+
                            s(Inflow.Q.m3,bs="cc",k=12)+
                            s(EcoZone3.f, bs = "re"), 
                          data = wq.hydro.sum.mean,method="GCV.Cp",
                          family= tw())
Chla_HGAM_mod_sum3 |> AIC()
Chla_HGAM_mod_sum3 |> summary()
Chla_HGAM_mod_sum3 |> plot(pages=1)
layout(matrix(1:4,2,2));gam.check(Chla_HGAM_mod_sum3,pch=21,bg="grey");abline(0,1,col="Red")
testResiduals(simulateResiduals(Chla_HGAM_mod_sum3))
LOOCV_fun(Chla_HGAM_mod_sum3); plot(pred~actual,LOOCV_fun(Chla_HGAM_mod_sum3));abline(0,1,col="red")

Chla.mod.test <- predict(Chla_HGAM_mod_sum3,wq.hydro.sum.mean.test,type="response")
plot(wq.hydro.sum.mean.test$Chla~Chla.mod.test);abline(0,1,col="red")
model_fit_params(wq.hydro.sum.mean.test$Chla,Chla.mod.test)


### TP Model ----------------------------------------------------------------
TP_HGAM_mod_sum2 <- gam(TP ~ 
                          s(mean.vol.km3,by=EcoZone3.f,k=10,bs=c("tp"))+
                          s(Outflow.Q.m3,k=10) +
                          s(EcoZone3.f, bs = "re"), 
                        data = wq.hydro.sum.mean,method="GCV.Cp",
                        family= tw())#
summary(TP_HGAM_mod_sum2)
LOOCV_fun(TP_HGAM_mod_sum2); plot(pred~actual,LOOCV_fun(TP_HGAM_mod_sum2));abline(0,1,col="red")

summary(TP_HGAM_mod_sum2)
plot(TP_HGAM_mod_sum2,page=1)
layout(matrix(1:4,2,2));gam.check(TP_HGAM_mod_sum2,pch=21,bg="grey");abline(0,1,col="red")
testResiduals(simulateResiduals(TP_HGAM_mod_sum2))

TP.mod.test <- predict(TP_HGAM_mod_sum2,wq.hydro.sum.mean.test,type="response")
plot(wq.hydro.sum.mean.test$TP~TP.mod.test);abline(0,1)
model_fit_params(wq.hydro.sum.mean.test$TP,TP.mod.test)


### DIN Model ---------------------------------------------------------------
DIN.gam.fam2 <- tw() # Gamma(link="log") #   
DIN_HGAM_mod_sum2 <- gam(DIN ~ Outflow.Q.m3 + mean.vol.km3 +
                           s(WRT.yr,by=EcoZone3.f,k=10)+
                           # s(mean.depth.m,bs="cc",k=19)+
                           # s(mean.vol.km3,by=EcoZone3.f,k=10)+
                           s(EcoZone3.f, bs = "re"), 
                         data = wq.hydro.sum.mean,method="GCV.Cp", #"REML"
                         family=DIN.gam.fam2)
summary(DIN_HGAM_mod_sum2)
LOOCV_fun(DIN_HGAM_mod_sum2); plot(pred~actual,LOOCV_fun(DIN_HGAM_mod_sum2));abline(0,1,col="red")

concurvity(DIN_HGAM_mod_sum2)
plot(DIN_HGAM_mod_sum2,page=1)
layout(matrix(1:4,2,2));gam.check(DIN_HGAM_mod_sum2,pch=21,bg="grey");abline(0,1,col="red")
testResiduals(simulateResiduals(DIN_HGAM_mod_sum2))

DIN.mod.test <- predict(DIN_HGAM_mod_sum2,wq.hydro.sum.mean.test,type="response")
plot(wq.hydro.sum.mean.test$DIN~DIN.mod.test);abline(0,1)
model_fit_params(wq.hydro.sum.mean.test$DIN,DIN.mod.test)

### LOOCV -------------------------------------------------------------------
# Leave-One-Out Cross-Validation

Chla_HGAM3.LOOCV <- LOOCV_fun(Chla_HGAM_mod_sum3)
Chla_HGAM3.LOOCV
plot(Chla_HGAM3.LOOCV$actual, Chla_HGAM3.LOOCV$pred);abline(0,1,col="red")

TP_HGAM2.LOOCV <- LOOCV_fun(TP_HGAM_mod_sum2)
TP_HGAM2.LOOCV
plot(TP_HGAM2.LOOCV$actual, TP_HGAM2.LOOCV$pred);abline(0,1,col="red")

DIN_HGAM2.LOOCV <- LOOCV_fun(DIN_HGAM_mod_sum2);DIN_HGAM2.LOOCV
plot(DIN_HGAM2.LOOCV$actual, DIN_HGAM2.LOOCV$pred);abline(0,1,col="red")

### TSI exploration ---------------------------------------------------------
lake.wq.mon.xtab2$TSI_TP <- with(lake.wq.mon.xtab2,14.42*log(TP*1000)+4.14)
lake.wq.mon.xtab2$TSI_TN <- with(lake.wq.mon.xtab2,54.45+14.43*log(TN))
lake.wq.mon.xtab2$TSI_SD <- with(lake.wq.mon.xtab2,60-14.41*log(SD))
lake.wq.mon.xtab2$TSI_Chla <- with(lake.wq.mon.xtab2,9.81*log(Chla+30.6))

lake.wq.mon.xtab2$TSI_TNTP_diff <- with(lake.wq.mon.xtab2,TSI_TN-TSI_TP)
lake.wq.mon.xtab2$TSI_ChlTP_diff <- with(lake.wq.mon.xtab2,TSI_Chla-TSI_TP)
lake.wq.mon.xtab2$TSI_ChlSD_diff <- with(lake.wq.mon.xtab2,TSI_Chla-TSI_SD)
head(lake.wq.mon.xtab2)

ggplot(subset(lake.wq.mon.xtab2,Season=="Summer"),
       aes(x=WY,y=TSI_TNTP_diff,color=EcoZone3.f))+
  geom_line()+
  geom_hline(yintercept=0,color="red")+
  labs(y=expression("TSI"[" TN"] * " - TSI"[" TP"]),
       subtitle = "Summer mean values",
       title = "Lake Okeechobee Ecological Zones",
       color="Zones")+
  facet_wrap(~EcoZone3.f)


## RSM Data ----------------------------------------------------------------
LOSOM.alts  <-  c("NA25f","PA25")
LOCAR.alts <-  c("PA_FWOLL","LCR1","ECB23L")

### Data aggregation and calcs ----------------------------------------------
min.stg <- 11.5
all.alts <- c(LOSOM.alts,LOCAR.alts)

LOK.RSM.stg <- read.csv(paste0(data.path,"RSMBN_LOK.csv"))|>
  mutate(Date=date.fun(Date),
         CY=as.numeric(format(Date,"%Y")),
         month=as.numeric(format(Date,"%m")),
         DOY=as.numeric(format(Date,"%j")),
         WY=WY(Date,"Fed"),
         FLWY=WY(Date),
         hydro.season=FL.Hydroseason(Date),
         delta.min.stg=pmax(STAGE-min.stg,0,na.rm=T),
         Alt = factor(Alt,levels=all.alts),
         STAGE.m = ft.to.m(STAGE),
         volume.acft = get_LOK_volume(STAGE),
         area.acres = get_LOK_area(STAGE),
         depth = get_LOK_depth(STAGE))

LOK.RSM.stg <-  LOK.RSM.stg[order(LOK.RSM.stg$Date,LOK.RSM.stg$Alt),]

RSM.mon.vol <- ddply(LOK.RSM.stg,c("Alt","CY","month","WY"),summarise,
                  mean.vol.km3=mean(acft.to.m3(volume.acft)/1e9,na.rm=T),
                  mean.area.km2=mean(acres.to.m2(area.acres)/1e6,na.rm=T),
                  mean.depth.m = mean(ft.to.m(depth),na.rm=T),
                  mean.STG = mean(STAGE,na.rm=T))
RSM.mon.vol$Date.monCY <- with(RSM.mon.vol,date.fun(paste(CY,month,"01",sep="-")))

vars = c("Date.monCY","Alt", "mean.vol.km3", "mean.area.km2","mean.depth.m", "mean.STG" )
melt(RSM.mon.vol[,vars],id.vars = vars[1:2])|>
  ggplot(aes(x=Date.monCY,y=value,color=Alt))+
  geom_line(linewidth = 0.75)+
  facet_wrap(~variable,scales = "free_y",ncol=1)+
  labs(
    x = "Date",
    y = "Value"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Customize panel labels
    legend.position = "bottom"  # Adjust legend position
  )


### RSM Inflow/Outflow  
## LOSOM
Q.inflow.sites <- c('S65E','FEC','TOTAL_ISTOK','S77BF','S4BP','S3','S2',
                    'C12ABP','C12BP',"C10BP",'C4ABP','S236','P5WPS','S308BF',
                    'TCNSQ','S154','S135'); 
Q.outflow.sites <-c('NELKSH_WS_QWS','NLKSH_WS_QWS','S77','S4_WS','S271','C12A','C12','C10',
                    'S352',"S351",'C4A','C3','S354','BRIGHTON_WS','S308')

## LOCAR - adjusted for changes in Indian Prairie Basin Refinement
Q.inflow.sites <- c(Q.inflow.sites,"S84","C40C41LSCOUT")
Q.outflow.sites <-c(Q.outflow.sites,'G207G208')

LOK.wb.sites <- rbind(
  data.frame(SITE = Q.inflow.sites,dir = "inflow"),
  data.frame(SITE = Q.outflow.sites,dir = "outflow")
)

LOK.RSMBN.Q <- read.csv(paste0(data.path,"RSMBN_LOK_Q.csv"))
LOK.RSMBN.Q <- LOK.RSMBN.Q|>
  mutate(Date=date.fun(Date),
         CY=as.numeric(format(Date,"%Y")),
         month=as.numeric(format(Date,"%m")),
         DOY=as.numeric(format(Date,"%j")),
         WY=WY(Date,"Fed"),
         FLWY=WY(Date),
         hydro.season=FL.Hydroseason(Date),
         Alt = factor(Alt,levels=all.alts))|>
  merge(month.season,"month")

unique(LOK.RSMBN.Q$proj)
unique(LOK.RSMBN.Q$Alt)

RSM.mon.Q <- dcast(LOK.RSMBN.Q,CY+month+WY+Alt~dir,value.var = "FLOW",
                   fun.aggregate = function(x) sum(cfs.to.m3d(x),na.rm=T))
RSM.mon.Q$Date.monCY <- with(RSM.mon.Q,date.fun(paste(CY,month,"01",sep="-")))

names(RSM.mon.Q)

RSM.Qdat.names <- names(RSM.mon.Q)# c(names(RSM.mon.Q)[!(names(RSM.mon.Q)%in%c("inflow","outflow"))],paste(c("Inflow","Outflow"),"Q.m3",sep="."))
RSM.Qdat.names <- gsub("inflow","Inflow.Q.m3",RSM.Qdat.names)
RSM.Qdat.names <- gsub("outflow","Outflow.Q.m3",RSM.Qdat.names)
colnames(RSM.mon.Q) <- RSM.Qdat.names
head(RSM.mon.Q)

RSM.mon.Q.vol <- merge(RSM.mon.Q,RSM.mon.vol,c("CY", "month","WY","Date.monCY",'Alt'),all.x=T)|>
  mutate(WRT.yr = mean.vol.km3/(Outflow.Q.m3/1e9)/12)|>
  merge(month.season,"month",all.x=T)
unique(lake.wq.mon.xtab2$EcoZone3.f)

### Stage Duration Curves ---------------------------------------------------
compare <- function(x, y) {
  n <- length(x); m <- length(y)
  w <- c(x, y)
  o <- order(w)
  z <- cumsum(ifelse(o <= n, m, -n))
  i <- which.max(abs(z))
  w[o[i]]
}

NA25.ecdf.dat  <-  ecdf_fun(subset(LOK.RSM.stg,Alt==all.alts[1])$STAGE.m)|>mutate(proportion = 1-proportion)
PA25.ecdf.dat <-  ecdf_fun(subset(LOK.RSM.stg,Alt==all.alts[2])$STAGE.m)|>mutate(proportion = 1-proportion)
FWOLL.ecdf.dat <-  ecdf_fun(subset(LOK.RSM.stg,Alt==all.alts[3])$STAGE.m)|>mutate(proportion = 1-proportion)
LCR1.ecdf.dat <-  ecdf_fun(subset(LOK.RSM.stg,Alt==all.alts[4])$STAGE.m)|>mutate(proportion = 1-proportion)

LOSOM.comp  <- compare(subset(LOK.RSM.stg,Alt==all.alts[1])$STAGE.m,
                     subset(LOK.RSM.stg,Alt==all.alts[2])$STAGE.m)
LOCAR.comp  <- compare(subset(LOK.RSM.stg,Alt==all.alts[3])$STAGE.m,
                     subset(LOK.RSM.stg,Alt==all.alts[4])$STAGE.m)

den.NA25 <- density(x = subset(LOK.RSM.stg,Alt==all.alts[1])$STAGE.m)
den.PA25 <- density(x = subset(LOK.RSM.stg,Alt==all.alts[2])$STAGE.m)
den.FWOLL <- density(x = subset(LOK.RSM.stg,Alt==all.alts[3])$STAGE.m)
den.LCR1 <- density(x = subset(LOK.RSM.stg,Alt==all.alts[4])$STAGE.m)


LOK.RSM.stg.sea <- ddply(LOK.RSM.stg,c("WY","Season","Alt","proj"),summarise,
                      mean.STG=mean(STAGE,na.rm=T),
                      mean.delta.min=mean(delta.min.stg,na.rm=T),
                      median.STG=median(STAGE,na.rm=T),
                      median.delta.min=median(delta.min.stg,na.rm=T),
                      N.val=N.obs(STAGE))|>
  mutate(Alt=factor(Alt,levels = c(LOSOM.alts,LOCAR.alts)),
         proj = as.factor(proj))

LOK.RSM.stg.sea.su <-  subset(LOK.RSM.stg.sea,Season == "Summer")

### LMM Model predictions -------------------------------------------------------
exp.dat  <-  expand.grid(
  WY = unique(LOK.RSM.stg.sea.su$WY),
  Alt = unique(LOK.RSM.stg.sea.su$Alt),
  EcoZone3.f = as.factor(c("Littoral_North", "Littoral_South", "Littoral_West", "nearshore", "pelagic"))
)
LOK.RSM.stg.sea.su  <-  merge(LOK.RSM.stg.sea.su,exp.dat,c("WY","Alt"))          
LOK.RSM.stg.sea.su<- LOK.RSM.stg.sea.su[order(LOK.RSM.stg.sea.su$Alt,LOK.RSM.stg.sea.su$EcoZone3.f,LOK.RSM.stg.sea.su$WY),]

LOK.RSM.stg.sea.su.pred<- LOK.RSM.stg.sea.su|>
  mutate(
    chla = exp(predict(mod_site.seamean1.stg,newdata = LOK.RSM.stg.sea.su)),
    f20 = predict(mod_site.sea.f20,newdata = LOK.RSM.stg.sea.su,type="response"),
    f40 = predict(mod_site.sea.f40,newdata = LOK.RSM.stg.sea.su,type="response"),
    chla.fxf = exp(predict(mod_site.seamean1.stg,newdata = LOK.RSM.stg.sea.su,
                           exclude=c("s(EcoZone3.f)", "s(EcoZone3.f,mean.delta.min)"))),
    f20.fxf = predict(mod_site.sea.f20,newdata = LOK.RSM.stg.sea.su,type="response",
                      exclude=c("s(EcoZone3.f)", "s(EcoZone3.f,mean.delta.min)")),
    f40.fxf = predict(mod_site.sea.f40,newdata = LOK.RSM.stg.sea.su,type="response",
                      exclude=c("s(EcoZone3.f)", "s(EcoZone3.f,mean.delta.min)")),
    Alt = factor(Alt, levels= c("NA25f", "PA25", "PA_FWOLL", "LCR1")) 
  )

fxf.vals  <-  ddply(LOK.RSM.stg.sea.su.pred,c("WY","Alt"),summarise,
                 chla.fxf = mean(chla.fxf),
                 f20.fxf = mean(f20.fxf),
                 f40.fxf = mean(f40.fxf))

range(LOK.RSM.stg.sea.su.pred$chla)
range(fxf.vals$chla.fxf)

range(LOK.RSM.stg.sea.su.pred$f20)
range(fxf.vals$f20.fxf)

range(LOK.RSM.stg.sea.su.pred$f40)
range(fxf.vals$f40.fxf)


### Comparisons -------------------------------------------------------------
alt.combo <- combn(all.alts, 2)
wc.pairwise.chla <- data.frame()
for(i in 1:ncol(alt.combo)){
  tmp <- wilcox.test(chla.fxf~Alt,
                    # alternative ="greater",
                    subset(fxf.vals,Alt%in%alt.combo[,i]),
                    paired = T)
  
  rslt <- data.frame(alt1 = alt.combo[1,i], alt2 = alt.combo[2,i],
                    comparison = paste(alt.combo[1,i],alt.combo[2,i],sep="-"),
                    stat = as.numeric(tmp$statistic),
                    pval = as.numeric(tmp$p.value)
  )
  wc.pairwise.chla <- rbind(wc.pairwise.chla,rslt)
}
wc.pairwise.chla$pval.adj <- p.adjust(wc.pairwise.chla$pval,"holm",n=6)
wc.pairwise.chla$param <- "chla"
wc.pairwise.chla_lts <- rcompanion::cldList(pval.adj ~ comparison,data=wc.pairwise.chla,threshold = 0.05)


wc.pairwise.f20 <- data.frame()
for(i in 1:ncol(alt.combo)){
  tmp <- wilcox.test(f20.fxf~Alt,
                    subset(fxf.vals,Alt%in%alt.combo[,i]),
                    paired = T)
  
  rslt <- data.frame(alt1 = alt.combo[1,i], alt2 = alt.combo[2,i],
                    comparison = paste(alt.combo[1,i],alt.combo[2,i],sep="-"),
                    stat = as.numeric(tmp$statistic),
                    pval = as.numeric(tmp$p.value)
  )
  wc.pairwise.f20 <- rbind(wc.pairwise.f20,rslt)
}
wc.pairwise.f20$pval.adj <- p.adjust(wc.pairwise.f20$pval,"holm",n=6)
wc.pairwise.f20$param <- "f20"
wc.pairwise.f20_lts <- rcompanion::cldList(pval.adj ~ comparison,data=wc.pairwise.f20,threshold = 0.05)

wc.pairwise.f40 <- data.frame()
for(i in 1:ncol(alt.combo)){
  tmp <- wilcox.test(f40.fxf~Alt,
                    subset(fxf.vals,Alt%in%alt.combo[,i]),
                    paired = T)
  
  rslt <- data.frame(alt1 = alt.combo[1,i], alt2 = alt.combo[2,i],
                    comparison = paste(alt.combo[1,i],alt.combo[2,i],sep="-"),
                    stat = as.numeric(tmp$statistic),
                    pval = as.numeric(tmp$p.value)
  )
  wc.pairwise.f40 <- rbind(wc.pairwise.f40,rslt)
}
wc.pairwise.f40$pval.adj <- p.adjust(wc.pairwise.f40$pval,"holm",n=6)
wc.pairwise.f40$param <- "f40"
wc.pairwise.f40_lts <- rcompanion::cldList(pval.adj ~ comparison,data=wc.pairwise.f40,threshold = 0.05)

## Percent Difference
perdiff.dat = rbind(
  dcast(LOK.RSM.stg.sea.su.pred,EcoZone3.f~Alt,value.var = "chla",mean)|>
    mutate(param = "chla"),
  dcast(LOK.RSM.stg.sea.su.pred,EcoZone3.f~Alt,value.var = "f20",mean)|>
    mutate(param = "f20"),
  dcast(LOK.RSM.stg.sea.su.pred,EcoZone3.f~Alt,value.var = "f40",mean)|>
    mutate(param = "f40")
)|>
  mutate(PerDiff_NA25PA25 =((PA25-NA25f)/NA25f)*100,
         PerDiff_PA25FWOLL =((PA_FWOLL-PA25)/PA25)*100,
         PerDiff_FWOLLLCR1 =((LCR1-PA_FWOLL)/PA_FWOLL)*100,
         PerDiff_PA25LCR1 =((LCR1-PA25)/LCR1)*100)

### BGChem Models (LOK HABAM) -------------------------------------------------

# Calculate summer mean values
vars <- c("month", "CY", "WY", "Date.monCY", "Season", "Alt", "Inflow.Q.m3", "Outflow.Q.m3", 
          "mean.vol.km3", "mean.area.km2", "mean.depth.m", "mean.STG", 
          "WRT.yr")
RSM.sum.Q.vol <- melt(subset(RSM.mon.Q.vol,Season=="Summer")[,vars],id.vars = vars[1:6])|>
  dcast(WY+Season+Alt~variable,
        value.var = "value",mean,na.rm=T)

eco.region.ls <- c("Littoral_North", "Littoral_South", "Littoral_West", "nearshore", "pelagic")
RSM.sum.Q.vol <- lapply(eco.region.ls,function(var){
  RSM.sum.Q.vol$EcoZone3.f <- var
  RSM.sum.Q.vol
})
head(RSM.sum.Q.vol)
RSM.sum.Q.vol <- do.call(rbind,RSM.sum.Q.vol)
RSM.sum.Q.vol$EcoZone3.f <- as.factor(RSM.sum.Q.vol$EcoZone3.f)
RSM.sum.Q.vol$Alt <- factor(RSM.sum.Q.vol$Alt,levels=c("NA25f","PA25","PA_FWOLL","LCR1"))

RSM.sum.Q.vol2 <- RSM.sum.Q.vol; # "clean" version for later

# TP 
TP.rsm.pred <- predict(TP_HGAM_mod_sum2,newdata=RSM.sum.Q.vol,type="response",se.fit=T)
RSM.sum.Q.vol$TP.fit.pred <- TP.rsm.pred$fit|>as.numeric()
RSM.sum.Q.vol$TP.fit.se <- TP.rsm.pred$se.fit|>as.numeric()

# DIN
DIN.rsm.pred <- predict(DIN_HGAM_mod_sum2,newdata=RSM.sum.Q.vol,type="response",se.fit=T)
RSM.sum.Q.vol$DIN.fit.pred <- DIN.rsm.pred$fit|>as.numeric()
RSM.sum.Q.vol$DIN.fit.se <- DIN.rsm.pred$se.fit|>as.numeric()

# Chla
Chla.rsm.pred <- predict(Chla_HGAM_mod_sum3,newdata=RSM.sum.Q.vol,type="response",se.fit=T)
RSM.sum.Q.vol$Chla.fit.pred <- RSM.sum.Q.vol$Chla <- Chla.rsm.pred$fit|>as.numeric()
RSM.sum.Q.vol$Chla.fit.se <- Chla.rsm.pred$se.fit|>as.numeric()

### WQ Scenario results --------------------------------------------------------
RSM.sum.Q.vol2.FWOLL <- subset(RSM.sum.Q.vol2,Alt=="PA_FWOLL") 

TP.rsm.pred <- predict(TP_HGAM_mod_sum2,newdata=RSM.sum.Q.vol2.FWOLL,type="response",se.fit=T)
RSM.sum.Q.vol2.FWOLL$TP.fit.pred <- TP.rsm.pred$fit|>as.numeric()

DIN.rsm.pred <- predict(DIN_HGAM_mod_sum2,newdata=RSM.sum.Q.vol2.FWOLL,type="response",se.fit=T)
RSM.sum.Q.vol2.FWOLL$DIN.fit.pred <- DIN.rsm.pred$fit|>as.numeric()

# Base
RSM.sum.Q.vol2.FWOLL$TP <- RSM.sum.Q.vol2.FWOLL$TP.fit.pred
RSM.sum.Q.vol2.FWOLL$DIN <- RSM.sum.Q.vol2.FWOLL$DIN.fit.pred

Chla.rsm.pred <- predict(Chla_HGAM_mod_sum3,newdata=RSM.sum.Q.vol2.FWOLL,type="response",se.fit=T)
RSM.sum.Q.vol2.FWOLL$Chla.fit.pred <- Chla.rsm.pred$fit|>as.numeric()

# tmp <- aggregate(Chla.fit.pred~EcoZone3.f,RSM.sum.Q.vol2.FWOLL,mean,na.rm=T)|>
#   mutate(scenario="base")

tmp <- ddply(RSM.sum.Q.vol2.FWOLL,c("EcoZone3.f"),summarise,
             mean.val = mean(Chla.fit.pred,na.rm=T),
             SE.val = SE(Chla.fit.pred))|>
  mutate(scenario="base")

# P10_N10
RSM.sum.Q.vol2.FWOLL$TP <- RSM.sum.Q.vol2.FWOLL$TP.fit.pred*0.9
RSM.sum.Q.vol2.FWOLL$DIN <- RSM.sum.Q.vol2.FWOLL$DIN.fit.pred*0.9

Chla.rsm.pred2 <- predict(Chla_HGAM_mod_sum3,newdata=RSM.sum.Q.vol2.FWOLL,type="response",se.fit=T)
RSM.sum.Q.vol2.FWOLL$Chla.fit.pred <- Chla.rsm.pred2$fit|>as.numeric()

tmp <- rbind(tmp,
             ddply(RSM.sum.Q.vol2.FWOLL,c("EcoZone3.f"),summarise,
                   mean.val = mean(Chla.fit.pred,na.rm=T),
                   SE.val = SE(Chla.fit.pred))|>
               mutate(scenario="P10_N10")
)

# P20_N10
RSM.sum.Q.vol2.FWOLL$TP <- RSM.sum.Q.vol2.FWOLL$TP.fit.pred*0.8
RSM.sum.Q.vol2.FWOLL$DIN <- RSM.sum.Q.vol2.FWOLL$DIN.fit.pred*0.9

Chla.rsm.pred2 <- predict(Chla_HGAM_mod_sum3,newdata=RSM.sum.Q.vol2.FWOLL,type="response",se.fit=T)
RSM.sum.Q.vol2.FWOLL$Chla.fit.pred <- Chla.rsm.pred2$fit|>as.numeric()

tmp <- rbind(tmp,
             ddply(RSM.sum.Q.vol2.FWOLL,c("EcoZone3.f"),summarise,
                   mean.val = mean(Chla.fit.pred,na.rm=T),
                   SE.val = SE(Chla.fit.pred))|>
               mutate(scenario="P20_N10")
)

# P10_N20
RSM.sum.Q.vol2.FWOLL$TP <- RSM.sum.Q.vol2.FWOLL$TP.fit.pred*0.9
RSM.sum.Q.vol2.FWOLL$DIN <- RSM.sum.Q.vol2.FWOLL$DIN.fit.pred*0.8

Chla.rsm.pred2 <- predict(Chla_HGAM_mod_sum3,newdata=RSM.sum.Q.vol2.FWOLL,type="response",se.fit=T)
RSM.sum.Q.vol2.FWOLL$Chla.fit.pred <- Chla.rsm.pred2$fit|>as.numeric()

tmp <- rbind(tmp,
             ddply(RSM.sum.Q.vol2.FWOLL,c("EcoZone3.f"),summarise,
                   mean.val = mean(Chla.fit.pred,na.rm=T),
                   SE.val = SE(Chla.fit.pred))|>
               mutate(scenario="P10_N20")
)

# P30_N30
RSM.sum.Q.vol2.FWOLL$TP <- RSM.sum.Q.vol2.FWOLL$TP.fit.pred*0.7
RSM.sum.Q.vol2.FWOLL$DIN <- RSM.sum.Q.vol2.FWOLL$DIN.fit.pred*0.7

Chla.rsm.pred2 <- predict(Chla_HGAM_mod_sum3,newdata=RSM.sum.Q.vol2.FWOLL,type="response",se.fit=T)
RSM.sum.Q.vol2.FWOLL$Chla.fit.pred <- Chla.rsm.pred2$fit|>as.numeric()

tmp <- rbind(tmp,
             ddply(RSM.sum.Q.vol2.FWOLL,c("EcoZone3.f"),summarise,
                   mean.val = mean(Chla.fit.pred,na.rm=T),
                   SE.val = SE(Chla.fit.pred))|>
               mutate(scenario="P30_N30")
)

tmp$scenario <- factor(tmp$scenario,levels=c("base", "P10_N10", "P20_N10", "P10_N20", "P30_N30"))

ggplot(tmp,aes(x=scenario,y=mean.val,fill=scenario))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = format(round(mean.val,1))),vjust = -1.5)+
  geom_errorbar(aes(ymin = mean.val - SE.val,ymax = mean.val+SE.val),
                width = 0.2,color="black")+
  scale_fill_manual(values = wesanderson::wes_palette("Zissou1",5,"continuous")) +
  facet_wrap(~EcoZone3.f)+
  ylim(c(0,35))+
  labs(title = "Alternative: FWOLL (LOSOM + EAA Res)",
       subtitle = "POS Mean \u00B1 SE",
       y = expression("Chl-a ("*mu*"g L"^" -1"*")"),
       x = "Scenarios")

HABAM_WQScenario <- tmp
HABAM_WQScenario$scenario <- factor(HABAM_WQScenario$scenario,
                                    levels =  c("base", "P10_N10", "P20_N10", 
                                                "P10_N20","P30_N30"))

# END ---------------------------------------------------------------------
# save.image(file=paste0(data.path,"LOK_AlgaeEval_analysis.RData"))
