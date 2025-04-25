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



# END ---------------------------------------------------------------------
# save.image(file=paste0(data.path,"LOK_AlgaeEval_analysis.RData"))
