# Source "download.R" to get the data downloaded and saved as an .rdata file
library(dplyr)
library(data.table)
library(deldir)
load("removal_bias_data.rdata")

# remove non-States (keep DC)
non.states <- c(72, 78, 80, 66)
all88101 <- subset(all88101, !(State.Code %in% non.states))
all88502 <- subset(all88502, !(State.Code %in% non.states))
allozone <- subset(allozone, !(State.Code %in% non.states))

# function for processing data so that POCs are averaged
processPOCs <- function(data){
  # PMdata <- all88101
  data <- group_by(data, State.Code, County.Code, Site.ID, Date)
  as.data.frame(summarize(data, Max.Value = mean(Max.Value)))
}

all88101 <- processPOCs(all88101)
all88502 <- processPOCs(all88502)
allozone <- processPOCs(allozone)

# add key column to sites data frame
sites$sites.key <- 1:dim(sites)[1]

# add sites key to PM and ozone data frames
State.Code <- as.integer(substr(sites$aqsid_dash, 1, 2))
County.Code <- as.integer(substr(sites$aqsid_dash, 4, 6))
Site.ID <- as.integer(substr(sites$aqsid_dash, 8, 11))
sites.key.df <- data.frame(State.Code, County.Code, Site.ID, 
                           sites.key = sites$sites.key, aqsid_dash = sites$aqsid_dash,
                           stringsAsFactors = FALSE)

# function for merging sites.key.df with daily values
mergeSiteKey <- function(pol.data, key.data, paste.POC = FALSE){
  # pol.data
  pol.data <- merge(pol.data, key.data, all.x = TRUE)
  if(paste.POC){
    pol.data$aqsid_dash <- paste(pol.data$aqsid_dash, pol.data$POC, sep = "-")
  }
  pol.data
}

all88101 <- mergeSiteKey(all88101, sites.key.df)
all88502 <- mergeSiteKey(all88502, sites.key.df)
allozone <- mergeSiteKey(allozone, sites.key.df)

# split data into .rdata files for quick access

# create directories for each pollutant
lapply(c("all88101", "all88502", "allozone"), dir.create)

# function for splitting data into .rdata files for each site and saving in 
# appropriate directory
splitData <- function(data, folder){
  data <- as.data.table(data)
  setkey(data, aqsid_dash)
  keys <- unique(data$aqsid_dash)
  for(i in keys){
    sub.data <- subset(data, aqsid_dash == i)
    file <- paste0(folder, "/", sub.data$aqsid_dash[1], ".rdata")
    sub.data <- sub.data[, c("Date", "Max.Value", "sites.key"), with = FALSE]
    save(sub.data, file = file)
  }
}

splitData(allozone, "allozone")
splitData(all88101, "all88101")
splitData(all88502, "all88502")

# create sites data frame for ozone and PM
ozone.sites <- merge(sites, data.frame(sites.key = unique(allozone$sites.key)))
pm.88101.sites <- merge(sites, data.frame(sites.key = unique(all88101$sites.key)))
pm.88502.sites <- merge(sites, data.frame(sites.key = unique(all88502$sites.key)))

# function for distance in kilometers between two long/lat positions (from "fossil" package)
earth.dist <- function (long1, lat1, long2, lat2) 
{
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.145
  d <- R * c
  return(d)
}


# create function that gets Dirichlet (aka Voronoi) tessellation for the ozone and
# PM monitors and calculates distances between nearest neighbors 
getDistances <- function(pol.sites){
  # pol.sites <- ozone.sites
  sites.deldir <- deldir(pol.sites$Longitude, pol.sites$Latitude)
  combos <- sites.deldir$delsgs
  
  aqsid_1 <- sapply(1:dim(combos)[1], function(x){
    ind1 <- combos[x, "ind1"]
    aqsid_1 <- pol.sites[ind1, "aqsid_dash"]
    aqsid_1
  })
  
  combos$ind1 <- aqsid_1
  
  aqsid_2 <- sapply(1:dim(combos)[1], function(x){
    ind2 <- combos[x, "ind2"]
    aqsid_2 <- pol.sites[ind2, "aqsid_dash"]
    aqsid_2
  })
  
  combos$ind2 <- aqsid_2
  
  dist <- mapply(FUN = earth.dist, long1 = combos[, 1],
                 lat1 = combos[, 2], long2 = combos[, 3],
                 lat2 = combos[, 4])
  
  combos$dist <- dist
  combos
}

ozone.combos <- getDistances(ozone.sites)
pm.88101.combos <- getDistances(pm.88101.sites)
pm.88502.combos <- getDistances(pm.88502.sites)



calcBias <- function(monitor.file, pollutant.directory, distances, sites){
  # pollutant.directory <- "allozone"
  # monitor.file <- "01-003-0010.rdata"
  # distances <- ozone.combos
  # sites <- ozone.sites
  
  # load the data for the monitor of interest. All of the .rdata files have a data.table object 
  # called sub.data
  load(paste0(pollutant.directory, "/", monitor.file))
  
  # get monitor id
  monitor <- strsplit(monitor.file, ".", fixed = TRUE)[[1]][1]
  
  # remove rows with NAs
  sub.data <- sub.data[!is.na(sub.data$Max.Value), ]
  
  # subset distances data.frame down to rows with the sites key for the monitor of interest
  distances <- subset(distances, ind1 == monitor  | ind2 == monitor)
  
  if(dim(distances)[1] > 0){
    # get the ids for the nearest neighbors
    neighbor.ids <- rbind(data.frame(id = distances$ind1, dist = distances$dist),
                          data.frame(id = distances$ind2, dist = distances$dist))
    neighbor.ids <- neighbor.ids[neighbor.ids$id != monitor, ]
    
    # create a data table for collecting the daily values
    monitors.dt <- sub.data[, c("Date", "Max.Value"), with = FALSE]
    setnames(monitors.dt, "Max.Value", paste0("site_", monitor))
    rm(sub.data)
    
    # load the monitor data for the neighbors and add to monitors.dt
    for(i in neighbor.ids$id){
      # i <- neighbor.ids$id[1]
      # i <- neighbor.ids$id[2]
      # i <- neighbor.ids$id[3]
      # i <- neighbor.ids$id[4]
      
      # load data for site 
      load(paste0(pollutant.directory, "/", i, ".rdata"))
      # subset the data table to just Date and Max.Value columns
      sub.data <- sub.data[, c("Date", "Max.Value"), with = FALSE]
      # rename the value column
      setnames(sub.data, "Max.Value", paste0("site_", i))
      
      # merge with monitors.dt
      monitors.dt <- merge(monitors.dt, sub.data, by = "Date", all.x = TRUE)
      
      # remove sub.data
      rm(sub.data)
      
    }
    
    # make matrix of monitor values
    values <- as.matrix(monitors.dt[, 3:dim(monitors.dt)[2], with = FALSE]) 
    
    # make matrix of weights with same dimension as value matrix 
    # (weights are inverse distances squared)
    weights <- matrix(1/(rep(distances$dist, dim(monitors.dt)[1])^2),
                      nrow = dim(monitors.dt)[1], byrow = TRUE)
    
    # replace weights with 0 where there is an NA in the values matrix
    weights[is.na(values)] <- 0
    
    # replace values with 0 where there is an NA 
    values[is.na(values)] <- 0
    
    # multiply the values and weights matrices and calculate inner product using 
    # a vector of ones to get the sums for each row 
    summed <- (values * weights) %*% rep(1, dim(values)[2])
    
    # calculate the sum of each row in the 
    denom <- weights %*% rep(1, dim(values)[2])
    
    # if all neighboring monitors are missing a value for a day, remove that day from
    # the analysis
    if(sum(denom == 0) > 0){
      monitors.dt <- as.data.table(as.data.frame(monitors.dt)[!denom == 0, ])
      summed <- summed[!denom == 0, ]
      denom <- denom[!denom == 0, ]
    }
    
    # if the denom vector has zeros, remove that index from denom and summed
    denom <- denom[denom != 0]
    summed <- summed[denom != 0]
    
    # calculate inverse distance squared weighted average for each day
    weighted.avg <- summed/denom
    
    # get the daily values for the monitor of interest as a vector
    daily <- as.data.frame(subset(monitors.dt, select = 2))[, 1]
    
    # calculate difference between each interpolated value and the actual
    # value for the monitor
    diff <-  weighted.avg - daily 
    
    # if there are zeros in the daily vector, get rid of that index
    rel.daily <- daily[daily != 0]
    rel.diff <- diff[daily != 0]
    
    # calculate relative bias
    rel.diff <- rel.diff/rel.daily
    
    # create data table with mean, min, max, standard deviation, and number of days
    data.table(siteID = monitor, bias_mean = mean(diff), bias_min = min(diff), 
               bias_max = max(diff), bias_sd = sd(diff), bias_n = length(diff),
               relbias_mean = mean(rel.diff), relbias_min = min(rel.diff), 
               relbias_max = max(rel.diff))
  } else {
    data.table(siteID = monitor, bias_mean = NA, bias_min = NA, bias_max = NA,
               bias_sd = NA, bis_n = NA, relbias_mean = NA, relbias_min = NA, 
               relbias_max = NA)
  }
  
  
}

# get list of files for each pollutant
ozone.files <- list.files("allozone")
pm.88101.files <- list.files("all88101")
pm.88502.files <- list.files("all88502")

# get biases for ozone
library(pbapply)
options("pbapply.pb"="txt")
then <- Sys.time()
ozone.bias.ls <- pblapply(ozone.files, calcBias, pollutant.directory = "allozone",
                          distances = ozone.combos, sites = ozone.sites)
Sys.time() - then 
ozone.bias.dt <- rbindlist(ozone.bias.ls)
write.csv(ozone.bias.dt, file = "ozoneBias.csv")

# get biases for 88101
then <- Sys.time()
pm.88101.bias.ls <- pblapply(pm.88101.files, calcBias, pollutant.directory = "all88101",
                          distances = pm.88101.combos, sites = pm.88101.sites)
Sys.time() - then 
pm.88101.bias.dt <- rbindlist(pm.88101.bias.ls)
write.csv(pm.88101.bias.dt, file = "pm88101Bias.csv")

# get biases for 88502
then <- Sys.time()
pm.88502.bias.ls <- pblapply(pm.88502.files, calcBias, pollutant.directory = "all88502",
                             distances = pm.88502.combos, sites = pm.88502.sites)
Sys.time() - then 
pm.88502.bias.dt <- rbindlist(pm.88502.bias.ls)
write.csv(pm.88502.bias.dt, file = "pm88502Bias.csv")

# make a sites table with mean bias for each pollutant
ozone.bias.dt <- ozone.bias.dt[, c("siteID", "bias_mean", "relbias_mean"), with = FALSE] 
setnames(ozone.bias.dt, c("bias_mean", "relbias_mean"), c("bias_ozone", "relbias_ozone"))
pm.88101.bias.dt <- pm.88101.bias.dt[, c("siteID", "bias_mean", "relbias_mean"), with = FALSE] 
setnames(pm.88101.bias.dt, c("bias_mean", "relbias_mean"), c("bias_pm88101", "relbias_pm88101"))
bias.dt <- merge(ozone.bias.dt, pm.88101.bias.dt, by = "siteID", all = TRUE)
pm.88502.bias.dt <- pm.88502.bias.dt[, c("siteID", "bias_mean", "relbias_mean"), with = FALSE] 
setnames(pm.88502.bias.dt, c("bias_mean", "relbias_mean"), c("bias_pm88502", "relbias_pm88502"))
bias.dt <- merge(bias.dt, pm.88502.bias.dt, by = "siteID", all = TRUE)
write.csv(bias.dt, file = "allSitesBias.csv")
