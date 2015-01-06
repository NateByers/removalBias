############################################################################
#                                                                          #
# Beginning of code for daily space/time interpolation                     # 
#                                                                          #
############################################################################
library(dplyr)
library(spacetime)
library(xts)
library(data.table)
# Source "download.R" to get the data downloaded and saved as an .rdata file
load("removal_bias_data.rdata")

# function for making an STSDF object (spacetime package)
makeSTSDF <- function(pollution, locations){
  # pollution <- allozone
  # locations <- sites
  
  # order pollution data by date then site id
  pollution <- arrange(pollution, Date, aqsid_dash)
  
  # order site locations by site id
  locations <- arrange(locations, aqsid_dash)
  
  # merge (inner join) site locations with pollution monitor locations
  ids <- unique(pollution$aqsid_dash)
  locations <- merge(locations[, -1], data.frame(aqsid_dash = ids))
  locations$index <- 1:dim(locations)[1]
  
  # create time data frame 
  date <- unique(pollution$date)
  date.df <- data.frame(date, index = 1:length(date))
  
  # create index for space and time
  pol.dt <- as.data.table(pollution[, c("aqsid_dash", "date")])
  sp.index <- merge(pol.dt, locations[, c("aqsid_dash", "index")], by = "aqsid_dash")
  time.index <- merge(pol.dt, date.df, by = "date")
  ind <- as.matrix(cbind(sp.index$index, time.index$index))
  
  # create Date object
  days <- as.Date(date.df$date, format = "%Y-%m-%d")
  
  # create sp object
  sp <- SpatialPoints(as.matrix(locations[, c("Longitude", "Latitude")]),
                      CRS("+proj=longlat +datum=WGS84"))
  
  # create STSDF object
  STSDF(sp = sp, time = days, data = pollution[, c("Max.Value", "aqsid_dash", "date")],
        index = ind)
  
  
}