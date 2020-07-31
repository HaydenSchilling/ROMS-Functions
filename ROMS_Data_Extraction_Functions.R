# R Functions to extract data from a 3d ROMS model
# Hayden Schilling (SIMS) 
# 30/7/2020
# ROMS has a non-linear grid so this function works by first aligning the lat/lon grid with the index values
# Then it uses index values to find the correct grid cells

# To use you need to provide a dataframe containing Longitude, Latitude and Date as well as the variable of interest.

#########################################################################################################
### Function to get a depth profile of data for each location (data for multiple depths at each location)
#########################################################################################################


Extract_Profile_Data_ROMS <- function(data, variable){
  library(angstroms)
  library(anglr)
  library(tidyverse)
  
  # load test data and recognise date column
  my_locations <- data %>% select(Longitude, Latitude, Date) # grabs relevant columns from input data
  my_locations$Date <- lubridate::dmy(my_locations$Date) # parses date
  my_locations_column_number <-
    ncol(my_locations) #this will help with assigning output later
  
  # Progress Bar
  pb <- txtProgressBar(min = 0, max = length(data$Latitude), style = 3)
  # loop through each row of input data
  for (d in 1:nrow(my_locations)) {
    # extract date components from the specific row to feed into the thredds url
    yr <- lubridate::year(my_locations$Date[d])
    mth <-
      str_pad(lubridate::month(my_locations$Date[d]), 2, "left", pad = "0") # makes sure month has two digits
    dy <-
      str_pad(lubridate::day(my_locations$Date[d]), 2, "left", pad = "0") # makes sure day has two digits
    
    # creates the url for each specific day
    u <-
      paste0(
        "http://dapds00.nci.org.au/thredds/dodsC/rr6/eReefs/4kmReanalyses/ocean_avg_",
        yr,
        mth,
        dy,
        ".nc"
      )
    #vars$name # can check variables in each file
    
    # get the grid of lat/lon values to use when matching the grid indexes to locations
    coords <- romscoords(u, c("lon_u", "lat_u"))
    
    # check grid shape/location
    #plot(values(coords), pch =".")
    #maps::map(add=TRUE)
    
    xy <-
      my_locations[d, 1:2] # Long/Lat, selects only a single row of the data (could be improved to filter by day in future)
    n <- nrow(xy)
    
    for (i in 1:30) {
      # iterate through all 30 depth levels in ROMS model
      # first extract a temperature layer
      temp <- roms_xy(u, varname = variable, slice = c(i, 1)) # the 1 represents the single time point in the daily .nc file
      
      # now extract your specific point (or set of points) using the coords grabbed earlier
      temp0 <-
        raster::extract(temp,
                        romsmap(
                          sp::SpatialPointsDataFrame(sp::SpatialPoints(xy), data.frame(n = 1:n)),
                          coords
                        ),
                        method = "bilinear")
      
      # assign the extracted data to new columns, 1 column per depth level
      my_locations[d, my_locations_column_number + i] <-  temp0 
      
    }
    setTxtProgressBar(pb, d)
  }
  # Check output
  #my_locations
  
  # Take temperature output and make it long format
  my_locations_long <- my_locations %>%
    pivot_longer(
      cols = starts_with("V"),
      names_to = "Layer",
      names_prefix = "V",
      values_to = variable
    ) %>%
    mutate(Layer = 31 - (as.numeric(Layer) - my_locations_column_number))
  
  # Check output
  #my_locations_long
  
  ### now need to get the depth data
  coordh <- romshcoords(u) # extracts mesh grid with depth data
  
  # Select lat/long from input data
  xy <- my_locations[, 1:2]
  n <- nrow(xy)
  
  # extract depth bin info from all layers ar once for each 
  depth0 <-
    raster::extract(coordh, romsmap(SpatialPointsDataFrame(SpatialPoints(xy), data.frame(n = 1:n)),
                                    coords),
                    method = "bilinear")
  
  # combine depth and lat/lon info
  location_depths <- cbind(xy, depth0)
  
  # Make depth data long format
  location_depth_long <- location_depths %>%
    pivot_longer(
      cols = starts_with("layer"),
      names_to = "Layer",
      names_prefix = "layer.",
      values_to = "Depth_m"
    ) %>%
    mutate(Layer = as.numeric(Layer))
  
  # Check output
  #location_depth_long
  
  # Combine depth and temperature data
  full_dat <- full_join(my_locations_long, location_depth_long)
  return(full_dat)
  
}

#########################################################################################################
### Function to get a single point value for each location (variable which has no depth dimension)
#########################################################################################################

Extract_2D_Point_Data_ROMS <- function(data, variable) {
  # From Code to extract temperatures at all depth from a ROMS model for specific dates and locations
  # Hayden Schilling
  data <-  read.csv("Test points.csv")
  library(angstroms)
  library(anglr)
  library(tidyverse)
  
  # load test data and recognise date column
  my_locations <- data %>% select(Longitude, Latitude, Date) # Grabs columns from input data
  my_locations$Date <- lubridate::dmy(my_locations$Date) # recognise date column
  my_locations_column_number <-
    ncol(my_locations) #this will help with assigning output later
  
  # Progress Bar
  pb <-
    txtProgressBar(min = 0,
                   max = length(data$Latitude),
                   style = 3)
  # loop through each row of input data
  for (d in 1:nrow(my_locations)) {
    # extract date components from the specific row to feed into the thredds url
    yr <- lubridate::year(my_locations$Date[d])
    mth <-
      str_pad(lubridate::month(my_locations$Date[d]), 2, "left", pad = "0") # makes sure month has two digits
    dy <-
      str_pad(lubridate::day(my_locations$Date[d]), 2, "left", pad = "0") # makes sure day has two digits
    
    # creates the url for each specific day
    u <-
      paste0(
        "http://dapds00.nci.org.au/thredds/dodsC/rr6/eReefs/4kmReanalyses/ocean_avg_",
        yr,
        mth,
        dy,
        ".nc"
      )
    #vars$name # can check variables in each file
    
    # get the grid of lat/lon values to use when matching the grid indexes to locations
    coords <- romscoords(u, c("lon_u", "lat_u"))
    
    # check grid shape/location
    #plot(values(coords), pch =".")
    #maps::map(add=TRUE)
    
    xy <-
      my_locations[d, 1:2] # selects only a single row of the data (could be improved to filter by day in future)
    n <- nrow(xy)
    
    # iterate through all 30 depth levels in ROMS model
    # first extract a temperature layer
    temp <- raster(u, varname = variable)
    
    # now extract your specific point (or set of points) using the grid index we extracted earlier
    temp0 <-
      raster::extract(temp,
                      romsmap(
                        sp::SpatialPointsDataFrame(sp::SpatialPoints(xy), data.frame(n = 1:n)),
                        coords
                      ),
                      method = "bilinear")
    
    # assign the extracted data to new columns, 1 column per depth level
    my_locations[d, my_locations_column_number + 1] <-  temp0
    setTxtProgressBar(pb, d)
  }
  
  # Check output
  #my_locations
  full_data <- my_locations %>% rename("MLD" = "V4")
  
  #my_locations
  
  return(full_data)
  
}