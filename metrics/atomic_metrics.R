library(glue)
library(lubridate)
library(terra)


##### ATOMIC EXTENT #####

extent <- function(subdomain_layer, q_threshold) {
  # creating a reclassification matrix to determine HAB incidence
  # 0:q_threshold -> 0 (no bloom)
  # q_threshold:200 -> 1 (bloom)
  reclass_m <- matrix(c(0, q_threshold, 0, q_threshold, 200, 1), ncol = 3, byrow = TRUE)
  
  # casting the subdomain_layer as a SpatRaster (terra package) so that we can use the classify() method
  # we should be able to convert the SpatRaster back into a matrix (2D array) once we're done reclassifying
  subdomain_layer_rast <- rast(subdomain_layer)
  
  # reclasify the SpatRaster
  # right=TRUE means that the classification rules are right-closed and left-open
  # right-closed and left open is (0,1] = {x | 0 < x <= 1}
  # i.e. right-closed means that the upper limit is inclusive
  # ***so chla values equal to q_threshold will NOT be classified as incident. run ?classify for more info
  subdomain_layer_reclassified <- classify(subdomain_layer_rast, reclass_m, right=TRUE)
  
  # this counts the total pixels/cells not including NA cells
  # Takis uses this line in his total_pixels for loop, so this will be our preferred method
  total_pixels <- sum(!is.na(values(subdomain_layer_reclassified)))
  
  # similarly, we now count the total number of incident pixels/cells
  incident_pixels <- sum(values(subdomain_layer_reclassified) == 1, na.rm = TRUE)
  
  # now we get the fraction of incident cells compared to the total number of cells
  frac_incident_cells <- incident_pixels/total_pixels
  
  return(c(frac_incident_cells, incident_pixels))
}

##### EXTENT WRAPPER #####

extent_wrap <- function(subdomain_arr, q_threshold) {
  return(apply(subdomain_arr, MARGIN = 3, FUN = extent, q_threshold))
}

##### ATOMIC INCIDENCE #####

incidence <- function(subdomain_layer, q_threshold, areal_threshold){
  # the extent function returns the fraction of cells that pass the q_threshold (incident cells) in the given subdomain_layer
  extent_fractions <- extent(subdomain_layer, q_threshold)[1]
  
  # if the fraction of incident cells is greater than or equal to the spatial threshold, then return 1 (incidence)
  # if not, return 0 (no incident)
  if (extent_fractions >= areal_threshold) {
    return(1)
  } else {return(0)}
}

##### INCIDENCE WRAPPER ####

incidence_wrap <- function(subdomain_arr, q_threshold, areal_threshold){
  return(apply(subdomain_arr, MARGIN = 3, FUN = incidence, q_threshold, areal_threshold))
}

##### INCIDENCE_FROM_EXTENT #####

incidence_from_extent <- function(extent_fractions, areal_threshold){
  # extent_fractions is the first vector returned by extent_wrap()
  # classifying the vector of incident pixel/non-incident pixel fractions is super quick and easy
  return(as.integer(extent_fractions >= areal_threshold))
}

##### DURATION_FROM_INCIDENCE #####

duration_from_incidence <- function(incidence_vect, time_index, window_size, step_size) {
  # ensure the incidence and time index vectors have the same length
  stopifnot(length(incidence_vect) == length(time_index))
  
  # these indices correspond to the left-edge of each rolling window
  window_start_indices <- seq(1, length(incidence_vect) - window_size+1, by = step_size)
  
  # creating an empty datetime vector to store the datetimes at the center of each window
  window_end_timestamps <- as.POSIXct(rep(NA, length(window_start_indices)), tz = "UTC")
  # empty vector to store number of incidences counted per window
  incidence_counts  <- numeric(length(window_start_indices))
  
  
  for (i in seq_along(window_start_indices)) {
    # calculating the indices of the left and right edges of the rolling window
    left_idx <- window_start_indices[i]
    # minus 1 because of R's 1-based indexing (1 + 15 = 16, so we have to add -1: 1+15-1=15)
    right_idx <- left_idx + window_size - 1
    
    # storing the right-end timestamp of each rolling window (inclusive in each window)
    window_end_timestamps[i] <- time_index[right_idx]
    # storing the incidence count for each window
    # note that indexing in R is inclusive at both bounds
    incidence_counts[i] <- sum(incidence_vect[left_idx:right_idx])
  }
  # for now we are assuming evenly-spaced timestamps, so we can just use the first 
  # always return duration in days
  timeseries_interval <- as.duration(time_index[2]-time_index[1])
  # apply that time interval to the incident counts to get duration in terms of time rather than index-based incident counts
  duration_times <- timeseries_interval*incidence_counts
  
  # window_end_timestamps: the last time stamp included in each respective window 
  # duration_times: 'Duration' object from the 'lubridate' class. Default representation of duration times is in seconds
  # To easily convert duration times to a more relevant unit, use: as.numeric(duration_times, 'days') or whatever time unit you want
  return(list(window_end_timestamps, duration_times))
}

##### DURATION #####

duration <- function(subdomain_rast, q_threshold, areal_threshold, window_size, step_size) {
  # first we have to calculate incidence
  inc_vect <- incidence_wrap(as.array(subdomain_rast), q_threshold, areal_threshold)
  
  # now we can calculate duration from incidence
  dur_list <- duration_from_incidence(inc_vect, time(subdomain_rast), window_size, step_size)
  
  return(dur_list)
}

#### CREATE SUBDOMAIN #####

# masks and clips a raster according to a shapefile
# the raster can be a single layer or a stack. If a stack, it will mask+crop every layer in the stack
create_subdomain <- function(r, shapefile) {
  # returns a new raster masked and cropped to the shapefile
  print(glue("Masking and cropping raster for lake segment: {shapefile$SEG_NAME}"))
  r_masked <- mask(r, shapefile)
  r_cropped <- crop(r_masked, shapefile)
  return(r_cropped)
}
