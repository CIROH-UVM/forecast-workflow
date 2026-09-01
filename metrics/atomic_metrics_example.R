source('atomic_metrics.R')

# set the random seed if you want reproducible mock data
# set.seed(42)

##### mock_raster_stack() #####
# create some mock data representing a raster stack / 3D array
mock_raster_stack <- function(nlayers, xy_dims, mean_range, stdev_range) {
  
  # create random means and stdevs so that each layer of the stack has a slightly different data
  rand_means <- runif(nlayers, mean_range[1], mean_range[2])
  rand_stdevs <- runif(nlayers, stdev_range[1], stdev_range[2])
  n_values <- prod(xy_dims) * nlayers
  
  # using a log-normal distribution because that seems to be a reasonable distribution for chla data
  # convert mean/sd to log-space parameters for all layers at once (vector math)
  sigma_sq <- log(1 + (rand_stdevs^2 / rand_means^2))
  meanlog  <- log(rand_means) - (sigma_sq / 2)
  sdlog    <- sqrt(sigma_sq)
  
  # repeat parameters for all cels in each layer
  # we're essentially repeating each value in meanlog and stlog x*y times
  # since meanlog and sdlong each have length = nlayers, meanlog_rep and sdlog_rep will each have lengths = nlayers*x*y
  meanlog_rep <- rep(meanlog, each = prod(xy_dims))
  sdlog_rep   <- rep(sdlog,   each = prod(xy_dims))
  
  # generate all log-normal values in one go
  values <- rlnorm(n_values, meanlog = meanlog_rep, sdlog = sdlog_rep)
  
  # reshape into 3D array and return
  return(array(values, dim = c(xy_dims, nlayers)))
}

nlayers <- 28
xy_dims <- c(248, 317)
mean_range <- c(5,15)
stdev_range <- c(1,5)

mockarr <- mock_raster_stack(nlayers, xy_dims, mean_range, stdev_range)

##### Mocking a 3D SpatRaster with a time index #####

# casting a 3D array to a SpatRaster is very simple:
mockrast <- rast(mockarr)

# now we want to create an index of datetimes (POSIXct)
# The POSIXct class stores date/time values as the number of seconds since January 1, 1970
first_date <- as.POSIXct("2019-05-28", tz = 'UTC')

# now we create the datetime vector
time_sequence <- seq(from = first_date, by = "day", length.out = nlyr(mockrast))

# and assign it to the SpatRaster
time(mockrast) <- time_sequence

##### Create Lake Subdomain #####

working_dir <- "/netfiles/ciroh/30dayHABsHindcast/results/Preprocessing"

# loading the shapefile that defines the St.Albans Bay subdomain
StAlbans <- vect(file.path(working_dir, "Detailed_Lake_Segments/StAlbans_prj.shp"))

# easiest way to test the create_subdomain function is to use some real lake data
# defining the baseline dir
baseline_file <- file.path(working_dir, "Monthly_FEE_Results7/CYANO_Baseline/CYANO_DailyAvg_MaySep_2019.nc")
# creating a raster
baseline_rast <- rast(baseline_file)

# create the st. albans subdomain
sta_subdomain <- create_subdomain(baseline_rast, StAlbans)

# Visually check the results to ensure the masking+cropping worked
plot(baseline_rast$CYANO_1)
plot(sta_subdomain$CYANO_1)

##### Extent #####

extents <- extent_wrap(mockarr, q_threshold = 15)

# array of extent fractions
extents[1,]

# array of incident pixel counts
extents[2,]

##### Incidence ######

# calling the atomic function on a single layer
incidence(subdomain_layer = mockarr[,,1], q_threshold = 15, areal_threshold = 0.20)

# calculating incidence from our extent results...
# slightly more efficient than calling incidence() or incidence_wrap(), since those functions have to calculate extent internally
incidences_ext <- incidence_from_extent(extent_fractions = extents[1,], areal_threshold = 0.20)

# calling the wrapper on the entire 3D array
incidences_wrp <- incidence_wrap(mockarr, q_threshold = 15, areal_threshold = 0.20)

# Can also work with a SpatRaster once you cast it as an array
incidences_wrp_rast <- incidence_wrap(as.array(mockrast), q_threshold = 15, areal_threshold = 0.20)

# gives us a vector with an incidence for each layer of the 3D array
incidences_wrp_rast

# check the ensure the results are the same from each incidence function
all(incidences_ext == incidences_wrp)
all(incidences_wrp == incidences_wrp_rast)

##### Duration  ######

# calculating duration from our incidence results which we already calculated
inc_dur <- duration_from_incidence(incidence_vect = incidences_ext, time_index = time(mockrast), window_size = 7, step_size = 7)


# calling just duration (computes incidence vector internally), which makes it slightly less efficient that duration_from_incidence()
# duration() function works with raster object, assuming the raster has a time index
dur <- duration(subdomain_rast = mockrast, q_threshold = 15, areal_threshold = 0.20, window_size = 7, step_size = 7)

# vector of the dates at the center of each rolling window
window_centers <- dur[[1]]
window_centers

# corresponding vector of the bloom duration time for each window
duration_time <- dur[[2]]
duration_time

# duration_time is a 'Duration' object from the lubridate library.
# to get a numeric vector, just cast as such and specify the desired time unit
as.numeric(duration_time, 'days')

