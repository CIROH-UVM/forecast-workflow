import netCDF4 as nc
import matplotlib.pyplot as plt
import os
import numpy as np

# Need to add outfiles and then filename to load...
run1_dir = '/data/forecastRuns/7dayforecast-20231115/aem3d-run/'
run2_dir = '/data/forecastRuns/7dayforecast-20231115-newProdPat/aem3d-run/'

output_file = 'outfiles/nc/sheet_surface_forecast.nc'

run1 = nc.Dataset(os.path.join(run1_dir, output_file))
run2 = nc.Dataset(os.path.join(run2_dir, output_file))

fig, axs = plt.subplots(2)
axs[0].imshow(run1.variables['TCHLA'][-1,:])
axs[0].set_title('Production Run')
axs[1].imshow(run2.variables['TCHLA'][-1,:])
axs[1].set_title('Test Run')
plt.show(block=False)

# Now, do a profile over time
output_file = 'outfiles/nc/profile51.nc'

run1 = nc.Dataset(os.path.join(run1_dir, output_file))
run2 = nc.Dataset(os.path.join(run2_dir, output_file))

plt.figure()
plt.plot(run1.variables['Ordinal_Dates'][:], np.mean(run1.variables['TEMPERATURE'][:,0,:], axis=1))
plt.plot(run2.variables['Ordinal_Dates'][:], np.mean(run2.variables['TEMPERATURE'][:,0,:], axis=1))
plt.title('Temperature')
plt.show(block=False)

plt.figure()
plt.plot(run1.variables['Ordinal_Dates'][:], np.mean(run1.variables['TCHLA'][:,0,:], axis=1))
plt.plot(run2.variables['Ordinal_Dates'][:], np.mean(run2.variables['TCHLA'][:,0,:], axis=1))
plt.title('TCHLA')
plt.show(block=False)

input("Press key to exit")

# to grab different indexes from a 2D array
'''
https://stackoverflow.com/questions/23435782/numpy-selecting-specific-column-index-per-row-by-using-a-list-of-indexes
>>> data = ds.variables['TEMPERATURE'][[0,-1],0,:]
>>> data[np.arange(data.shape[0]), [32,33]]
'''

# profile variables
'''
   dimensions(sizes): T(7776), STRING_LENGTH(32), Static2DVariables(6), Dynamic2DVariables(7), Dynamic3DVariables(26), Z(68), Station(1)
    variables(dimensions): float64 T(T), float64 Ordinal_Dates(T), int32 Year(T), int32 Month(T), int32 Day(T), int32 Hour(T), int32 Minute(T), float64 Second(T), |S1 Static2DVariables(Static2DVariables, STRING_LENGTH), |S1 Dynamic2DVariables(Dynamic2DVariables, STRING_LENGTH), |S1 Dynamic3DVariables(Dynamic3DVariables, STRING_LENGTH), float32 Z(Z), float32 DZ(Z), float32 Station(Station), float32 Iway_pts(Station), float32 Jway_pts(Station), float32 BATH(Station), float32 LATITUDE(Station), float32 LONGITUDE(Station), float32 GRID_X(Station), float32 GRID_Y(Station), float32 UTM_ZONE(Station), float32 HEIGHT(T, Station), float32 BOTTOM(T, Station), float32 WIND_SPEED(T, Station), float32 WIND_DIR(T, Station), float32 BLUEICE_DEPTH(T, Station), float32 WHITEICE_DEPTH(T, Station), float32 SNOW_DEPTH(T, Station), float32 TEMPERATURE(T, Station, Z), float32 DENSITY(T, Station, Z), float32 U_VELOCITY(T, Station, Z), float32 V_VELOCITY(T, Station, Z), float32 RETENTION_T(T, Station, Z), float32 TRACER_1(T, Station, Z), float32 TRACER_2(T, Station, Z), float32 TRACER_3(T, Station, Z), float32 PAR_EXTINCTION(T, Station, Z), float32 FDIAT(T, Station, Z), float32 CYANO(T, Station, Z), float32 NO3(T, Station, Z), float32 NH4(T, Station, Z), float32 PO4(T, Station, Z), float32 TN(T, Station, Z), float32 TP(T, Station, Z), float32 SiO2(T, Station, Z), float32 DO(T, Station, Z), float32 TOC(T, Station, Z), float32 SSOL1(T, Station, Z), float32 TCHLA(T, Station, Z), float32 GPP(T, Station, Z), float32 CYANO_IN(T, Station, Z), float32 CYANO_IP(T, Station, Z), float32 FDIAT_IN(T, Station, Z), float32 FDIAT_IP(T, Station, Z)
'''


# sheet_surface variables
'''
    maxI: 194
    maxJ: 75
    sheet_value: 1.0
    dimensions(sizes): T(168), STRING_LENGTH(32), Static2DVariables(6), Dynamic2DVariables(7), Dynamic3DVariables(26), X(194), Y(75)
    variables(dimensions): float64 T(T), float64 Ordinal_Dates(T), int32 Year(T), int32 Month(T), int32 Day(T), int32 Hour(T), int32 Minute(T), float64 Second(T), |S1 Static2DVariables(Static2DVariables, STRING_LENGTH), |S1 Dynamic2DVariables(Dynamic2DVariables, STRING_LENGTH), |S1 Dynamic3DVariables(Dynamic3DVariables, STRING_LENGTH), float32 X(X), float32 Y(Y), float32 DX_i(X), float32 DY_j(Y), float32 BATH(Y, X), float32 LATITUDE(Y, X), float32 LONGITUDE(Y, X), float32 GRID_X(Y, X), float32 GRID_Y(Y, X), float32 UTM_ZONE(Y, X), float32 WIND_SPEED(T, Y, X), float32 WIND_DIR(T, Y, X), float32 HEIGHT(T, Y, X), float32 BOTTOM(T, Y, X), float32 BLUEICE_DEPTH(T, Y, X), float32 WHITEICE_DEPTH(T, Y, X), float32 SNOW_DEPTH(T, Y, X), float32 TEMPERATURE(T, Y, X), float32 U_VELOCITY(T, Y, X), float32 V_VELOCITY(T, Y, X), float32 RETENTION_T(T, Y, X), float32 DENSITY(T, Y, X), float32 TRACER_1(T, Y, X), float32 TRACER_2(T, Y, X), float32 TRACER_3(T, Y, X), float32 PAR_EXTINCTION(T, Y, X), float32 FDIAT(T, Y, X), float32 CYANO(T, Y, X), float32 NO3(T, Y, X), float32 NH4(T, Y, X), float32 PO4(T, Y, X), float32 TN(T, Y, X), float32 TP(T, Y, X), float32 SiO2(T, Y, X), float32 DO(T, Y, X), float32 TOC(T, Y, X), float32 SSOL1(T, Y, X), float32 TCHLA(T, Y, X), float32 GPP(T, Y, X), float32 CYANO_IN(T, Y, X), float32 CYANO_IP(T, Y, X), float32 FDIAT_IN(T, Y, X), float32 FDIAT_IP(T, Y, X)
'''