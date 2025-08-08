import datetime as dt
from .utils import parse_to_datetime, add_units
import requests
import pandas as pd
from urllib.error import HTTPError
import subprocess
from pathlib import Path
import rasterio
import numpy as  np
from osgeo import gdal, osr
import time
import shutil
import os

import glob
import xarray as xr
import rioxarray as riox

import contextily as cx
import matplotlib.pyplot as plt
import matplotlib.colors as mc

class CIFilesDownloadProcess:
    def __init__ (self, app_key, start_date, end_date, output_dir, cropbox, 
                 areaid, convert_to_ci=False, remove_temp=True):
        """
        
        Initializes the CIFilesDownloadProcess class.

        Args:
            app_key (str): the app key for authentication.
            start_date (str): Start date for downloading Cyanobacterial Index (CI) files in 'YYYY-MM-DD' format.
            end_date (str)                : End date for downloading CI files in 'YYYY-MM-DD' format.
            output_dir (str)              : Directory where output files will be saved.
            cropbox (lat,lon,lat,lon)       : corners of the area to crop the geotiff with
            areaid(str)                     : string designating the tileid/areaid for the tiles to download, i.e. "8_2" for Champlain Valley
            convert_to_ci (bool, optional): Flag indicating whether to convert Digital Number (DN) values to CI values. Defaults to False.
            remove_temp (bool, optional)  : Flag indicating whether to remove temporary files after processing. Defaults to True.

        Attributes:
            app_key (str): The app key for authentication.
            start_date (str): Start date for downloading CI files.
            end_date (str): End date for downloading CI files.
            output_dir (Path): Path object for the output directory.
            ci_python_script_path (str): Path to the Python script for downloading CI files.
            temp_dir (Path): Path object for the temporary directory.
            urls_file (Path): Path object for the file storing the list of URLs to download.
            geotiff_path (Path): Path object for the directory storing downloaded GeoTIFF files.
            projected_path (Path): Path object for the directory storing reprojected GeoTIFF files.
            crop_dir_path (Path): Path object for the directory storing cropped DN GeoTIFF files.
            ci_output_path (Path, optional): Path object for the directory storing converted CI GeoTIFF files if `convert_to_ci` is True.
            convert_to_ci (bool): Indicates whether to convert DN to CI values.
            remove_temp (bool): Indicates whether to remove temporary files after processing.
        
        """
        
        self.app_key = app_key
            
        self.start_date = start_date
        self.end_date = end_date
        
        self.output_dir = Path(output_dir)
            
        self.cropbox = cropbox
        
        self.areaid = areaid
     
        self.temp_dir = self.output_dir / 'temp'
        self.temp_dir.mkdir(parents=True, exist_ok=True) # parents=True, means create all the 
                                                            # parents directories too if they dont exist
        urls_path = self.output_dir / 'Urls_File'
        urls_path.mkdir(parents=True, exist_ok=True) 
        self.urls_file = urls_path / 'urls.txt'

        geotiff_path = self.temp_dir / 'GeoTiffs' # to store downloaded GeoTiffs
        geotiff_path.mkdir(parents=True, exist_ok=True)
        self.geotiff_path = geotiff_path

        projected_path = self.temp_dir / 'Projected_Geotiffs'  # To store projected GeoTiffs into EPSG:4326
        projected_path.mkdir(parents=True, exist_ok=True) # later to use crop Lake Champlain
        self.projected_path = projected_path


        crop_dir_path = self.output_dir / 'DN_Files'
        crop_dir_path.mkdir(parents=True, exist_ok=True)
        self.crop_dir_path = crop_dir_path

        self.convert_to_ci = convert_to_ci
        if self.convert_to_ci:
            ci_output_path = self.output_dir / 'CI_Files'
            ci_output_path.mkdir(parents=True, exist_ok=True)
            self.ci_output_path = ci_output_path

        self.remove_temp = remove_temp

    def download_urls(self):
        """
        To Fetch URLs For CI Files and save them in a text file.
        See: https://oceancolor.gsfc.nasa.gov/about/projects/cyan/
        
        """
        url = "https://oceandata.sci.gsfc.nasa.gov/api/cyan_file_search"
        # end_date = dt.now().strftime('%Y-%m-%d')
        data = {
            "sdate": self.start_date, # "%Y-%m-%d"
            "edate": self.end_date, # "%Y-%m-%d"
            "period": "2", # either 14 days (1), daily (2), weekly (3)
            "product": "1", # Cyanobacterial index (1), True Color Image (2)
            "results_as_file": "1",
            "addurl": "1",
            "region": "1",
            "areaids": '["13"]', # area id or tile number
            "resolution": "1", # 300m (1), 500m (2) - CONUS=300m, Alaska = 500m
            "cyan_search": "cyan_search"
        }
        response = requests.post(url, data=data)
        urls = [url.strip() for url in response.text.split('\n') if url.strip()]
        
        with open(self.urls_file, 'w') as f:
            f.write('\n'.join(urls))



    def download_tiles(self):
        """
        A Method to Download Complete Tiles and store in the temp directory. 
        
        """
        # find directory script is running from 
        downloadfct = os.path.dirname(os.path.realpath(__file__))+'/ci_data_download/ci-download.py'
        print("prepping to run ",downloadfct)
        subprocess.run([
            
            "python", downloadfct,
            "--filelist", str(self.urls_file),
            "--odir", str(self.geotiff_path),
            "--appkey", self.app_key])
        

    def convert_projection(self):
        """
        Convertes files into WGS 84 Projections
        
        """
    
        # we will loop through each file and project it and store it into the projected path
        for file in self.geotiff_path.glob("*.tif"):
            output_file_path = self.projected_path / file.name
            self.reproject_to_wgs84(file, output_file_path)
    
    @staticmethod
    def reproject_to_wgs84(input_file, output_file):
        """
        Converts projection to WGS 84.
        
        """
        dst_src = osr.SpatialReference()
        dst_src.ImportFromEPSG(4326)
        gdal.Warp(str(output_file), str(input_file), dstSRS=dst_src, format='GTiff')

    def crop_aoi(self):
        """
        Crops files to the area of interest.
        """
        
        for file in self.projected_path.glob('*.tif'):
            output_file = self.crop_dir_path / file.name
            self.crop_using_gdal(file, output_file)

    def crop_using_gdal(self, input_file, output_file):
        """
        
        Crops a GeoTiffs to the Specified AOI.
        
        """
 
        #coords = self.json_data['features'][0]['geometry']['coordinates'][0]
        #min_long, min_lat = min(coords, key=lambda c: (c[0], c[1]))
        #max_long, max_lat = max(coords, key=lambda c: (c[0], c[1]))
        translate_options = gdal.TranslateOptions(
            format='GTiff', projWin=self.cropbox
        )
        gdal.Translate(str(output_file), str(input_file), options=translate_options)


    def convert_dn_to_ci(self):
        """
        Converts Digital Number values to Cyanobacterial Index.
        Problematic results... not currently recommended 2025-07-30
        
        """
        
        for file in self.crop_dir_path.glob('*.tif'):
            output_file = self.ci_output_path / file.name
            self.dn_to_ci(file, output_file)

    @staticmethod
    def dn_to_ci(input_file, output_file):
        """
        Transforms DN values to CI values.

        For DN to CI Formula See - https://oceancolor.gsfc.nasa.gov/about/projects/cyan/
        
        """
        with rasterio.open(input_file) as src:
            dn_img = src.read(1)
            ci_img = np.power(10, dn_img * 0.011714 - 4.1870866)

            with rasterio.open(
                output_file, 'w', driver='GTiff', height=ci_img.shape[0],
                width=ci_img.shape[1], count=1, dtype=ci_img.dtype,
                crs='+proj=latlong', transform=src.transform
            ) as dst:
                dst.write(ci_img, 1)
    
    def remove_temp_dir(self):
        """
        Removes the temporary directory and its contents.
        
        """
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
            
    def run_pipeline(self):
        """
        A Method to Execute the entire pipeline.

        
        """
        urls_file = self.download_urls()
        self.download_tiles()

        self.convert_projection()
        self.crop_aoi()

        if self.convert_to_ci: # if flag is true, call convert_dn_to_ci method.
            self.convert_dn_to_ci()

        if self.remove_temp: # remove the redundant data
            self.remove_temp_dir()
            
    def time_index_from_filenames(self, filenames, string_slice=slice(0, 10)):
        '''
        Helper function to generate a Pandas datetimeindex object
            from dates contained in a file path string
        '''
    
        date_strings = [os.path.basename(i)[string_slice] for i in filenames]
    
        return pd.to_datetime(date_strings,format='%Y%j')


    def tif_to_ds(self):
        # Get file paths and obtain list of dates from file
        print("Load Downloaded TIF Files into Dataset")
        mosaic_files = sorted(list(self.crop_dir_path.glob('*.tif')))
        mosaic_list = [str(path) for path in mosaic_files]
        print(mosaic_list)

        # Import data and create xarray dask array labelled by timestamps from files
        time_var = xr.Variable('time', self.time_index_from_filenames(mosaic_list, string_slice=slice(1, 8)))
        #chunks = {'x': 1452, 'y': 1839, 'band': 1}
        concat_arrays = xr.concat([riox.open_rasterio(i) for i in mosaic_list], dim=time_var)
        
        # Convert to dataset and set band names
        concat_ds = concat_arrays.to_dataset(dim='band')
        concat_ds = concat_ds.rename({1: 'DN'})
        print("Set Coordinate Reference System to Epsg4326")
        concat_ds.rio.write_crs('epsg:4326', inplace=True)
        concat_ds = concat_ds.sortby('time')
        return concat_ds

def get_data(start_date,
			end_date,
			appkey,
			cropbox,
			output_dir,
            locations={'8_2'},
            remove_temp = True,
            convert_to_ci = False, 
            areaid = "8_2",
			variables={'streamflow':'Débit (m³/s)'},
			service='iv'):
    """
    A function to download and process Sentinel 3 observational cyano index image data

    Args:
    -- start_date (str, date, or datetime) [req]: the start date for which to grab Canadian Instantaneous data
    -- end_date (str, date, or datetime) [req]: the end date for which to grab Canadian Instantaneous data
	-- service (str) [opt]: what frequency of data to get. Default is 'iv', or instantaneous data (15-min frequency). Other option is 'dv', which returns daily data. For more info, see https://www.cehq.gouv.qc.ca/hydrometrie/historique_donnees/fiche_station.asp?NoStation=030425
    -- app_key (str): string containing the app key for authentication.
    -- output_dir (str)              : Directory where output files will be saved.
    -- cropbox (lat,lon,lat,lon)       : corners of the area to crop the geotiff with
    -- areaid(str)                     : string designating the tileid/areaid for the tiles to download, i.e. "8_2" for Champlain Valley
    -- convert_to_ci (bool, optional): Flag indicating whether to convert Digital Number (DN) values to CI values. Defaults to False.
    -- remove_temp (bool, optional)  : Flag indicating whether to remove temporary files after processing. Defaults to True.

        Attributes:
            app_key (str): The app key for authentication.
            start_date (str): Start date for downloading CI files.
            end_date (str): End date for downloading CI files.
            output_dir (Path): Path object for the output directory.
            temp_dir (Path): Path object for the temporary directory.
            urls_file (Path): Path object for the file storing the list of URLs to download.
            geotiff_path (Path): Path object for the directory storing downloaded GeoTIFF files.
            projected_path (Path): Path object for the directory storing reprojected GeoTIFF files.
            crop_dir_path (Path): Path object for the directory storing cropped DN GeoTIFF files.
            ci_output_path (Path, optional): Path object for the directory storing converted CI GeoTIFF files if `convert_to_ci` is True.
            convert_to_ci (bool): Indicates whether to convert DN to CI values.
            remove_temp (bool): Indicates whether to remove temporary files after processing.	
 
	Returns:
	Sentinel 3 satellite  observed data for the specified Tile ID in an xarray dataset where time slices represent each downloaded Tile for variable "DN"
	"""
    start_date = parse_to_datetime(start_date)
    end_date = parse_to_datetime(end_date)
    yearlist = list(range(start_date.year, end_date.year+1)) # the list of all years data is requested for
    print('sentinel3_ob.get_data years to process:')
    print(yearlist)
    print('sentinel3_ob.get_data() locations')
    print(locations)

	# instantiate the data object subsequent methods will use
    sentiles =  CIFilesDownloadProcess(
        app_key=appkey,
        start_date=start_date,
        end_date=end_date,
        output_dir=output_dir,
        cropbox= cropbox,
        areaid=areaid,
        convert_to_ci=convert_to_ci,
        remove_temp=remove_temp )
 
    urls_file = sentiles.download_urls()
    sentiles.download_tiles()

    sentiles.convert_projection()
    sentiles.crop_aoi()

    if sentiles.convert_to_ci: # if flag is true, call convert_dn_to_ci method.
        print("Generating CI Files")
        sentiles.convert_dn_to_ci()

    sentinel3_data = sentiles.tif_to_ds()

    if remove_temp: # remove the redundant data
        sentiles.remove_temp_dir()
        
        
    return sentinel3_data

def plot (da, cmap=None, base=True, title=None, **kwargs) :
        """ Plot a Sentinel3 DN for a time slice (2 dimensional array)
                If no colormap (cmap) passed, build a custom color map for DN encoding
                If base=True, display a basemap using CRS from passed dataset
                If title=None, let plot function build default title from data context
        """
            
        if (cmap == None) :
            # build custom colormap for DN encoding
            new_colors = plt.colormaps['coolwarm'](np.linspace(0, 1, plt.colormaps['coolwarm'].N)) 
            DN_cmap = mc.ListedColormap(new_colors) # new listed color colormap
            DN_cmap.colors[254] = (1,1,1,1) # ground value grey
            DN_cmap.colors[255] = (0,0,0,1) # no-data/cloud value black
        else : 
            DN_cmap = cmap
        
        fig, axs = plt.subplots()
        
        # Plot without including the data identified as ground (254)
        # It does include pixels identified as no-data/clouds (255)
        da.where(da.values!=254).plot(ax=axs,cmap=DN_cmap)
            
        if (base) :
            cx.add_basemap(axs, crs=da.rio.crs)
        if (title != None) :   
            plt.title(title)
        return(axs)