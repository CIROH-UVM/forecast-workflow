

import requests
import json
import subprocess
from pathlib import Path
from datetime import datetime
import rasterio
import numpy as  np
from osgeo import gdal, osr
import time
import shutil


class CIFilesDownloadProcess:
    def __init__ (self, app_key_path, start_date, end_date, output_dir, aoi_json_path, 
                 ci_python_script_path, convert_to_ci=False, remove_temp=True):
        """
        
        Initializes the CIFilesDownloadProcess class.

        Args:
            app_key_path (str): Path to the file containing the app key for authentication.
            start_date (str): Start date for downloading Cyanobacterial Index (CI) files in 'YYYY-MM-DD' format.
            end_date (str)                : End date for downloading CI files in 'YYYY-MM-DD' format.
            output_dir (str)              : Directory where output files will be saved.
            aoi_json_path (str)           : Path to the JSON file defining the area of interest (AOI) for cropping GeoTIFF files.
            ci_python_script_path (str)   : Path to the Python script used for downloading CI files.
            convert_to_ci (bool, optional): Flag indicating whether to convert Digital Number (DN) values to CI values. Defaults to False.
            remove_temp (bool, optional)  : Flag indicating whether to remove temporary files after processing. Defaults to True.

        Attributes:
            app_key (str): The app key for authentication.
            start_date (str): Start date for downloading CI files.
            end_date (str): End date for downloading CI files.
            output_dir (Path): Path object for the output directory.
            json_data (dict): Parsed JSON data defining the AOI.
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
        with open(app_key_path, 'r') as f:
            self.app_key = f.read().strip()
            
        self.start_date = start_date
        self.end_date = end_date
        
        self.output_dir = Path(output_dir)
        with open(aoi_json_path, 'r') as j_file:
            self.json_data = json.load(j_file)

        self.ci_python_script_path = ci_python_script_path
        
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


        crop_dir_path = self.output_dir / 'Lake_Champlain_DN_Files'
        crop_dir_path.mkdir(parents=True, exist_ok=True)
        self.crop_dir_path = crop_dir_path

        self.convert_to_ci = convert_to_ci
        if self.convert_to_ci:
            ci_output_path = self.output_dir / 'Lake_Champlain_CI_Files'
            ci_output_path.mkdir(parents=True, exist_ok=True)
            self.ci_output_path = ci_output_path

        self.remove_temp = remove_temp

    def download_urls(self):
        """
        To Fetch URLs For CI Files and save them in a text file.
        See: https://oceancolor.gsfc.nasa.gov/about/projects/cyan/
        
        """
        url = "https://oceandata.sci.gsfc.nasa.gov/api/cyan_file_search"
        # end_date = datetime.now().strftime('%Y-%m-%d')
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
        subprocess.run([
            
            "python", self.ci_python_script_path,
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
 
        coords = self.json_data['features'][0]['geometry']['coordinates'][0]
        min_long, min_lat = min(coords, key=lambda c: (c[0], c[1]))
        max_long, max_lat = max(coords, key=lambda c: (c[0], c[1]))
        translate_options = gdal.TranslateOptions(
            format='GTiff', projWin=[min_long, max_lat, max_long, min_lat]
        )
        gdal.Translate(str(output_file), str(input_file), options=translate_options)


    def convert_dn_to_ci(self):
        """
        Converts Digital Number values to Cyanobacterial Index.
        
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

       
     


