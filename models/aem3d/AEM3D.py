from default_config import defaults
import argparse
import configparser

## Check on naming convention for classes and functions and variables
class Model:
  
  required_inputs = {}

  file_list = {}

  def Model(path, starttime, endtime, config={}):
    # All this from current IAM aem3d worker and prep worker -- why divided?
    
    # Copy base files
    # Mod setup
    # generate_input_files()
    pass

  def generate_input_files():
    # Use starttime / endtime from model creation
    
    #generate_temp()
    # and so on... See current AEM3D
    pass

  def input_from_calibration_file(filename, newRange, origRange = None):
    #If origRange = None, start at first mm/dd same as newRange
    #Else Check origRange == new Range

    #Load calibration file
    #Read calibration file until arrive at starting mm/dd
    #Then, build list within list of values using split(), and converting split[0] to new datetime
    
    #Create Pandas dataframe from lists within list
    #return dataframe
    pass
  
  def input_from_NWS_v21_hindcast():
    #read in hindcast to dataFrame
    pass

class InputFile:
  def write(data, headers = {}):
    pass
  
  def read():
    # Return {{Headers}, pandas DF}
    pass


## Emulate Python datetime... just add conversion capabilities to/from AEM3D datetime

class Datetime:
  
  datetime = None
  aem3d_datetime = None  ## Or, just format on the fly each time

  def from_AEM3D_datetime(aem3d_datetime):
    pass
  
  def to_AEM3D_datetime():
    # return aem3d_datetime
    pass

  def to_datetime():
    # return datetime
    pass

def parse_config():
  pass

def parse_args():
  parser = argparse.ArgumentParser(description="command-line arguments for running the AEM3D-based HABs Forecast",
                                   epilog="more details and documentation to come soon")
  pass

def get_settings():
  # check if a config file is passed
  # check if there are command line args passed
  # Use default settings for any that were not passed in either of the methods above
  # return settings dictionary
  pass
  
def main():
  # settings = parse_args(sys.argv[1:])
  pass

  #Create model

def test():
  SETTINGS_KEYS = ['forecast_start',
                    'forecast_end',
                    'spinup_start',
                    'blending_variable',
                    'blending_ratio',
                    'weather_dataset_observed',
                    'weather_dataset_forecast',
                    'hydrology_dataset_observed',
                    'hydrology_dataset_forecast']
  config_fpath = '/data/users/n/b/nbeckage/forecast-workflow/test_settings.cfg'
  config_file = configparser.ConfigParser()
  # read the configuration file
  config_file.read(config_fpath)
  print(config_file.sections())
  # the settings section should be the first and only section
  settings_section = config_file.sections()[0]
  print(settings_section)
  # pull the settings section and make it a dictionary
  config_settings = dict(config_file[settings_section])
  print(config_settings)



if __name__ == '__main__':
    # main()
    # test()
    parse_args()




