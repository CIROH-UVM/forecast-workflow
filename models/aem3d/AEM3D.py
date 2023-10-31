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
  pass
  
def main():
  # settings = parse_args(sys.argv[1:])
  pass

  #Create model

def test():
  print('hello world')

if __name__ == '__main__':
    # main()
    test()




