import pandas as pd
import datetime as dt
import re
import os
from .AEM3D import ordinalToDatetime

# AEM3D Input File class to read formatted streamflow data
class AEM3DInputFile:
    header = ''
    data = None
    float_format = None

    def __init__(self, filename):
        with open(filename) as infile:
            num_header_rows = 0
            regex = r'^[0-9]{7}\.[0-9]{4}\s'
            for line in infile:
                if re.search(regex, line) != None or line == '':
                    break
                self.header = self.header + line
                num_header_rows = num_header_rows + 1                
            if re.search(regex, line) != None and line != '':
                decimal_places = len(re.split(r'\s', line)[1].split('.')[1])
                self.float_format = f'%.{decimal_places}f'
        self.data = pd.read_csv(filename, 
                           delimiter=r'\s+',
                           index_col=False,
                           skiprows=num_header_rows-1,
                           dtype={'TIME': 'string'})
        self.data = self.data.set_index(self.data['TIME'].apply(ordinalToDatetime)).rename_axis(['DATETIME'])

    def write(self, filename, overwrite = False):
        if os.path.exists(filename) and not overwrite:
            print("File exists, not writing to file")
            return        
        with open(filename, 'w') as outfile:
            outfile.write(self.header)
            self.data.to_csv(outfile, sep='\t', float_format=self.float_format, header=False, index=False)

    def replace_years(self, year_dict):
        for year, orig_year in year_dict.items():
            new_data = self.data[self.data.index.year == orig_year]
            new_data = new_data.set_index(new_data.index + pd.DateOffset(years = (year - orig_year)))
            new_data['TIME'] = new_data['TIME'].replace(f'^{orig_year}', f'{year}', regex=True)
            self.data = pd.concat([self.data, new_data])

        #print(self.data)
        # Remove any original data points for the new years
        self.data = self.data[~self.data.index.duplicated(keep = 'last')]
        #print(self.data)
        
        # Copy last data point to new year, first day
        last_row = self.data.tail(1).copy()
        last_row.loc[:, 'TIME'] = last_row['TIME'].apply(lambda x: str(int(x[:4]) + 1) + "001.0000").astype('string')
        last_row = last_row.set_index(last_row['TIME'].apply(ordinalToDatetime))
        self.data = pd.concat([self.data, last_row])
        #print(self.data)