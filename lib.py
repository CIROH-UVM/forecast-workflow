import os


class IAMBAY:
    '''
    IAMBAY Class contains the most common information used to describe a lake bay.

        bayid - Which Bay
        sourcelist - list of flow source IDs
        hydromodel - which hydrology model feeds the bay
    '''

    def __init__(self, bayid='ILS', year=2000, **kwargs):
        '''
        Init.
        '''
        self.bayid = bayid
        self.year = year
        self.FirstDate = str(year)+'001.000'  # default first date of data
        self.LastDate =  str(year)+'365.000'  # default last date of data
        self.wrfgridxy = (1,1)      # a default grid selection
        self.infile_dir = 'AEM3D-inputs/infiles'            # directory for boundary condition inputs
        self.run_dir = 'AEM3D-inputs'     # directory to contain everything the lake model needs to run
        self.bayfiles = []                           # initial (empty) list of boundary condition files for bay modeling
        self.flowdf =  None       # eventually to contain bay input flow time series dataframe
        self.tempdf = None        # will contain wtr_temp time series dataframe



        # list of water source IDs for multiple bays
        sourcelist_MAP = {
            'MB': {'201', '202', '203', '204', '21', '22'},
            'STA': {'17', '19'},
            'ILS': {'201', '202', '203', '204', '21', '22', '17', '19' }
        }

        # which WRF grid cell corresponds to lake data
        wrfgrid_MAP = {
            'MB': (32,45),
            'STA': (31,38),
            'ILS': (31,38)
        }
        
        climateZones_MAP = {
            'MB': [{}],
            'STA': [{}],
            'ILS': [{'vars': ['T2', 'AEMLW', 'SWDOWN', 'U10', 'V10'],
                     'zones': [{
                        'name': '401',
                        'desc': 'MB',
                        'lat': 45.035531,
                        'lon': -73.127665},
                        {
                        'name': '402',
                        'desc': 'SAB',
                        'lat': 44.790926,
                        'lon': -73.155695},
                        {
                        'name': '403',
                        'desc': 'IS',
                        'lat': 44.778645,
                        'lon': -73.17273}]},
                    {'vars': ['RH2', 'RAIN'],
                     'zones': [{
                         'name': '0',
                         'desc': 'Unified',
                         'lat': 44.778645,
                         'lon': -73.17273}]}
                     ]}
                     
        #
        #  bay source list structure
        #      map of bay names
        #           bayname: map of source IDs to watershed and the flow proportion of watershed for that sourceID
        #               Watershed designations:
        #                   MS - Missisquoi (rhessys or swat)
        #                   ML - Mill (swat)
        #                   JS - Jewett Stevens (swat)
        sourcelist_STRUCT = {
            'MB': {
                '201': {'name': 'MissisquoiRiverDC', 'wshed': 'msflow', 'prop': 0.57},      # proportions of Miss River Flow
                '202': {'name': 'MissisquoiRiverNE', 'wshed': 'msflow', 'prop': 0.15},
                '203': {'name': 'MissisquoiRiverCE', 'wshed': 'msflow', 'prop': 0.14},
                '204': {'name': 'MissisquoiRiverNW', 'wshed': 'msflow', 'prop': 0.14},
                '21':  {'name': 'RockRiver',         'wshed': 'msflow', 'prop': 0.038},   # .03 of Total Miss Bay Inflow
                '22':  {'name': 'PikeRiver',         'wshed': 'msflow', 'prop': 0.228}    # .18 of Total Miss Bay Inflow
            },
            'STA': {
                '17': {'name': 'MillRiver', 'wshed': 'mlflow', 'prop': 1.0},               # millriver flow
                '19': {'name': 'JewettStevens', 'wshed': 'jsflow', 'prop': 1.0}             # JewettStevens flow
            },
            'ILS': {
                '201': {'name': 'MissisquoiRiverDC', 'wshed': 'msflow', 'prop': 0.57},      # proportions of Miss River Flow
                '202': {'name': 'MissisquoiRiverNE', 'wshed': 'msflow', 'prop': 0.15},
                '203': {'name': 'MissisquoiRiverCE', 'wshed': 'msflow', 'prop': 0.14},
                '204': {'name': 'MissisquoiRiverNW', 'wshed': 'msflow', 'prop': 0.14},
                '21':  {'name': 'RockRiver',         'wshed': 'msflow', 'prop': 0.038},   # .03 of Total Miss Bay Inflow
                '22':  {'name': 'PikeRiver',         'wshed': 'msflow', 'prop': 0.228},    # .18 of Total Miss Bay Inflow
                '17':  {'name': 'MillRiver',         'wshed': 'mlflow', 'prop': 1.0},    #  about 1/3 volume of Rock
                '19':  {'name': 'JewettStevens',     'wshed': 'jsflow', 'prop': 1.0}     #  about the volume of Rock
            }
        }



        # hydrology model specifications by bay ID
        hydromodel_MAP = {
            'MB': 'rhessys', # rhessys model
            'STA': 'swat',  # swat model
            'ILS': 'swat' # swat and rhessys models needed
        }


        if(bayid == 'MB' or bayid == 'STA' or bayid == 'ILS'):
            self.sourcelist = sourcelist_MAP.get(bayid)     # map of flow source IDs for specified bay
            self.sourcemap = sourcelist_STRUCT.get(bayid)   # structure giving name and proportion for sources
            self.hydromodel = hydromodel_MAP.get(bayid)     # the hydrology model associated with bay sources
            self.wrfgridxy = wrfgrid_MAP.get(bayid)         # which wrf climate grid coordinates
            self.climateZones = climateZones_MAP.get(bayid) # climate zones for weather inputs
        else:
            raise Exception('Bay ID is neither "MB", "STA", nor "ILS" ')


        self.validate()

    def validate(self):
        if True == False:
            raise Exception('The Bay is Busted')

    #
    #   method to keep track of files used to model the bay
    #       fname : unpathed name of file
    #       ftype : descriptor of file type (default = "boundary_condition_file")
    def addfile(self, fname, ftype = 'boundary_condition_file'):

        self.bayfiles.append((fname, ftype))
    ##
    #       End of IAMBAY Class
    ##


@contextmanager
def cd(path):
    '''
    A context manager for changing directories.

    Will change the current working directory to path, and then once the
    context exits, change back to the original working directory, optionally
    re-raising an exception that was generated along the way.
    '''
    origin = os.getcwd()
    os.chdir(path)

    try:
        yield
    except Exception as e:
        os.chdir(origin)
        raise e

    os.chdir(origin)

class IAMLogger(logging.Logger):

    logging_initialized = False

    @classmethod
    def setup_logging(cls, name=None, default_level=logging.INFO):
        log_name = 'workers.lib'
        calling_package = get_calling_package()
        if(calling_package is not None):
            log_name = calling_package.split('.')[-1]
        logging_config = {
            "version": 1,
            "disable_existing_loggers": False,
            "formatters": {
                "simple": {
                    "format": "%(asctime)s - %(levelname)s - %(name)s - %(message)s"
                }
            },
            "handlers": {
                "console": {
                    "class": "logging.StreamHandler",
                    "level": "DEBUG",
                    "formatter": "simple",
                    "stream": "ext://sys.stdout"
                },
                "info_file_handler": {
                    "class": "logging.handlers.RotatingFileHandler",
                    "level": "INFO",
                    "formatter": "simple",
                    "filename": f"{log_name}.log",
                    "maxBytes": 10485760,
                    "backupCount": 20,
                    "encoding": "utf8"
                }
            },
            "root": {
                "level": "INFO",
                "handlers": ["console", "info_file_handler"]
            }
        }
        logging.config.dictConfig(logging_config)

        logger = logging.getLogger(log_name)
        logger.info('IAMlogging initialized')
        cls.logging_initialized = True

    @classmethod
    def getLogger(cls, name):
        if(not cls.logging_initialized):
            cls.setup_logging(name)

        return logging.getLogger(name)


IAMLogger.setup_logging()
logger = logging.getLogger(get_calling_package())
