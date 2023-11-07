#   Worker for preparing the AEM3D Lake Model Inputs
#       Initial coding branches from the prep_RCA_worker content
#
#   Determine runtime specifics -
#   - the bay of interest (Missisquoi, St. Albans)
#   - the hydrology model (SWAT, RHESSys)

from lib import cd, logger, IAMBAY
from .AEM3D_prep_IAM import *
from sh import cp
import os
import sys
import datetime
import traceback


def is_num(value):
    try:
        int_ed = int(value)

        return True
    except:
        return


# def get_output_file_name(year):
#     global SCENARIO

#     return '{}-AEM3D-inputs.tar.gz'.format(SCENARIO.id)

def main():
    prep_path = 'aem3d-run'

    # for i in range(len(sys.argv)):
    #     if sys.argv[i] == '--aem3d-dir':
    #         logger.info(f'Setting path to {sys.argv[i+1]}')
    #         prep_path=sys.argv[i+1]

    # Create prep_path if doesn't exist
    # if not os.path.exists(prep_path):
    #     logger.info(f'Creating prep_path at {prep_path}')
    #     os.makedirs(prep_path)
    
    #THEBAY = IAMBAY(settings['whichbay'])   # Create Bay Object for Bay specified
    THEBAY = IAMBAY(bayid='ILS')   # Create Bay Object for Bay specified

    # Make today and today at midnight
    today = datetime.date.today()
    #today = datetime.date(2023,9,13)
    todayMidnight = datetime.datetime.combine(today, datetime.datetime.min.time())

    THEBAY.FirstDate = datetimeToOrdinal(datetime.datetime.combine(datetime.date(2023,1,2), datetime.datetime.min.time()))
    THEBAY.LastDate = datetimeToOrdinal(todayMidnight + datetime.timedelta(days=7))

    ## Need dataframes for hydrology from Missisquoi, Mill, JewittStevens
    
    # cp('p_reductions.csv', 'AEM3D-file-prep/')

    # template_files = [
    #     file
    #     for file in pkg_resources.resource_listdir(
    #         'workers.prep_aem3d_worker', 'resources')
    #     if file.endswith('.txt')
    # ]

    # for f in template_files:
    #     with open(os.path.join('AEM3D-file-prep','TEMPLATES',f), 'w') as output_file:
    #         output_file.write(
    #             pkg_resources.resource_string('workers.prep_aem3d_worker.resources', f).decode('utf-8')
    #         )

    THEBAY.infile_dir = os.path.join(prep_path, 'infiles')
    THEBAY.template_dir = os.path.join(prep_path, 'TEMPLATES')
    THEBAY.run_dir = prep_path

    # Probably don't need this with the cp anymore...
    # if not os.path.exists(prep_path):
    #     os.makedirs(prep_path)
    
    # Copy current aem3d run template
    cp('-R', '/netfiles/ciroh/models/aem3d/current/AEM3D-inputs', prep_path)

    with cd('.'):

        try:
            # Create Dir for infiles
            if not os.path.exists(THEBAY.infile_dir):
                os.makedirs(THEBAY.infile_dir)

            # source the python file prep script
            preprc = AEM3D_prep_IAM(forecastDate=today, theBay=THEBAY)

        except Exception as e:
            logger.info('AEM3D_prep_IAM.py failed. Exiting.')
            logger.info(traceback.print_exc())
            sys.exit(1)

        # # build this tar with everything under the run directory, but without the run directory itself
        # tar(
        #     'czf',
        #     '%s-AEM3D-inputs.tar.gz' % SCENARIO.id,
        #     '-C'
        #     'AEM3D-inputs',
        #     '.'
        # )
        # mv('%s-AEM3D-inputs.tar.gz' % SCENARIO.id, '../')

    logger.info('Fin')

if __name__ == '__main__':
    main()