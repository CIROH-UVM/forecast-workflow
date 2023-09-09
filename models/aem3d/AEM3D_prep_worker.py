#   Worker for preparing the AEM3D Lake Model Inputs
#       Initial coding branches from the prep_RCA_worker content
#
#   Determine runtime specifics -
#   - the bay of interest (Missisquoi, St. Albans)
#   - the hydrology model (SWAT, RHESSys)

from ...lib import cd, logger, IAMBAY
from .AEM3D_prep_IAM import *
from sh import cp
import os
import sys
from string import Template
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
    #THEBAY = IAMBAY(settings['whichbay'])   # Create Bay Object for Bay specified
    THEBAY = IAMBAY(bayid='ILS')   # Create Bay Object for Bay specified

    today = datetime.date.today()

    THEBAY.FirstDate = datetimeToOrdinal(today - datetime.timedelta(days=90))
    THEBAY.LastDate = datetimeToOrdinal(today + datetime.timedelta(days=7))

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

    with cd('AEM3D-file-prep'):

        try:
            # source the python file prep script
            preprc = AEM3D_prep_IAM(forecastDate=today, whichbay = THEBAY)

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
