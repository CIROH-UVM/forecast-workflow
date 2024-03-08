#   Worker for preparing the AEM3D Lake Model Inputs
#       Initial coding branches from the prep_RCA_worker content
#
#   Determine runtime specifics -
#   - the bay of interest (Missisquoi, St. Albans)
#   - the hydrology model (SWAT, RHESSys)

from lib import cd, logger, IAMBAY
from .AEM3D_prep_IAM import *
from .get_args import get_args
from sh import cp
import os
import sys
import datetime as dt
import traceback


def is_num(value):
    try:
        int_ed = int(value)

        return True
    except:
        return

def main():
    prep_path = 'aem3d-run'


    SETTINGS = get_args()

    #THEBAY = IAMBAY(settings['whichbay'])   # Create Bay Object for Bay specified
    THEBAY = IAMBAY(bayid='ILS')   # Create Bay Object for Bay specified


    # These two settings return datetime objs set to midnight already
    THEBAY.FirstDate = datetimeToOrdinal(SETTINGS['spinup_date'])
    THEBAY.LastDate = datetimeToOrdinal(SETTINGS['forecast_end'])


    ## Need dataframes for hydrology from Missisquoi, Mill, JewittStevens

    THEBAY.infile_dir = os.path.join(prep_path, 'infiles')
    THEBAY.template_dir = os.path.join(prep_path, 'TEMPLATES')
    THEBAY.run_dir = prep_path

    # Copy current aem3d run template
    cp('-R', SETTINGS["aem3d_input_dir"], prep_path)

    with cd('.'):

        try:
            # Create Dir for infiles
            if not os.path.exists(THEBAY.infile_dir):
                os.makedirs(THEBAY.infile_dir)

            # source the python file prep script
            preprc = AEM3D_prep_IAM(theBay=THEBAY, settings=SETTINGS)

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