from ...lib import (
    cd,
    logger,
)
from sh import Command, rm


# def activate_restart_file(aem3d_cntl_file_name):
#     '''
#     Activate AEM3D's restart file system by altering a line in the control file.
#     Regular expression finds the irestart line and flips the "0" to a "1"

#     '''

#     if not os.path.exists(aem3d_cntl_file_name):
#         logger.critical(
#             'File {} does not exist. Exiting.'.format(
#                 aem3d_cntl_file_name
#             )
#         )
#         exit_with_code_based_on_host()


#     cp(aem3d_cntl_file_name, aem3d_cntl_file_name+'.template')     ## Copy original file to a template.

#     with open(aem3d_cntl_file_name+'.template') as inFile:

#         with open(aem3d_cntl_file_name, 'w') as outFile:

#             for line in inFile:                                                                                                 ## I think, when Morgan uses perl, he reads the entire file into a string.  You can do that too.
#                 line = re.sub('[0](\s*irestart)', f'1\\1', line)  # flip the 0 to a 1 in irestart line
#                 outFile.write(line)

#     '''
#     replacement = '1                                                       irestart'
#     new_content = ''
#     with open(aem3d_cntl_file_name) as input_file:
#         new_content = input_file.read()

#     if replacement not in new_content:
#         logger.critical(
#             'Restart file functionality is not enabled after a perl pie of the control file. Exiting.'
#         )

#         exit_with_code_based_on_host()
#     '''


# def run_model():

#     args = [
#      ]

#     uuid_log(args)
#     aem3d = Command(AEM3D_BIN)

#     try:
#         for line in aem3d(*args, _iter=True):
#             uuid_log(line)
#     except Exception as e:
#         uuid_log('A run of AEM3D has failed. Exiting.')
#         uuid_log(e)
#         raise Exception('AEM3D failure.')  # Kill pool with an exception

#     uuid_log('AEM3D complete.')


# def setup(settings):
#     if settings['restart']:
#         logger.info('Enabling Lake Model Restart')

#         # Get the IAM id for the preceeding year
#         previous_aem3d_args = scenario.__dict__.copy()
#         previous_aem3d_args['decade'] = previous_aem3d_args['decade'] - 1
#         previous_aem3d = '%s-aem3d.tar.gz' % IAMScenario(**previous_aem3d_args).id
#         tar('xzf', '../'+previous_aem3d,'outfiles/unf/restart_file.unf')

#         #previous_aem3d_dir = previous_aem3d.replace('.tar.gz', '')
#         cp(
#             'outfiles/unf/restart_file.unf',
#             'infiles'
#         )
#         restartnextyear('infiles/restart_file.unf', scenario.decade)   # modify restart file to begin on Jan 1 of current year
#         activate_restart_file('run_aem3d.dat')

#     else:
#         logger.info('Running without lake model restart')


# def pick_aem3d_bin(settings):
#     '''
#     Select the correct aem3d_openmh binary for the executing system.
#     '''
#     global AEM3D_BIN

#     fqdn = socket.getfqdn()
#     if fqdn == 'epscor-pascal.uvm.edu':
#         AEM3D_BIN = '/usr/local/bin/aem3d_openmp'
#     elif fqdn == 'epscor.uvm.edu':
#         AEM3D_BIN = '/usr/local/bin/aem3d_openmp'
#     elif 'ucar' in fqdn:
#         AEM3D_BIN = 'aem3d_openmp'
#     else:
#         raise NotImplementedError('No {} AEM3D binary is available for {0:s}'.format(
#             settings['aem3d_ver'], fqdn
#         ))



def main():
    # settings = parse_args(sys.argv[1:])
    # scenario = IAMScenario.from_id(settings['scenario'])

    # pick_aem3d_bin(settings)
    aem3d = Command('/usr/local/bin/aem3d_openmp')

    with cd('today'):        # the working dir for the run

        # setup(settings)

        try:
            for line in aem3d(*args, _iter=True):
                logger.info(line)
        except Exception as e:
            logger.info('A run of AEM3D has failed. Exiting.')
            logger.info(e)
            raise Exception('AEM3D failure.')  # Kill pool with an exception

        logger.info('AEM3D complete.')

        # Remove large 3D file???
        #rm('-rf', 'outfiles/nc/All3D.nc')

    logger.info('Fin')


if __name__ == '__main__':
    main()
