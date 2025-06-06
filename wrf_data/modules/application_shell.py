import os
def make_final_application_shell(rttovCoef, varShape, satChannels, rttov_install_path):
    with open('modules/run_wrf_example_fwd.sh', 'r') as file:
        runScript = file.readlines()
    file.close()
    for index, line in enumerate(runScript):
        line = line.lstrip() #remove indentation spaces
        if line.startswith('COEF_FILENAME='):
            runScript[index] = 'COEF_FILENAME='+rttovCoef+'\n'
        if line.startswith('NLEVELS='):
            runScript[index] = 'NLEVELS='+varShape+'\n'
        if line.startswith('CHAN_LIST='):
            runScript[index] = 'CHAN_LIST="'+satChannels+'"\n'
        if line.startswith('TEST_DIR='):
            runScript[index] = 'TEST_DIR='+rttov_install_path+'/rttov_test/test_example.1\n'
        if line.startswith('homeDir='):
            runScript[index] = 'homeDir='+os.getcwd()+'/\n'
    with open('run_wrf_example_fwd.sh', 'w') as rttovRunFile:
        rttovRunFile.writelines(runScript)
        os.chmod('run_wrf_example_fwd.sh', 0o755)
    rttovRunFile.close()

def make_final_dust_application_shell(rttovCoef, varShape, satChannels, rttov_install_path, aerosol_coeff_path):
    with open('modules/run_wrfchem_dust_example_fwd.sh', 'r') as file:
        runScript = file.readlines()
    file.close()
    for index, line in enumerate(runScript):
        line = line.lstrip() #remove indentation spaces
        if line.startswith('COEF_FILENAME='):
            runScript[index] = 'COEF_FILENAME='+rttovCoef+'\n'
        if line.startswith('NLEVELS='):
            runScript[index] = 'NLEVELS='+varShape+'\n'
        if line.startswith('CHAN_LIST='):
            runScript[index] = 'CHAN_LIST="'+satChannels+'"\n'
        if line.startswith('TEST_DIR='):
            runScript[index] = 'TEST_DIR='+rttov_install_path+'/rttov_test/test_example.1\n'
        if line.startswith('homeDir='):
            runScript[index] = 'homeDir='+os.getcwd()+'/\n'
        if line.startswith('AER_COEF_FILENAME='):
            runScript[index] = 'AER_COEF_FILENAME="'+aerosol_coeff_path+'"\n'
    with open('run_wrfchem_dust_example_fwd.sh', 'w') as rttovRunFile:
        rttovRunFile.writelines(runScript)
        os.chmod('run_wrfchem_dust_example_fwd.sh', 0o755)
    rttovRunFile.close()
