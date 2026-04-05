import os
def make_final_application_shell(rttovCoef, varShape, satChannels, rttov_install_path):
    with open('modules14plus/run_era5_example_fwd.sh', 'r') as file:
        runScript = file.readlines()
    file.close()
    print("Creating final application shell script (run_era5_example_fwd.sh) for RTTOV forward simulation with ERA5 input data...")
    for index, line in enumerate(runScript):
        if line.strip().startswith('COEF_FILENAME='):
            runScript[index] = 'COEF_FILENAME='+rttovCoef+'\n'
        if line.strip().startswith('NLEVELS='):
            runScript[index] = 'NLEVELS='+varShape+'\n'
        if line.strip().startswith('CHAN_LIST='):
            runScript[index] = 'CHAN_LIST="'+satChannels+'"\n'
        if line.strip().startswith('TEST_DIR='):
            runScript[index] = 'TEST_DIR='+rttov_install_path+'/rttov_test/test_example.1\n'
    with open('run_era5_example_fwd.sh', 'w') as rttovRunFile:
        rttovRunFile.writelines(runScript)
        os.chmod('run_era5_example_fwd.sh', 0o755)
    rttovRunFile.close()