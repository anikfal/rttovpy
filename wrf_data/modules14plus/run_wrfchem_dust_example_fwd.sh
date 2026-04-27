#!/bin/sh
# Script to run the example_fwd.F90 example program
# 

# Set BIN directory if supplied
BIN=$(perl -e 'for(@ARGV){m/^BIN=(\S+)$/o && print "$1";}' $*)
if [ "x$BIN" = "x" ]
then
  BIN=bin
fi

homeDir=/home/anikfal/WRFDA/atmospheric_science/rttov/wrf_data/
TEST_DIR=/home/anikfal/WRFDA/rttov12/rttov_test/test_example.1
wrffilename=`basename $(sed -n '/wrf_file_path/ p' namelist_wrf.yaml | awk '{print $2}')`
inputDirSuffix=`sed -n '/rttov_inputdata_directory_suffix/ p' namelist_wrf.yaml | awk '{print $2}'`
rttov_output_dir_name=`sed -n '/rttov_outputdata_directory_suffix/ p' namelist_wrf.yaml | awk '{print $2}'`
profile_directory=$homeDir$wrffilename"_"$inputDirSuffix
outputDir=$wrffilename"_"$rttov_output_dir_name
if [ -d "$outputDir" ]; then
    rm -rf "$outputDir"
fi
mkdir $outputDir

for myprofile in $profile_directory/prof-*.dat; do
  AER_PROF_FILENAME="aer_prof-"`basename $myprofile | cut -d "-" -f 2`
  ARG_ARCH=$(perl -e 'for(@ARGV){m/^ARCH=(\S+)$/o && print "$1";}' $*)
  if [ ! "x$ARG_ARCH" = "x" ]; then
    ARCH=$ARG_ARCH
  fi
  if [ "x$ARCH" = "x" ];
  then
    echo 'Please supply ARCH and run again.'
    echo 'For example, if RTTOV is compiled with gfortran, run it as below:'
    echo "$0 ARCH=gfortran"
    exit
  fi

  filePath=`realpath $myprofile`
  profileName=`basename $myprofile`
  echo Simulating based on $filePath
  ln -sf $filePath $TEST_DIR
  ln -sf $profile_directory/$AER_PROF_FILENAME $TEST_DIR

  # Test case input data
  COEF_FILENAME=/home/anikfal/WRFDA/rttov12/rtcoef_rttov12/rttov7pred54L/rtcoef_noaa_16_avhrr.dat
  AER_COEF_FILENAME="scaercoef_msg_3_seviri.dat"
  if [ ! -f $COEF_FILENAME ]; then
    echo "Coef file $COEF_FILENAME not found, aborting..."
    exit 1
  fi

  # PROF_FILENAME="prof.dat"                  # Input profile(s), usually found in $TEST_DIR set below
  PROF_FILENAME=$profileName                  # Input profile(s), usually found in $TEST_DIR set below
  NPROF=1                                   # Number of profiles defined in prof.dat
  NLEVELS=37
  DO_SOLAR=0                                # 0 = solar off / 1 = solar on
  NCHAN=3                                   # Number of channels to simulate for each profile
  CHAN_LIST="1 2 3"
  NTHREADS=1                                # Number of threads to use (compile RTTOV with OpenMP to exploit this)

  # Paths relative to the rttov_test/${TEST_DIR} directory:
  BIN_DIR=../../$BIN                           # BIN directory (may be set with BIN= argument)

  CWD=$(pwd)
  cd $TEST_DIR

  echo " "
  echo " "
  echo " Test forward "
  echo " "

# $BIN_DIR/example_fwd.exe << EOF
# "${COEF_FILENAME}", Coefficient filename
# ${PROF_FILENAME},   Input profile filename
# ${NPROF}        ,   Number of profiles
# ${NLEVELS}      ,   Number of levels
# ${DO_SOLAR}     ,   Turn solar radiation on/off
# ${NCHAN}        ,   Number of channels
# ${CHAN_LIST}    ,   Channel numbers
# ${NTHREADS}     ,   Number of threads
# EOF

$BIN_DIR/example_aer_file_fwd.exe << EOF
"${COEF_FILENAME}"    , Coefficient filename
"${AER_COEF_FILENAME}", Aerosol coefficient filename
${PROF_FILENAME}      , Input profile filename
${AER_PROF_FILENAME}  , Input aerosol profile filename
${NPROF}              , Number of profiles
${NLEVELS}            , Number of levels
${DO_SOLAR}           , Turn solar radiation on/off
${NCHAN}              , Number of channels
${CHAN_LIST}          , Channel numbers
${NTHREADS}           , Number of threads
EOF

if [ $? -eq 0 ]; then
  cd $CWD
  mv $TEST_DIR/output_example_aer_file_fwd.dat $outputDir/output_example_fwd.dat_$profileName
else
  echo "Profile data " $profileName " has an issue (zenith angle > 75, etc)."
  echo "Skipping this profile, and filling with missing value."
  echo "missing_value" > $CWD/$outputDir/output_example_fwd.dat_$profileName
fi
cd $CWD
rm $TEST_DIR/$profileName

done
