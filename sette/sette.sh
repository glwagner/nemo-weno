#!/bin/bash
# ======================================================================
#                       ***  SCRIPT  sette.sh  ***
#  SETTE: SETTE test-run setup
# ======================================================================
#
# ----------------------------------------------------------------------
# NEMO/SETTE 5.1.a, NEMO Consortium (2026)
# Software governed by the CeCILL license (see ./LICENSE.txt)
# ----------------------------------------------------------------------

set +x

# initialise user dependent variable
export cmd=$0 ; export cmdargs=$@
SETTE_DIR=$(cd $(dirname "$0"); pwd)
MAIN_DIR=$(dirname $SETTE_DIR)
export CMPL_CORES=8            # Number of threads to use for compiling
export SETTE_STG="_ST"         # Base suffix to append to configuration name
NEMO_DEBUG=""
dry_run=0
SETTE_REPORT=0
WAIT_SETTE=0
export DOTENV_FILE=""

# source sette functions
. ./all_functions.sh

# Controls that determine the SETTE build and run-time-configuration variant
#
# Controls for common compile-time keys:
export USING_QCO='yes'         # Default: yes => add key_qco            ; use -q to delete key_qco
export USING_SI3_1D='no'       # Default: no  => add key_si3_1D         ; hard-coded option
export USING_XIOS='yes'        # Default: yes => add key_xios           ; use -X to delete key_xios
                               #    Note: changing USING_XIOS may require a change in arch file
# Control of compilation options and associated run-length adjustment:
export USING_DEBUG='no'        # Default: no  => default compil. options; use -b to select debug compilation options
                               # and shortened test runs
# Controls for common run-time options:
export USING_TIMING='yes'      # Default: yes => set ln_timing=.true.   ; use -T to disable
export USING_ICEBERGS='yes'    # Default: yes => set ln_icebergs=.true. ; use -i to disable
export USING_ABL='no'          # Default: no  => set ln_abl=.false.     ; use -a to set ln_abl=.true.
export USING_EXTRA_HALO='no'   # Default: no  => set nn_hls=2           ; use -e to set nn_hls=3
export USING_COLLECTIVES='yes' # Default: yes => set nn_comm=2          ; use -C to set nn_comm=1
export USING_NOGATHER='yes'    # Default: yes => set ln_nnogather=.true.; use -N to set ln_nnogather=.false.
export USING_TILING='yes'      # Default: yes => set ln_tile=.true.     ; use -t to disable
export USING_MPMD='yes'        # Default: yes => run with detached XIOS servers ; use -A to run in attached (SPMD) mode
                               #    Note: 'yes' also ensures key_xios but -A will not remove it

# Controls that affect SETTE behaviour
export USER_INPUT='yes'        # Default: yes => request user input on decisions. For example:
                               #                 1. regarding mismatched options
                               #                 2. regardin incompatible options
                               #                 3. regarding creation of directories

# Avaliable source-code transformation options: currently only 'passthrough'
SCTRANSFORMS=(passthrough)
# Default source-code transformation option: none
TRANSFORM=""

# Parse command-line arguments
if [ $# -gt 0 ]; then
  while getopts n:x:g:cybrshTqQteiACFNXuawm:z:p:d: option; do

     case $option in
        p) export SETTE_TEST_PHASES=($OPTARG)
           echo "-p: SETTE (sub)phase(s) ${SETTE_TEST_PHASES[@]} selected"
           echo "";;
        x) export SETTE_TEST_TYPES=(${OPTARG})
           echo "-x: SETTE test types ${SETTE_TEST_TYPES[@]} selected"
           echo "";;
        c) export SETTE_CLEAN_CONFIGS='yes'
           export SETTE_SYNC_CONFIGS='yes'
           echo "-c: Configuration(s) ${SETTE_TEST_CONFIGS[@]} will be cleaned; this option enforces also synchronisation"
           echo "";;
        y) dry_run=1
           echo "";;
        r) SETTE_REPORT=1
           echo "-r: Sette report will be printed once jobs are finished"
           WAIT_SETTE=1
           echo "";;
        w) WAIT_SETTE=1
           echo "-w: Sette will wait for jobs to finish"
           echo "";;
        s) export SETTE_SYNC_CONFIGS='yes'
           echo "-s: MY_SRC and EXP00 in ${SETTE_TEST_CONFIGS[@]} will be synchronised with the MY_SRC and EXPREF from the reference configuration"
           echo "";;
        n) OPTSTR="$OPTARG"
           SETTE_TEST_CONFIGS=(${OPTSTR})
           echo "-n: Configuration(s) ${SETTE_TEST_CONFIGS[@]} will be tested if they are available"
           echo "";;
        g) case $OPTARG in
             [0-9,a-z,A-Z] ) echo "-g: Using ${SETTE_STG}${OPTARG} as the configuration suffix";;
             * ) echo "-g only accepts a single, alphanumeric character. Processing halted"; exit 42;;
           esac
           export SETTE_STG=${SETTE_STG}${OPTARG}
           echo "";;
        u) export USER_INPUT='no'
           echo "-u: sette.sh will not expect any user interaction == no safety net!" 
           echo "";;
        m) export SETTE_COMPILER=($OPTARG)
           echo "-m: $SETTE_COMPILER archfile and batch template will be used"
           echo "";;
        z) for sct in ${SCTRANSFORMS[@]}; do [[ "${OPTARG}" == "${sct}" ]] && TRANSFORM="${sct}"; done
	   if [[ -n ${TRANSFORM} ]]; then
             echo "-z: source-code transformation option ${TRANSFORM} selected"
           else
             echo "-z: source-code transformation unavailable"; exit 42
	   fi
	   echo "";;
	# SETTE controls
        q) export USING_QCO='no'
           echo "-q: key_qco and key_linssh will NOT be activated"
           echo "";;
        X) export USING_XIOS='no'
           echo "-X: key_xios will not be activated"
           echo "";;
        b) export USING_DEBUG='yes'
           echo "-b: NEMO will be compiled with DEBUG options if available in ARCH file and test runs will be shortened"
           echo "";;
        T) export USING_TIMING='no'
           echo "-T: ln_timing will be set to false"
           echo "";;
        i) export USING_ICEBERGS='no'
           echo "-i: ln_icebergs will be set to false"
           echo "";;
        a) export USING_ABL='yes'
           echo "-a: ln_abl will be set to true"
           echo "";;
        e) export USING_EXTRA_HALO='yes'
           echo "-e: nn_hls will be set to 3"
           echo "";;
        C) export USING_COLLECTIVES='no'
           echo "-C: nn_comm will be set to 1"
           echo "";;
        N) export USING_NOGATHER='no'
           echo "-N: ln_nnogather will be set to false"
           echo "";;
        t) export USING_TILING='no'
           echo "-t: ln_tile will be set to false"
           echo "";;
        A) export USING_MPMD='no'
           echo "-A: Tasks will be run in attached (SPMD) mode"
           echo "";;
        d) export DOTENV_FILE="${OPTARG}"
           echo "-d: sette will use ${DOTENV_FILE} .env file instead of .git directory to define related variables"
           echo "";;
        # Usage message
        h | *) echo 'sette.sh with no arguments (in this case all configuration will be tested with default options)'
               echo '-p space-separated list of SETTE test phases (if omitted, the available'
               echo '   phases, COMPILE and RUN, are selected); each phase comprises two'
               echo '   subphases which can be selected explicitly (COMPILE_BASELINE,'
               echo '   COMPILE_VARIANTS, RUN_REFERENCE, and RUN_TESTS for the compilation of'
               echo '   baseline test configurations, the compilation of test-configuration'
               echo '   variants, the execution of reference test runs, and the execution of'
               echo '   additional test runs as required for the selected test types,'
               echo '   respectively'
               echo '-x space-separated list of SETTE test types (if omitted, all available types,'
               echo '   RESTART, REPRO, PHYOPTS, ROTSYM, COUPLING, and VARIANTS, are selected)'
               echo '-T to set ln_timing false for configurations (default: true)'
               echo '-t set ln_tile false in all tests that support it (default: true)'
               echo '-e set nn_hls=3 (default: nn_hls=2)'
               echo '-i set ln_icebergs false (default: true)'
               echo '-a set ln_abl true (default: false)'
               echo '-C set nn_comm=1 (default: nn_comm=2 ==> use MPI3 collective comms)'
               echo '-N set ln_nnogather false for ORCA2 configurations (default: true)'
               echo '-q to remove the key_qco key (default: added)'
               echo '-X to remove the key_xios key (default: added)'
               echo '-A to run tests in attached (SPMD) mode (default: MPMD with key_xios)'
               echo '-n "CFG1_to_test CFG2_to_test ..." to test some specific configurations'
               echo '-g "group_suffix" single character suffix to be appended to the standard _ST suffix used'
               echo '                  for SETTE-built configurations (needed if sette.sh invocations may overlap)'
               echo '-r to execute without waiting to run sette_rpt.sh at the end (useful for chaining sette.sh invocations)'
               echo '-y to perform a dryrun to simply report what settings will be used'
               echo '-b to compile NEMO with debug options (only if %DEBUG_FCFLAGS if defined in your arch file)'
               echo '-c to clean each configuration'
               echo '-s to synchronise the sette MY_SRC and EXP00 with the reference MY_SRC and EXPREF'
               echo '-w to wait for Sette jobs to finish'
               echo '-r to print Sette report after Sette jobs completion'
               echo '-d to use .env file instead of .git directory to define related variables'
               echo '-u to run sette.sh without any user interaction. This means no checks on creating'
               echo '          directories etc. i.e. no safety net!'
               echo '-m to select computing architecture file (arch_<machine_name>.fcm) and batch template to run SETTE'
               echo '-z to select a source-code transformation option' ; exit 42 ;;
     esac
  done
  shift $((OPTIND - 1))
fi

# Get SETTE parameters
. ./param.default
if [ -f ./param.cfg ]; then
   . ./param.cfg
else
   echo "warning: \"param.cfg\" file not found; SETTE will use default paramaters from \"param.default\" file"
fi

# Define sette variables from local .git repository (default) or .env file ("-d" option)
if [ -z "${DOTENV_FILE}" ]; then
  set_git_var
else
  set_dotenv_var
fi

# Retrieval of the source-code revision: it will be used as the primary
# identifier for accessing the SETTE validation database; a '+' suffix
# indicates local changes with regard to the upstream source-code revision
# indicated by the hash-function value
[ ${localchanges} -gt 0 ] && VALID_REV="${VALID_REV}+"

# Check the validity of the compilation-environment name, the
# source-code-transformation name, and the source-code-revision identifier:
#
#   compilation-environment name
[[ ! "${COMPILER}" =~ ^[[:alnum:].+_-]{1,64}$ ]] && echo "Error: incompatible compilation-environment name" && exit 1
found=$( find ${MAIN_DIR}/arch/ -name "arch-${COMPILER}.fcm" | wc -l )
[[ ! "${COMPILER}" == "auto"* ]] && [[ "${found}" == "0" ]] && echo "Error: compilation-environment is not available" && exit 1
#   source-code-transformation identifier
if [[ -n "${TRANSFORM}" ]] && [[ ! "${TRANSFORM}" =~ ^[[:alnum:]]{1,64}$ ]]; then
  echo "Error: incompatible name of the source-code-transformation option"
  exit 1
fi
#   source-code-revision identifier
VALID_REV=$(echo ${VALID_REV} | tr '[:upper:]' '[:lower:]' | tr -d -c '[:xdigit:][+]')
if [[ ${#VALID_REV} -lt 8 ]] || [[ ${#VALID_REV} -gt 40 ]]; then
  echo "Error: incompatible source-code-revision identifier"
  exit 1
fi

# If requested, select NEMO debug build options
[[ "${USING_DEBUG}" == "yes" ]] && NEMO_DEBUG="-b"
# Also select the SETTE debug-compilation option if the name of the architecture-configuration file suggests the use of debug
# compilation options
COMPILER_L=$(echo ${COMPILER} | tr '[:upper:]' '[:lower:]')
if [[ ${COMPILER_L} =~ ("debug"|"dbg") ]]; then
  echo "warning: the name of the architecture-configuration file suggests the use of"
  echo "         debug compilation options, therefore SETTE test runs will be shortened"
  export USING_DEBUG='yes'
fi

# Set the common compile keys to add or delete based on command-line arguments:
#
export ADD_KEYS="" ; export DEL_KEYS=""
if [ ${USING_XIOS} == "yes" ] ; then export ADD_KEYS="${ADD_KEYS}key_xios " ; fi
if [ ${USING_XIOS} == "no" ]  ; then export DEL_KEYS="${DEL_KEYS}key_xios " ; fi
#
if [ ${USING_QCO} == "yes" ] ; then export ADD_KEYS="${ADD_KEYS}key_qco " ; fi
if [ ${USING_QCO} == "no" ]  ; then export DEL_KEYS="${DEL_KEYS}key_qco key_linssh " ; fi
#
if [ ${USING_SI3_1D} == "yes" ] ; then export ADD_KEYS="${ADD_KEYS}key_si3_1D " ; fi
if [ ${USING_SI3_1D} == "no" ]  ; then export DEL_KEYS="${DEL_KEYS}key_si3_1D " ; fi

#
# Set validation record sub-directories (if required)
#
if [ ! -d $NEMO_VALIDATION_DIR ] ; then
 if [ ${dry_run} -eq 0 ] ; then
  if [ ${USER_INPUT} == "yes" ] ; then
   while true; do
       read -p "$NEMO_VALIDATION_DIR does not exist. Do you wish to create it? " yn
       case $yn in
           [Yy]* ) echo "Ok, creating $NEMO_VALIDATION_DIR"; mkdir $NEMO_VALIDATION_DIR; break;;
           [Nn]* ) echo "Ok, exiting instead"; exit 42;;
           * ) echo "Please answer yes or no.";;
       esac
   done
  else
       # Without user input, carry on regardless
       echo "$NEMO_VALIDATION_DIR does not exist. It will be created"
       mkdir $NEMO_VALIDATION_DIR
  fi
 else
  echo "$NEMO_VALIDATION_DIR does not exist"
  echo "but this is a dry run so it will not be created"
 fi
fi
export NEMO_VALIDATION_DIR=${NEMO_VALIDATION_DIR}
if [ ${dry_run} -eq 0 ] && [ ! -d ${NEMO_VALIDATION_DIR}/${VALID_REV}/ ]; then
  [ -a ${NEMO_VALIDATION_DIR}/${VALID_REV} ] && echo "Error: spurious file found in SETTE validation database" && exit 1
  mkdir ${NEMO_VALIDATION_DIR}/${VALID_REV}
fi

# Encode the control variables that determine the SETTE build and
# run-time-configuration variant in a string; a hash-function value of this
# string will be used as SETTE variant identifier (variable 'VALID_VAR') for
# grouping SETTE-test output in the validation database at
#   '${NEMO_VALIDATION_DIR}/${VALID_REV}/${VALID_VAR}/',
# and the original string will be recoreded at
#   '${NEMO_VALIDATION_DIR}/${VALID_REV}/${VALID_VAR}.lookup'
# for lookup purposes.
VAR=""
VAR="\"COMPILER\",\"${COMPILER}\";\"TRANSFORM\",\"${TRANSFORM}\""
VAR="${VAR};\"USING_QCO\",\"${USING_QCO}\""
VAR="${VAR};\"USING_SI3_1D\",\"${USING_SI3_1D}\""
VAR="${VAR};\"USING_XIOS\",\"${USING_XIOS}\""
VAR="${VAR};\"USING_DEBUG\",\"${USING_DEBUG}\""
VAR="${VAR};\"USING_TIMING\",\"${USING_TIMING}\""
VAR="${VAR};\"USING_ICEBERGS\",\"${USING_ICEBERGS}\""
VAR="${VAR};\"USING_ABL\",\"${USING_ABL}\""
VAR="${VAR};\"USING_EXTRA_HALO\",\"${USING_EXTRA_HALO}\""
VAR="${VAR};\"USING_COLLECTIVES\",\"${USING_COLLECTIVES}\""
VAR="${VAR};\"USING_NOGATHER\",\"${USING_NOGATHER}\""
VAR="${VAR};\"USING_TILING\",\"${USING_TILING}\""
VAR="${VAR};\"USING_MPMD\",\"${USING_MPMD}\""
VALID_VAR=$(echo ${VAR} | md5sum | cut -f 1 -d ' ' | tr '[:upper:]' '[:lower:]' | tr -d -c '[:xdigit:]')
if [[ ${#VALID_VAR} == 32 ]]; then
  echo ${VAR} > ${NEMO_VALIDATION_DIR}/${VALID_REV}/${VALID_VAR}.lookup
else
  echo "Error: SETTE-identifier computation failed"
  exit
fi

if [ ${#SETTE_TEST_CONFIGS[@]} -eq 0 ]; then
   echo "=================================="
   echo "Configurations $TEST_CONFIGS will be tested if they are available"
fi
echo "Test phases                       : ${TEST_PHASES[@]}"
echo "Test types                        : ${TEST_TYPES[@]}"
echo "requested by the command          : "$cmd $cmdargs
echo "on revision                       : "${VALID_REV}
VAR2="${VAR#\"}"
VAR2="${VAR2%\"}"
for v in ${VAR2//\";\"/ }; do printf "%-33s : %s\n" ${v/\",\"/ }; done
printf "%-33s : %s\n" "USER_INPUT" "${USER_INPUT}"
printf "%-33s : %s\n" "Common compile keys to be added" "$ADD_KEYS"
printf "%-33s : %s\n" "Common compile keys to be deleted" "$DEL_KEYS"
echo "Validation records to appear under: "$NEMO_VALIDATION_DIR
echo "=================================="
echo ""
#
# Option compatibility tests
#
if [ ${USING_MPMD} == "yes" ] && [ ${USING_XIOS} == "no" ] ; then echo "Incompatible choices. MPMD mode requires the XIOS server" ; exit ; fi

# Set the external SETTE build directory if an external SETTE build location
# has been requested
if [ -n "${CUSTOM_DIR}" ]; then
  export CMP_DIR=${CUSTOM_DIR}/${VALID_REV}_${VALID_VAR}
else
  export CMP_DIR=''
fi

# End of dry run
if [ ${dry_run} -eq 1 ] ; then echo "dryrun only: no tests performed" ; exit ; fi

# run sette on reference configuration
. ./sette_configuration_test.sh
if [[ $? != 0 ]]; then
   echo ""
   echo "--------------------------------------------------------------"
   echo "./sette_configuration_test.sh didn't finish properly, need investigations"
   echo "--------------------------------------------------------------"
   echo ""
   exit 42
fi

# wait for sette jobs to finish
if [[ ${WAIT_SETTE} -eq 1 && "${TEST_TYPES[@]}" =~ (RESTART|REPRO|PHYOPTS|ROTSYM|COUPLING|VARIANTS) ]]; then
   echo ""
   echo "-------------------------------------------------------------"
   echo "wait for sette jobs to finish..."
   echo "-------------------------------------------------------------"
   echo ""
   NRUN=999
   NIT=0
   while [[ $NRUN -ne 0 && $nit -le 1080 ]]; do
      nit=$((nit+1))
      if [[ "${BATCH_STAT}" == "squeue" && -n "${BATCH_LST[@]}" ]]; then
        BATCH_STAT="squeue -j $(echo ${BATCH_LST[@]} | tr ' ' ',') -h -o %j"
        NRUN=$( ${BATCH_STAT} | wc -l )
        echo "currently running jobs: "$(${BATCH_STAT})
      else
        NRUN=$( ${BATCH_STAT} | grep ${BATCH_NAME} | wc -l )
      fi
      if [[ $NRUN -ne 0 ]]; then
         echo $NRUN "sette jobs still in queue or running ..."
         sleep 10
      else
         echo "all sette jobs completed"
         break
      fi
   done
   # check jobs exit codes
   if [[ "${BATCH_STAT}" == "squeue"* && -n "${BATCH_LST[@]}" ]]; then
     BATCH_EXITCODE=( $(sacct --job=$(tr ' ' ',' <<< ${BATCH_LST[@]}) --format=ExitCode --noheader --allocations | cut -d":" -f 1) )
     if [ $(IFS=+; echo "$((${BATCH_EXITCODE[*]}))") -gt 0 ]; then
       echo ""
       echo "---------------------------------------------------------------"
       echo "sette jobs failed: please check jobs log and ocean.output files"
       echo "---------------------------------------------------------------"
       echo ""
       exit 42
     fi
   fi
fi

# run sette report
if [ ${SETTE_REPORT} -eq 1 ]; then
   echo ""
   echo "-------------------------------------------------------------"
   echo "./sette_rpt.sh"
   echo "-------------------------------------------------------------"
   echo ""
   ./sette_rpt.sh ${NEMO_DEBUG} -n "${TEST_CONFIGS[*]}"
   if [[ $? != 0 ]]; then
      echo ""
      echo "-----------------------------------------------------------------"
      echo "./sette_rpt.sh didn't finish properly, need investigations"
      echo "-----------------------------------------------------------------"
      echo ""
      exit 42
   fi
fi

exit 0
