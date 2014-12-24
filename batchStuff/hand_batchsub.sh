#!/bin/bash
#SBATCH -A belle        # this is the account to which time is charged --
#                       # use `belle' for belle analysis 
#SBATCH -t 12:00:00     # time limit for job in HH:MM:SS
#SBATCH -N 1            # number of CPU cores requested
#
# Redirect the job's stdout
#SBATCH -o /pic-disk/projects/belle/voss771/ptSpectOut/jobId_%j.out
#
# Redirect the job's stderr
#SBATCH -e /pic-disk/projects/belle/voss771/ptSpectOut/jobId_%j.err
#
# Set the job name to be `ypipi-e000053r000077-b20090127_0910.mdst'
#SBATCH -J ptSpect

# `.brofile' is my equivalent to `.bashrc'
echo Zuhause: $HOME
 source $HOME/.cshrc

# This sets up your environment to run code in the Belle paradigm
source /pic/projects/belle/scripts/general_env.sh

# Add my modifications to the standard Belle environment
source $HOME/custom_setup.sh

# This is necessary to get code that fills HBOOK ntuples running comfortably
# on Olympus
shopt -s expand_aliases
belleset_hbook32768_libs

echo ANSELM_BASF_START `date`
export USE_GRAND_REPROCESS_DATA=1
export BELLE_MESSAGE_LEVEL=INFO
export LD_LIBRARY_PATH=/people/voss771/ptSpect/:`root-config --libdir`:$LD_LIBRARY_PATH
export BASF_MODULE_DIR=/people/voss771/ptSpect/:$BASF_MODULE_DIR
#setenv BASF_NPROCESS 0

echo    $LD_LIBRARY_PATH
echo    $BASF_MODULE_DIR
echo    $PANTHER_TABLE_DIR
MODULE=ptSpect
HBKFILE=/dev/null
# I put the BASF script inline like so
basf <<EOF
path create main
path create analysis
path add_module main fix_mdst
path add_condition main >:0:analysis
path add_condition main <=:0:KILL
path add_module analysis ${MODULE}
module put_parameter ptSpect rfname\/pic/projects/belle/voss771/mcEx55New.root
module put_parameter fix_mdst Make_pi0_option\2
module put_parameter fix_mdst Make_pi0_lower_limit\-5.0
module put_parameter fix_mdst Make_pi0_upper_limit\5.0
module put_parameter fix_mdst Correct_ecl_option\1
initialize
nprocess set 0
histogram define ${HBKFILE}
process_event /pic-disk/projects/belle/ZMDST/MC/MC_4S/MC_4S_EXP55/continuum/uds/evtgen-uds-00-all-e000055r001635-b20090127_0910.mdst 0
output close

terminate

EOF

# Grab the exit code of BASF
BASFRET=$?
echo VOSS771_BASF_FINISH `date`


# Mark the file as bad if BASF returned something other than 0
if [ $BASFRET -ne 0 ]; then
  for p in /pic/projects/belle/voss771/ptSpectOut/*.h5; do
    mv $p $p.badret
  done
fi

