
# `.brofile' is my equivalent to `.bashrc'
echo Zuhause: $HOME
#source $HOME/.bashrc

# This sets up your environment to run code in the Belle paradigm
#source /pic/projects/belle/scripts/general_env.sh

# Add my modifications to the standard Belle environment
#source $HOME/custom_env.sh

module purge
source $HOME/.bashrc

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
