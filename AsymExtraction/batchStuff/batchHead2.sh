
# `.brofile' is my equivalent to `.bashrc'
echo Zuhause: $HOME
#source $HOME/.bashrc

# This sets up your environment to run code in the Belle paradigm
#source /pic/projects/belle/scripts/general_env.sh

# Add my modifications to the standard Belle environment
#source $HOME/custom_env.sh

source $HOME/.bashrc
source /sw/belle/local/etc/bashrc_general


# This is necessary to get code that fills HBOOK ntuples running comfortably
# on Olympus
shopt -s expand_aliases

echo ANSELM_BASF_START `date`
export USE_GRAND_REPROCESS_DATA=1
export BELLE_MESSAGE_LEVEL=INFO
export LD_LIBRARY_PATH=/home/belle/vossen/myProjects/ptSpect/AsymExtraction/:`root-config --libdir`:$LD_LIBRARY_PATH

#setenv BASF_NPROCESS 0

echo    $LD_LIBRARY_PATH

