#!/bin/sh
export USE_GRAND_REPROCESS_DATA=1
export BELLE_MESSAGE_LEVEL=INFO
export LD_LIBRARY_PATH=.:`root-config --libdir`:$LD_LIBRARY_PATH
export BASF_Module_Dir=.:$basf_Module_DIR
#setenv BASF_NPROCESS 0

echo    $LD_LIBRARY_PATH
echo    $BASF_MODULE_DIR
echo    $PANTHER_TABLE_DIR


MODULE=ptSpect
HBKFILE=/dev/null

basf <<EOF
path create main
path create analysis
path add_module main fix_mdst
path add_condition main >:0:analysis
path add_condition main <=:0:KILL
path add_module analysis ${MODULE}
module put_parameter ptSpect rfname\/home/belle/vossen/myProjects/ptSpect/mcEx55.root
#module put_parameter ptSpect onlyGen\7
module put_parameter fix_mdst Make_pi0_option\2
module put_parameter fix_mdst Make_pi0_lower_limit\-20.0
module put_parameter fix_mdst Make_pi0_upper_limit\25.0
#module put_parameter fix_mdst Correct_ecl_option\1
initialize
nprocess set 0
histogram define ${HBKFILE}

process_event /group/belle/nagoya/cp201910/chai/sgt2/bfactory/belle_data/caseB/mcprod/b20090127_0910.4/eecc.e55.3/continuum/08/eecc-e00005\
5r000805-b20090127_0910-3.mdst 0
#process_event /group/belle/nagoya/cp201910/chai/sgt2/bfactory/belle_data/caseB/mcprod/b20090127_0910.2/tautau.e55.0/continuum/07/tautau-e000055r000793-b20090127_0910-0.mdst 0
#process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/04/all/0127/continuum/15/evtgen-uds-04-all-e000055r001585-b20090127_0910.mdst 0

# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001602-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001603-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001604-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001605-b20090127_0910.mdst 0 
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001606-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001607-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001609-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001610-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001611-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001619-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001620-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001621-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001626-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001627-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001631-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001632-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001633-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001634-b20090127_0910.mdst 0
# process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/05/all/0127/continuum/16/evtgen-uds-05-all-e000055r001635-b20090127_0910.mdst 0
#debug run
#process_event /group/belle/bdata_b/dstprod/dat/e000055/HadronBJ/0127/continuum/16/HadronBJ-e000055r001621-b20090127_0910.mdst 0

#job 25
#process_event /group/belle/bdata_b/dstprod/dat/e000055/HadronBJ/0127/continuum/16/HadronBJ-e000055r001660-b20090127_0910.mdst 1000000
#process_event /group/belle/bdata_b/dstprod/dat/e000055/HadronBJ/0127/continuum/16/HadronBJ-e000055r001621-b20090127_0910.mdst 0
#process_event /group/belle/bdata_b/dstprod/dat/e000055/HadronBJ/0127/continuum/08/HadronBJ-e000055r000821-b20090127_0910.mdst 0

#process_event /group/belle/bdata_b/dstprod/dat/e000055/HadronBJ/0127/continuum/16/HadronBJ-e000055r001602-b20090127_0910.mdst 100000


#process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/uds/00/all/0127/continuum/07/evtgen-uds-00-all-e000055r000793-b20090127_0910.mdst 0
#process_event /group/belle/bdata_b/dstprod/dat/e000055/HadronBJ/0127/continuum/08/HadronBJ-e000055r000820-b20090127_0910.mdst 0
#process_event /group/belle/bdata_b/dstprod/dat/e000055/HadronBJ/0127/continuum/16/HadronBJ-e000055r001621-b20090127_0910.mdst 1000

output close
terminate
EOF
