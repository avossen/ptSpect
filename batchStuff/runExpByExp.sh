#!/bin/bash

#for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
for ex in 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
#for res in on_resonance continuum
for res in on_resonance
do
#    for d in `ls subData_ex$ex\_$res/*.sh`;do echo $d; done
for d in `ls subData_ex$ex\_$res/*.sh`; do bsub -q s $PWD/$d; done
    num_jobs=`bjobs | wc -l`
    while [[ ${num_jobs} -ge 1 ]] ;  do 
	echo " waiting another 2 mins for $num_jobs jobs"
	sleep 120;
	num_jobs=`bjobs | wc -l`
    done
# all jobs are done, run analysis code
cd /home/belle/vossen/myProjects/ptSpect/AsymExtraction

for pt in kT qT 
do
    mkdir mDstOut/subData_ex$ex\_$res\_$pt/
    ./TwoHadAsymsCMod  ../groupData/ptSpect/subData_ex$ex\_$res/ mDstOut/subData_ex$ex\_$res\_$pt/ $pt
done
#delete large output
rm ../groupData/ptSpect/subData_ex$ex\_$res/*.root
rm ~/.lsf/*.out
cd /home/belle/vossen/myProjects/ptSpect/batchStuff
done
done
