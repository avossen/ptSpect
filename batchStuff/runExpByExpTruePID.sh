#!/bin/bash

#for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
#for ex in 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
#no 55 below since we already ran that
for ex in 31 33 35 37 39 41 43 45 47 49 51 53 61 63 65 67 69 71 73
#for ex in 55
#for ex in 61 63 65 67 69 71 73
#for ex in 33 35 37 39 41 43 45 47 49 51 53 
do
#for res in on_resonance continuum
for res in on_resonance 
do
#truePID only needed for mixed charged, since that will be subtracted
for spec in charged mixed
do
for stream in 00 01 02 03 04 05 06 07 08 09 
do
#    for d in `ls subData_ex$ex\_$res/*.sh`;do echo $d; done

echo looking at subMC_ex$ex\_$res\_$spec\_$stream/ 
for d in `ls subMC_ex$ex\_$res\_$spec\_$stream/*.sh`
do
 bsub -q s $PWD/$d
done
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
    mkdir mDstOut/subMCTruePID_ex$ex\_$res\_$spec\_$stream\_$pt/
    ./TwoHadAsymsCMod  ../groupData/ptSpect/subMC_ex$ex\_$res\_$spec\_$stream/ mDstOut/subMCTruePID_ex$ex\_$res\_$spec\_$stream\_$pt/ $pt truePID
done
#delete large output
rm ../groupData/ptSpect/subMC_ex$ex\_$res\_$spec\_$stream/*.root
rm ~/.lsf/*.out
cd /home/belle/vossen/myProjects/ptSpect/batchStuff
done
done
done
done
