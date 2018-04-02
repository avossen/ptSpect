#!/bin/bash                                                                                                                                    

#echo $d " counter is " $counter;                                                                                                              
counter=0;
subCounter=0;
dateString=`date +%d%b%Y`

mkdir /group/belle/users/vossen/ptSpect/AsymExtraction/PlotCompOutMC

for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
for res in on_resonance continuum
do
for spec in uds charm
do


myDir=subMC_ex$ex\_$res\_$spec
#myOutDir=/group/belle/users/vossen/ptSpect/AsymExtraction/PlotCompOutMC/$myDir
myOutDir=/home/belle/vossen/myProjects/ptSpect/AsymExtraction/PlotCompOutMC/
dataDir=/group/belle/users/vossen/ptSpect/$myDir


mkdir $myOutDir

#echo " dir: $myDir " ;
#echo " out dir: $myOutDir " ;

targetShFile=job_MCEx$ex\_$res\_$spec.sh
#cp batchHead.sh $targetShFile
cp batchHead1.sh $targetShFile

cat batchHead2.sh >> $targetShFile 
echo "/home/belle/vossen/myProjects/ptSpect/AsymExtraction/TwoHadAsymsCMod $dataDir $myOutDir mc" >>$targetShFile
cat batchEnd.sh >> $targetShFile
#fi
#fi
done 
done
done







