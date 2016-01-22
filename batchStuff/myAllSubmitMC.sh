#!/bin/bash                                                                                                                                    

#echo $d " counter is " $counter;                                                                                                              
counter=0;
subCounter=0;
dateString=`date +%d%b%Y`



for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
  for res in  continuum on_resonance
do
for spec in uds charm 
#for spec in mixed charged 
do

myDir=subMC_ex$ex\_$res\_$spec
mkdir $myDir
echo " dir: $myDir " ;
mkdir /pic/projects/belle/voss771/ptSpect/$myDir
mkdir /pic/projects/belle/voss771/ptSpectOut/$myDir


#for d in `cat mc$1_onRes.list`
for d in `cat lists/mc$ex\_$res\_$spec.list`
do
let "counter+=1"
#if [ "$counter" -le "$5" ]
#then
#if [ "$counter" -ge "$4" ]
#then
let "subCounter+=1"
targetShFile=$myDir/job_$ex\_$res\_$spec\_$subCounter.sh
#cp batchHead.sh $targetShFile
cp batchHead1.sh $targetShFile
echo "#SBATCH -o /pic/projects/belle/voss771/ptSpectOut/$myDir/jobId_$subCounter.out" >> $targetShFile
echo "#SBATCH -e /pic/projects/belle/voss771/ptSpectOut/$myDir/jobId_$subCounter.err" >> $targetShFile
echo "#SBATCH -J handedness_$subCounter"  >> $targetShFile 
cat batchHead2.sh >> $targetShFile 
echo "module put_parameter ptSpect rfname\\/pic/projects/belle/voss771/ptSpect/$myDir/job_$subCounter.root" >> $targetShFile
cat batchMiddle.sh >> $targetShFile
echo "process_event $d 0" >> $targetShFile
cat batchEnd.sh >> $targetShFile
#fi
#fi

done
done
done 
done








