#!/bin/bash                                                                                                                                    

#echo $d " counter is " $counter;                                                                                                              
counter=0;

dateString=`date +%d%b%Y`

for tune in 13
#for tune in 00 01 10 11 12 13 14
do
for yn in yes no
do
for spec in uds charm 
do
subCounter=0;
myDir=genMC_Tune$tune\_ISR$yn\_$spec
mkdir $myDir
echo " dir: $myDir " ;
mkdir /group/belle/users/vossen/ptSpect/ISRStudies/$myDir
mkdir /group/belle/users/vossen/ptSpectOut/ISRStudies/$myDir

# so this gives the zipped files
for d in `cat lists/genTune$tune\_ISR$yn\_$spec.list`
do
let "counter+=1"
#unpacked=unpackedTune$tune\_ISR$yn\_$spec\_$counter.pgen
#unzip


#if [ "$counter" -le "$5" ]
#then
#if [ "$counter" -ge "$4" ]
#then
let "subCounter+=1"
targetShFile=$myDir/job_Tune$tune\_ISR$yn\_$spec\_$subCounter.sh

#cp batchHead.sh $targetShFile
cp batchHead1.sh $targetShFile
#echo "gunzip -c $d > /btmp/$unpacked" >> $targetShFile


echo "cp $d /btmp/" >> $targetShFile
#get file w/o path
fWOPath=$(basename $d)
fWOExt=$(basename "${d%.*}")
echo "gunzip /btmp/$fWOPath" >> $targetShFile


echo "#BSUB -o  /group/belle/users/vossen/ptSpectOut/ISRStudies/$myDir/jobId_$subCounter.out" >> $targetShFile
echo "#BSUB -e  /group/belle/users/vossen/ptSpectOut/ISRStudies/$myDir/jobId_$subCounter.err" >> $targetShFile
echo "#BSUB -J ISRStudy_$subCounter"  >> $targetShFile 
cat batchHead2.sh >> $targetShFile 
echo "module put_parameter ptSpect rfname\\/group/belle/users/vossen/ptSpect/ISRStudies/$myDir/job_$subCounter.root" >> $targetShFile
cat batchMiddle.sh >> $targetShFile
echo "process_event /btmp/$fWOExt 0" >> $targetShFile
cat batchEnd.sh >> $targetShFile
echo "rm /btmp/$fWOExt" >> $targetShFile
#fi
#fi
chmod a+x $targetShFile
done
done
done 
done








