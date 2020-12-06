#!/bin/bash                                                                                                                                    

#echo $d " counter is " $counter;                                                                                                              
counter=0;

dateString=`date +%d%b%Y`

#for tune in 00
#for tune in 12 13 14
for tune in 00 01 10 11 12 13 14
do
#for yn in yes 
for yn in yes no
do
#for spec in uds 
for spec in uds charm 
do
subCounter=0;
myDir=genMC_red_Tune$tune\_ISR$yn\_$spec
mkdir $myDir
myFailedDir=$myDir\_Failed
mkdir $myFailedDir
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
targetShFileNoDir=job_Tune$tune\_ISR$yn\_$spec\_$subCounter.sh
#cp batchHead.sh $targetShFile
cp batchHead1.sh $targetShFile
#echo "gunzip -c $d > /btmp/$unpacked" >> $targetShFile


echo "cp $d /tmp/" >> $targetShFile
#get file w/o path
fWOPath=$(basename $d)
fWOExt=$(basename "${d%.*}")
echo "gunzip /tmp/$fWOPath" >> $targetShFile


echo "#BSUB -o  /group/belle/users/vossen/ptSpectOut/ISRStudies/$myDir/jobId_$subCounter.out" >> $targetShFile
echo "#BSUB -e  /group/belle/users/vossen/ptSpectOut/ISRStudies/$myDir/jobId_$subCounter.err" >> $targetShFile
echo "#BSUB -J ISRStudy_$subCounter"  >> $targetShFile 
cat batchHead2.sh >> $targetShFile 
echo "module put_parameter ptSpect onlyGen\\7" >> $targetShFile
echo "module put_parameter ptSpect rfname\\/group/belle/users/vossen/ptSpect/ISRStudies/$myDir/job_$subCounter.root" >> $targetShFile
cat batchMiddle.sh >> $targetShFile
echo "process_event /tmp/$fWOExt 0" >> $targetShFile
#have the same now in this file
#cat batchEnd.sh >> $targetShFile
#fi
#fi

echo "output close">> $targetShFile

echo "terminate" >>$targetShFile

echo "EOF" >>$targetShFile
echo "rm /tmp/$fWOExt" >> $targetShFile
echo "rm /tmp/$fWOPath" >> $targetShFile

echo 'dateString=`date +%d%b%Y`'>> $targetShFile

# Grab the exit code of BASF
echo 'BASFRET=$?' >> $targetShFile
echo 'echo VOSSEN_BASF_FINISH `date`' >> $targetShFile


# Mark the file as bad if BASF returned something other than 0
echo 'if [ $BASFRET -ne 0 ]; then' >> $targetShFile
echo   "cp $PWD/$targetShFile $PWD/$myFailedDir/$targetShFile" >> $targetShFile
echo "fi"  >> $targetShFile

chmod a+x $targetShFile
done
done
done 
done





