#!/bin/bash                                                                                                                                    
#need to split up in batches 

#echo $d " counter is " $counter;                                                                          
HOME=/home/belle/vossen/myProjects/ptSpect/AsymExtraction


#for tune in 00 01 10 11 12 13 14
for tune in 13 
do
#for yn in yes 
for yn in yes no
do
#for spec in uds 
for spec in uds charm 
do
myDir=genMC_red_Tune$tune\_ISR$yn\_$spec
myDirOut=ISROut_Tune$tune\_ISR$yn\_$spec
mkdir $HOME/$myDirOut
for batch in 1 2 3
do


myDirBatch=genMC_red_Tune$tune\_ISR$yn\_$spec\_$batch
myDirBatchOut=ISROut_Tune$tune\_ISR$yn\_$spec\_$batch



echo " dir: $myDirBatch " ;
mkdir $HOME/$myDirBatchOut
echo running $HOME/TwoHadAsymsCMod /group/belle/users/vossen/ptSpect/ISRStudies/$myDirBatch/ $HOME/$myDirBatchOut/ mc
$HOME/TwoHadAsymsCMod /group/belle/users/vossen/ptSpect/ISRStudies/$myDirBatch/ $HOME/$myDirBatchOut/ mc
cp $HOME/$myDirBatchOut/NormalWoA_*.root $HOME/$myDirOut/NormalWoA_$batch.root
done
done
done 
done








