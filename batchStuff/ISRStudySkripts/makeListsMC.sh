#!/bin/bash


for tune in 00 01 10 11 12 13 14
do
for yn in yes no
do
for spec in uds charm 
do
#find /pic/projects/Belle/ZMDST/MC/MC_4S/MC_4S_EXP$ex/MC_4S_EXP$ex/$res/$spec/ > lists/mc$ex\_$res\_$spec.list
find  /ghi/fs01/belle/bdata2/users/rseidl/mcgen/tune$yn\isr_weak$tune\_$spec/ -wholename "*.pgen.gz" > lists/genTune$tune\_ISR$yn\_$spec.list
done
done
done

#/pic/projects/Belle/ZMDST/MC/MC_4S/MC_4S_EXP55/MC_4S_EXP55/continuum/charm/