#for ex in 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
for ex in 55
do
for res in on_resonance continuum
do
for fl in uds charm
do
for pt in kT qT
do
for num in 00 01 02 03 04 05
do
myDir=subMCWoFixMdst_ex$ex\_$res\_$fl\_$num\_$pt
for d in `ls mDstOut/$myDir/NormalWoA_*.root`
do
newFName=`basename "$d" .root`
newFName=$newFName\_$fl\_$num.root
echo looking at $d new name $newFName
# $d has the full path, but newFName not
cp $d mDstOut/all_MCWoA_$res\_$pt/$newFName
done
for d in `ls mDstOut/$myDir/smearingPlotterRaw_*.root`
do
newFName=`basename "$d" .root`
newFName=$newFName\_$fl\_$num.root
echo looking at $d new name $newFName
# $d has the full path, but newFName not
cp $d mDstOut/all_MCWoA_$res\_$pt/$newFName
done
done
done
done
done
done
