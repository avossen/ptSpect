for ex in 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
for stream in 00 01 02 03 04 05
do
for res in on_resonance 
do
for fl in mixed charged
do
for pt in kT qT
do
myDir=subMCTruePID_ex$ex\_on_resonance_$fl\_$stream\_$pt
for d in `ls mDstOut/$myDir/Normal_*.root`
do
newFName=`basename "$d" .root`
newFName=$newFName\TruePID_$fl\_$stream\_$pt.root
echo looking at $d new name $newFName
# $d has the full path, but newFName not
cp $d mDstOut/all_$fl/$newFName
done
done
done
done
done
done
