#!/bin/bash
energy="3p0"
preidx="Test1."

for num in {0..0}
do
	#mkdir -p cutSet${num}
	#mv cutSet${num}InvM.txt cutSet${num}
	#mv cutSet${num}.txt cutSet${num}
	#./generateHists ${energy}Tree.list ${energy}cutSet${preidx}${num}/cutSet${preidx}${num}InvM.txt ${energy}cutSet${preidx}${num}/${energy}cutSet${preidx}${num}InvM.root
	NSigma=`cat 3p2cutSet3.1/cutSet3.1InvM.txt | grep "NSigmaMass" | sed "s/.*NSigmaMass //g"`
	root -l -b -q ../script/Kshort0PurityScan.C\(\"${energy}cutSet${preidx}${num}/${energy}cutSet${preidx}${num}InvM.root\",\"${energy}cutSet${preidx}${num}\",\"$NSigma\"\);
	#./generateHists ${energy}Tree.list ${energy}cutSet${preidx}${num}/cutSet${preidx}${num}.txt ${energy}cutSet${preidx}${num}/${energy}cutSet${preidx}${num}.root
done
