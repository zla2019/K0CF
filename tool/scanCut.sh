#!/bin/bash
energy="3p2"
preidx="Test3."

for num in {1..1}
do
	#mkdir -p cutSet${num}
	#mv cutSet${num}InvM.txt cutSet${num}
	#mv cutSet${num}.txt cutSet${num}
	./../K0CF/analysis/generateHists ${energy}Tree.list ${energy}cutSet${preidx}${num}/cutSet${preidx}${num}InvM.txt ${energy}cutSet${preidx}${num}/${energy}cutSet${preidx}${num}InvM.root
	NSigma=`cat ${energy}cutSet${preidx}${num}/cutSet${preidx}${num}InvM.txt | grep "NSigmaMass" | sed "s/.*NSigmaMass //g"`
	root -l -b -q ../K0CF/script/Kshort0PurityScan.C\(\"${energy}cutSet${preidx}${num}/${energy}cutSet${preidx}${num}InvM.root\",\"${energy}cutSet${preidx}${num}\",\"$NSigma\"\);
	./../K0CF/analysis/generateHists ${energy}Tree.list ${energy}cutSet${preidx}${num}/cutSet${preidx}${num}.txt ${energy}cutSet${preidx}${num}/${energy}cutSet${preidx}${num}.root
done
