#!/bin/bash
energy="3p5"
preidx="Acc3_3.1."

for num in {15..18}
do
	#mkdir -p cutSet${num}
	#mv cutSet${num}InvM.txt cutSet${num}
	#mv cutSet${num}.txt cutSet${num}
	#./../K0CF/analysis/generateHists ${energy}Tree.list ${energy}cutSet${preidx}${num}/cutSet${preidx}${num}InvM.txt ${energy}cutSet${preidx}${num}/${energy}cutSet${preidx}${num}InvM2.root
	#NSigma=`cat ${energy}cutSet${preidx}${num}/cutSet${preidx}${num}InvM.txt | grep "NSigmaMass" | sed "s/.*NSigmaMass //g"`
	#root -l -b -q ../K0CF/script/Kshort0PurityScan.C\(\"${energy}cutSet${preidx}${num}/${energy}cutSet${preidx}${num}InvM2.root\",\"${energy}cutSet${preidx}${num}\",\"$NSigma\"\);
	./../K0CF/analysis/generateHists ${energy}Tree.list ${energy}cutSet${preidx}${num}/cutSet${preidx}${num}.txt ${energy}cutSet${preidx}${num}/${energy}cutSet${preidx}${num}.root
done
