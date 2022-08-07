#!/bin/bash
energy="3p2"
preidx="3."

for num in {0..5}
do
	#mkdir -p cutSet${num}
	#mv cutSet${num}InvM.txt cutSet${num}
	#mv cutSet${num}.txt cutSet${num}
	./generateHists ${energy}Tree.list ${energy}cutSet${preidx}${num}/cutSet${preidx}${num}InvM.txt ${energy}cutSet${preidx}${num}/${energy}cutSet${preidx}${num}InvM.root
	root -l -b -q Kshort0PurityScan.C\(\"${energy}cutSet${preidx}${num}/${energy}cutSet${preidx}${num}InvM.root\",\"${energy}cutSet${preidx}${num}\"\);
	./generateHists ${energy}Tree.list ${energy}cutSet${preidx}${num}/cutSet${preidx}${num}.txt ${energy}cutSet${preidx}${num}/${energy}cutSet${preidx}${num}.root
done
