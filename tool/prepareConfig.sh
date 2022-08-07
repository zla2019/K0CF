#!/bin/bash
energy="3p2"
preidx="2."
newpreidx="3."

for num in {0..5}
do
	cp -r ${energy}cutSet${preidx}${num} ${energy}cutSet${newpreidx}${num}
	rm ${energy}cutSet${newpreidx}${num}/*.png
	rm ${energy}cutSet${newpreidx}${num}/*.root
	mv ${energy}cutSet${newpreidx}${num}/cutSet${preidx}${num}.txt ${energy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}.txt
	mv ${energy}cutSet${newpreidx}${num}/cutSet${preidx}${num}InvM.txt ${energy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}InvM.txt

	sed -i 's/Cut: Chi2Topo -1e7 4/Cut: Chi2Topo -1e7 3/g' ${energy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}.txt
	sed -i 's/Cut: Chi2NDF -1e7 4/Cut: Chi2NDF -1e7 3/g' ${energy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}.txt
	sed -i "s|SET: PurityPath ${energy}cutSet${preidx}${num}|SET: PurityPath ${energy}cutSet${newpreidx}${num}|g" ${energy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}.txt

	sed -i 's/Cut: Chi2Topo -1e7 4/Cut: Chi2Topo -1e7 3/g' ${energy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}InvM.txt
	sed -i 's/Cut: Chi2NDF -1e7 4/Cut: Chi2NDF -1e7 3/g' ${energy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}InvM.txt
done
