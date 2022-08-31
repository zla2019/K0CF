#!/bin/bash
energy="3p5"
newenergy="3p5"
preidx="Acc4_3.1."
newpreidx="Acc2_3.1."

for num in {0..14}
do
	cp -r ${energy}cutSet${preidx}${num} ${newenergy}cutSet${newpreidx}${num}
	rm ${newenergy}cutSet${newpreidx}${num}/*.png
	rm ${newenergy}cutSet${newpreidx}${num}/*.root
	mv ${newenergy}cutSet${newpreidx}${num}/cutSet${preidx}${num}.txt ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}.txt
	mv ${newenergy}cutSet${newpreidx}${num}/cutSet${preidx}${num}InvM.txt ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}InvM.txt

	sed -i 's/Rap -1.0 0.4/Rap -1.0 0.0/g' ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}.txt
	#cat ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}.txt
	#sed -i 's/Energy 3.2/Energy 3.5/g' ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}.txt
	#sed -i 's/Energy 3.2/Energy 3.5/g' ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}InvM.txt
	#sed -i 's/3p2/3p5/g' ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}.txt
	#sed -i 's/_1.1./_3.1./g' ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}.txt

	#sed -i 's/Cut: Chi2Topo 0 5/Cut: Chi2Topo 0 3/g' ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}.txt
	#sed -i 's/Cut: Chi2NDF 0 5/Cut: Chi2NDF 0 3/g' ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}.txt
	#sed -i 's/Cut: Chi2NDF -1e7 4/Cut: Chi2NDF -1e7 3/g' ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}.txt
	#sed -i "s|SET: PurityPath ${energy}cutSet${preidx}${num}|SET: PurityPath ${newenergy}cutSet${newpreidx}${num}|g" ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}.txt

	#sed -i 's/Cut: Chi2Topo 0 5/Cut: Chi2Topo 0 3/g' ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}InvM.txt
	#sed -i 's/Cut: Chi2NDF 0 5/Cut: Chi2NDF 0 3/g' ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}InvM.txt
	#sed -i 's/Cut: Chi2NDF -1e7 4/Cut: Chi2NDF -1e7 3/g' ${newenergy}cutSet${newpreidx}${num}/cutSet${newpreidx}${num}InvM.txt
done
