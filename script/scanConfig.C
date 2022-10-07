#include <iostream>
#include "Kshort0PurityScan.C"

void scanConfig()
{
	//inital {{{
	std::string prefix;
	std::string energy;
	std::string debugMode;
	std::string nSigma = "2.0";
	int loopStart = 0, loopEnd = 0;
	//}}}

	//config {{{
	std::cout << "Energy(eg: \"3p0\" means 3.0 GeV): ";
	std::cin >> energy;

	std::cout << "prefix: ";
	std::cin >> prefix;

	std::cout << "loopStart index: ";
	std::cin >> loopStart;

	std::cout << "loopEnd index: ";
	std::cin >> loopEnd;

	std::cout << "debugging mode?(Y/N): ";
	std::cin >> debugMode;

	//}}}

	//processing {{{
	for(int icfg = loopStart; icfg <= loopEnd; ++icfg) {
		std::ifstream config(Form("%scutSet%s%d/cutSet%s%dInvM.txt", energy.c_str(), prefix.c_str(), icfg, prefix.c_str(), icfg));
		while(getline(config, nSigma)) {
			std::string::size_type pos = nSigma.find("NSigmaMass");
			if(pos == std::string::npos) {
				continue;
			} else {
				nSigma = nSigma.substr(pos + 11);
				std::cout << "nSigma: " << nSigma << std::endl;
			}
		}

		if(debugMode == "N") {
			gSystem->Exec(Form("./../K0CF/analysis/generateHists %sTree.list %scutSet%s%d/cutSet%s%dInvM.txt %scutSet%s%d/%scutSet%s%dInvM.root", 
						energy.c_str(), energy.c_str(), prefix.c_str(), icfg, prefix.c_str(), icfg, energy.c_str(), prefix.c_str(), icfg, energy.c_str(), prefix.c_str(), icfg));

			Kshort0PurityScan(Form("%scutSet%s%d/%scutSet%s%dInvM.root", energy.c_str(), prefix.c_str(), icfg, energy.c_str(), prefix.c_str(), icfg), Form("%scutSet%s%d", energy.c_str(), prefix.c_str(), icfg), nSigma);

			gSystem->Exec(Form("./../K0CF/analysis/generateHists %sTree.list %scutSet%s%d/cutSet%s%d.txt %scutSet%s%d/%scutSet%s%d.root", 
						energy.c_str(), energy.c_str(), prefix.c_str(), icfg, prefix.c_str(), icfg, energy.c_str(), prefix.c_str(), icfg, energy.c_str(), prefix.c_str(), icfg));
		} else {
			std::cout << Form("./../K0CF/analysis/generateHists %sTree.list %scutSet%s%d/cutSet%s%dInvM.txt %scutSet%s%d/%scutSet%s%dInvM.root\n", 
						energy.c_str(), energy.c_str(), prefix.c_str(), icfg, prefix.c_str(), icfg, energy.c_str(), prefix.c_str(), icfg, energy.c_str(), prefix.c_str(), icfg);
			std::cout << "Kshort0PurityScan(" << Form("\"%scutSet%s%d/%scutSet%s%dInvM.root\"", energy.c_str(), prefix.c_str(), icfg, energy.c_str(), prefix.c_str(), icfg) << ", " << Form("\"%scutSet%s%d\"", energy.c_str(), prefix.c_str(), icfg) << ", nSigma)\n";
			std::cout << Form("./../K0CF/analysis/generateHists %sTree.list %scutSet%s%d/cutSet%s%d.txt %scutSet%s%d/%scutSet%s%d.root\n", 
						energy.c_str(), energy.c_str(), prefix.c_str(), icfg, prefix.c_str(), icfg, energy.c_str(), prefix.c_str(), icfg, energy.c_str(), prefix.c_str(), icfg);
		}
	}
	//}}}
}
