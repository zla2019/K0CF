#include <iostream>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include "MyTree.h"
#include "Config.h"
#include <unordered_map>
#include <TCanvas.h>
#include <TLatex.h>
#include <fstream>
#include <ctime>
#include <TRandom.h>
#include "LLFunc.h"
#include "Hist.h"

//CRAB header file
#define double_complex Complex
#include "source_files/volya_complex.h"
#define REDUCED_MOM
long IDUM=-1234;
#define NBMAX 5
#define NPHASEMAX 1000000
#define ABSOLUTE_MAXMOM 100

#include "source_files/crab.h"
//#include "interactions/crab_interaction_kpluskplus.cpp"
#include "interactions/crab_interaction_k0k0.cpp"
//#include "interactions/crab_interaction_pp.cpp"
#include "binnings/crab_bindefs_qinv.cpp"
#include "source_files/crab_main.cpp"
#include "source_files/crab_prinput.cpp"
#ifdef STRONG_INTERACTION
#include "source_files/crab_partwaveinit.cpp"
#endif

#ifdef COULOMB
#include "source_files/crab_coulomb.cpp"
#endif
#include "filters/crab_filter_acceptall.cpp"
#include "source_files/crab_corrcalc.cpp"
#include "source_files/crab_misc.cpp"
#include "source_files/crab_random.cpp"


const int nPBins = 23;
bool passCut(double a, double lower, double upper);
int getMomBin(float mom, const float* MomLower, const float* MomUpper, const float NMom);
bool loadBadRun(std::unordered_map<int, bool>& badRunList, std::ifstream& ifBadRunList);
bool passCut(double a, Config& config, std::string cutName);
bool passCut(double a, double b, Config& config, std::string cutName);
int getCent(int refMult);
inline bool passAllCuts(MyTree::Particle& p, Config& config);
MyTree::Particle boostParticle(TLorentzVector& p_cm, MyTree::Particle& p);

										//centraility defined by -1.0 < eta < 1.0
const int centFull_7p7[9] = { 11, 22, 39, 65, 102, 152, 220, 314, 379 }; 	// 7.7 GeV
const int centFull_3p9[9] = { 7, 13, 24, 39, 60, 89, 128, 181, 217 }; 		// 3.9 GeV
const int centFull_3p5[9] = { 6, 12, 21, 35, 54, 79, 113, 159, 190 }; 		// 3.5 GeV
const int centFull_3p2[9] = { 6, 11, 20, 32, 49, 71, 101, 142, 169 }; 		// 3.2 GeV
const int centFull_3p0[9] = { 5, 10, 18, 29, 45, 65, 93, 130, 155 }; 		// 3.0 GeV
const int* centFull = nullptr;

int main(int argc, char **argv)
{
	time_t first, second;
	first = time(NULL);
	//initial{{{
	std::string ifListName, ifCutListName, ofName;
	if(argc < 4) {
		return 1;
	}
	ifListName = argv[1];
	ifCutListName = argv[2];
	ofName = argv[3];
	std::cout << "LOG: Reading file list: " << ifListName << std::endl;
	std::cout << "LOG: Reading Cut File: " << ifCutListName << std::endl;
	std::ifstream ifFilelist(ifListName);
	std::ifstream ifCutList(ifCutListName);
	Config config(&ifCutList);
	const bool MuteWarning = true;	//used for debug

	float beamRapidity;
	if(config.mSetList["Energy"] == "3.0") {
		beamRapidity = 1.045;
		centFull = centFull_3p0;
	} else if(config.mSetList["Energy"] == "3.2") {
		beamRapidity = 1.135;
		centFull = centFull_3p2;
	} else if(config.mSetList["Energy"] == "3.5") {
		beamRapidity = 1.24;
		centFull = centFull_3p5;
	} else if(config.mSetList["Energy"] == "3.9") {
		centFull = centFull_3p9;
		beamRapidity = 1.36;
	}

	bininit();
#ifdef STRONG_INTERACTION
	partwaveinit();
#endif
#ifdef COULOMB
	coulset();
#endif
	//}}}

	//prepare plots{{{
	Hist hist;
	hist.init();

	//SI weight function
	gRandom->SetSeed(2);
	float rMean = 0, rErr = 0;
	float lambdaMean = 0, lambdaErr = 0;

	if(config.mSetList["Energy"] == "3.0")  {
		rMean = 2.05, rErr = 0.61;
		lambdaMean = 0.41, lambdaErr = 0.26;
	} else if(config.mSetList["Energy"] == "3.2") {
		rMean = 2.15, rErr = 0.49;
		lambdaMean = 0.44, lambdaErr = 0.18;
	} else if(config.mSetList["Energy"] == "3.5") {
		rMean = 2.65, rErr = 0.60;
		lambdaMean = 0.64, lambdaErr = 0.22;
	} else if(config.mSetList["Energy"] == "3.9") {
		rMean = 2.36, rErr = 0.45;
		lambdaMean = 0.40, lambdaErr = 0.16;
	}
	//}}}

	//processing{{{
	std::string ifName;
	while(getline(ifFilelist, ifName)) {
		TFile* ifTree = TFile::Open(ifName.c_str());
		if(!ifTree) {
			std::cout << ifName.c_str() << " can not open" << std::endl;
			continue;
		}
		MyTree *myTree = new MyTree((TTree*)ifTree->Get("urqmd"));
		if(!myTree->mTree) {
			std::cout << "ERROR: Open tree fail, skipping to next file" << std::endl;
			continue;
		}
		Long64_t nEvent = myTree->getNEvent();

		std::cout << "LOG: Reading file: " << ifName << std::endl;
		std::cout << "LOG: nEvent: " << nEvent << std::endl;
		for(int ievt = 0; ievt < nEvent; ++ievt) {
			if(ievt % 10000 == 0) {
				std::cout << "LOG: ievent: " << ievt << std::endl;
			}
			myTree->getEntry(ievt);
			//pre select
                        std::vector<MyTree::Particle> vUrqmdK;
                        std::vector<MyTree::Particle> vLabK;
			int refMult = 0;
                        for(int itrk = 0; itrk < myTree->mMul; ++itrk) {
                                if(myTree->mFrt[itrk] == 0) continue;
                                if(fabs(myTree->mFrt[itrk]) > 50) continue;
                                float p = sqrt(myTree->mPx[itrk]*myTree->mPx[itrk] + myTree->mPy[itrk]*myTree->mPy[itrk] + myTree->mPz[itrk]*myTree->mPz[itrk]);
                                float eta = 0.5 * log((p + myTree->mPz[itrk]) / (p - myTree->mPz[itrk]));
                                MyTree::Particle urqmdK = myTree->getParticle(itrk); 
                                TLorentzVector p_cm(0, 0, 1.1686, 1.498);
                                MyTree::Particle labCurK = boostParticle(p_cm, urqmdK);
				//if(fabs(urqmdK.eta) > 1) continue;
                                
                                if(fabs(myTree->mPid[itrk]) == 321 || fabs(myTree->mPid[itrk]) == 2212 || fabs(myTree->mPid[itrk]) == 211) ++refMult;
                                
                                if(fabs(myTree->mPid[itrk] != 311)) continue;
                                labCurK.y = -(labCurK.y + 1.045);
				if(labCurK.y < -1.0 || labCurK.y > 0.0) continue;
				if(labCurK.pt < 0.2 || labCurK.pt > 1.8) continue;
				if(labCurK.eta < -2.0 || labCurK.eta > 0) continue;

                                vLabK.push_back(std::move(labCurK));
                                vUrqmdK.push_back(std::move(urqmdK));
                        }

			if(refMult == 0) continue;
			myTree->mBufferCent9 = getCent(refMult);
			int cent9 = (int)myTree->mBufferCent9;
			hist.hCent9->Fill(cent9);
			hist.hRefMult->Fill(refMult);
			if(cent9 < 0 || cent9 > 8) continue;
			for(int iK = 0; iK < vLabK.size(); ++iK) {
				MyTree::Particle urqmdK = vUrqmdK[iK];
				MyTree::Particle labCurK = vLabK[iK];
				hist.Fill(labCurK);
				for(int iK2 = iK + 1; iK2 < vLabK.size(); ++iK2) {
					MyTree::Particle urqmdK2 = vUrqmdK[iK2];
					MyTree::Particle labCurK2 = vLabK[iK2];

					//declear lorentz vector for two particle
                                        TLorentzVector k1_v4, k2_v4;
                                        k1_v4.SetXYZT(labCurK.px, labCurK.py, labCurK.pz, labCurK.energy);
                                        k2_v4.SetXYZT(labCurK2.px, labCurK2.py, labCurK2.pz, labCurK2.energy);
                                        TLorentzVector kDiff_v4 = k1_v4 - k2_v4;
					float qinv = fabs(kDiff_v4.Mag());

					//CRAB correction
					float rr = labCurK2.frt - labCurK.frt;
					float pp = labCurK2.energy + labCurK.energy;
					float pdotr = rr * pp;
					float kdotr = (labCurK2.energy - labCurK.energy) * rr;
					float ptot2 = pp*pp;
					float r = -rr*rr;

					float rx = labCurK2.frx - labCurK.frx;
					float px = labCurK2.px + labCurK.px;
					float ry = labCurK2.fry - labCurK.fry;
					float py = labCurK2.py + labCurK.py;
					float rz = labCurK2.frz - labCurK.frz;
					float pz = labCurK2.pz + labCurK.pz;

					pdotr = pdotr - px * rx;
					pdotr = pdotr - py * ry;
					pdotr = pdotr - pz * rz;
					kdotr = kdotr - (labCurK2.px - labCurK.px) * rx;
					kdotr = kdotr - (labCurK2.py - labCurK.py) * ry;
					kdotr = kdotr - (labCurK2.pz - labCurK.pz) * rz;
					ptot2 = ptot2 - px*px;
					ptot2 = ptot2 - py*py;
					ptot2 = ptot2 - pz*pz;
					r = r + rx*rx;
					r = r + ry*ry;
					r = r + rz*rz;
					float qdotr = kdotr;
					r = sqrt(r + pdotr*pdotr / ptot2);
					double corr = corrcalc(fabs(kDiff_v4.Mag()) * 500., qdotr * 500., r);

					TF1 fLLWeight = TF1("fLLTmp", CFLLSI, 0, 0.4, 2);
					TF1 fGausR = TF1("fGausR", "gaus", rMean - rErr, rMean + rErr);
					TF1 fGausLambda = TF1("fGausLambda", "gaus", lambdaMean - lambdaErr, lambdaMean + lambdaErr);
					fGausR.SetParameters(1, rMean, rErr);
					fGausLambda.SetParameters(1, lambdaMean, lambdaErr);

					float lambdaTmp = fGausLambda.GetRandom();
					float rTmp = fGausR.GetRandom();
					fLLWeight.SetParameter(0, lambdaTmp);
					fLLWeight.SetParameter(1, rTmp);

					double llWeight = fabs(kDiff_v4.Mag()) > 0.4 ? 0 : fLLWeight.Eval(fabs(kDiff_v4.Mag()));

					hist.hQinvCorr->Fill(fabs(kDiff_v4.Mag()), corr);
					hist.hQdotrCorr->Fill(qdotr, corr);
					hist.hRCorr->Fill(r, corr);

					hist.hLLWeight->Fill(fabs(kDiff_v4.Mag()), llWeight);
					hist.hTotWeight->Fill(fabs(kDiff_v4.Mag()), llWeight + corr);

					// fill histo
					hist.hSameKqinvSI[cent9]->Fill(fabs(kDiff_v4.Mag()), 1 + llWeight);
					hist.hSameKqinvQS[cent9]->Fill(fabs(kDiff_v4.Mag()), corr);
					hist.hSameKqinv[cent9]->Fill(fabs(kDiff_v4.Mag()), corr + llWeight);
					hist.hSameKqinvWoCrab[cent9]->Fill(fabs(kDiff_v4.Mag()));
					//hist.FillSame(fabs(kDiff_v4.Mag()), cent9);
				}
				for(int imixevt = 0; imixevt < myTree->mMaxMixEvent[cent9] + 1; ++imixevt) {
                                        unsigned int nK = myTree->mMixBuffer[(int)cent9][imixevt].mBufferNTrack;
                                        for(int imixK = 0; imixK < nK; ++imixK) {
                                                MyTree::Particle labMixK = myTree->mMixBuffer[(int)cent9][imixevt].labParticle[imixK];
                                                MyTree::Particle urqmdMixK = myTree->mMixBuffer[(int)cent9][imixevt].urqmdParticle[imixK];

                                                TLorentzVector k1_v4, k2_v4;
                                                k1_v4.SetXYZT(labCurK.px, labCurK.py, labCurK.pz, labCurK.energy);
                                                k2_v4.SetXYZT(labMixK.px, labMixK.py, labMixK.pz, labMixK.energy);
                                                TLorentzVector kDiff_v4 = k1_v4 - k2_v4;
                                                hist.FillMix(fabs(kDiff_v4.Mag()), cent9);
                                        }     
				}
			}
			myTree->copyToBuffer(vLabK, vUrqmdK);
		}
		delete myTree;
		ifTree->Close();
		delete ifTree;
		myTree = nullptr;
		ifTree = nullptr;
	}
	//}}}

	//prepare output{{{
	TFile* ofPlots = new TFile(ofName.c_str(), "RECREATE");
	hist.Write(ofPlots);
	ofPlots->Close();
	delete ofPlots;
	//}}}
	second = time(NULL);
	std::cout << "processing time: " << difftime(second, first) << std::endl;
	return 0;
}

//my tools{{{
bool loadBadRun(std::unordered_map<int, bool>& badRunList, std::ifstream& ifBadRunList)
{
	if(!ifBadRunList.is_open()) {
		std::cout << "ERROR: Bad run list file not open" << std::endl;
	}
	while(!ifBadRunList.eof()) {
		std::string badRunStr = "";
		ifBadRunList >> badRunStr;
		if(badRunStr == "") {
			break;
		}
		int badRun = std::stoi(badRunStr);
		badRunList[badRun] = true;
	}
	std::cout << "LOG: " << badRunList.size() << " Bad Runs Loaded" << std::endl;
	return true;
}

bool passCut(double a, double lower, double upper)
{
	if(a < lower) {
		return false;
	} else if(a > upper) {
		return false;
	}
	return true;
}

int getMomBin(float mom, const float* MomLower, const float* MomUpper, const float NMom)
{
	int momBin = -1;
	for(int imom = 0; imom < NMom; ++imom) {
		if(mom > MomLower[imom] && mom < MomUpper[imom]) {
			momBin = imom;
			break;
		}
	}
	return momBin;
}

bool passCut(double a, Config& config, std::string cutName)
{
	if(a < config.mCutList[cutName].first) {
		return false;
	} else if(a > config.mCutList[cutName].second) {
		return false;
	}
	return true;
}
bool passCut(double a, double b, Config& config, std::string cutName)
{
	bool hasTHS = false;
	float ths = 0;
	bool pattern = 0;
	if(config.mThsList.find(cutName) == config.mThsList.end()) {
		hasTHS = false;
	} else {
		hasTHS = true;
		ths = config.mThsList[cutName].first;
		pattern = config.mThsList[cutName].second;
	}
	if(hasTHS) {
		bool passTHS = false;
		if(!pattern) {
			passTHS = passCut(b, -1e7, ths);
			return !(passTHS && !passCut(a, config, cutName));
		} else {
			passTHS = passCut(b, ths, 1e7);
			return !(passTHS && !passCut(a, config, cutName));
		}
	} else {
		return passCut(a, config, cutName);
	}

	return false;
}

inline bool passAllCuts(MyTree::Particle& p, Config& config)
{
	return true;
}

MyTree::Particle boostParticle(TLorentzVector& p_cm, MyTree::Particle& p)
{
        MyTree::Particle p2 = p;
        float mass = 0;
        if(fabs(p.pid) == 321) mass = 0.49368;
        else if(fabs(p.pid) == 2212) mass = 0.93827;
        else if(fabs(p.pid) == 211) mass = 0.13957;
        else if(fabs(p.pid) == 311) mass = 0.49761;
        TVector3 beta = -1 * p_cm.BoostVector();
        TLorentzVector p_v4(p2.px, p2.py, p2.pz, p2.energy);
        TLorentzVector r_v4(p2.frx, p2.fry, p2.frz, p2.frt);
        p_v4.Boost(beta);
        r_v4.Boost(beta);
        p2.px = p_v4.Px(), p2.py = p_v4.Py(), p2.pz = p_v4.Pz(), p2.energy = p_v4.E();
        p2.frx = r_v4.X(), p2.fry = r_v4.Y(), p2.frz = r_v4.Z(), p2.frt = r_v4.E();
        p2.pt = p_v4.Pt();
        p2.p = p_v4.P();
        p2.eta = p_v4.Eta();
        p2.energy = p_v4.E();
        p2.y = p_v4.Rapidity();
        p2.phi = p_v4.Phi();
        //p2.theta = p_v4.Theta();
	return p2;
}

int getCent(int refMult)
{
        int centrality;

        if      (refMult>=centFull[8]) centrality=8;
        else if (refMult>=centFull[7]) centrality=7;
        else if (refMult>=centFull[6]) centrality=6;
        else if (refMult>=centFull[5]) centrality=5;
        else if (refMult>=centFull[4]) centrality=4;
        else if (refMult>=centFull[3]) centrality=3;
        else if (refMult>=centFull[2]) centrality=2;
        else if (refMult>=centFull[1]) centrality=1;
        else if (refMult>=centFull[0]) centrality=0;
        else centrality = -1;

        return centrality;
}

//}}}
