#include <iostream>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include "MyTree.h"
#include "Config.h"
#include <unordered_map>
#include <TCanvas.h>
#include <TLatex.h>
#include <fstream>
#include <ctime>
#include "Hist.h"

const int nPBins = 23;
//PID calibration
const float p_low[nPBins]  = {0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8};
const float p_high[nPBins] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 10.0};
float pidCalib_pion_3p9[nPBins]  = {1.33, 1.50, 1.54, 1.55, 1.55, 1.56, 1.57, 1.58, 1.60, 1.62, 1.63, 1.64, 1.66, 1.67, 1.69, 1.70, 1.72, 1.73, 1.75, 1.78, 1.80, 1.82, 1.84};
float pidCalib_pion_3p5[nPBins]  = {1.30, 1.49, 1.53, 1.54, 1.54, 1.55, 1.56, 1.57, 1.59, 1.61, 1.62, 1.63, 1.65, 1.66, 1.68, 1.69, 1.71, 1.73, 1.75, 1.77, 1.79, 1.82, 1.83};
float pidCalib_pion_3p2[nPBins]  = {1.39, 1.57, 1.61, 1.62, 1.64, 1.65, 1.66, 1.68, 1.69, 1.71, 1.73, 1.74, 1.75, 1.77, 1.79, 1.80, 1.82, 1.84, 1.86, 1.88, 1.91, 1.94, 1.94};
float pidCalib_pion_3p0[nPBins]  = {0.354, 0.354, 0.388, 0.396, 0.404, 0.414, 0.425, 0.437, 0.450, 0.461, 0.472, 0.482, 0.493, 0.499, 0.509, 0.516, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526};
float* pidCalib_pion = nullptr;

bool passCut(double a, double lower, double upper);
int getMomBin(float mom, const float* MomLower, const float* MomUpper, const float NMom);
bool loadBadRun(std::unordered_map<int, bool>& badRunList, std::ifstream& ifBadRunList);
bool passCut(double a, Config& config, std::string cutName);
bool passCut(double a, double b, Config& config, std::string cutName);
inline bool passAllCuts(MyTree::Particle& p, Config& config);
void loopVect(std::vector<MyTree::Particle>& v1, std::vector<MyTree::Particle>& v2, Hist& hist, int cent9, std::string opt, TH2F* hPurity[] = { 0 });

int main(int argc, char **argv)
{
	time_t first, second;
	first = time(NULL);
	//initial{{{
	std::string ifListName, ifRotListName, ifCutListName, ofName;
	if(argc < 4) {
		return 1;
	}
	ifListName = argv[1];
	ifCutListName = argv[2];
	ofName = argv[3];
	std::cout << "LOG: Reading file list: " << ifListName << std::endl;
	std::cout << "LOG: Reading rotation file list: " << ifRotListName << std::endl;
	std::cout << "LOG: Reading Cut File: " << ifCutListName << std::endl;
	std::ifstream ifFilelist(ifListName);
	std::ifstream ifRotFilelist(ifRotListName);
	std::ifstream ifCutList(ifCutListName);
	Config config(&ifCutList);
	std::unordered_map<int, bool> badRunList;
	const float Kmass = 0.497611;
	//mass window info
	//const float Mean[3][4] = { { 0.4979, 0.4981, 0.4979, 0.4949 }, { 0.4978, 0.4981, 0.4979, 0.4979 }, { 0.4979, 0.4981, 0.4980, 0.4980 } };
	//const float Sigma[3][4] = { { 0.0027, 0.0035, 0.0033, 0.0033 }, { 0.0027, 0.0034, 0.0036, 0.0036 }, { 0.0027, 0.0035, 0.0039, 0.0039 } };
	const float Mean[3][4] = { { 0.4981, 0.4981, 0.4981, 0.4981 }, { 0.4981, 0.4981, 0.4981, 0.4981 }, { 0.4981, 0.4981, 0.4981, 0.4981 } };
	const float Sigma[3][4] = { { 0.0035, 0.0035, 0.0035, 0.0035 }, { 0.0035, 0.0035, 0.0035, 0.0035 }, { 0.0035, 0.0035, 0.0035, 0.0035 } };
	const float NMassSigma = std::stof(config.mSetList["NSigmaMass"]);
	const bool MuteWarning = true;	//used for debug

	float beamRapidity;
	if(config.mSetList["Energy"] == "3.0") beamRapidity = 1.045;
	else if(config.mSetList["Energy"] == "3.2") beamRapidity = 1.135;
	else if(config.mSetList["Energy"] == "3.5") beamRapidity = 1.24;
	else if(config.mSetList["Energy"] == "3.9") beamRapidity = 1.36;

	if(config.mSetList["Energy"] == "3.0") pidCalib_pion = pidCalib_pion_3p0;
	else if(config.mSetList["Energy"] == "3.2") pidCalib_pion = pidCalib_pion_3p2;
	else if(config.mSetList["Energy"] == "3.5") pidCalib_pion = pidCalib_pion_3p5;
	else if(config.mSetList["Energy"] == "3.9") pidCalib_pion = pidCalib_pion_3p9;

	if(config.mSetList["Energy"] == "3.2") {
		std::ifstream ifBadRun("/home/zla/CF/K0CF/analysis/3p2badRun.list");
		loadBadRun(badRunList, ifBadRun);
	}
	//}}}

	//Save cut info{{{
	TCanvas* caCut = new TCanvas("cuts", "cuts", 1920, 1080);
	caCut->cd();
	TLatex* ltx = new TLatex(0.0, 0.9, "Cuts:");
	//}}}

	//prepare plots{{{
	Hist hist;
	hist.init();
	TFile* ifPurity;
	TH2F* hPurity[3];
	if(config.mSwitchList["OpenPairPurity"]) {
		ifPurity = TFile::Open(config.mSetList["PurityPath"].c_str());
		if(!ifPurity->IsOpen()) {
			std::cout << "ERROR: Purity file is not open" << std::endl;
		}
		for(int icent = 0; icent < 3; ++icent) {
			hPurity[icent] = (TH2F*)ifPurity->Get(Form("hKPurity_cent%i", icent));
		}
	}
	//}}}

	//processing{{{
	std::string ifName;
	std::string ifRotName;
	while(getline(ifFilelist, ifName)) {
		ifName.replace(0, ifName.find_last_of("/"), "/gpfs/mnt/gpfs01/star/pwg/lazhang/CF/correlation_test/out");
		ifRotName = ifName;
		ifRotName.replace(ifRotName.find_last_of("_"), 9, "_rot_ana.root");
		std::cout << "LOG: open signal file: " << ifName << std::endl;
		std::cout << "LOG: open rotsig file: " << ifRotName << std::endl;
		TFile* ifTree = TFile::Open(ifName.c_str());
		TFile* ifRotTree = TFile::Open(ifRotName.c_str());
		if(!ifTree || !ifRotTree) {
			std::cout << ifName.c_str() << " can not open" << std::endl;
			continue;
		}
		MyTree *myTree = new MyTree((TTree*)ifTree->Get("tree"));
		MyTree *myRotTree = new MyTree((TTree*)ifRotTree->Get("tree"));
		if(!myTree->mTree || !myRotTree->mTree) {
			std::cout << "ERROR: Open tree fail, skipping to next file" << std::endl;
			continue;
		}
		Long64_t nEvent = myTree->getNEvent();
		Long64_t nRotEvent = myRotTree->getNEvent();
		if(nEvent != nRotEvent) {
			std::cout << "ERROR: event not match" << std::endl;
			continue;
		}
		for(int ievt = 0; ievt < nEvent; ++ievt) {
			if(ievt % 10000 == 0) std::cout << "LOG: ievent: " << ievt << std::endl;
			myTree->getEntry(ievt);
			myRotTree->getEntry(ievt);
			float vx = myTree->mBufferVx;
			float vy = myTree->mBufferVy;
			float vz = myTree->mBufferVz;
			float vr = sqrt(vx*vx + vy*vy);
			int cent9 = (int)myTree->mBufferCent9;
			unsigned int nK = myTree->mBufferNTrack;
			unsigned int nRotK = myRotTree->mBufferNTrack;
			if(cent9 < 0 || cent9 > 8) continue;
			hist.hVz->Fill(vz);
			hist.hVr->Fill(vx, vy);
			hist.hCent9->Fill(cent9);
			std::vector<MyTree::Particle> vKaon, vRotKaon;
			for(int iK = 0; iK < nK; ++iK) {
				MyTree::Particle kaon = myTree->getParticle(iK, beamRapidity);
				if(!passAllCuts(kaon, config)) continue;
				hist.hSameKPtRapMass[cent9]->Fill(kaon.mass, kaon.rap, kaon.pt);
				if(!passCut(kaon.rap, config, "Rap")) continue;
				if(!passCut(kaon.pt, config, "Pt")) continue;
				if(!passCut(kaon.mass, config, "Mass")) continue;
				vKaon.push_back(std::move(kaon));
			}
			for(int iK = 0; iK < nRotK; ++iK) {
				MyTree::Particle kaon = myRotTree->getParticle(iK, beamRapidity);
				if(!passAllCuts(kaon, config)) continue;
				hist.hRotKPtRapMass[cent9]->Fill(kaon.mass, kaon.rap, kaon.pt);
				if(!passCut(kaon.rap, config, "Rap")) continue;
				if(!passCut(kaon.pt, config, "Pt")) continue;
				if(!passCut(kaon.mass, config, "Mass")) continue;
				vRotKaon.push_back(std::move(kaon));
			}

			for(int iK = 0; iK < vKaon.size(); ++iK) {
				hist.FillAll(vKaon[iK], cent9);
				hist.Fill(vKaon[iK]);
				hist.FillCut(vKaon[iK]);
			}
			loopVect(vKaon, vKaon, hist, cent9, "SS");
			loopVect(vKaon, vRotKaon, hist, cent9, "SR");
			loopVect(vRotKaon, vKaon, hist, cent9, "RS");
			loopVect(vRotKaon, vRotKaon, hist, cent9, "RR");

			std::list<std::vector<MyTree::Particle>>::iterator iter;
			for(iter = myTree->mMixBuffer[cent9].begin(); iter != myTree->mMixBuffer[cent9].end(); ++iter) {
				loopVect(vKaon, *iter, hist, cent9, "SMSS");
				loopVect(vKaon, *iter, hist, cent9, "SMW", hPurity);
				loopVect(vRotKaon, *iter, hist, cent9, "SMRS");
			}
			for(iter = myRotTree->mMixBuffer[cent9].begin(); iter != myRotTree->mMixBuffer[cent9].end(); ++iter) {
				loopVect(vKaon, *iter, hist, cent9, "SMSR");
				loopVect(vRotKaon, *iter, hist, cent9, "SMRR");
			}

			myTree->copyToBuffer(vKaon);
			myRotTree->copyToBuffer(vRotKaon);
		}

		delete myTree, myRotTree;
		ifTree->Close(), ifRotTree->Close();
		delete ifTree, ifRotTree;
		myTree = nullptr, myRotTree = nullptr;
		ifTree = nullptr, ifRotTree = nullptr;
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
	TVector3 mom(p.px, p.py, p.pz);
	TVector3 pos(p.bx, p.by, p.bz);
	TVector3 vectPmPos = pos - mom;
	float cosTheta = vectPmPos.CosTheta();
	int momBinA = getMomBin(p.pA, p_low, p_high, nPBins);
	int momBinB = getMomBin(p.pB, p_low, p_high, nPBins);

	float nsigmaA = p.nSigmaA;
	float nsigmaB = p.nSigmaB;
	if(config.mSwitchList["OpenNSigmaShift"]) {
		nsigmaA -= pidCalib_pion[momBinA];
		nsigmaB -= pidCalib_pion[momBinB];
	}

	if(!passCut(p.chi2Topo, config, "Chi2Topo")) return false;
	if(!passCut(p.chi2NDF, config, "Chi2NDF")) return false;
	if(!passCut(p.chi2PrimPip, config, "Chi2PrimPip")) return false;
	if(!passCut(p.chi2PrimPim, config, "Chi2PrimPim")) return false;
	//if(!passCut(p.rap, config, "Rap")) return false;
	//if(!passCut(p.pt, config, "Pt")) return false;
	if(!passCut(p.nHitsA, config, "NHitsA")) return false;
	if(!passCut(p.nHitsB, config, "NHitsB")) return false;
	if(!passCut(p.dcaA, config, "DCAA")) return false;
	if(!passCut(p.dcaB, config, "DCAB")) return false;
	if(!passCut(p.m2A, p.pA, config, "Mass2Pip")) return false;
	if(!passCut(p.m2B, p.pB, config, "Mass2Pim")) return false;
	if(!passCut(p.dca, config, "KS0DCA")) return false;
	//if(!passCut(p.decayLength, config, "DecayLength")) return false;
	//if(!passCut(p.mass, config, "Mass")) return false;
	if(!passCut(nsigmaA, config, "NSigmaPi")) return false;
	if(!passCut(nsigmaB, config, "NSigmaPi")) return false;
	//if(p.rap > 0 && p.pt < 0.3) return false;
	if(p.etaA < -2.0 || p.etaA > 0) return false;
	if(p.etaB < -2.0 || p.etaB > 0) return false;
	return true;
}

void loopVect(std::vector<MyTree::Particle>& v1, std::vector<MyTree::Particle>& v2, Hist& hist, int cent9, std::string opt, TH2F* hPurity[])
{
	for(int iK1 = 0; iK1 < v1.size(); ++iK1) {
		int startIdx = ((opt == "SS") || (opt == "RR")) ? iK1 + 1 : 0;
		for(int iK2 = startIdx; iK2 < v2.size(); ++iK2) {
			if(v1[iK1].ptPip == v2[iK2].ptPip && v1[iK1].rapPip == v2[iK2].rapPip && v1[iK1].phiA == v2[iK2].phiA) continue;
			if(v1[iK1].ptPim == v2[iK2].ptPim && v1[iK1].rapPim == v2[iK2].rapPim && v1[iK1].phiB == v2[iK2].phiB) continue;
			TLorentzVector k1_v4, k2_v4;
			k1_v4.SetXYZT(v1[iK1].px, v1[iK1].py, v1[iK1].pz, v1[iK1].energy);
			k2_v4.SetXYZT(v2[iK2].px, v2[iK2].py, v2[iK2].pz, v2[iK2].energy);
			TLorentzVector kDiff_v4 = (k1_v4 - k2_v4);
			float qinv = fabs(kDiff_v4.Mag());
			if(opt == "SS") {
				hist.hSameKqinv[cent9]->Fill(qinv);
			} else if(opt == "SR") {
				hist.hSRQinv[cent9]->Fill(qinv);
			} else if(opt == "RS") {
				hist.hRSQinv[cent9]->Fill(qinv);
			} else if(opt == "RR") {
				hist.hRRQinv[cent9]->Fill(qinv);
			} else if(opt == "SMW") {
				int binX1 = hPurity[2]->GetXaxis()->FindBin(v1[iK1].rap);
				int binY1 = hPurity[2]->GetYaxis()->FindBin(v1[iK1].pt);
				int binX2 = hPurity[2]->GetXaxis()->FindBin(v2[iK2].rap);
				int binY2 = hPurity[2]->GetYaxis()->FindBin(v2[iK2].pt);
				float purity1 = hPurity[2]->GetBinContent(binX1, binY1) / 100.;
				float purity2 = hPurity[2]->GetBinContent(binX2, binY2) / 100.;
				float ppSS = purity1 * purity2;
				float ppSR = (1 - purity2);
				float ppRS = (1 - purity1);
				float ppRR = -1 * (1 - purity1) * (1 - purity2);
				hist.hMixKqinvWeight[cent9][0]->Fill(qinv, ppSS);
				hist.hMixKqinvWeight[cent9][1]->Fill(qinv, ppSR);
				hist.hMixKqinvWeight[cent9][2]->Fill(qinv, ppRS);
				hist.hMixKqinvWeight[cent9][3]->Fill(qinv, ppRR);
			} else if(opt == "SMSS") {
				hist.hMixKqinv[cent9]->Fill(qinv);
			} else if(opt == "SMSR") {
				hist.hMixSRQinv[cent9]->Fill(qinv);
			} else if(opt == "SMRS") {
				hist.hMixRSQinv[cent9]->Fill(qinv);
			} else if(opt == "SMRR") {
				hist.hMixRRQinv[cent9]->Fill(qinv);
			}
		}
	}
}
//}}}
