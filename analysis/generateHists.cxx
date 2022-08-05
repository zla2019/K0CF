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

const int nPBins = 23;
//PID calibration
const float p_low[nPBins]  = {0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8};
const float p_high[nPBins] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 10.0};
const float pidCalib_pion[nPBins]  = {1.30, 1.49, 1.53, 1.54, 1.54, 1.55, 1.56, 1.57, 1.59, 1.61, 1.62, 1.63, 1.65, 1.66, 1.68, 1.69, 1.71, 1.73, 1.75, 1.77, 1.79, 1.82, 1.83};

bool passCut(double a, double lower, double upper);
int getMomBin(float mom, const float* MomLower, const float* MomUpper, const float NMom);
bool loadBadRun(std::unordered_map<int, bool>* badRunList, std::ifstream* ifBadRunList);
bool passCut(double a, Config& config, std::string cutName);
bool passCut(double a, double b, Config& config, std::string cutName);
inline bool passAllCuts(MyTree::Particle& p, Config& config);

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
	//std::ifstream ifBadRunList("/star/u/lazhang/data01/CF/correlation_test/new_badruns_3.list");
	std::unordered_map<int, bool> badRunList;
	//loadBadRun(&badRunList, &ifBadRunList);
	const float Kmass = 0.497611;
	//mass window info
	//const float Mean[3][4] = { { 0.4979, 0.4981, 0.4979, 0.4949 }, { 0.4978, 0.4981, 0.4979, 0.4979 }, { 0.4979, 0.4981, 0.4980, 0.4980 } };
	//const float Sigma[3][4] = { { 0.0027, 0.0035, 0.0033, 0.0033 }, { 0.0027, 0.0034, 0.0036, 0.0036 }, { 0.0027, 0.0035, 0.0039, 0.0039 } };
	const float Mean[3][4] = { { 0.4981, 0.4981, 0.4981, 0.4981 }, { 0.4981, 0.4981, 0.4981, 0.4981 }, { 0.4981, 0.4981, 0.4981, 0.4981 } };
	const float Sigma[3][4] = { { 0.0035, 0.0035, 0.0035, 0.0035 }, { 0.0035, 0.0035, 0.0035, 0.0035 }, { 0.0035, 0.0035, 0.0035, 0.0035 } };
	const int NMassSigma = std::stoi(config.mSetList["NSigmaMass"]);
	const bool MuteWarning = true;	//used for debug

	float beamRapidity;
	if(config.mSetList["Energy"] == "3.0") beamRapidity = 1.045;
	else if(config.mSetList["Energy"] == "3.2") beamRapidity = 1.135;
	else if(config.mSetList["Energy"] == "3.5") beamRapidity = 1.24;
	//}}}

	//Save cut info{{{
	TCanvas* caCut = new TCanvas("cuts", "cuts", 1920, 1080);
	caCut->cd();
	TLatex* ltx = new TLatex(0.0, 0.9, "Cuts:");
	//}}}

	//prepare plots{{{
	TH1F* hSameKPDG = new TH1F("hSameKPDG", "Kaon pdg distribution", 5, 280, 340);
	TH1F* hVz = new TH1F("hVz", "V_{z} distribution;cm;cnts", 100, 198, 202);
	TH2F* hVr = new TH2F("hVr", "V_{r} distribution;cm;cm", 100, -2, 2, 100, -4, 0);
	TH1F* hCent9 = new TH1F("hCent9", "Cent9 dist.", 10, -1, 9);

	TH1F* hCosTheta = new TH1F("hCosTheta", "cos(#theta) distribution", 50, 0.95, 1.0);
	TH1F* hDecayLength = new TH1F("hDecayLength", "K^{0}_{s} decay length", 100, 0, 20);
	TH1F* hDgDCA = new TH1F("hDgDCA", "K_{s}^{0} daughter DCA", 100, 0, 2);
	TH1F* hDCA = new TH1F("hDCA", "K_{s}^{0} DCA", 100, 0, 2);

	TH2F* hDedx = new TH2F("hDedx", "dEdx vs p*q;p*q;dEdx", 1000, -5, 5, 500, 1.5, 6.5);
	TH2F* hMass2 = new TH2F("hMass2", "m^{2} vs p*q;p*q;m^{2}", 1000, -5, 5, 500, -2, 3);

	TH2F* hSameKRapPt = new TH2F("hSameKRapPt", "K_{s}^{0} Acc.", 200, -1, 1, 300, 0, 3);
	TH1F* hSameKMass = new TH1F("hSameKMass", "K_{s}^{0} mass distribution;M_{inv};cnts", 160, 0.48, 0.52);
	TH1F* hSameKPhi = new TH1F("hSameKPhi", "K_{s}^{0} #phi distribution;#phi;cnts", 100, -TMath::Pi(), TMath::Pi());
	TH2F* hSameKPipRapPt = new TH2F("hSameKPipRapPt", "K_{s}^{0} daughter #pi^{+} Acc.", 200, -1, 1, 500, 0, 5);
	TH2F* hSameKPimRapPt = new TH2F("hSameKPimRapPt", "K_{s}^{0} daughter #pi^{-} Acc.", 200, -1, 1, 500, 0, 5);
	TH1F* hDaughterPipDCA = new TH1F("hDaughterPipDCA", "K_{s}^{0} daughter #pi^{+} DCA distribution", 50, 0, 10);
	TH1F* hDaughterPimDCA = new TH1F("hDaughterPimDCA", "K_{s}^{0} daughter #pi^{-} DCA distribution", 50, 0, 10);

	TH2F* hLeftSideSameKRapPt = new TH2F("hLeftSideSameKRapPt", "K_{s}^{0} Acc.", 200, -1, 1, 300, 0, 3);
	TH1F* hLeftSideSameKMass = new TH1F("hLeftSideSameKMass", "K_{s}^{0} mass distribution;M_{inv};cnts", 1200, 0.4, 0.70);
	TH1F* hLeftSideSameKPhi = new TH1F("hLeftSideSameKPhi", "K_{s}^{0} #phi distribution;#phi;cnts", 100, -TMath::Pi(), TMath::Pi());
	TH2F* hLeftSideSameKPipRapPt = new TH2F("hLeftSideSameKPipRapPt", "K_{s}^{0} daughter #pi^{+} Acc.", 200, -1, 1, 500, 0, 5);
	TH2F* hLeftSideSameKPimRapPt = new TH2F("hLeftSideSameKPimRapPt", "K_{s}^{0} daughter #pi^{-} Acc.", 200, -1, 1, 500, 0, 5);

	TH2F* hRightSideSameKRapPt = new TH2F("hRightSideSameKRapPt", "K_{s}^{0} Acc.", 200, -1, 1, 300, 0, 3);
	TH1F* hRightSideSameKMass = new TH1F("hRightSideSameKMass", "K_{s}^{0} mass distribution;M_{inv};cnts", 1200, 0.4, 0.70);
	TH1F* hRightSideSameKPhi = new TH1F("hRightSideSameKPhi", "K_{s}^{0} #phi distribution;#phi;cnts", 100, -TMath::Pi(), TMath::Pi());
	TH2F* hRightSideSameKPipRapPt = new TH2F("hRightSideSameKPipRapPt", "K_{s}^{0} daughter #pi^{+} Acc.", 200, -1, 1, 500, 0, 5);
	TH2F* hRightSideSameKPimRapPt = new TH2F("hRightSideSameKPimRapPt", "K_{s}^{0} daughter #pi^{-} Acc.", 200, -1, 1, 500, 0, 5);

	TH1F* hPipNSigma = new TH1F("hPipNSigma", "#pi^{+} n#sigma distribution", 100, -5, 5);
	TH1F* hPimNSigma = new TH1F("hPimNSigma", "#pi^{-} n#sigma distribution", 100, -5, 5);
	TH1F* hPipNSigma2 = new TH1F("hPipNSigma2", "#pi^{+} n#sigma distribution", 100, -5, 5);
	TH1F* hPimNSigma2 = new TH1F("hPimNSigma2", "#pi^{-} n#sigma distribution", 100, -5, 5);

	TH1F* hDTheta = new TH1F("hDTheta", "K_{s}^{0}-K_{s}^{0} #Delta#theta dist.", 314, 0, TMath::Pi());
	TH1F* hDPhi = new TH1F("hDPhi", "K_{s}^{0}-K_{s}^{0} #Delta#phi dist.", 314, 0, TMath::Pi());
	TH2F* hDThetaDPhi = new TH2F("hDThetaDPhi", "K_{s}^{0}-K_{s}^{0} #DeltaTheta vs. #DeltaPhi", 314, 0, TMath::Pi(), 314, 0, TMath::Pi());

	TH3F* hSameKPtRapMass[9];
	TH1F* hSameKqinv[9];
	TH1F* hMixKqinv[9];
	TH1F* hMixKqinvWeight[9][4];		//cent, rap, case;	//case0: peak*side, case1: side*peak, case2: side*side, case3: peak*peak
	TH1F* hMixKqinvLeftWeight[9][4];	//cent, rap, case;	//case0: peak*side, case1: side*peak, case2: side*side, case3: peak*peak
	TH1F* hMixKqinvRightWeight[9][4];	//cent, rap, case;	//case0: peak*side, case1: side*peak, case2: side*side, case3: peak*peak
	TH1F* hSameKqlong[9];
	TH1F* hMixKqlong[9];
	TH1F* hSameKqout[9];
	TH1F* hMixKqout[9];
	TH1F* hSameKqside[9];
	TH1F* hMixKqside[9];

	TH1F* hSameLeftSideKqinv[9][4];	//cent, rap, case. where case mean how we make a side pair 
	TH1F* hMixLeftSideKqinv[9][4];
	TH1F* hSameRightSideKqinv[9][4];	//cent, rap, case. where case mean how we make a side pair 
	TH1F* hMixRightSideKqinv[9][4];
	for(int icent = 0; icent < 9; ++icent) {
		hSameKPtRapMass[icent] = new TH3F(Form("hMassCascadeRapidityvsPt_cent%d", icent),Form("hMassCascadeRapidityvsPt_cent%d", icent),160, 0.42, 0.58, 20, -1, 1, 30, 0, 3);
		hSameKqinv[icent] = new TH1F(Form("hSameKqinv_cent%i", icent), Form("K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 500, 0, 1);
		hMixKqinv[icent] = new TH1F(Form("hMixKqinv_cent%i", icent), Form("Mix K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 500, 0, 1);
		hSameKqlong[icent] = new TH1F(Form("hSameKqlong_cent%i", icent), Form("K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 2000, -2, 2);
		hMixKqlong[icent] = new TH1F(Form("hMixKqlong_cent%i", icent), Form("Mix K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 2000, -2, 2);
		hSameKqout[icent] = new TH1F(Form("hSameKqout_cent%i", icent), Form("K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 2000, -2, 2);
		hMixKqout[icent] = new TH1F(Form("hMixKqout_cent%i", icent), Form("Mix K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 2000, -2, 2);
		hSameKqside[icent] = new TH1F(Form("hSameKqside_cent%i", icent), Form("K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 2000, -2, 2);
		hMixKqside[icent] = new TH1F(Form("hMixKqside_cent%i", icent), Form("Mix K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 2000, -2, 2);
		for(int icase = 0; icase < 4; ++icase) {
			hSameLeftSideKqinv[icent][icase] = new TH1F(Form("hSameLeftSideKqinv_cent%i_case%i", icent, icase), Form("LeftSide Band K_{s}^{0} q_{inv}(k^{*}) @cent%i, case%i", icent, icase), 500, 0, 1);
			hMixLeftSideKqinv[icent][icase] = new TH1F(Form("hMixLeftSideKqinv_cent%i_case%i", icent, icase), Form("LeftSide Band Mix K_{s}^{0} q_{inv}(k^{*}) @cent%i, case%i", icent, icase), 500, 0, 1);
			hSameRightSideKqinv[icent][icase] = new TH1F(Form("hSameRightSideKqinv_cent%i_case%i", icent, icase), Form("RightSide Band K_{s}^{0} q_{inv}(k^{*}) @cent%i, case%i", icent, icase), 500, 0, 1);
			hMixRightSideKqinv[icent][icase] = new TH1F(Form("hMixRightSideKqinv_cent%i_case%i", icent, icase), Form("RightSide Band Mix K_{s}^{0} q_{inv}(k^{*}) @cent%i, case%i", icent, icase), 500, 0, 1);
			hMixKqinvWeight[icent][icase] = new TH1F(Form("hMixKqinvWeight_cent%i_case%i", icent, icase), Form("Weighted Mix K_{s}^{0} q_{inv} for SideBand @cent%i, case%i", icent, icase), 500, 0, 1);
			hMixKqinvLeftWeight[icent][icase] = new TH1F(Form("hMixKqinvLeftWeight_cent%i_case%i", icent, icase), Form("LeftWeighted Mix K_{s}^{0} q_{inv} for SideBand @cent%i, case%i", icent, icase), 500, 0, 1);
			hMixKqinvRightWeight[icent][icase] = new TH1F(Form("hMixKqinvRightWeight_cent%i_case%i", icent, icase), Form("RightWeighted Mix K_{s}^{0} q_{inv} for SideBand @cent%i, case%i", icent, icase), 500, 0, 1);
		}
	}

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
	while(getline(ifFilelist, ifName)) {
		TFile* ifTree = TFile::Open(ifName.c_str());
		if(!ifTree) {
			std::cout << ifName.c_str() << " can not open" << std::endl;
			continue;
		}
		MyTree *myTree = new MyTree((TTree*)ifTree->Get("tree"));
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
			if(badRunList[myTree->mBufferRunId] == true) {
				if(!MuteWarning) {
					std::cout << "WARNING: Bad run " << myTree->mBufferRunId << " event found, will continue this event" << std::endl;
				}
				continue;
			}
			float vx = myTree->mBufferVx;
			float vy = myTree->mBufferVy;
			float vz = myTree->mBufferVz;
			float vr = sqrt(vx*vx + vy*vy);
			float cent9 = myTree->mBufferCent9;
			int centForPurity;
			if(cent9 >= 7) {
				centForPurity = 0;
			} else if(cent9 >= 2) {
				centForPurity = 1;
			} else {
				centForPurity = 2;
			}
			unsigned int nK = myTree->mBufferNTrack;
			if(cent9 < 0 || cent9 > 8) {
				continue;
			}

			hVz->Fill(vz);
			hVr->Fill(vx, vy);
			hCent9->Fill((int)cent9);
			//pre select
			std::vector<int> idxK;
			for(int icurK = 0; icurK < nK; ++icurK) {
				MyTree::Particle curKTmp = myTree->getParticle(icurK, beamRapidity);
				if(!passAllCuts(curKTmp, config)) continue;
				idxK.push_back(icurK);
			}
			int nPassCutK = idxK.size();
			for(int iK1 = 0; iK1 < nPassCutK; ++iK1) {
				int icurK = idxK[iK1];
				MyTree::Particle curK = myTree->getParticle(icurK, beamRapidity);
				TVector3 p(curK.px, curK.py, curK.pz);
				TVector3 pos(curK.bx, curK.by, curK.bz);
				TVector3 vectPmPos = pos - p;
				float cosTheta = vectPmPos.CosTheta();
				int momBinA = getMomBin(curK.pA, p_low, p_high, nPBins);
				int momBinB = getMomBin(curK.pB, p_low, p_high, nPBins);

				float nsigmaA = curK.nSigmaA;
				float nsigmaB = curK.nSigmaB;
				if(config.mSwitchList["OpenNSigmaShift"]) {
					nsigmaA -= pidCalib_pion[momBinA];
					nsigmaB -= pidCalib_pion[momBinB];
				}

				float purity = 1;
				if(config.mSwitchList["OpenPairPurity"]) {
					purity = hPurity[centForPurity]->GetBinContent(hPurity[centForPurity]->GetXaxis()->FindBin(curK.rap), hPurity[centForPurity]->GetYaxis()->FindBin(curK.pt)) / 100.;
				}
				int isSideBand = -1;	//defualt: -1; peak region: 0; left side: 1; right side: 2;

				hDaughterPipDCA->Fill(curK.dcaA);
				hDaughterPimDCA->Fill(curK.dcaB);
				hCosTheta->Fill(cosTheta);
				hDecayLength->Fill(curK.decayLength);
				hDgDCA->Fill(curK.dgDCA);
				hDCA->Fill(curK.dca);
				hSameKPtRapMass[(int)cent9]->Fill(curK.mass, curK.rap, curK.pt);

				int rapBin = 1;
				if(passCut(curK.mass, Mean[2][rapBin] - (NMassSigma * Sigma[2][rapBin]), Mean[2][rapBin] + (NMassSigma * Sigma[2][rapBin]))) {
					//if(passCut(curK.mass, 0.48, 0.51)) {
					hSameKRapPt->Fill(curK.rap, curK.pt);
					isSideBand = 0;
					//} else if(passCut(mass, SideBand2Lower, SideBand2Upper)) {
				} else if(passCut(curK.mass, config, "SideBand2")) {
					isSideBand = 1;
				} else if(passCut(curK.mass, config, "SideBand")) {
					isSideBand = 2;
				} else {
					continue;
				}

				//if(curK.rap < -0.8 || curK.rap > 0.4) continue;
				if(isSideBand == 0) {
					hSameKPDG->Fill(curK.pdg);
					hSameKMass->Fill(curK.mass);
					hSameKPipRapPt->Fill(curK.rapPip, curK.ptPip);
					hSameKPimRapPt->Fill(curK.rapPim, curK.ptPim);
					hSameKPhi->Fill(curK.phi);
					hDedx->Fill(curK.pA, curK.dEdxA);
					hDedx->Fill(-curK.pB, curK.dEdxB);
					hMass2->Fill(curK.pA, curK.m2A);
					hMass2->Fill(-curK.pB, curK.m2B);
					hPipNSigma->Fill(nsigmaA);
					hPimNSigma->Fill(nsigmaB);
					hPipNSigma2->Fill(curK.nSigmaA);
					hPimNSigma2->Fill(curK.nSigmaB);
				} else if(isSideBand == 1) {
					hLeftSideSameKMass->Fill(curK.mass);
					hLeftSideSameKPipRapPt->Fill(curK.rapPip, curK.ptPip);
					hLeftSideSameKPimRapPt->Fill(curK.rapPim, curK.ptPim);
					hLeftSideSameKRapPt->Fill(curK.rap, curK.pt);
					hLeftSideSameKPhi->Fill(curK.phi);
				} else if(isSideBand == 2) {
					hRightSideSameKMass->Fill(curK.mass);
					hRightSideSameKPipRapPt->Fill(curK.rapPip, curK.ptPip);
					hRightSideSameKPimRapPt->Fill(curK.rapPim, curK.ptPim);
					hRightSideSameKRapPt->Fill(curK.rap, curK.pt);
					hRightSideSameKPhi->Fill(curK.phi);
				}
				if(!config.mSwitchList["OpenCF"]) continue;
				//same pair{{{
				for(int iK2 = iK1 + 1; iK2 < nPassCutK; ++iK2) {
					int icurK2 = idxK[iK2];
					MyTree::Particle curK2 = myTree->getParticle(icurK2, beamRapidity);
					TVector3 p(curK2.px, curK2.py, curK2.pz);
					TVector3 pos(curK2.bx, curK2.by, curK2.bz);
					TVector3 vectPmPos = pos - p;
					float cosTheta2 = vectPmPos.CosTheta();
					int momBinA = getMomBin(curK2.pA, p_low, p_high, nPBins);
					int momBinB = getMomBin(curK2.pB, p_low, p_high, nPBins);
					float purity2 = 1;
					if(config.mSwitchList["OpenPairPurity"]) {
						purity2 = hPurity[centForPurity]->GetBinContent(hPurity[centForPurity]->GetXaxis()->FindBin(curK2.rap), hPurity[centForPurity]->GetYaxis()->FindBin(curK2.pt)) / 100.;
					}

					int isSideBand2 = -1;
					int rapBin2 = 1;
					//if(curK2.rap < -0.8 || curK2.rap > 0.4) continue;

					if(passCut(curK2.mass, Mean[2][rapBin2] - (NMassSigma * Sigma[2][rapBin2]), Mean[2][rapBin2] + (NMassSigma * Sigma[2][rapBin2]))) {
						//if(passCut(curK2.mass, 0.48, 0.51)) {
						isSideBand2 = 0;
					} else if(passCut(curK2.mass, config, "SideBand")) {
						isSideBand2 = 2;
					} else if(passCut(curK2.mass, config, "SideBand2")) {
						isSideBand2 = 1;
					} else {
						continue;
					}

					if(icurK == icurK2) {
						continue;
					}
					if(curK.ptPip == curK2.ptPip && curK.rapPip == curK2.rapPip/* && curK.trkIdA == curK2.trkIdA*/) {
						continue;
					}
					if(curK.ptPim == curK2.ptPim && curK.rapPim == curK2.rapPim/* && curK.trkIdB == curK2.trkIdB*/) {
						continue;
					}

					TLorentzVector k1_v4, k2_v4;
					//k1_v4.SetXYZT(curK.px, curK.py, curK.pz, sqrt(curK.px*curK.px + curK.py*curK.py + curK.pz*curK.pz + 0.497611*0.497611));
					//k2_v4.SetXYZT(curK2.px, curK2.py, curK2.pz, sqrt(curK2.px*curK2.px + curK2.py*curK2.py + curK2.pz*curK2.pz + 0.497611*0.497611));
					k1_v4.SetXYZT(curK.px, curK.py, curK.pz, curK.energy);
					k2_v4.SetXYZT(curK2.px, curK2.py, curK2.pz, curK2.energy);
					TLorentzVector kDiff_v4 = (k1_v4 - k2_v4);

					if(isSideBand == 0 && isSideBand2 == 0) {
						hSameKqinv[(int)cent9]->Fill(fabs(kDiff_v4.Mag()));
						hSameLeftSideKqinv[(int)cent9][3]->Fill(fabs(kDiff_v4.Mag()));
						hSameRightSideKqinv[(int)cent9][3]->Fill(fabs(kDiff_v4.Mag()));
					} else if(isSideBand == 0 && isSideBand2 == 1) {
						hSameLeftSideKqinv[(int)cent9][0]->Fill(fabs(kDiff_v4.Mag()));
					} else if(isSideBand == 0 && isSideBand2 == 2) {
						hSameRightSideKqinv[(int)cent9][0]->Fill(fabs(kDiff_v4.Mag()));
					} else if(isSideBand == 1 && isSideBand2 == 0) {
						hSameLeftSideKqinv[(int)cent9][1]->Fill(fabs(kDiff_v4.Mag()));
					} else if(isSideBand == 2 && isSideBand2 == 0) {
						hSameRightSideKqinv[(int)cent9][1]->Fill(fabs(kDiff_v4.Mag()));
					} else if(isSideBand == 1 && isSideBand2 == 1) {
						hSameLeftSideKqinv[(int)cent9][2]->Fill(fabs(kDiff_v4.Mag()));
					} else if(isSideBand == 2 && isSideBand2 == 2) {
						hSameRightSideKqinv[(int)cent9][2]->Fill(fabs(kDiff_v4.Mag()));
					}
				}
				//}}}
				//mix pair{{{
				for(int imixevt = 0; imixevt < myTree->mMaxMixEvent[(int)cent9] + 1; ++imixevt) {
					unsigned int nK = myTree->mMixBuffer[(int)cent9][imixevt].mBufferNTrack;
					for(int imixK = 0; imixK < nK; ++imixK) {
						MyTree::Particle mixK = myTree->getMixParticle((int)cent9, imixevt, imixK, beamRapidity);
						int momBinA = getMomBin(mixK.pA, p_low, p_high, nPBins);
						int momBinB = getMomBin(mixK.pB, p_low, p_high, nPBins);
						float purity2 = 1;
						if(config.mSwitchList["OpenPairPurity"]) {
							purity2 = hPurity[centForPurity]->GetBinContent(hPurity[centForPurity]->GetXaxis()->FindBin(mixK.rap), hPurity[centForPurity]->GetYaxis()->FindBin(mixK.pt)) / 100.;
						}
						int isSideBand2 = -1;
						int rapBin2 = 1;
						//if(mixK.rap < -0.8 || mixK.rap > 0.4) continue;

						if(passCut(mixK.mass, Mean[2][rapBin2] - (NMassSigma * Sigma[2][rapBin2]), Mean[2][rapBin2] + (NMassSigma * Sigma[2][rapBin2]))) {
							//if(passCut(mixK.mass, 0.48, 0.51)) {
							isSideBand2 = 0;
						} else if(passCut(mixK.mass, config, "SideBand")){
							isSideBand2 = 2;
						} else if(passCut(mixK.mass, config, "SideBand2")) {
							isSideBand2 = 1;
						} else {
							continue;
						}

						float sideBandWeight[4] = { 0 };
						if(config.mSwitchList["OpenPairPurity"]) {
							sideBandWeight[0] = purity * (1 - purity2);
							sideBandWeight[1] = (1 - purity) * purity2;
							sideBandWeight[2] = (1 - purity) * (1 - purity2);
							sideBandWeight[3] = purity * purity2;
						}
						TLorentzVector k1_v4, k2_v4;
						k1_v4.SetXYZT(curK.px, curK.py, curK.pz, curK.energy);
						k2_v4.SetXYZT(mixK.px, mixK.py, mixK.pz, mixK.energy);
						//k1_v4.SetXYZT(curK.px, curK.py, curK.pz, sqrt(curK.px*curK.px + curK.py*curK.py + curK.pz*curK.pz + 0.497611*0.497611));
						//k2_v4.SetXYZT(mixK.px, mixK.py, mixK.pz, sqrt(mixK.px*mixK.px + mixK.py*mixK.py + mixK.pz*mixK.pz + 0.497611*0.497611));
						TLorentzVector kDiff_v4 = (k1_v4 - k2_v4);

						if(isSideBand == 0 && isSideBand2 == 0) {
							hMixKqinv[(int)cent9]->Fill(fabs(kDiff_v4.Mag()));
							hMixLeftSideKqinv[(int)cent9][3]->Fill(fabs(kDiff_v4.Mag()));
							hMixRightSideKqinv[(int)cent9][3]->Fill(fabs(kDiff_v4.Mag()));
							hMixKqinvWeight[(int)cent9][0]->Fill(fabs(kDiff_v4.Mag()), sideBandWeight[0]);
							hMixKqinvWeight[(int)cent9][1]->Fill(fabs(kDiff_v4.Mag()), sideBandWeight[1]);
							hMixKqinvWeight[(int)cent9][2]->Fill(fabs(kDiff_v4.Mag()), sideBandWeight[2]);
							hMixKqinvWeight[(int)cent9][3]->Fill(fabs(kDiff_v4.Mag()), sideBandWeight[3]);
							hMixKqinvLeftWeight[(int)cent9][3]->Fill(fabs(kDiff_v4.Mag()), sideBandWeight[3]);
							hMixKqinvRightWeight[(int)cent9][3]->Fill(fabs(kDiff_v4.Mag()), sideBandWeight[3]);
						} else if(isSideBand == 0 && isSideBand2 == 1) {
							hMixLeftSideKqinv[(int)cent9][0]->Fill(fabs(kDiff_v4.Mag()));
							hMixKqinvLeftWeight[(int)cent9][0]->Fill(fabs(kDiff_v4.Mag()), sideBandWeight[0]);
						} else if(isSideBand == 0 && isSideBand2 == 2) {
							hMixRightSideKqinv[(int)cent9][0]->Fill(fabs(kDiff_v4.Mag()));
							hMixKqinvRightWeight[(int)cent9][0]->Fill(fabs(kDiff_v4.Mag()), sideBandWeight[0]);
						} else if(isSideBand == 1 && isSideBand2 == 0) {
							hMixLeftSideKqinv[(int)cent9][1]->Fill(fabs(kDiff_v4.Mag()));
							hMixKqinvLeftWeight[(int)cent9][1]->Fill(fabs(kDiff_v4.Mag()), sideBandWeight[1]);
						} else if(isSideBand == 2 && isSideBand2 == 0) {
							hMixRightSideKqinv[(int)cent9][1]->Fill(fabs(kDiff_v4.Mag()));
							hMixKqinvRightWeight[(int)cent9][1]->Fill(fabs(kDiff_v4.Mag()), sideBandWeight[1]);
						} else if(isSideBand == 1 && isSideBand2 == 1) {
							hMixLeftSideKqinv[(int)cent9][2]->Fill(fabs(kDiff_v4.Mag()));
							hMixKqinvLeftWeight[(int)cent9][2]->Fill(fabs(kDiff_v4.Mag()), sideBandWeight[2]);
						} else if(isSideBand == 2 && isSideBand2 == 2) {
							hMixRightSideKqinv[(int)cent9][2]->Fill(fabs(kDiff_v4.Mag()));
							hMixKqinvRightWeight[(int)cent9][2]->Fill(fabs(kDiff_v4.Mag()), sideBandWeight[2]);
						}
					}
				}
				//}}}
			}
			myTree->copyToBuffer(idxK);
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
	ofPlots->cd();
	caCut->Write();
	hVz->Write();
	hVr->Write();
	hCent9->Write();
	hSameKPDG->Write();
	hSameKMass->Write();
	hSameKPhi->Write();
	hSameKRapPt->Write();
	hSameKPipRapPt->Write();
	hSameKPimRapPt->Write();

	hLeftSideSameKMass->Write();
	hLeftSideSameKPhi->Write();
	hLeftSideSameKRapPt->Write();
	hLeftSideSameKPipRapPt->Write();
	hLeftSideSameKPimRapPt->Write();
	hRightSideSameKMass->Write();
	hRightSideSameKPhi->Write();
	hRightSideSameKRapPt->Write();
	hRightSideSameKPipRapPt->Write();
	hRightSideSameKPimRapPt->Write();

	hDaughterPipDCA->Write();
	hDaughterPimDCA->Write();
	hCosTheta->Write();
	hDecayLength->Write();
	hDgDCA->Write();
	hDCA->Write();
	hDedx->Write();
	hMass2->Write();
	hPipNSigma->Write();
	hPimNSigma->Write();
	hPipNSigma2->Write();
	hPimNSigma2->Write();

	hDTheta->Write();
	hDPhi->Write();
	hDThetaDPhi->Write();
	for(int icent = 0; icent < 9; ++icent) {
		hSameKPtRapMass[icent]->Write();
		hSameKqinv[icent]->Write();
		hMixKqinv[icent]->Write();
		if(config.mSwitchList["Open3DCF"]) {
			hSameKqlong[icent]->Write();
			hMixKqlong[icent]->Write();
			hSameKqout[icent]->Write();
			hMixKqout[icent]->Write();
			hSameKqside[icent]->Write();
			hMixKqside[icent]->Write();
		}
		for(int icase = 0; icase < 4; ++icase) {
			hSameLeftSideKqinv[icent][icase]->Write();
			hMixLeftSideKqinv[icent][icase]->Write();
			hSameRightSideKqinv[icent][icase]->Write();
			hMixRightSideKqinv[icent][icase]->Write();
			hMixKqinvWeight[icent][icase]->Write();
			hMixKqinvLeftWeight[icent][icase]->Write();
			hMixKqinvRightWeight[icent][icase]->Write();
		}
	}
	delete ofPlots;
	//}}}
	second = time(NULL);
	std::cout << "processing time: " << difftime(second, first) << std::endl;
	return 0;
}

//my tools{{{
bool loadBadRun(std::unordered_map<int, bool>* badRunList, std::ifstream* ifBadRunList)
{
	while(!ifBadRunList->eof()) {
		std::string badRunStr = "";
		*ifBadRunList >> badRunStr;
		if(badRunStr == "") {
			continue;
		}
		int badRun = std::stoi(badRunStr);
		(*(badRunList))[badRun] = true;
	}
	std::cout << "LOG: " << badRunList->size() << " Bad Runs Loaded" << std::endl;
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
	if(!passCut(p.rap, config, "Rap")) return false;
	if(!passCut(p.pt, config, "Pt")) return false;
	if(!passCut(p.nHitsA, config, "NHitsA")) return false;
	if(!passCut(p.nHitsB, config, "NHitsB")) return false;
	if(!passCut(p.dcaA, config, "DCAA")) return false;
	if(!passCut(p.dcaB, config, "DCAB")) return false;
	if(!passCut(p.m2A, p.pA, config, "Mass2Pip")) return false;
	if(!passCut(p.m2B, p.pB, config, "Mass2Pim")) return false;
	if(!passCut(p.dca, config, "KS0DCA")) return false;
	//if(!passCut(p.decayLength, config, "DecayLength")) return false;
	if(!passCut(nsigmaA, config, "NSigmaPi")) return false;
	if(!passCut(nsigmaB, config, "NSigmaPi")) return false;
	if(p.rap > 0 && p.pt < 0.3) return false;
	return true;
}
//}}}
