#ifndef HIST_H
#define HIST_H
#include <TTree.h>
#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include "MyTree.h"

class Hist
{
public:
	Hist() {};
	~Hist() {};
	void init();
	void Fill(MyTree::Particle& p);
	void FillLeft(MyTree::Particle& p);
	void FillRight(MyTree::Particle& p);
	void FillAll(MyTree::Particle& p, int cent9);

	void FillCut(MyTree::Particle& p);

	void Write(TFile* of);

	//hist list
	TH1F* hSameKPDG;
	TH1F* hVz;
	TH2F* hVr;
	TH1F* hCent9;
	TH1F* hPsi;
	TH1F* hRefMult;

	TH1F* hCosTheta;
	TH1F* hDecayLength;
	TH1F* hDgDCA;
	TH1F* hDCA;
	TH2F* hDedx;
	TH2F* hMass2;
	TH2F* hAllKRapPt;
	TH2F* hSameKRapPt;
	TH1F* hSameKMass;
	TH1F* hSameKPhi;
	TH2F* hSameKPipRapPt;
	TH2F* hSameKPimRapPt;
	TH1F* hDaughterPipDCA;
	TH1F* hDaughterPimDCA;
	TH2F* hKstarCheck;

	TH2F* hLeftSideSameKRapPt;
	TH1F* hLeftSideSameKMass;
	TH1F* hLeftSideSameKPhi;
	TH2F* hLeftSideSameKPipRapPt;
	TH2F* hLeftSideSameKPimRapPt;
	TH2F* hRightSideSameKRapPt;
	TH1F* hRightSideSameKMass;
	TH1F* hRightSideSameKPhi;
	TH2F* hRightSideSameKPipRapPt;
	TH2F* hRightSideSameKPimRapPt;

	TH1F* hPipNSigma;
	TH1F* hPimNSigma;
	TH1F* hPipNSigma2;
	TH1F* hPimNSigma2;
	TH3F* hSameKPtRapMass[9];

	TH1F* hDTheta;
	TH1F* hDPhi;
	TH2F* hDThetaDPhi;
	TH1F* hSameKqinv[9];
	TH1F* hMixKqinv[9];
	TH1F* hMixKqinvWeight[9][4];

	TH1F* hMix0000Qinv[9];
	TH1F* hMix0001Qinv[9];
	TH1F* hMix0011Qinv[9];
	TH1F* hMix0023Qinv[9];
	TH1F* hMix0101Qinv[9];
	TH1F* hMix0102Qinv[9];
	TH1F* hMix0103Qinv[9];
	TH1F* hMix0123Qinv[9];

	TH1F* hMix000_Qinv[9];
	TH1F* hMix001_Qinv[9];
	TH1F* hMix010_Qinv[9];
	TH1F* hMix011_Qinv[9];
	TH1F* hMix012_Qinv[9];
	TH1F* hMix00_0Qinv[9];
	TH1F* hMix00_1Qinv[9];
	TH1F* hMix01_0Qinv[9];
	TH1F* hMix01_1Qinv[9];
	TH1F* hMix01_2Qinv[9];

	TH1F* hMixBkgMass;

	//cut plots
	TH1F* hCutChi2Topo;
	TH1F* hCutChi2NDF;
	TH1F* hCutChi2PrimA;
	TH1F* hCutChi2PrimB;
	TH1F* hCutNHitsA;
	TH1F* hCutNHitsB;
	TH1F* hCutDCAA;
	TH1F* hCutDCAB;
	TH1F* hCutDecayLength;
	TH2F* hEtaPt;
};
#endif
