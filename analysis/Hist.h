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

	void Fill(float qinv, int cent9, int isSideBand1, int isSideBand2);
	void FillMix(float qinv, int cent9, int isSideBand1, int isSideBand2, float sideBandWeight[]);

	void Write(TFile* of);

	//hist list
	TH1F* hSameKPDG;
	TH1F* hVz;
	TH2F* hVr;
	TH1F* hCent9;

	TH1F* hCosTheta;
	TH1F* hDecayLength;
	TH1F* hDgDCA;
	TH1F* hDCA;
	TH2F* hDedx;
	TH2F* hMass2;
	TH2F* hSameKRapPt;
	TH1F* hSameKMass;
	TH1F* hSameKPhi;
	TH2F* hSameKPipRapPt;
	TH2F* hSameKPimRapPt;
	TH1F* hDaughterPipDCA;
	TH1F* hDaughterPimDCA;
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
	TH1F* hMixKqinvLeftWeight[9][4];
	TH1F* hMixKqinvRightWeight[9][4];
	TH1F* hSameKqlong[9];
	TH1F* hMixKqlong[9];
	TH1F* hSameKqout[9];
	TH1F* hMixKqout[9];
	TH1F* hSameKqside[9];
	TH1F* hMixKqside[9];
	TH1F* hSameLeftSideKqinv[9][4];
	TH1F* hMixLeftSideKqinv[9][4];
	TH1F* hSameRightSideKqinv[9][4];
	TH1F* hMixRightSideKqinv[9][4];
};
#endif
