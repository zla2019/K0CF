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
	void FillSame(float qinv, int cent9);
	void FillMix(float qinv, int cent9);
	void Write(TFile* of);

	//hist list
	TH1F* hCent9;
	TH1F* hRefMult;

	TH2F* hAllKRapPt;
	TH2F* hKRapPt;
	TH2F* hKPipRapPt;
	TH2F* hKPimRapPt;
	TH1F* hFrt;
	TH1F* hFrx;
	TH1F* hFry;
	TH1F* hFrz;
	TH2F* hFrxFry;

	TH1F* hSameKqinv[9];
	TH1F* hSameKqinvSI[9];
	TH1F* hSameKqinvQS[9];
	TH1F* hSameKqinvWoCrab[9];
	TH1F* hMixKqinv[9];

	//crab test plots
	TH2F* hQinvCorr;
	TH2F* hQdotrCorr;
	TH2F* hRCorr;
	TH2F* hLLWeight;
	TH2F* hTotWeight;
};
#endif
