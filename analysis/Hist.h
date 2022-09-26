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

	TH1F* hSameKqinv[9];
	TH1F* hMixKqinv[9];
};
#endif
