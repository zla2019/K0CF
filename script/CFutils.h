#ifndef _CFUTILS_H_
#define _CFUTILS_H_
#include <iostream>
#include "utils.h"

#ifndef SETTING
#define SETTING
const int NCent9 = 9;
const int NCent = 3;
const int NCase = 3;
const int Rebin = 20;
const float PtEdge[2] = { 0.2, 1.8 };
const float RapEdge[2] = { 1.0, 0.0 };
const int Cent9To3[NCent9] = { -1, -1, 1, 1, 1, 1, 0, 0, 0 };
const int Cent9To1[NCent9] = { -1, -1, 2, 2, 2, 2, 2, 2, 2 };
const std::string CentName[NCent] = { "0-20%", "20-60%", "0-60%" };
const std::string CaseName[NCase] = { "K_{s}^{0}#tilde{K_{s}^{0}}", "#tilde{K_{s}^{0}}K_{s}^{0}", "#tilde{K_{s}^{0}}#tilde{K_{s}^{0}}" };
float weightLeft, weightRight;
TPad* thisPad = 0;
#endif

void RawCF(TFile* ifPlots, TH1F* hCFRaw[], float normLower = 0.3, float normUpper = 1.0);
void MisidCF(TFile* ifPlots, TH1F* hCFMisid[], float normLower = 0.3, float normUpper = 1.0);
float SBWeight(TFile* ifPlots, float leftL = 0.42, float leftU = 0.48, float rightL = 0.52, float rightU = 0.58);
void PureCF(TH1F* hCFRaw[], TH1F* hCFMisid[], TH1F* hCFPure[], TH1F* hPairPurity[]);
void PairPurity(TFile* ifPlots, TH1F* hPairPurity[]);

void RawCF(TFile* ifPlots, TH1F* hCFRaw[], float normLower, float normUpper)
{
	//plots list {{{
	TH1F* hSameQinvCent9[NCent9] = { 0 };
	TH1F* hMixQinvCent9[NCent9] = { 0 };

	TH1F* hSameQinv[NCent] = { 0 };
	TH1F* hMixQinv[NCent] = { 0 };
	//}}}

	//get plots {{{
	for(int icent = 0; icent < NCent9; ++icent) {
		hSameQinvCent9[icent] = (TH1F*)getCopy(ifPlots, Form("hSameKqinv_cent%d", icent));
		hMixQinvCent9[icent] = (TH1F*)getCopy(ifPlots, Form("hMixKqinv_cent%d", icent));
		hSameQinvCent9[icent]->Sumw2();
		hMixQinvCent9[icent]->Sumw2();
		hSameQinvCent9[icent]->Rebin(Rebin);
		hMixQinvCent9[icent]->Rebin(Rebin);
	}
	for(int icent = 2; icent < NCent9; ++icent) {
		if(Cent9To3[icent] < 0 || Cent9To1[icent] < 0) {
			continue;
		}
		addHist(hSameQinv[Cent9To3[icent]], hSameQinvCent9[icent], Form("hSameQinv_cent%d", Cent9To3[icent]));
		addHist(hMixQinv[Cent9To3[icent]], hMixQinvCent9[icent], Form("hMixQinv_cent%d", Cent9To3[icent]));

		addHist(hSameQinv[Cent9To1[icent]], hSameQinvCent9[icent], Form("hSameQinv_cent%d", Cent9To1[icent]));
		addHist(hMixQinv[Cent9To1[icent]], hMixQinvCent9[icent], Form("hMixQinv_cent%d", Cent9To1[icent]));
	}
	//}}}

	//calculate CF {{{
	for(int icent = 0; icent < NCent; ++icent) {
		hCFRaw[icent] = (TH1F*)hSameQinv[icent]->Clone(Form("hCFRaw_cent%d", icent));
		hCFRaw[icent]->Divide(hSameQinv[icent], hMixQinv[icent], 1 / getIntegral(hSameQinv[icent], normLower, normUpper), 1 / getIntegral(hMixQinv[icent], normLower, normUpper));
	}
	//}}}

	//delete {{{
	for(int icent = 0; icent < NCent9; ++icent) {
		delete hSameQinvCent9[icent];
		delete hMixQinvCent9[icent];
	}
	for(int icent = 0; icent < NCent; ++icent) {
		delete hSameQinv[icent];
		delete hMixQinv[icent];
	}
	//}}}
}

void MisidCF(TFile* ifPlots, TH1F* hCFMisid[], float normLower, float normUpper)
{
	//inital {{{
	weightLeft = SBWeight(ifPlots);
	weightRight = 1 - weightLeft;
	//}}}

	//plots list {{{
	TH1F* hSameLeftSideKqinvCent9[NCent9][NCase + 1] = { 0 };
	TH1F* hMixLeftSideKqinvCent9[NCent9][NCase + 1] = { 0 };
	TH1F* hMixKqinvLeftWeightCent9[NCent9][NCase + 1] = { 0 };
	TH1F* hSameRightSideKqinvCent9[NCent9][NCase + 1] = { 0 };
	TH1F* hMixRightSideKqinvCent9[NCent9][NCase + 1] = { 0 };
	TH1F* hMixKqinvRightWeightCent9[NCent9][NCase + 1] = { 0 };

	TH1F* hSameLeftSideKqinv[NCent][NCase + 1] = { 0 };
	TH1F* hMixLeftSideKqinv[NCent][NCase + 1] = { 0 };
	TH1F* hMixKqinvLeftWeight[NCent][NCase + 1] = { 0 };
	TH1F* hSameRightSideKqinv[NCent][NCase + 1] = { 0 };
	TH1F* hMixRightSideKqinv[NCent][NCase + 1] = { 0 };
	TH1F* hMixKqinvRightWeight[NCent][NCase + 1] = { 0 };

	TH1F* hLeftWeight[NCent][NCase + 1] = { 0 };
	TH1F* hRightWeight[NCent][NCase + 1] = { 0 };
	TH1F* hLeftSideCF[NCent][NCase + 1] = { 0 };
	TH1F* hRightSideCF[NCent][NCase + 1] = { 0 };

	TH1F* hLeftSideCFTot[NCent] = { 0 };
	TH1F* hRightSideCFTot[NCent] = { 0 };

	TH1F* hLeftSideCFWeight[NCent] = { 0 };
	TH1F* hRightSideCFWeight[NCent] = { 0 };
	//}}}

	//get plots {{{
	for(int icent = 0; icent < NCent9; ++icent) {
		for(int icase = 0; icase < NCase + 1; ++icase) {
			hSameLeftSideKqinvCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("hSameLeftSideKqinv_cent%d_case%d", icent, icase));
			hMixLeftSideKqinvCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("hMixLeftSideKqinv_cent%d_case%d", icent, icase));
			hMixKqinvLeftWeightCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("hMixKqinvLeftWeight_cent%d_case%d", icent, icase));
			hSameRightSideKqinvCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("hSameRightSideKqinv_cent%d_case%d", icent, icase));
			hMixRightSideKqinvCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("hMixRightSideKqinv_cent%d_case%d", icent, icase));
			hMixKqinvRightWeightCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("hMixKqinvRightWeight_cent%d_case%d", icent, icase));

			hSameLeftSideKqinvCent9[icent][icase]->Rebin(Rebin);
			hMixLeftSideKqinvCent9[icent][icase]->Rebin(Rebin);
			hMixKqinvLeftWeightCent9[icent][icase]->Rebin(Rebin);
			hSameRightSideKqinvCent9[icent][icase]->Rebin(Rebin);
			hMixRightSideKqinvCent9[icent][icase]->Rebin(Rebin);
			hMixKqinvRightWeightCent9[icent][icase]->Rebin(Rebin);
		}
	}

	for(int icase = 0; icase < NCase + 1; ++icase) {
		for(int icent = 2; icent < NCent9; ++icent) {
			if(Cent9To3[icent] < 0 || Cent9To1[icent] < 0) {
				continue;
			}
			addHist(hSameLeftSideKqinv[Cent9To3[icent]][icase], hSameLeftSideKqinvCent9[icent][icase], Form("hSameLeftSideKqinv_cent%d_case%d", Cent9To3[icent], icase));
			addHist(hMixLeftSideKqinv[Cent9To3[icent]][icase], hMixLeftSideKqinvCent9[icent][icase], Form("hMixLeftSideKqinv_cent%d_case%d", Cent9To3[icent], icase));
			addHist(hMixKqinvLeftWeight[Cent9To3[icent]][icase], hMixKqinvLeftWeightCent9[icent][icase], Form("hMixKqinvLeftWeight_cent%d_case%d", Cent9To3[icent], icase));
			addHist(hSameRightSideKqinv[Cent9To3[icent]][icase], hSameRightSideKqinvCent9[icent][icase], Form("hSameRightSideKqinv_cent%d_case%d", Cent9To3[icent], icase));
			addHist(hMixRightSideKqinv[Cent9To3[icent]][icase], hMixRightSideKqinvCent9[icent][icase], Form("hMixRightSideKqinv_cent%d_case%d", Cent9To3[icent], icase));
			addHist(hMixKqinvRightWeight[Cent9To3[icent]][icase], hMixKqinvRightWeightCent9[icent][icase], Form("hMixKqinvRightWeight_cent%d_case%d", Cent9To3[icent], icase));

			addHist(hSameLeftSideKqinv[Cent9To1[icent]][icase], hSameLeftSideKqinvCent9[icent][icase], Form("hSameLeftSideKqinv_cent%d_case%d", Cent9To1[icent], icase));
			addHist(hMixLeftSideKqinv[Cent9To1[icent]][icase], hMixLeftSideKqinvCent9[icent][icase], Form("hMixLeftSideKqinv_cent%d_case%d", Cent9To1[icent], icase));
			addHist(hMixKqinvLeftWeight[Cent9To1[icent]][icase], hMixKqinvLeftWeightCent9[icent][icase], Form("hMixKqinvLeftWeight_cent%d_case%d", Cent9To1[icent], icase));
			addHist(hSameRightSideKqinv[Cent9To1[icent]][icase], hSameRightSideKqinvCent9[icent][icase], Form("hSameRightSideKqinv_cent%d_case%d", Cent9To1[icent], icase));
			addHist(hMixRightSideKqinv[Cent9To1[icent]][icase], hMixRightSideKqinvCent9[icent][icase], Form("hMixRightSideKqinv_cent%d_case%d", Cent9To1[icent], icase));
			addHist(hMixKqinvRightWeight[Cent9To1[icent]][icase], hMixKqinvRightWeightCent9[icent][icase], Form("hMixKqinvRightWeight_cent%d_case%d", Cent9To1[icent], icase));
		}
	}
	//}}}

	//get weight vs. qinv{{{
	for(int icent = 0; icent < NCent; ++icent) {
		for(int icase = 0; icase < NCase + 1; ++icase) {
			addHist(hLeftWeight[icent][icase], hMixKqinvLeftWeight[icent][icase], Form("hLeftWeight_cent%d_case%d", icent, icase));
			addHist(hRightWeight[icent][icase], hMixKqinvRightWeight[icent][icase], Form("hRightWeight_cent%d_case%d", icent, icase));
			hLeftWeight[icent][icase]->Divide(hLeftWeight[icent][icase], hMixLeftSideKqinv[icent][icase]);
			hRightWeight[icent][icase]->Divide(hRightWeight[icent][icase], hMixRightSideKqinv[icent][icase]);
		}
	}
	//}}}

	//calculate Side Band CF{{{
	for(int icent = 0; icent < NCent; ++icent) {
		for(int icase = 0; icase < NCase; ++icase) {
			addHist(hLeftSideCF[icent][icase], hSameLeftSideKqinv[icent][icase], Form("hLeftSideCF_cent%d_case%d", icent, icase));
			addHist(hRightSideCF[icent][icase], hSameRightSideKqinv[icent][icase], Form("hRightSideCF_cent%d_case%d", icent, icase));
			hLeftSideCF[icent][icase]->Divide(hSameLeftSideKqinv[icent][icase], hMixLeftSideKqinv[icent][icase], 1 / getIntegral(hSameLeftSideKqinv[icent][icase], normLower, normUpper), 1 / getIntegral(hMixLeftSideKqinv[icent][icase], normLower, normUpper));
			hRightSideCF[icent][icase]->Divide(hSameRightSideKqinv[icent][icase], hMixRightSideKqinv[icent][icase], 1 / getIntegral(hSameRightSideKqinv[icent][icase], normLower, normUpper), 1 / getIntegral(hMixRightSideKqinv[icent][icase], normLower, normUpper));
		}
	}

	for(int icent = 0; icent < NCent; ++icent) {
		for(int icase = 0; icase < 1; ++icase) {
			addHist(hLeftSideCFTot[icent], hLeftSideCF[icent][icase], Form("hLeftSideCFTot_cent%d", icent), hLeftWeight[icent][icase], 1);
			addHist(hRightSideCFTot[icent], hRightSideCF[icent][icase], Form("hRightSideCFTot_cent%d", icent), hRightWeight[icent][icase], 1);
		}
	}

	for(int icent = 0; icent < NCent; ++icent) {
		hLeftSideCFWeight[icent] = (TH1F*)hLeftSideCFTot[icent]->Clone();
		hLeftSideCFWeight[icent]->SetName(Form("hLeftSideCFWeight_cent%d", icent));
		scale(hLeftSideCFWeight[icent], weightLeft, 1);
		hRightSideCFWeight[icent] = (TH1F*)hRightSideCFTot[icent]->Clone();
		hRightSideCFWeight[icent]->SetName(Form("hRightSideCFWeight_cent%d", icent));
		scale(hRightSideCFWeight[icent], (1 - weightLeft), 1);

		addHist(hCFMisid[icent], hLeftSideCFWeight[icent], Form("hCFMisid_cent%d", icent), 1, 1);
		addHist(hCFMisid[icent], hRightSideCFWeight[icent], Form("hCFMisid_cent%d", icent), 1, 1);
	}
	//}}}

	//delete {{{
	for(int icent = 0; icent < NCent9; ++icent) {
		for(int icase = 0; icase < NCase + 1; ++icase) {
			delete hSameLeftSideKqinvCent9[icent][icase];
			delete hMixLeftSideKqinvCent9[icent][icase];
			delete hMixKqinvLeftWeightCent9[icent][icase];
			delete hSameRightSideKqinvCent9[icent][icase];
			delete hMixRightSideKqinvCent9[icent][icase];
			delete hMixKqinvRightWeightCent9[icent][icase];
		}
	}

	for(int icent = 0; icent < NCent; ++icent) {
		delete hLeftSideCFTot[icent];
		delete hRightSideCFTot[icent];
		delete hLeftSideCFWeight[icent];
		delete hRightSideCFWeight[icent];
		for(int icase = 0; icase < NCase + 1; ++icase) {
			delete hSameLeftSideKqinv[icent][icase];
			delete hMixLeftSideKqinv[icent][icase];
			delete hMixKqinvLeftWeight[icent][icase];
			delete hSameRightSideKqinv[icent][icase];
			delete hMixRightSideKqinv[icent][icase];
			delete hMixKqinvRightWeight[icent][icase];

			delete hLeftWeight[icent][icase];
			delete hRightWeight[icent][icase];
			delete hLeftSideCF[icent][icase];
			delete hRightSideCF[icent][icase];
		}
	}
	//}}}
}

void PureCF(TH1F* hCFRaw[], TH1F* hCFMisid[], TH1F* hCFPure[], TH1F* hPairPurity[])
{
	//extract hCFPure {{{
	for(int icent = 0; icent < NCent; ++icent) {
		addHist(hCFPure[icent], hCFRaw[icent], Form("hCFGenuine_cent%d, icent", icent), 1, 1);
		addHist(hCFPure[icent], hCFMisid[icent], Form("hCFGenuine_cent%d, icent", icent), -1, 1);
		divide(hCFPure[icent], hPairPurity[icent], 1);
	}
	//}}}
}

float SBWeight(TFile* ifPlots, float leftL, float leftU, float rightL, float rightU)
{
	//plots list {{{
	TH3F* hMassRapPtCent9[NCent9] = { 0 };
	TH3F* hMassRapPt[NCent] = { 0 };
	TH1F* hMass[NCent] = { 0 };
	//}}}

	//get plots{{{
	for(int icent = 0; icent < NCent9; ++icent) {
		hMassRapPtCent9[icent] = (TH3F*)getCopy(ifPlots, Form("hMassCascadeRapidityvsPt_cent%d", icent));
		if(Cent9To3[icent] == -1 || Cent9To1[icent] == -1) {
			continue;
		}
		addHist(hMassRapPt[Cent9To3[icent]], hMassRapPtCent9[icent], Form("hMassRapPt_cent%d", Cent9To3[icent]));
		addHist(hMassRapPt[Cent9To1[icent]], hMassRapPtCent9[icent], Form("hMassRapPt_cent%d", Cent9To1[icent]));
	}
	for(int icent = 0; icent < NCent; ++icent) {
		hMass[icent] = projectionX(hMassRapPt[icent], RapEdge[0], RapEdge[1], PtEdge[0], PtEdge[1], Form("hMass_cent%d", icent));
	}
	//}}}

	//calculate weight {{{
	for(int icent = 0; icent < NCent; ++icent) {
		weightLeft = getIntegral(hMass[icent], leftL, leftU);
		weightRight = getIntegral(hMass[icent], rightL, rightU);
		float tot = weightLeft + weightRight;
		weightLeft = weightLeft / tot;
		weightRight = 1 - weightLeft;
	}
	//}}}

	//delete {{{
	for(int icent = 0; icent < NCent9; ++icent) {
		delete hMassRapPtCent9[icent];
	}
	for(int icent = 0; icent < NCent; ++icent) {
		delete hMassRapPt[icent];
		delete hMass[icent];
	}
	//}}}
	return weightLeft;
}

void PairPurity(TFile* ifPlots, TH1F* hPairPurity[])
{
	//plots list {{{
	TH1F* hMixKqinvWeightCent9[NCent9] = { 0 };
	TH1F* hMixKqinvCent9[NCent9] = { 0 };
	TH1F* hMixKqinvWeight[NCent] = { 0 };
	TH1F* hMixKqinv[NCent] = { 0 };
	//}}}

	//get plots {{{
	for(int icent = 0; icent < NCent9; ++icent) {
		hMixKqinvWeightCent9[icent] = (TH1F*)getCopy(ifPlots, Form("hMixKqinvLeftWeight_cent%d_case%d", icent, 3));
		hMixKqinvCent9[icent] = (TH1F*)getCopy(ifPlots, Form("hMixLeftSideKqinv_cent%d_case%d", icent, 3));
		hMixKqinvWeightCent9[icent]->Rebin(Rebin);
		hMixKqinvCent9[icent]->Rebin(Rebin);
		if(Cent9To3[icent] < 0 || Cent9To1[icent] < 0) {
			continue;
		}
		addHist(hMixKqinvWeight[Cent9To3[icent]], hMixKqinvWeightCent9[icent], Form("hMixQinvWeight_cent%d", Cent9To3[icent]));
		addHist(hMixKqinv[Cent9To3[icent]], hMixKqinvCent9[icent], Form("hMixQinv_cent%d", Cent9To3[icent]));

		addHist(hMixKqinvWeight[Cent9To1[icent]], hMixKqinvWeightCent9[icent], Form("hMixQinvWeight_cent%d", Cent9To1[icent]));
		addHist(hMixKqinv[Cent9To1[icent]], hMixKqinvCent9[icent], Form("hMixQinv_cent%d", Cent9To1[icent]));
	}
	//}}}

	//calculate the pair purity {{{
	for(int icent = 0; icent < NCent; ++icent) {
		addHist(hPairPurity[icent], hMixKqinvWeight[icent], Form("hPairPurity_cent%d", icent));
		hPairPurity[icent]->Divide(hMixKqinvWeight[icent], hMixKqinv[icent]);
	}
	//}}}

	//delete {{{
	for(int icent = 0; icent < NCent9; ++icent) {
		delete hMixKqinvWeightCent9[icent];
		delete hMixKqinvCent9[icent];
	}
	for(int icent = 0; icent < NCent; ++icent) {
		delete hMixKqinvWeight[icent];
		delete hMixKqinv[icent];
	}
	//}}}
}
#endif
