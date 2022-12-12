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
void RotCF(TFile* ifPlots, TH1F* hCFMisid[], TH1F* hRotCF[][NCase + 1], float normLower = 0.3, float normUpper = 1.0);
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
	const std::string SameHistName[4] = { "hSameKqinv", "hSRQinv", "hRSQinv", "hRRQinv" };
	const std::string MixHistName[4] = { "hMixKqinv", "hMixSRQinv", "hMixRSQinv", "hMixRRQinv" };
	//}}}

	//plots list {{{
	TH1F* hRotSameKqinvCent9[NCent9][NCase + 1] = { 0 };
	TH1F* hRotMixKqinvCent9[NCent9][NCase + 1] = { 0 };
	TH1F* hRotMixKqinvWeightCent9[NCent9][NCase + 1] = { 0 };

	TH1F* hRotSameKqinv[NCent][NCase + 1] = { 0 };
	TH1F* hRotMixKqinv[NCent][NCase + 1] = { 0 };
	TH1F* hRotMixKqinvWeight[NCent][NCase + 1] = { 0 };

	TH1F* hRotWeight[NCent][NCase + 1] = { 0 };
	TH1F* hRotCF[NCent][NCase + 1] = { 0 };

	TH1F* hRotCFTot[NCent] = { 0 };

	TH1F* hRotCFWeight[NCent] = { 0 };
	//}}}

	//get plots {{{
	for(int icent = 0; icent < NCent9; ++icent) {
		for(int icase = 0; icase < NCase + 1; ++icase) {
			hRotSameKqinvCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("%s_cent%d", SameHistName[icase].c_str(), icent));
			hRotMixKqinvCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("%s_cent%d", MixHistName[icase].c_str(), icent));
			hRotMixKqinvWeightCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("hMixKqinvWeight_cent%d_case%d", icent, icase));

			hRotSameKqinvCent9[icent][icase]->Rebin(Rebin);
			hRotMixKqinvCent9[icent][icase]->Rebin(Rebin);
			hRotMixKqinvWeightCent9[icent][icase]->Rebin(Rebin);
		}
	}

	for(int icase = 0; icase < NCase + 1; ++icase) {
		for(int icent = 2; icent < NCent9; ++icent) {
			if(Cent9To3[icent] < 0 || Cent9To1[icent] < 0) {
				continue;
			}
			addHist(hRotSameKqinv[Cent9To3[icent]][icase], hRotSameKqinvCent9[icent][icase], Form("hRotSameKqinv_cent%d_case%d", Cent9To3[icent], icase));
			addHist(hRotMixKqinv[Cent9To3[icent]][icase], hRotMixKqinvCent9[icent][icase], Form("hRotMixKqinv_cent%d_case%d", Cent9To3[icent], icase));
			addHist(hRotMixKqinvWeight[Cent9To3[icent]][icase], hRotMixKqinvWeightCent9[icent][icase], Form("hRotMixKqinvWeight_cent%d_case%d", Cent9To3[icent], icase));

			addHist(hRotSameKqinv[Cent9To1[icent]][icase], hRotSameKqinvCent9[icent][icase], Form("hRotSameKqinv_cent%d_case%d", Cent9To1[icent], icase));
			addHist(hRotMixKqinv[Cent9To1[icent]][icase], hRotMixKqinvCent9[icent][icase], Form("hRotMixKqinv_cent%d_case%d", Cent9To1[icent], icase));
			addHist(hRotMixKqinvWeight[Cent9To1[icent]][icase], hRotMixKqinvWeightCent9[icent][icase], Form("hRotMixKqinvWeight_cent%d_case%d", Cent9To1[icent], icase));
		}
	}
	//}}}

	//get weight vs. qinv{{{
	for(int icent = 0; icent < NCent; ++icent) {
		for(int icase = 0; icase < NCase + 1; ++icase) {
			addHist(hRotWeight[icent][icase], hRotMixKqinvWeight[icent][icase], Form("hWeight_cent%d_case%d", icent, icase));
			hRotWeight[icent][icase]->Divide(hRotWeight[icent][icase], hRotMixKqinv[icent][0]);
		}
	}
	//}}}

	//calculate Side Band CF{{{
	for(int icent = 0; icent < NCent; ++icent) {
		for(int icase = 0; icase < NCase + 1; ++icase) {
			addHist(hRotCF[icent][icase], hRotSameKqinv[icent][icase], Form("hRotCF_cent%d_case%d", icent, icase));
			hRotCF[icent][icase]->Divide(hRotSameKqinv[icent][icase], hRotMixKqinv[icent][icase], 1 / getIntegral(hRotSameKqinv[icent][icase], normLower, normUpper), 1 / getIntegral(hRotMixKqinv[icent][icase], normLower, normUpper));
		}
	}

	for(int icent = 0; icent < NCent; ++icent) {
		for(int icase = 1; icase < NCase + 1; ++icase) {
			addHist(hRotCFTot[icent], hRotCF[icent][icase], Form("hRotCFTot_cent%d", icent), hRotWeight[icent][icase], 1);
		}
	}

	for(int icent = 0; icent < NCent; ++icent) {
		hRotCFWeight[icent] = (TH1F*)hRotCFTot[icent]->Clone();
		hRotCFWeight[icent]->SetName(Form("hRotCFWeight_cent%d", icent));
		scale(hRotCFWeight[icent], 1, 1);

		addHist(hCFMisid[icent], hRotCFWeight[icent], Form("hCFMisid_cent%d", icent), 1, 1);
	}
	//}}}

	//delete {{{
	for(int icent = 0; icent < NCent9; ++icent) {
		for(int icase = 0; icase < NCase + 1; ++icase) {
			delete hRotSameKqinvCent9[icent][icase];
			delete hRotMixKqinvCent9[icent][icase];
			delete hRotMixKqinvWeightCent9[icent][icase];
		}
	}

	for(int icent = 0; icent < NCent; ++icent) {
		delete hRotCFTot[icent];
		delete hRotCFWeight[icent];
		for(int icase = 0; icase < NCase + 1; ++icase) {
			delete hRotSameKqinv[icent][icase];
			delete hRotMixKqinv[icent][icase];
			delete hRotMixKqinvWeight[icent][icase];

			delete hRotWeight[icent][icase];
			delete hRotCF[icent][icase];
		}
	}
	//}}}
}

void RotCF(TFile* ifPlots, TH1F* hCFMisid[], TH1F* hRotCF[][NCase + 1], float normLower, float normUpper)
{
	//inital {{{
	const std::string SameHistName[4] = { "hSameKqinv", "hSRQinv", "hRSQinv", "hRRQinv" };
	const std::string MixHistName[4] = { "hMixKqinv", "hMixSRQinv", "hMixRSQinv", "hMixRRQinv" };
	//}}}

	//plots list {{{
	TH1F* hRotSameKqinvCent9[NCent9][NCase + 1] = { 0 };
	TH1F* hRotMixKqinvCent9[NCent9][NCase + 1] = { 0 };
	TH1F* hRotMixKqinvWeightCent9[NCent9][NCase + 1] = { 0 };

	TH1F* hRotSameKqinv[NCent][NCase + 1] = { 0 };
	TH1F* hRotMixKqinv[NCent][NCase + 1] = { 0 };
	TH1F* hRotMixKqinvWeight[NCent][NCase + 1] = { 0 };

	TH1F* hRotWeight[NCent][NCase + 1] = { 0 };

	TH1F* hRotCFTot[NCent] = { 0 };

	TH1F* hRotCFWeight[NCent] = { 0 };
	//}}}

	//get plots {{{
	for(int icent = 0; icent < NCent9; ++icent) {
		for(int icase = 0; icase < NCase + 1; ++icase) {
			hRotSameKqinvCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("%s_cent%d", SameHistName[icase].c_str(), icent));
			hRotMixKqinvCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("%s_cent%d", MixHistName[icase].c_str(), icent));
			hRotMixKqinvWeightCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("hMixKqinvWeight_cent%d_case%d", icent, icase));

			hRotSameKqinvCent9[icent][icase]->Rebin(Rebin);
			hRotMixKqinvCent9[icent][icase]->Rebin(Rebin);
			hRotMixKqinvWeightCent9[icent][icase]->Rebin(Rebin);
		}
	}

	for(int icase = 0; icase < NCase + 1; ++icase) {
		for(int icent = 2; icent < NCent9; ++icent) {
			if(Cent9To3[icent] < 0 || Cent9To1[icent] < 0) {
				continue;
			}
			addHist(hRotSameKqinv[Cent9To3[icent]][icase], hRotSameKqinvCent9[icent][icase], Form("hRotSameKqinv_cent%d_case%d", Cent9To3[icent], icase));
			addHist(hRotMixKqinv[Cent9To3[icent]][icase], hRotMixKqinvCent9[icent][icase], Form("hRotMixKqinv_cent%d_case%d", Cent9To3[icent], icase));
			addHist(hRotMixKqinvWeight[Cent9To3[icent]][icase], hRotMixKqinvWeightCent9[icent][icase], Form("hRotMixKqinvWeight_cent%d_case%d", Cent9To3[icent], icase));

			addHist(hRotSameKqinv[Cent9To1[icent]][icase], hRotSameKqinvCent9[icent][icase], Form("hRotSameKqinv_cent%d_case%d", Cent9To1[icent], icase));
			addHist(hRotMixKqinv[Cent9To1[icent]][icase], hRotMixKqinvCent9[icent][icase], Form("hRotMixKqinv_cent%d_case%d", Cent9To1[icent], icase));
			addHist(hRotMixKqinvWeight[Cent9To1[icent]][icase], hRotMixKqinvWeightCent9[icent][icase], Form("hRotMixKqinvWeight_cent%d_case%d", Cent9To1[icent], icase));
		}
	}
	//}}}

	//get weight vs. qinv{{{
	for(int icent = 0; icent < NCent; ++icent) {
		for(int icase = 0; icase < NCase + 1; ++icase) {
			addHist(hRotWeight[icent][icase], hRotMixKqinvWeight[icent][icase], Form("hWeight_cent%d_case%d", icent, icase));
			hRotWeight[icent][icase]->Divide(hRotWeight[icent][icase], hRotMixKqinv[icent][0]);
		}
	}
	//}}}

	//calculate Side Band CF{{{
	for(int icent = 0; icent < NCent; ++icent) {
		for(int icase = 0; icase < NCase + 1; ++icase) {
			addHist(hRotCF[icent][icase], hRotSameKqinv[icent][icase], Form("hRotCF_cent%d_case%d", icent, icase));
			hRotCF[icent][icase]->Divide(hRotSameKqinv[icent][icase], hRotMixKqinv[icent][icase], 1 / getIntegral(hRotSameKqinv[icent][icase], normLower, normUpper), 1 / getIntegral(hRotMixKqinv[icent][icase], normLower, normUpper));
		}
	}

	for(int icent = 0; icent < NCent; ++icent) {
		for(int icase = 1; icase < NCase + 1; ++icase) {
			addHist(hRotCFTot[icent], hRotCF[icent][icase], Form("hRotCFTot_cent%d", icent), hRotWeight[icent][icase], 1);
		}
	}

	for(int icent = 0; icent < NCent; ++icent) {
		hRotCFWeight[icent] = (TH1F*)hRotCFTot[icent]->Clone();
		hRotCFWeight[icent]->SetName(Form("hRotCFWeight_cent%d", icent));
		scale(hRotCFWeight[icent], 1, 1);

		addHist(hCFMisid[icent], hRotCFWeight[icent], Form("hCFMisid_cent%d", icent), 1, 1);
	}
	//}}}

	//delete {{{
	for(int icent = 0; icent < NCent9; ++icent) {
		for(int icase = 0; icase < NCase + 1; ++icase) {
			delete hRotSameKqinvCent9[icent][icase];
			delete hRotMixKqinvCent9[icent][icase];
			delete hRotMixKqinvWeightCent9[icent][icase];
		}
	}

	for(int icent = 0; icent < NCent; ++icent) {
		delete hRotCFTot[icent];
		delete hRotCFWeight[icent];
		for(int icase = 0; icase < NCase + 1; ++icase) {
			delete hRotSameKqinv[icent][icase];
			delete hRotMixKqinv[icent][icase];
			delete hRotMixKqinvWeight[icent][icase];

			delete hRotWeight[icent][icase];
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
		hMixKqinvWeightCent9[icent] = (TH1F*)getCopy(ifPlots, Form("hMixKqinvWeight_cent%d_case%d", icent, 0));
		hMixKqinvCent9[icent] = (TH1F*)getCopy(ifPlots, Form("hMixKqinv_cent%d", icent));
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
