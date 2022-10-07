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
const int Cent9To3[NCent9] = { -1, -1, 1, 1, 1, 1, 1, 0, 0 };
const int Cent9To1[NCent9] = { -1, -1, 2, 2, 2, 2, 2, 2, 2 };
const std::string CentName[NCent] = { "0-20%", "20-60%", "0-60%" };
float weightLeft, weightRight;
TPad* thisPad = 0;
#endif

void getUrqmdCF(TFile* ifPlots, std::string prefix, TH1F* hCF[], float normLower = 0.3, float normUpper = 1.0);

void getUrqmdCF(TFile* ifPlots, std::string prefix, TH1F* hCF[], float normLower, float normUpper)
{
	//plots list {{{
	TH1F* hSameKqinv[NCent9] = { 0 };
	TH1F* hMixKqinv[NCent9] = { 0 };

	TH1F* hSameKqinv_Cent[NCent] = { 0 };
	TH1F* hMixKqinv_Cent[NCent] = { 0 };
	//}}}

	//get plots{{{
	for(int icent = 2; icent < NCent9; ++icent) {
		hSameKqinv[icent] = (TH1F*)getCopy(ifPlots, Form("hSameKqinv%s_cent%i", prefix.c_str(), icent));
		hMixKqinv[icent] = (TH1F*)getCopy(ifPlots, Form("hMixKqinv_cent%i", icent));
		hSameKqinv[icent]->Rebin(Rebin);
		hMixKqinv[icent]->Rebin(Rebin);
		//split cent9 to cent3{{{
		if(icent == 7) {
			hSameKqinv_Cent[0] = (TH1F*)hSameKqinv[icent]->Clone();
			hSameKqinv_Cent[0]->SetName(Form("h%sSameKqinv_cent30", prefix.c_str()));
			hMixKqinv_Cent[0] = (TH1F*)hMixKqinv[icent]->Clone();
			hMixKqinv_Cent[0]->SetName(Form("h%sMixKqinv_cent30", prefix.c_str()));
		} else if(icent <= 8 && icent > 7) {
			hSameKqinv_Cent[0]->Add(hSameKqinv[icent]);
			hMixKqinv_Cent[0]->Add(hMixKqinv[icent]);
		}

		if(icent == 2) {
			hSameKqinv_Cent[1] = (TH1F*)hSameKqinv[icent]->Clone();
			hSameKqinv_Cent[1]->SetName(Form("h%sSameKqinv_cent31", prefix.c_str()));
			hMixKqinv_Cent[1] = (TH1F*)hMixKqinv[icent]->Clone();
			hMixKqinv_Cent[1]->SetName(Form("h%sMixKqinv_cent31", prefix.c_str()));
		} else if(icent > 2 && icent < 7) {
			hSameKqinv_Cent[1]->Add(hSameKqinv[icent]);
			hMixKqinv_Cent[1]->Add(hMixKqinv[icent]);
		}

		if(icent == 2) {
			hSameKqinv_Cent[2] = (TH1F*)hSameKqinv[icent]->Clone();
			hSameKqinv_Cent[2]->SetName(Form("h%sSameKqinv_cent31", prefix.c_str()));
			hMixKqinv_Cent[2] = (TH1F*)hMixKqinv[icent]->Clone();
			hMixKqinv_Cent[2]->SetName(Form("h%sMixKqinv_cent31", prefix.c_str()));
		} else if(icent <= 8 && icent > 2) {
			hSameKqinv_Cent[2]->Add(hSameKqinv[icent]);
			hMixKqinv_Cent[2]->Add(hMixKqinv[icent]);
		}
		//}}}
	}
	//}}}

	//calculate CF{{{
	for(int icent = 0; icent < NCent; ++icent) {
		hCF[icent] = (TH1F*)hSameKqinv_Cent[icent]->Clone();
		TH1F* tmpMixhist = (TH1F*)hMixKqinv_Cent[icent]->Clone();
		float scale = getIntegral(hSameKqinv_Cent[icent], normLower, normUpper) / getIntegral(tmpMixhist, normLower, normUpper);
		tmpMixhist->Scale(scale);
		hCF[icent]->SetName(Form("h%sCF_cent%i", prefix.c_str(), icent));
		hCF[icent]->Divide(tmpMixhist);
		delete tmpMixhist;
	}
	//}}}

	//release memory {{{
	for(int icent = 0; icent < NCent9; ++icent) {
		delete hSameKqinv[icent];
		delete hMixKqinv[icent];
	}
	for(int icent = 0; icent < NCent; ++icent) {
		delete hSameKqinv_Cent[icent];
		delete hMixKqinv_Cent[icent];
	}
	//}}}
}
#endif
