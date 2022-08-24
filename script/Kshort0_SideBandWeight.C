#include <iostream>
#include "utils.h"
#include "func.h"

void Kshort0_SideBandWeight()
{
	//initial{{{
	//TFile* ifPlots = TFile::Open("/home/zla/CF/data/check/3p0WensongAcc_2sigma.root");
	TFile* ifPlots = TFile::Open("/home/zla/CF/analysis/3p5WensongAcc_2sigma.root");
	//TFile* ifPlots = TFile::Open("/home/zla/CF/analysis_3p2/3p2WensongAcc_2sigma.root");
	//TFile* ifPlots = TFile::Open("~/CF/analysis_3p2/3p2OldAcc.root");
	//TFile* ifPlots = TFile::Open("~/CF/analysis/3p5OldAcc.root");
	const int NCent9 = 9;
	const int NCent = 3;
	const int NRap = 3;
	const float RapEdge[NRap + 1] = { -1.0, -0.5, 0.0, 0.5 };
	const std::string CentName[NCent] = { "0-10%", "10-60%", "0-60%" };
	const int Cent9To3[NCent9] = { -1, -1, 1, 1, 1, 1, 1, 0, 0 };
	const int Cent9To1[NCent9] = { -1, -1, 2, 2, 2, 2, 2, 2, 2 };
	const float SideBandLower = 0.42, SideBandUpper = 0.48;
	const float SideBand2Lower = 0.52, SideBand2Upper = 0.58;
	float weightLeft, weightRight;
	//}}}

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
		hMass[icent] = projectionX(hMassRapPt[icent], RapEdge[1], RapEdge[2], 0.2, 1.2, Form("hMass_cent%d", icent));
		//TH1F* hTmp = projectionX(hMassRapPt[icent], RapEdge[1], 0.0, 0.2, 0.3, "hTmp");
		//hMass[icent]->Add(hTmp);
	}
	//}}}

	//calculate weight {{{
	for(int icent = 0; icent < NCent; ++icent) {
		weightLeft = getIntegral(hMass[icent], SideBandLower, SideBandUpper);
		weightRight = getIntegral(hMass[icent], SideBand2Lower, SideBand2Upper);
		float tot = weightLeft + weightRight;
		weightLeft = weightLeft / tot;
		weightRight = 1 - weightLeft;
	}
	//}}}

	//fitting{{{
	TF1* fitFuncTmp = new TF1("fitFuncTmp", fitFuncGaus, 0.45, 0.55, /*7*/6);
	fitFuncTmp->SetParLimits(1, 0.497, 0.499);
	fitFuncTmp->SetParLimits(0, 0., 100000);
	fitFuncTmp->SetParLimits(2, 0.002, 0.012);
	fitFuncTmp->SetParLimits(3, 0., 10000);
	hMass[2]->Fit(fitFuncTmp, "IMNR", "", 0.45, 0.55);
	const double* parErr = fitFuncTmp->GetParErrors();
	TF1* fDG = new TF1("fDG", fitFuncDoubleGaus, 0.45, 0.55, 8);
	fDG->SetParLimits(1, 0.497, 0.499);
	fDG->SetParLimits(2, fitFuncTmp->GetParameter(2) - 100.*parErr[2] > 0 ? fitFuncTmp->GetParameter(2) - 100.*parErr[2] : 0, fitFuncTmp->GetParameter(2) - 60.*parErr[2] > 0 ? fitFuncTmp->GetParameter(2) - 60.*parErr[2] : fitFuncTmp->GetParameter(2));
	fDG->SetParLimits(3, 0.33*fitFuncTmp->GetParameter(0), 0.63*fitFuncTmp->GetParameter(0));
	fDG->SetParLimits(4, fitFuncTmp->GetParameter(2)*1.6, fitFuncTmp->GetParameter(2)*1.9);
	fDG->SetParLimits(0, fitFuncTmp->GetParameter(0) * 0.72, fitFuncTmp->GetParameter(0) * 0.8);
	hMass[2]->Fit(fDG, "MR", "", 0.45, 0.55);
	//}}}

	//plotting {{{
	TCanvas* caMass = new TCanvas("caMass", "caMass", 1080, 1080);
	TPad* thisPad = (TPad*)caMass->cd();
	thisPad->SetGrid();
	//caMass->Divide(3, 1);
	hMass[2]->SetTitle("");
	hMass[2]->SetMarkerStyle(24);
	hMass[2]->SetStats(0);
	setXYTitle(hMass[2], "M_{#pi^{+}#pi^{-}}", "Counts");
	hMass[2]->Draw("ep");
	drawXRegion(0.42, 0.48, thisPad, kBlue, 2, 4);
	drawXRegion(0.52, 0.58, thisPad, kBlue, 2, 4);
	drawText(0.7, 0.7, 0.035, Form("#omega_{left}: %.2f", weightLeft), Form("#omega_{right}: %.2f", weightRight));
	drawText(0.15, 0.7, 0.035, "K_{s}^{0} candidates:", "Centrality: 0~60%", "y: -0.5~0.0", "p_{T}:0.2~1.2");
	//}}}
}

float Kshort0_SideBandWeight(TFile* ifPlots, float SideBandLower, float SideBandUpper, float SideBand2Lower, float SideBand2Upper)
{
	//initial{{{
	//TFile* ifPlots = TFile::Open("/home/zla/CF/data/check/production_202206012249.root");
	const int NCent9 = 9;
	const int NCent = 3;
	const int NRap = 3;
	const float RapEdge[NRap + 1] = { -1.0, -0.5, 0.0, 0.5 };
	const std::string CentName[NCent] = { "0-10%", "10-60%", "0-60%" };
	const int Cent9To3[NCent9] = { -1, -1, 1, 1, 1, 1, 1, 0, 0 };
	const int Cent9To1[NCent9] = { -1, -1, 2, 2, 2, 2, 2, 2, 2 };
	//const float SideBandLower = 0.43, SideBandUpper = 0.47;
	//const float SideBand2Lower = 0.53, SideBand2Upper = 0.57;
	float weightLeft, weightRight;
	//}}}

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
		hMass[icent] = projectionX(hMassRapPt[icent], RapEdge[1], RapEdge[2], 0.2, 1.5, Form("hMass_cent%d", icent));
		//TH1F* hTmp = projectionX(hMassRapPt[icent], RapEdge[1], 0.0, 0.2, 0.3, "hTmp");
		//hMass[icent]->Add(hTmp);
	}
	//}}}

	//calculate weight {{{
	for(int icent = 0; icent < NCent; ++icent) {
		weightLeft = getIntegral(hMass[icent], SideBandLower, SideBandUpper);
		weightRight = getIntegral(hMass[icent], SideBand2Lower, SideBand2Upper);
		float tot = weightLeft + weightRight;
		weightLeft = weightLeft / tot;
		weightRight = 1 - weightLeft;
	}
	//}}}
	return weightLeft;
}
