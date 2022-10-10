#include <iostream>
#include "CFLLFunc.C"
#define SETTING
const int NCent9 = 9;
const int NCent = 3;
const int NCase = 3;
const int Rebin = 5;
const int NEnergy = 4;
const float Energy[NEnergy] = { 3.0, 3.2, 3.5, 3.9 };
const float EnergyWidth[NEnergy] = { 0, 0, 0, 0 };
const float PtEdge[2] = { 0.2, 1.8 };
const float RapEdge[2] = { 1.0, 0.0 };
const int Cent9To3[NCent9] = { -1, -1, 1, 1, 1, 1, 1, 0, 0 };
const int Cent9To1[NCent9] = { -1, -1, 2, 2, 2, 2, 2, 2, 2 };
const std::string CentName[NCent] = { "0-20%", "20-60%", "0-60%" };
float weightLeft, weightRight;
TPad* thisPad = 0;
#include "CFutils.h"

void EnergyDependence()
{
	//inital {{{
	TFile* ifPlots[NEnergy] = { 0 };
	ifPlots[0] = TFile::Open("/home/zla/CF/data/urqmd/3p0/result_202210040411.root");	//3.0
	ifPlots[1] = TFile::Open("/home/zla/CF/data/urqmd/3p2/result_202210060518.root");	//3.2
	ifPlots[2] = TFile::Open("/home/zla/CF/data/urqmd/3p5/result_202210050157.root");	//3.5
	ifPlots[3] = TFile::Open("/home/zla/CF/data/urqmd/3p9/result_202210040728.root");	//3.9
	//ifPlots[4] = TFile::Open("/home/zla/CF/data/tmp/result_202209292213.root");		//7.7
	//}}}

	//plots list {{{
	TH1F* hCF[NEnergy][NCent] = { 0 };
	TH1F* hCFQS[NEnergy][NCent] = { 0 };
	TH1F* hCFSI[NEnergy][NCent] = { 0 };
	TH1F* hCFWoCrab[NEnergy][NCent] = { 0 };

	TF1* fLL[NEnergy][NCent] = { 0 };
	TF1* fGaus[NEnergy][NCent] = { 0 };

	TGraphErrors* grSS[NCent] = { 0 };
	TGraphErrors* grLambda[NCent] = { 0 };
	TGraphErrors* grGausSS[NCent] = { 0 };
	TGraphErrors* grGausLambda[NCent] = { 0 };
	//}}}

	//processing {{{
	for(int ie = 0; ie < NEnergy; ++ie) {
		getUrqmdCF(ifPlots[ie], "", hCF[ie]);
		getUrqmdCF(ifPlots[ie], "QS", hCFQS[ie]);
		getUrqmdCF(ifPlots[ie], "SI", hCFSI[ie]);
		getUrqmdCF(ifPlots[ie], "WoCrab", hCFWoCrab[ie]);
	}
	//}}}

	//fitting {{{
	for(int ie = 0; ie < NEnergy; ++ie) {
		for(int icent = 0; icent < NCent; ++icent) {
			fLL[ie][icent] = new TF1(Form("fLL_e%d_cent%d", ie, icent), CFLL, 0, 0.4, 2);
			fGaus[ie][icent] = new TF1(Form("fGaus_e%d_cent%d", ie, icent), "1 + [0] * exp((-1 * [1]^2 * (x)^2) / 0.038937929230)", 0, 1);

			fLL[ie][icent]->SetParameters(0.7, 3.0);
			fLL[ie][icent]->SetParLimits(0, 0.2, 1.0);
			fLL[ie][icent]->SetParLimits(1, 1, 7);
			fGaus[ie][icent]->SetParameters(0.7, 3.0);
			fGaus[ie][icent]->SetParLimits(0, 0.2, 1.0);
			fGaus[ie][icent]->SetParLimits(1, 1, 7);

			hCF[ie][icent]->Fit(fLL[ie][icent]);
			hCFQS[ie][icent]->Fit(fGaus[ie][icent]);
		}
	}
	//}}}

	//plotting {{{
	float lambda[NCent][NEnergy] = { 0 };	//the index is different with histo, the reason is we are going to draw an energy dependence of parameters, the energy must be last index
	float ss[NCent][NEnergy] = { 0 };
	float lambdaErr[NCent][NEnergy] = { 0 };
	float ssErr[NCent][NEnergy] = { 0 };

	float lambdaGaus[NCent][NEnergy] = { 0 };
	float ssGaus[NCent][NEnergy] = { 0 };
	float lambdaGausErr[NCent][NEnergy] = { 0 };
	float ssGausErr[NCent][NEnergy] = { 0 };

	for(int ie = 0; ie < NEnergy; ++ie) {
		for(int icent = 0; icent < NCent; ++icent) {
			lambda[icent][ie] = fLL[ie][icent]->GetParameter(0);
			ss[icent][ie] = fLL[ie][icent]->GetParameter(1);
			lambdaErr[icent][ie] = fLL[ie][icent]->GetParError(0);
			ssErr[icent][ie] = fLL[ie][icent]->GetParError(1);

			lambdaGaus[icent][ie] = fGaus[ie][icent]->GetParameter(0);
			ssGaus[icent][ie] = fGaus[ie][icent]->GetParameter(1);
			lambdaGausErr[icent][ie] = fGaus[ie][icent]->GetParError(0);
			ssGausErr[icent][ie] = fGaus[ie][icent]->GetParError(1);
		}
	}
	for(int icent = 0; icent < NCent; ++icent) {
		grSS[icent] = new TGraphErrors(NEnergy, Energy, ss[icent], EnergyWidth, ssErr[icent]);
		grLambda[icent] = new TGraphErrors(NEnergy, Energy, lambda[icent], EnergyWidth, lambdaErr[icent]);

		grGausSS[icent] = new TGraphErrors(NEnergy, Energy, ssGaus[icent], EnergyWidth, ssGausErr[icent]);
		grGausLambda[icent] = new TGraphErrors(NEnergy, Energy, lambdaGaus[icent], EnergyWidth, lambdaGausErr[icent]);
	}

	TCanvas* caSS = new TCanvas("caSS", "caSS", 720, 720);
	TH2F* hCanvas = new TH2F("hCanvasSS", "Energy vs. Source Size;Energy;Source size", 10, 1.5, 15, 10, 1.5, 5.0);
	hCanvas->Draw();
	grSS[2]->SetMarkerStyle(53);
	grSS[2]->Draw("same p");
	grGausSS[2]->SetMarkerStyle(54);
	grGausSS[2]->SetMarkerColor(kRed);
	grGausSS[2]->Draw("same p");
	//}}}
}
