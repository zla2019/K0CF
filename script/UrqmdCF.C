#include <iostream>
#define SETTING
const int NCent9 = 9;
const int NCent = 3;
const int NCase = 3;
const int Rebin = 5;
const float PtEdge[2] = { 0.2, 1.8 };
const float RapEdge[2] = { 1.0, 0.0 };
const int Cent9To3[NCent9] = { -1, -1, 1, 1, 1, 1, 1, 0, 0 };
const int Cent9To1[NCent9] = { -1, -1, 2, 2, 2, 2, 2, 2, 2 };
const std::string CentName[NCent] = { "0-20%", "20-60%", "0-60%" };
float weightLeft, weightRight;
TPad* thisPad = 0;
#include "CFutils.h"

void UrqmdCF(double energy = 3.0)
{
	//inital {{{
	TFile* ifPlots = nullptr;
	if(energy == 3.0) {
		ifPlots = TFile::Open("/home/zla/CF/data/urqmd/3p0/result_202210040411.root");
	} else if(energy == 3.2) {
		ifPlots = TFile::Open("/home/zla/CF/data/urqmd/3p2/result_202210060518.root");
	} else if(energy == 3.5) {
		ifPlots = TFile::Open("/home/zla/CF/data/urqmd/3p5/result_202210050157.root");
	} else if(energy == 3.9) {
		ifPlots = TFile::Open("/home/zla/CF/data/urqmd/3p9/result_202210040728.root");
	} else if(energy == 7.7) {
		ifPlots = TFile::Open("/home/zla/CF/data/tmp/result_202209292213.root");
	} else {
		std::cout << "ERROR: energy does not set" << std::endl;
		return;
	}
	//}}}

	//plots list {{{
	TH1F* hCF[NCent] = { 0 };
	TH1F* hCFQS[NCent] = { 0 };
	TH1F* hCFSI[NCent] = { 0 };
	TH1F* hCFWoCrab[NCent] = { 0 };
	//}}}

	//processing {{{
	getUrqmdCF(ifPlots, "", hCF);
	//getUrqmdCF(ifPlots, "QS", hCFQS);
	//getUrqmdCF(ifPlots, "SI", hCFSI);
	//getUrqmdCF(ifPlots, "WoCrab", hCFWoCrab);
	//}}}

	//Plotting {{{
	TCanvas* caCF = new TCanvas("", "", 900, 900);
	thisPad = (TPad*)caCF->cd();
	thisPad->SetMargin(0.12, 0.02, 0.12, 0.02);
	thisPad->SetGridx();
	thisPad->SetGridy();
	TLegend* legCF = new TLegend(0.5, 0.5, 0.85, 0.7);
	legCF->SetLineWidth(0);
	hCF[2]->SetTitle("");
	hCF[2]->GetYaxis()->SetRangeUser(0.6, 2.0);
	hCF[2]->SetStats(0);
	setXYTitle(hCF[2], "q_{inv}(GeV/c)", "CF(q_{inv})");
	hCF[2]->SetTitleSize(0.05, "XY");

	setStyle(hCF[2], 53, 1.4, kRed, 2, kRed);
	setStyle(hCFQS[2], 53, 1.4, kBlack, 2, kBlack);
	setStyle(hCFSI[2], 53, 1.4, kGreen - 5, 2, kGreen - 5);
	setStyle(hCFWoCrab[2], 53, 1.4, kBlue, 2, kBlue);
	legCF->AddEntry(hCF[2], "UrQMD with QS + SI", "lpe");
	legCF->AddEntry(hCFQS[2], "UrQMD with QS", "lpe");
	legCF->AddEntry(hCFSI[2], "UrQMD with SI", "lpe");
	legCF->AddEntry(hCFWoCrab[2], "UrQMD", "lpe");

	hCF[2]->GetXaxis()->SetRangeUser(0, 0.6);
	hCF[2]->Draw("ep");
	hCFQS[2]->Draw("same ep");
	hCFSI[2]->Draw("same ep");
	hCFWoCrab[2]->Draw("same ep");
	drawYBaseLine(1, thisPad, kBlue, 2, 2);
	drawText(0.3, 0.85, 0.035, "CF(q_{inv}) @0~60%, y:-1.0~0.0, p_{T}:0.2~1.8", Form("~9M %.1f GeV UrQMD Events", energy));
	legCF->Draw("same");
	//}}}
}
