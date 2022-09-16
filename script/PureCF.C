#include <iostream>
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
#include "CFutils.h"

void PureCF(double Energy)
{
	//inital {{{
	TFile* ifPlots;
	if(Energy == 3.0) {
		ifPlots = TFile::Open("~/CF/K0CF_output/3p0cutSetAcc2_1.1.0/3p0cutSetAcc2_1.1.0.root", "READ");
	} else if(Energy == 3.2) {
		ifPlots = TFile::Open("~/CF/K0CF_output/3p2cutSetAcc2_3.1.0/3p2cutSetAcc2_3.1.0.root", "READ");
	} else if(Energy == 3.5) {
		ifPlots = TFile::Open("~/CF/K0CF_output/3p5cutSetAcc2_3.1.0/3p5cutSetAcc2_3.1.0.root", "READ");
	} else if(Energy == 3.9) {
		ifPlots = TFile::Open("~/CF/K0CF_output/3p9cutSetAcc2_3.1.0/3p9cutSetAcc2_3.1.0.root", "READ");
	}
	//}}}

	//plots list {{{
	TH1F* hCFRaw[NCent] = { 0 };
	TH1F* hCFMisid[NCent] = { 0 };
	TH1F* hCFPure[NCent] = { 0 };
	TH1F* hPairPurity[NCent] = { 0 };
	//}}}

	//processing {{{
	RawCF(ifPlots, hCFRaw);
	MisidCF(ifPlots, hCFMisid);
	PairPurity(ifPlots, hPairPurity);

	PureCF(hCFRaw, hCFMisid, hCFPure, hPairPurity);
	//}}}

	//plotting {{{
	TCanvas* caPure = new TCanvas("caPure", "caPure", 1080, 1080);
	thisPad = (TPad*)caPure->cd();
	thisPad->SetGrid();
	thisPad->SetMargin(0.12, 0.02, 0.12, 0.02);
	setMarker(hCFPure[2], 20, 2, kBlack);
	hCFPure[2]->GetXaxis()->SetRangeUser(0, 0.6);
	hCFPure[2]->SetLineColor(kBlack);
	hCFPure[2]->SetLineWidth(2);
	hCFPure[2]->GetXaxis()->CenterTitle();
	hCFPure[2]->GetYaxis()->CenterTitle();
	hCFPure[2]->SetXTitle("q_{inv}(GeV/c)");
	hCFPure[2]->SetYTitle("C_{pure}(q_{inv})");
	hCFPure[2]->SetStats(0);
	hCFPure[2]->SetTitle("");
	hCFPure[2]->SetMaximum(hCFPure[2]->GetMaximum() * 1.3);
	hCFPure[2]->SetMinimum(hCFPure[2]->GetMinimum() * 0.7);
	hCFPure[2]->Draw("ep");
	drawYBaseLine(1, thisPad, kBlue, 2, 2);
	drawText(0.35, 0.9, 0.05, Form("%.1f GeV Au+Au collisions", Energy));
	drawText(0.35, 0.8, 0.03, "-0.8<y<0.4, 0.2<p_{T}<1.5(p_{T}>0.3 @ y>0) GeV/c");
	//}}}
}
