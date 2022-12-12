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
const std::string CaseName[NCase] = { "K_{s}^{0}(#tilde{K_{s}^{0}}+K_{s}^{0})", "(#tilde{K_{s}^{0}}+K_{s}^{0})K_{s}^{0}", "#tilde{K_{s}^{0}}#tilde{K_{s}^{0}}" };
float weightLeft, weightRight;
TPad* thisPad = 0;

#include "CFutils.h"

void RotCF(float Energy)
{
	//inital {{{
	TFile* ifPlots;
	if(Energy == 3.0) {
		//ifPlots = TFile::Open("~/CF/K0CF_output/3p0cutSetAcc2_1.1.0/3p0cutSetAcc2_1.1.0.root", "READ");
		//ifPlots = TFile::Open("/home/zla/CF/data/tmp/cutSetRotAcc2_1.1.0.root", "READ");
		ifPlots = TFile::Open("/home/zla/CF/data/tmp/newRot.root", "READ");
		//ifPlots = TFile::Open("/home/zla/CF/data/tmp/cutSetRotAcc10_1.1.0.root", "READ");
	} else if(Energy == 3.2) {
		ifPlots = TFile::Open("~/CF/K0CF_output/3p2cutSetAcc2_3.1.0/3p2cutSetAcc2_3.1.0.root", "READ");
	} else if(Energy == 3.5) {
		ifPlots = TFile::Open("~/CF/K0CF_output/3p5cutSetAcc2_3.1.0/3p5cutSetAcc2_3.1.0.root", "READ");
	} else if(Energy == 3.9) {
		ifPlots = TFile::Open("~/CF/K0CF_output/3p9cutSetAcc2_3.1.0/3p9cutSetAcc2_3.1.0.root", "READ");
	}
	//}}}

	//plots list {{{
	TH1F* hCFMisid[NCent] = { 0 };
	TH1F* hRotCF[NCent][NCase + 1] = { 0 };
	//}}}

	//calculate misid CF {{{
	RotCF(ifPlots, hCFMisid, hRotCF);
	//}}}

	//plotting {{{
	int col[NCase + 1] = { kBlack, kRed, kBlue, kGreen + 3 };

	TCanvas* caCFMisid = new TCanvas("caCFMisid", "caCFMisid", 1080, 1080);
	thisPad = (TPad*)caCFMisid->cd();
	thisPad->SetGrid();
	thisPad->SetMargin(0.15, 0.05, 0.15, 0.05);
	TLegend* leg = new TLegend(0.4, 0.6, 0.7, 0.8);
	leg->SetLineWidth(0);
	for(int icase = 1; icase < NCase + 1; ++icase) {
		hRotCF[2][icase]->SetMarkerStyle(24);
		hRotCF[2][icase]->SetMarkerColor(col[icase]);
		hRotCF[2][icase]->SetLineColor(col[icase]);
		hRotCF[2][icase]->SetLineWidth(2);

		hRotCF[2][icase]->SetMaximum(2.5);
		hRotCF[2][icase]->SetMinimum(0.7);
		hRotCF[2][icase]->SetStats(0);
		hRotCF[2][icase]->SetTitle(";q_{inv}(GeV/c);CF(q_{inv})");
		hRotCF[2][icase]->Draw("same");
		leg->AddEntry(hRotCF[2][icase], CaseName[icase - 1].c_str(), "lep");
	}
	leg->Draw("same");
	//}}}
}
