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

void MisidCF(float Energy)
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
	TH1F* hCFMisid[NCent] = { 0 };
	//}}}

	//calculate misid CF {{{
	MisidCF(ifPlots, hCFMisid);
	//}}}

	//plotting {{{
	TCanvas* caCFMisid = new TCanvas("caCFMisid", "caCFMisid", 1080, 1080);
	thisPad = (TPad*)caCFMisid->cd();
	thisPad->SetGrid();
	thisPad->SetMargin(0.15, 0.05, 0.15, 0.05);
	hCFMisid[2]->SetMaximum(2.2);
	hCFMisid[2]->SetMinimum(0.6);
	setMarker(hCFMisid[2], 20, 1.4, kBlue);
	hCFMisid[2]->SetLineColor(kBlue);
	hCFMisid[2]->SetLineWidth(3);
	hCFMisid[2]->GetXaxis()->SetRangeUser(0, 0.6);
	hCFMisid[2]->SetStats(0);
	hCFMisid[2]->SetTitle(";q_{inv} (GeV/c);C_{misid}(q_{inv})");
	hCFMisid[2]->Draw("ep");
	hCFMisid[2]->SetTitleSize(0.05, "XY");
	hCFMisid[2]->SetLabelSize(0.05, "XY");
	drawYBaseLine(1, thisPad, kBlue, 2, 2);
	drawText(0.25, 0.9, 0.04, Form("Au+Au collision @#sqrt{s_{NN}} = %.1f GeV", Energy), "Centrality: 0-60%", "y: -1.0~0.0", "p_{T}: 0.2~1.8");
	hCFMisid[2]->Draw("ep same");
	//}}}
}
