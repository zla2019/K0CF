#include <iostream>
//#include "CFLLFunc.C"
#include "/home/zla/CF/code/CFLLFunc.C"
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
		ifPlots = TFile::Open("~/CF/K0CF_output/20220913StarCollaboration/3p0cutSetAcc2_1.1.0/3p0cutSetAcc2_1.1.0.root", "READ");
		//ifPlots = TFile::Open("~/CF/K0CF_output/3p0cutSetAcc2_1.1.0/3p0cutSetAcc2_1.1.0.root", "READ");
	} else if(Energy == 3.2) {
		//ifPlots = TFile::Open("~/CF/K0CF_output/20220913StarCollaboration/3p2cutSetAcc2_3.1.0/3p2cutSetAcc2_3.1.0.root", "READ");
		ifPlots = TFile::Open("~/CF/K0CF_output/3p2cutSetAcc2_3.1.0/3p2cutSetAcc2_3.1.0.root", "READ");
	} else if(Energy == 3.5) {
		//ifPlots = TFile::Open("~/CF/K0CF_output/20220913StarCollaboration/3p5cutSetAcc2_3.1.0/3p5cutSetAcc2_3.1.0.root", "READ");
		ifPlots = TFile::Open("~/CF/K0CF_output/3p5cutSetAcc2_3.1.0/3p5cutSetAcc2_3.1.0.root", "READ");
	} else if(Energy == 3.9) {
		//ifPlots = TFile::Open("~/CF/K0CF_output/20220913StarCollaboration/3p9cutSetAcc2_3.1.0/3p9cutSetAcc2_3.1.0.root", "READ");
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

	//fitting {{{
	TF1* fCFGaus = new TF1("fCF", "1 + [0] * exp((-1 * [1]^2 * (x)^2) / 0.038937929230)", 0, 1);
	fCFGaus->SetLineColor(kBlack);
	fCFGaus->SetParameter(0, 1);
	fCFGaus->SetParameter(1, 3);
	float fitLower = hCFPure[2]->GetBinLowEdge(1);
	float fitUpper = hCFPure[2]->GetBinLowEdge(hCFPure[2]->GetNbinsX()) + hCFPure[2]->GetBinWidth(hCFPure[2]->GetNbinsX());
	hCFPure[2]->Fit(fCFGaus, "RN", "", fitLower, fitUpper);

	TF1* fCFLL = new TF1("fCF", CFLL, 0, 0.4, 2);
	fCFLL->SetParameters(fCFGaus->GetParameters());
	fCFLL->SetParLimits(1, 0.2, 7);
	fCFLL->SetParLimits(0, 0.2, 1);
	hCFPure[2]->Fit(fCFLL, "RN", "", fitLower, 0.4);
	//}}}

	//plotting {{{
	TCanvas* caPure = new TCanvas("caPure", "caPure", 900, 900);
	TF1* fCFDrawGaus = new TF1("fCFDrawGaus", "1 + [0] * exp((-1 * [1]^2 * (x)^2) / 0.038937929230)", 0, 1);
	TF1* fCFDrawLL = new TF1("fCFDrawLL", CFLL, 0, 0.4, 2);
	TF1* fCFDrawLLSI = new TF1("fCFDrawLLSI", CFLLSI, 0, 0.4, 2);
	fCFDrawGaus->SetLineColor(kBlack);
	fCFDrawGaus->SetParameters(fCFGaus->GetParameters());
	fCFDrawGaus->SetLineWidth(3);
	fCFDrawLL->SetLineColor(kRed);
	fCFDrawLL->SetParameters(fCFLL->GetParameters());
	fCFDrawLL->SetLineWidth(3);
	fCFDrawLLSI->SetLineColor(kGreen-5);
	fCFDrawLLSI->SetParameters(fCFLL->GetParameters());
	fCFDrawLLSI->SetLineWidth(3);

	thisPad = (TPad*)caPure->cd();
	thisPad->SetGrid();
	thisPad->SetMargin(0.12, 0.02, 0.12, 0.02);
	TLegend* leg = new TLegend(0.25, 0.8, 0.6, 0.6);
	leg->SetLineWidth(0);
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
	hCFPure[2]->SetMaximum(2.0);
	hCFPure[2]->SetMinimum(0.6);
	hCFPure[2]->Draw("ep");
	fCFDrawGaus->Draw("same");
	fCFDrawLL->Draw("same");
	fCFDrawLLSI->Draw("same");
	leg->AddEntry(fCFDrawGaus, "Gaus", "l");
	leg->AddEntry(fCFDrawLL, "LL", "l");
	leg->AddEntry(fCFDrawLLSI, "SI(= LL - Gaus.)", "l");
	drawYBaseLine(1, thisPad, kBlue, 2, 2);
	drawText(0.35, 0.9, 0.05, Form("%.1f GeV Au+Au collisions", Energy));
	drawText(0.35, 0.85, 0.03, "-1.0<y<0.0, 0.2<p_{T}<1.8 GeV/c");
	drawText(0.6, 0.75, 0.04, Form("R_{inv}: %.2f#pm%.2f", fCFGaus->GetParameter(1), fCFGaus->GetParError(1)), Form("#lambda: %.2f#pm%.2f", fCFGaus->GetParameter(0), fCFGaus->GetParError(0)));
	drawText(0.6, 0.57, 0.04, kRed, Form("R_{inv}: %.2f#pm%.2f", fCFLL->GetParameter(1), fCFLL->GetParError(1)), Form("#lambda: %.2f#pm%.2f", fCFLL->GetParameter(0), fCFLL->GetParError(0)));
	leg->Draw("same");
	//}}}
}
