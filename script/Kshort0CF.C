#include <iostream>
#include "utils.h"
#include "Kshort0_SideBandWeight.C"

void Kshort0CF_v1()
{
	//inital{{{
	//TFile* ifPlots = TFile::Open("~/CF/analysis_dev/3p5Default.root");	//all rapidity 3p5
	TFile* ifPlots = TFile::Open("~/CF/analysis_dev/cutSet3/3p5cutSet3.root");	//all rapidity 3p5
	TFile* ifRotPlots = TFile::Open("RotCF.root");
	const int NCent9 = 9;
	const int NCent = 3;
	const int NRap = 3;
	const int NCase = 3;
	const int Rebin = 20;
	const float NormalizeLower = 0.5, NormalizeUpper = 1.0;
	const float RapEdge[NRap + 1] = { -1.00, -0.8, 0.4, 0.5 };
	const std::string CentName[NCent] = { "0-10%", "10-60%", "0-60%" };
	const std::string CaseName[NCase] = { "K_{s}^{0}#tilde{K_{s}^{0}}", "#tilde{K_{s}^{0}}K_{s}^{0}", "#tilde{K_{s}^{0}}#tilde{K_{s}^{0}}" };
	const int Cent9To3[NCent9] = { -1, -1, 1, 1, 1, 1, 1, 0, 0 };
	const int Cent9To1[NCent9] = { -1, -1, 2, 2, 2, 2, 2, 2, 2 };
	const float LeftWeight = 0.61, RightWeight = 0.39;
	TPad* thisPad = 0;
	//}}}

	//plots list{{{
	TH1F* hSameQinvCent9[NCent9] = { 0 };
	TH1F* hMixQinvCent9[NCent9] = { 0 };
	TH1F* hSameLeftSideKqinvCent9[NCent9][NCase + 1] = { 0 };
	TH1F* hMixLeftSideKqinvCent9[NCent9][NCase + 1] = { 0 };
	TH1F* hMixKqinvLeftWeightCent9[NCent9][NCase + 1] = { 0 };
	TH1F* hSameRightSideKqinvCent9[NCent9][NCase + 1] = { 0 };
	TH1F* hMixRightSideKqinvCent9[NCent9][NCase + 1] = { 0 };
	TH1F* hMixKqinvRightWeightCent9[NCent9][NCase + 1] = { 0 };

	TH1F* hSameQinv[NCent] = { 0 };
	TH1F* hMixQinv[NCent] = { 0 };
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
	TH1F* hCFMisid[NCent] = { 0 };
	TH1F* hCFRaw[NCent] = { 0 };
	TH1F* hCFGenuine[NCent] = { 0 };
	//}}}

	//get plots{{{
	for(int icent = 0; icent < NCent9; ++icent) {
		hSameQinvCent9[icent] = (TH1F*)getCopy(ifPlots, Form("hSameKqinv_cent%d", icent));
		hMixQinvCent9[icent] = (TH1F*)getCopy(ifPlots, Form("hMixKqinv_cent%d", icent));
		hSameQinvCent9[icent]->Sumw2();
		hMixQinvCent9[icent]->Sumw2();
		hSameQinvCent9[icent]->Rebin(Rebin);
		hMixQinvCent9[icent]->Rebin(Rebin);
		for(int icase = 0; icase < NCase + 1; ++icase) {
			hSameLeftSideKqinvCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("hSameLeftSideKqinv_cent%d_case%d", icent, icase));
			hMixLeftSideKqinvCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("hMixLeftSideKqinv_cent%d_case%d", icent, icase));
			hMixKqinvLeftWeightCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("hMixKqinvLeftWeight_cent%d_case%d", icent, icase));
			hSameRightSideKqinvCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("hSameRightSideKqinv_cent%d_case%d", icent, icase));
			hMixRightSideKqinvCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("hMixRightSideKqinv_cent%d_case%d", icent, icase));
			hMixKqinvRightWeightCent9[icent][icase] = (TH1F*)getCopy(ifPlots, Form("hMixKqinvRightWeight_cent%d_case%d", icent, icase));

			hSameLeftSideKqinvCent9[icent][icase]->Sumw2();
			hMixLeftSideKqinvCent9[icent][icase]->Sumw2();
			hMixKqinvLeftWeightCent9[icent][icase]->Sumw2();
			hSameRightSideKqinvCent9[icent][icase]->Sumw2();
			hMixRightSideKqinvCent9[icent][icase]->Sumw2();
			hMixKqinvRightWeightCent9[icent][icase]->Sumw2();

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
			hLeftSideCF[icent][icase]->Divide(hSameLeftSideKqinv[icent][icase], hMixLeftSideKqinv[icent][icase], 1 / getIntegral(hSameLeftSideKqinv[icent][icase], NormalizeLower, NormalizeUpper), 1 / getIntegral(hMixLeftSideKqinv[icent][icase], NormalizeLower, NormalizeUpper));
			hRightSideCF[icent][icase]->Divide(hSameRightSideKqinv[icent][icase], hMixRightSideKqinv[icent][icase], 1 / getIntegral(hSameRightSideKqinv[icent][icase], NormalizeLower, NormalizeUpper), 1 / getIntegral(hMixRightSideKqinv[icent][icase], NormalizeLower, NormalizeUpper));
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
		scale(hLeftSideCFWeight[icent], LeftWeight, 1);
		hRightSideCFWeight[icent] = (TH1F*)hRightSideCFTot[icent]->Clone();
		hRightSideCFWeight[icent]->SetName(Form("hRightSideCFWeight_cent%d", icent));
		scale(hRightSideCFWeight[icent], (1 - LeftWeight), 1);

		addHist(hCFMisid[icent], hLeftSideCFWeight[icent], Form("hCFMisid_cent%d", icent), 1, 1);
		addHist(hCFMisid[icent], hRightSideCFWeight[icent], Form("hCFMisid_cent%d", icent), 1, 1);
	}

	for(int icent = 0; icent < NCent; ++icent) {
		hCFRaw[icent] = (TH1F*)hSameQinv[icent]->Clone();
		hCFRaw[icent]->SetName(Form("hCFRaw_cent%d", icent));
		hCFRaw[icent]->Divide(hSameQinv[icent], hMixQinv[icent], 1 / getIntegral(hSameQinv[icent], NormalizeLower, NormalizeUpper), 1 / getIntegral(hMixQinv[icent], NormalizeLower, NormalizeUpper));

		addHist(hCFGenuine[icent], hCFRaw[icent], Form("hCFGenuine_cent%d", icent), 1, 1);
		addHist(hCFGenuine[icent], hCFMisid[icent], Form("hCFGenuine_cent%d", icent), -1, 1);
		divide(hCFGenuine[icent], hLeftWeight[icent][3], 1);
	}
	//}}}

	//plotting{{{
	TCanvas* caWeight = new TCanvas("caWegiht", "caWeight", 1920, 1080);
	caWeight->Divide(4, 3);
	for(int icent = 0; icent < NCent; ++icent) {
		for(int icase = 0; icase < NCase + 1; ++icase) {
			caWeight->cd(icent * 4 + icase + 1);
			setMarker(hLeftWeight[icent][icase]);
			setMarker(hRightWeight[icent][icase], 20, 1, kRed);

			hLeftWeight[icent][icase]->SetMaximum(1.1);
			hLeftWeight[icent][icase]->SetMinimum(0);
			hLeftWeight[icent][icase]->Draw("ep");
			hRightWeight[icent][icase]->Draw("same ep");
		}
	}

	TCanvas* caPP = new TCanvas("caPP", "caPP", 640, 820);
	for(int icent = 2; icent < NCent; ++icent) {
		TPad* thisPad = (TPad*)caPP->cd();
		thisPad->SetGrid();
		setMarker(hLeftWeight[icent][3]);
		hLeftWeight[icent][3]->SetMaximum(1.1);
		hLeftWeight[icent][3]->SetMinimum(0);
		hLeftWeight[icent][3]->SetStats(0);
		hLeftWeight[icent][3]->SetTitle("");
		setXYTitle(hLeftWeight[icent][3], "q_{inv}(GeV/c)", "Pair Purity");
		hLeftWeight[icent][3]->Draw("ep");
		drawYBaseLine(1, thisPad, kBlue, 2, 2);
		drawText(0.3, 0.5, 0.04, Form("Centrality:%s", CentName[icent].c_str()));
		drawText(0.3, 0.45, 0.03, Form("average pair purity: %.0f%%", hLeftWeight[icent][3]->Integral(1, 25) / 0.25));
	}

	TCanvas* caQinv = new TCanvas("caQinv", "caQinv", 1920, 1080);
	caQinv->Divide(3, 3);
	for(int icent = 0; icent < NCent; ++icent) {
		for(int icase = 0; icase < NCase; ++icase) {
			caQinv->cd(icent * 3 + icase + 1)->SetLogy();
			hMixLeftSideKqinv[icent][icase]->GetYaxis()->SetRangeUser(1, hMixLeftSideKqinv[icent][icase]->GetMaximum() * 5);
			setMarker(hMixLeftSideKqinv[icent][icase]);
			setMarker(hSameLeftSideKqinv[icent][icase], 20, 1, kRed);
			setMarker(hMixRightSideKqinv[icent][icase], 21, 1, kBlack);
			setMarker(hSameRightSideKqinv[icent][icase], 21, 1, kRed);
			hMixLeftSideKqinv[icent][icase]->Draw("ep");
			hSameLeftSideKqinv[icent][icase]->Draw("same ep");
			hMixRightSideKqinv[icent][icase]->Draw("same ep");
			hSameRightSideKqinv[icent][icase]->Draw("same ep");
		}
	}

	TCanvas* caSideCF = new TCanvas("caSideCF", "caSideCF", 1920, 1080);
	caSideCF->Divide(3, 3);
	for(int icent = 0; icent < NCent; ++icent) {
		for(int icase = 0; icase < NCase; ++icase) {
			TPad* thisPad = (TPad*)caSideCF->cd(icent * 3 + icase + 1);
			TLegend* leg = new TLegend(0.4, 0.5, 0.7, 0.8);
			leg->SetLineWidth(0);
			setMarker(hLeftSideCF[icent][icase]);
			setMarker(hRightSideCF[icent][icase], 20, 1, kRed);
			hLeftSideCF[icent][icase]->SetMaximum(getMax(hLeftSideCF[icent][icase], hRightSideCF[icent][icase]) * 1.2);
			hLeftSideCF[icent][icase]->SetStats(0);
			hLeftSideCF[icent][icase]->SetTitle(Form("C_{%s,left/right}(q_{inv}) dist. @Cent:%s, y:-0.8~0.4", CaseName[icase].c_str(), CentName[icent].c_str()));
			setXYTitle(hLeftSideCF[icent][icase], "q_{inv}", "C_{misid}(q_{inv})");
			hLeftSideCF[icent][icase]->Draw("ep");
			hRightSideCF[icent][icase]->Draw("same ep");
			drawYBaseLine(1, thisPad, kBlue, 2, 2);
			leg->AddEntry(hLeftSideCF[icent][icase], Form("C_{%s,left}(q_{inv})", CaseName[icase].c_str()), "lp");
			leg->AddEntry(hRightSideCF[icent][icase], Form("C_{%s,right}(q_{inv})", CaseName[icase].c_str()), "lp");
			leg->Draw("same");
		}
	}

	TCanvas* caTotSideCF = new TCanvas("caTotSideCF", "caTotSideCF", 1920, 880);
	caTotSideCF->Divide(3, 1);
	for(int icent = 0; icent < NCent; ++icent) {
		TPad* thisPad = (TPad*)caTotSideCF->cd(icent + 1);
		thisPad->SetTopMargin(0.05);
		thisPad->SetRightMargin(0.05);
		TLegend* leg = new TLegend(0.3, 0.5, 0.8, 0.7);
		leg->SetLineWidth(0);
		setMarker(hLeftSideCFTot[icent]);
		setMarker(hRightSideCFTot[icent], 20, 1, kRed);
		hLeftSideCFTot[icent]->SetMaximum(getMax(hLeftSideCFTot[icent], hRightSideCFTot[icent]) * 1.2);
		hLeftSideCFTot[icent]->SetTitle("");
		setXYTitle(hLeftSideCFTot[icent], "q_{inv}", "C_{misid}(q_{inv})");
		hLeftSideCFTot[icent]->SetStats(0);
		hLeftSideCFTot[icent]->Draw("ep");
		hRightSideCFTot[icent]->Draw("same ep");
		leg->AddEntry(hLeftSideCFTot[icent], Form("C_{misid,left}(q_{inv}) @%s", CentName[icent].c_str()), "lp");
		leg->AddEntry(hRightSideCFTot[icent], Form("C_{misid,right}(q_{inv}) @%s", CentName[icent].c_str()), "lp");
		leg->Draw("same");
		drawYBaseLine(1, thisPad, kBlue, 2, 2);
	}

	//TCanvas* caCFMisid = new TCanvas("caCFMisid", "caCFMisid", 1640, 620);
	//TCanvas* caCFMisid = new TCanvas("caCFMisid", "caCFMisid", 720, 1080);
	TCanvas* caCFMisid = new TCanvas("caCFMisid", "caCFMisid", 1080, 1080);
	//caCFMisid->Divide(3, 2);
	//caCFMisid->Divide(2, 1);
	for(int icent = 2; icent < NCent; ++icent) {
		//TPad* thisPad = (TPad*)caCFMisid->cd(icent + 1);
		TPad* thisPad = (TPad*)caCFMisid->cd();
		thisPad->SetTopMargin(0.02);
		thisPad->SetRightMargin(0.02);
		thisPad->SetBottomMargin(0.12);
		thisPad->SetLeftMargin(0.12);
		TLegend* leg = new TLegend(0.4, 0.4, 0.8, 0.7);
		leg->SetLineWidth(0);
		setMarker(hCFMisid[icent], 20, 2, kBlue);
		setMarker(hLeftSideCFWeight[icent], 21, 2, kBlack);
		setMarker(hRightSideCFWeight[icent], 21, 2, kRed);
		setMarker(hCFRaw[icent], 20, 2, 8);
		hCFMisid[icent]->SetMaximum(getMax(hCFRaw[icent], hCFMisid[icent]) * 1.2);
		hCFMisid[icent]->SetTitle("");
		hCFMisid[icent]->SetStats(0);
		setXYTitle(hCFMisid[icent], "q_{inv}(GeV/c)", "C(q_{inv})");
		hCFMisid[icent]->GetXaxis()->SetTitleSize(0.055);
		hCFMisid[icent]->GetYaxis()->SetTitleSize(0.055);
		//hCFMisid[icent]->SetMinimum(0.5);
		hCFMisid[icent]->Draw("ep");
		hLeftSideCFWeight[icent]->Draw("same ep");
		hRightSideCFWeight[icent]->Draw("same ep");
		hCFRaw[icent]->Draw("same ep");
		drawYBaseLine(1, thisPad, kBlue, 2, 2);
		leg->AddEntry(hLeftSideCFWeight[icent], "#omega_{left}#lambda_{misid}C_{misid,left}(q_{inv})", "lp");
		leg->AddEntry(hRightSideCFWeight[icent], "#omega_{right}#lambda_{misid}C_{misid,right}(q_{inv})", "lp");
		leg->AddEntry(hCFMisid[icent], "#lambda_{misid}C_{misid}(q_{inv})", "lp");
		leg->AddEntry(hCFRaw[icent], "C_{raw}(q_{inv})", "lp");
		leg->Draw("same");
		drawText(0.3, 0.9, 0.035, "3.5 GeV Au+Au collisions", "-0.5<y<0.0, 0.2<p_{T}<1.2 GeV/c", Form("Centrality:%s", CentName[icent].c_str()));
		//drawText(0.3, 0.9, 0.035, "3.0 GeV Au+Au collisions", "-0.8<y<0.4, 0.2<p_{T}<1.5(p_{T}>0.3 @ y>0) GeV/c", Form("Centrality:%s", CentName[icent].c_str()));

		////thisPad = (TPad*)caCFMisid->cd(icent + 4);
		//thisPad = (TPad*)caCFMisid->cd();
		//thisPad->SetTopMargin(0.02);
		//thisPad->SetRightMargin(0.02);
		//thisPad->SetBottomMargin(0.12);
		//thisPad->SetLeftMargin(0.12);
		//TLegend* leg2 = new TLegend(0.4, 0.4, 0.8, 0.7);
		//leg2->SetLineWidth(0);
		//setMarker(hCFGenuine[icent], 20, 2, kBlack);
		//hCFGenuine[icent]->SetTitle("");
		//setXYTitle(hCFGenuine[icent], "q_{inv}(GeV/c)", "C(q_{inv})");
		//hCFGenuine[icent]->SetStats(0);
		//hCFGenuine[icent]->GetYaxis()->SetRangeUser(0.8, 2.0);
		//hCFGenuine[icent]->GetXaxis()->SetTitleSize(0.058);
		//hCFGenuine[icent]->GetYaxis()->SetTitleSize(0.058);
		//hCFGenuine[icent]->Draw("ep");
		//hCFRaw[icent]->Draw("same ep");
		//hCFMisid[icent]->Draw("same ep");
		//drawYBaseLine(1, thisPad, kBlue, 2, 2);
		//leg2->AddEntry(hCFRaw[icent], "C_{raw}(q_{inv})", "lp");
		//leg2->AddEntry(hCFGenuine[icent], "C_{pure(genuine)}(q_{inv})", "lp");
		//leg2->AddEntry(hCFMisid[icent], "C_{misid}(q_{inv})", "lp");
		//leg2->Draw("same ep");
		//hCFGenuine[icent]->Draw("same ep");
		////drawText(0.3, 0.9, 0.035, "3.2 GeV Au+Au collisions", "-0.8<y<0.4, 0.2<p_{T}<1.5(p_{T}>0.3 @ y>0) GeV/c", Form("Centrality:%s", CentName[icent].c_str()));
		//drawText(0.3, 0.9, 0.035, "3.5 GeV Au+Au collisions", "-0.5<y<0.0, 0.2<p_{T}<1.2 GeV/c", Form("Centrality:%s", CentName[icent].c_str()));
	}

	TCanvas* caCFGenuine = new TCanvas("caCFGenuine", "caCFGenuine", 1024, 820);
	thisPad = (TPad*)caCFGenuine->cd();
	thisPad->SetGridx();
	thisPad->SetGridy();
	thisPad->SetTopMargin(0.02);
	thisPad->SetRightMargin(0.02);
	thisPad->SetBottomMargin(0.12);
	thisPad->SetLeftMargin(0.12);
	int col[NCent] = { kBlue, kBlack, kRed };
	TLegend* legCFGenuine = new TLegend(0.42, 0.6, 0.65, 0.75);
	legCFGenuine->SetLineWidth(0);
	for(int icent = 0; icent < NCent; ++icent) {
		setMarker(hCFGenuine[icent], 20, 1, col[icent]);
		hCFGenuine[icent]->SetLineColor(col[icent]);
		hCFGenuine[icent]->SetLineWidth(1);
		hCFGenuine[icent]->GetXaxis()->CenterTitle();
		hCFGenuine[icent]->GetYaxis()->CenterTitle();
		legCFGenuine->AddEntry(hCFGenuine[icent], Form("%s", CentName[icent].c_str()), "lpe");
		if(icent == 0) {
			hCFGenuine[icent]->Draw("ep");
			hCFGenuine[icent]->GetYaxis()->SetTitle("C_{pure}(q_{inv})");
			hCFGenuine[icent]->SetTitle("");
			hCFGenuine[icent]->SetStats(0);
			hCFGenuine[icent]->SetMaximum(hCFGenuine[icent]->GetMaximum() * 1.3);
		} else {
			hCFGenuine[icent]->Draw("same ep");
		}
	}
	drawYBaseLine(1, thisPad, kBlue, 2, 2);
	drawText(0.35, 0.9, 0.05, "3 GeV Au+Au collisions");
	//TLatex* ltx = new TLatex(0.42, 0.8, "-0.8<y<0.4, 0.2<p_{T}<1.5(0.3<p_{T} if y>0)");
	//ltx->SetNDC();
	//ltx->SetTextAlign(23);
	//ltx->SetTextSize(0.03);
	drawText(0.35, 0.8, 0.03, "-0.8<y<0.4, 0.2<p_{T}<1.5(p_{T}>0.3 @ y>0) GeV/c");
	legCFGenuine->Draw("same");
	//ltx->Draw("same");
	//}}}

	//save cf pure{{{
	TFile* ofCFPure = new TFile("CFPure.root", "RECREATE");
	ofCFPure->cd();
	for(int icent = 0; icent < NCent; ++icent) {
		for(int irap = 0; irap < NRap; ++irap) {
			hCFGenuine[icent]->Write();
		}
	}
	//}}}
}


void Kshort0CF_v1(const char* ifName, TH1F***& hCFMisid, float SideBandLower, float SideBandUpper, float SideBand2Lower, float SideBand2Upper)
{
	//inital{{{
	TFile* ifPlots = TFile::Open(ifName);	//side band left and right
	const int NCent9 = 9;
	const int NCent = 3;
	const int NRap = 3;
	const int NCase = 3;
	const int Rebin = 20;
	const float NormalizeLower = 0.5, NormalizeUpper = 1.0;
	const float RapEdge[NRap + 1] = { -1.00, -0.8, 0.4, 0.5 };
	const std::string CentName[NCent] = { "0-10%", "10-60%", "0-60%" };
	const std::string CaseName[NCase] = { "K_{s}^{0}#tilde{K_{s}^{0}}", "#tilde{K_{s}^{0}}K_{s}^{0}", "#tilde{K_{s}^{0}}#tilde{K_{s}^{0}}" };
	const int Cent9To3[NCent9] = { -1, -1, 1, 1, 1, 1, 1, 1, 0 };
	const int Cent9To1[NCent9] = { -1, -1, 2, 2, 2, 2, 2, 2, 2 };
	float LeftWeight = 0.61, RightWeight = 0.39;
	//}}}

	//plots list{{{
	TH1F* hSameQinvCent9[NCent9][NRap] = { 0 };
	TH1F* hMixQinvCent9[NCent9][NRap] = { 0 };
	TH1F* hSameLeftSideKqinvCent9[NCent9][NRap][NCase + 1] = { 0 };
	TH1F* hMixLeftSideKqinvCent9[NCent9][NRap][NCase + 1] = { 0 };
	TH1F* hMixKqinvLeftWeightCent9[NCent9][NRap][NCase + 1] = { 0 };
	TH1F* hSameRightSideKqinvCent9[NCent9][NRap][NCase + 1] = { 0 };
	TH1F* hMixRightSideKqinvCent9[NCent9][NRap][NCase + 1] = { 0 };
	TH1F* hMixKqinvRightWeightCent9[NCent9][NRap][NCase + 1] = { 0 };

	TH1F* hSameQinv[NCent][NRap] = { 0 };
	TH1F* hMixQinv[NCent][NRap] = { 0 };
	TH1F* hSameLeftSideKqinv[NCent][NRap][NCase + 1] = { 0 };
	TH1F* hMixLeftSideKqinv[NCent][NRap][NCase + 1] = { 0 };
	TH1F* hMixKqinvLeftWeight[NCent][NRap][NCase + 1] = { 0 };
	TH1F* hSameRightSideKqinv[NCent][NRap][NCase + 1] = { 0 };
	TH1F* hMixRightSideKqinv[NCent][NRap][NCase + 1] = { 0 };
	TH1F* hMixKqinvRightWeight[NCent][NRap][NCase + 1] = { 0 };

	TH1F* hLeftWeight[NCent][NRap][NCase + 1] = { 0 };
	TH1F* hRightWeight[NCent][NRap][NCase + 1] = { 0 };
	TH1F* hLeftSideCF[NCent][NRap][NCase + 1] = { 0 };
	TH1F* hRightSideCF[NCent][NRap][NCase + 1] = { 0 };

	TH1F* hLeftSideCFTot[NCent][NRap] = { 0 };
	TH1F* hRightSideCFTot[NCent][NRap] = { 0 };

	TH1F* hLeftSideCFWeight[NCent][NRap] = { 0 };
	TH1F* hRightSideCFWeight[NCent][NRap] = { 0 };
	//TH1F* hCFMisid[NCent][NRap] = { 0 };
	TH1F* hCFRaw[NCent][NRap] = { 0 };
	TH1F* hCFGenuine[NCent][NRap] = { 0 };
	//}}}

	//get weightLeft and weightRight {{{
	LeftWeight = Kshort0_SideBandWeight(ifPlots, SideBandLower, SideBandUpper, SideBand2Lower, SideBand2Upper);
	RightWeight = 1 - LeftWeight;
	std::cout << ifName << " LeftWeight: " << LeftWeight << " RightWeight: " << RightWeight << std::endl;
	//}}}

	//get plots{{{
	for(int icent = 0; icent < NCent9; ++icent) {
		for(int irap = 0; irap < NRap; ++irap) {
			hSameQinvCent9[icent][irap] = (TH1F*)getCopy(ifPlots, Form("hSameKqinv_cent%d_rap%d", icent, irap));
			hMixQinvCent9[icent][irap] = (TH1F*)getCopy(ifPlots, Form("hMixKqinv_cent%d_rap%d", icent, irap));
			hSameQinvCent9[icent][irap]->Sumw2();
			hMixQinvCent9[icent][irap]->Sumw2();
			hSameQinvCent9[icent][irap]->Rebin(Rebin);
			hMixQinvCent9[icent][irap]->Rebin(Rebin);
			for(int icase = 0; icase < NCase + 1; ++icase) {
				hSameLeftSideKqinvCent9[icent][irap][icase] = (TH1F*)getCopy(ifPlots, Form("hSameLeftSideKqinv_cent%d_rap%d_case%d", icent, irap, icase));
				hMixLeftSideKqinvCent9[icent][irap][icase] = (TH1F*)getCopy(ifPlots, Form("hMixLeftSideKqinv_cent%d_rap%d_case%d", icent, irap, icase));
				hMixKqinvLeftWeightCent9[icent][irap][icase] = (TH1F*)getCopy(ifPlots, Form("hMixKqinvLeftWeight_cent%d_rap%d_case%d", icent, irap, icase));
				hSameRightSideKqinvCent9[icent][irap][icase] = (TH1F*)getCopy(ifPlots, Form("hSameRightSideKqinv_cent%d_rap%d_case%d", icent, irap, icase));
				hMixRightSideKqinvCent9[icent][irap][icase] = (TH1F*)getCopy(ifPlots, Form("hMixRightSideKqinv_cent%d_rap%d_case%d", icent, irap, icase));
				hMixKqinvRightWeightCent9[icent][irap][icase] = (TH1F*)getCopy(ifPlots, Form("hMixKqinvRightWeight_cent%d_rap%d_case%d", icent, irap, icase));

				hSameLeftSideKqinvCent9[icent][irap][icase]->Sumw2();
				hMixLeftSideKqinvCent9[icent][irap][icase]->Sumw2();
				hMixKqinvLeftWeightCent9[icent][irap][icase]->Sumw2();
				hSameRightSideKqinvCent9[icent][irap][icase]->Sumw2();
				hMixRightSideKqinvCent9[icent][irap][icase]->Sumw2();
				hMixKqinvRightWeightCent9[icent][irap][icase]->Sumw2();

				hSameLeftSideKqinvCent9[icent][irap][icase]->Rebin(Rebin);
				hMixLeftSideKqinvCent9[icent][irap][icase]->Rebin(Rebin);
				hMixKqinvLeftWeightCent9[icent][irap][icase]->Rebin(Rebin);
				hSameRightSideKqinvCent9[icent][irap][icase]->Rebin(Rebin);
				hMixRightSideKqinvCent9[icent][irap][icase]->Rebin(Rebin);
				hMixKqinvRightWeightCent9[icent][irap][icase]->Rebin(Rebin);
			}
		}
	}

	for(int irap = 0; irap < NRap; ++irap) {
		for(int icase = 0; icase < NCase + 1; ++icase) {
			for(int icent = 2; icent < NCent9; ++icent) {
				if(Cent9To3[icent] < 0 || Cent9To1[icent] < 0) {
					continue;
				}
				addHist(hSameLeftSideKqinv[Cent9To3[icent]][irap][icase], hSameLeftSideKqinvCent9[icent][irap][icase], Form("hSameLeftSideKqinv_cent%d_rap%d_case%d", Cent9To3[icent], irap, icase));
				addHist(hMixLeftSideKqinv[Cent9To3[icent]][irap][icase], hMixLeftSideKqinvCent9[icent][irap][icase], Form("hMixLeftSideKqinv_cent%d_rap%d_case%d", Cent9To3[icent], irap, icase));
				addHist(hMixKqinvLeftWeight[Cent9To3[icent]][irap][icase], hMixKqinvLeftWeightCent9[icent][irap][icase], Form("hMixKqinvLeftWeight_cent%d_rap%d_case%d", Cent9To3[icent], irap, icase));
				addHist(hSameRightSideKqinv[Cent9To3[icent]][irap][icase], hSameRightSideKqinvCent9[icent][irap][icase], Form("hSameRightSideKqinv_cent%d_rap%d_case%d", Cent9To3[icent], irap, icase));
				addHist(hMixRightSideKqinv[Cent9To3[icent]][irap][icase], hMixRightSideKqinvCent9[icent][irap][icase], Form("hMixRightSideKqinv_cent%d_rap%d_case%d", Cent9To3[icent], irap, icase));
				addHist(hMixKqinvRightWeight[Cent9To3[icent]][irap][icase], hMixKqinvRightWeightCent9[icent][irap][icase], Form("hMixKqinvRightWeight_cent%d_rap%d_case%d", Cent9To3[icent], irap, icase));

				addHist(hSameLeftSideKqinv[Cent9To1[icent]][irap][icase], hSameLeftSideKqinvCent9[icent][irap][icase], Form("hSameLeftSideKqinv_cent%d_rap%d_case%d", Cent9To1[icent], irap, icase));
				addHist(hMixLeftSideKqinv[Cent9To1[icent]][irap][icase], hMixLeftSideKqinvCent9[icent][irap][icase], Form("hMixLeftSideKqinv_cent%d_rap%d_case%d", Cent9To1[icent], irap, icase));
				addHist(hMixKqinvLeftWeight[Cent9To1[icent]][irap][icase], hMixKqinvLeftWeightCent9[icent][irap][icase], Form("hMixKqinvLeftWeight_cent%d_rap%d_case%d", Cent9To1[icent], irap, icase));
				addHist(hSameRightSideKqinv[Cent9To1[icent]][irap][icase], hSameRightSideKqinvCent9[icent][irap][icase], Form("hSameRightSideKqinv_cent%d_rap%d_case%d", Cent9To1[icent], irap, icase));
				addHist(hMixRightSideKqinv[Cent9To1[icent]][irap][icase], hMixRightSideKqinvCent9[icent][irap][icase], Form("hMixRightSideKqinv_cent%d_rap%d_case%d", Cent9To1[icent], irap, icase));
				addHist(hMixKqinvRightWeight[Cent9To1[icent]][irap][icase], hMixKqinvRightWeightCent9[icent][irap][icase], Form("hMixKqinvRightWeight_cent%d_rap%d_case%d", Cent9To1[icent], irap, icase));
			}
		}
		for(int icent = 2; icent < NCent9; ++icent) {
			if(Cent9To3[icent] < 0 || Cent9To1[icent] < 0) {
				continue;
			}
			addHist(hSameQinv[Cent9To3[icent]][irap], hSameQinvCent9[icent][irap], Form("hSameQinv_cent%d_rap%d", Cent9To3[icent], irap));
			addHist(hMixQinv[Cent9To3[icent]][irap], hMixQinvCent9[icent][irap], Form("hMixQinv_cent%d_rap%d", Cent9To3[icent], irap));

			addHist(hSameQinv[Cent9To1[icent]][irap], hSameQinvCent9[icent][irap], Form("hSameQinv_cent%d_rap%d", Cent9To1[icent], irap));
			addHist(hMixQinv[Cent9To1[icent]][irap], hMixQinvCent9[icent][irap], Form("hMixQinv_cent%d_rap%d", Cent9To1[icent], irap));
		}
	}
	//}}}

	//get weight vs. qinv{{{
	for(int icent = 0; icent < NCent; ++icent) {
		for(int irap = 0; irap < NRap; ++irap) {
			for(int icase = 0; icase < NCase + 1; ++icase) {
				addHist(hLeftWeight[icent][irap][icase], hMixKqinvLeftWeight[icent][irap][icase], Form("hLeftWeight_cent%d_rap%d_case%d", icent, irap, icase));
				addHist(hRightWeight[icent][irap][icase], hMixKqinvRightWeight[icent][irap][icase], Form("hRightWeight_cent%d_rap%d_case%d", icent, irap, icase));
				hLeftWeight[icent][irap][icase]->Divide(hLeftWeight[icent][irap][icase], hMixLeftSideKqinv[icent][irap][icase]);
				hRightWeight[icent][irap][icase]->Divide(hRightWeight[icent][irap][icase], hMixRightSideKqinv[icent][irap][icase]);
			}
		}
	}
	//}}}

	//calculate Side Band CF{{{
	for(int icent = 0; icent < NCent; ++icent) {
		for(int irap = 0; irap < NRap; ++irap) {
			for(int icase = 0; icase < NCase; ++icase) {
				addHist(hLeftSideCF[icent][irap][icase], hSameLeftSideKqinv[icent][irap][icase], Form("hLeftSideCF_cent%d_rap%d_case%d", icent, irap, icase));
				addHist(hRightSideCF[icent][irap][icase], hSameRightSideKqinv[icent][irap][icase], Form("hRightSideCF_cent%d_rap%d_case%d", icent, irap, icase));
				//hLeftSideCF[icent][irap][icase]->Divide(hSameLeftSideKqinv[icent][irap][icase], hMixLeftSideKqinv[icent][irap][icase], 1 / hSameLeftSideKqinv[icent][irap][icase]->Integral(), 1 / hMixLeftSideKqinv[icent][irap][icase]->Integral());
				//hRightSideCF[icent][irap][icase]->Divide(hSameRightSideKqinv[icent][irap][icase], hMixRightSideKqinv[icent][irap][icase], 1 / hSameRightSideKqinv[icent][irap][icase]->Integral(), 1 / hMixRightSideKqinv[icent][irap][icase]->Integral());
				hLeftSideCF[icent][irap][icase]->Divide(hSameLeftSideKqinv[icent][irap][icase], hMixLeftSideKqinv[icent][irap][icase], 1 / getIntegral(hSameLeftSideKqinv[icent][irap][icase], NormalizeLower, NormalizeUpper), 1 / getIntegral(hMixLeftSideKqinv[icent][irap][icase], NormalizeLower, NormalizeUpper));
				hRightSideCF[icent][irap][icase]->Divide(hSameRightSideKqinv[icent][irap][icase], hMixRightSideKqinv[icent][irap][icase], 1 / getIntegral(hSameRightSideKqinv[icent][irap][icase], NormalizeLower, NormalizeUpper), 1 / getIntegral(hMixRightSideKqinv[icent][irap][icase], NormalizeLower, NormalizeUpper));
			}
		}
	}

	for(int icent = 0; icent < NCent; ++icent) {
		for(int irap = 0; irap < NRap; ++irap) {
			for(int icase = 0; icase < 1; ++icase) {
				addHist(hLeftSideCFTot[icent][irap], hLeftSideCF[icent][irap][icase], Form("hLeftSideCFTot_cent%d_rap%d", icent, irap), hLeftWeight[icent][irap][icase], 1);
				addHist(hRightSideCFTot[icent][irap], hRightSideCF[icent][irap][icase], Form("hRightSideCFTot_cent%d_rap%d", icent, irap), hRightWeight[icent][irap][icase], 1);
			}
		}
	}

	for(int icent = 0; icent < NCent; ++icent) {
		for(int irap = 0; irap < NRap; ++irap) {
			hLeftSideCFWeight[icent][irap] = (TH1F*)hLeftSideCFTot[icent][irap]->Clone();
			hLeftSideCFWeight[icent][irap]->SetName(Form("hLeftSideCFWeight_cent%d_rap%d", icent, irap));
			scale(hLeftSideCFWeight[icent][irap], LeftWeight, 1);
			hRightSideCFWeight[icent][irap] = (TH1F*)hRightSideCFTot[icent][irap]->Clone();
			hRightSideCFWeight[icent][irap]->SetName(Form("hRightSideCFWeight_cent%d_rap%d", icent, irap));
			scale(hRightSideCFWeight[icent][irap], LeftWeight, 1);

			addHist(hCFMisid[icent][irap], hLeftSideCFWeight[icent][irap], Form("hCFMisid_cent%d_rap%d", icent, irap), 1, 1);
			addHist(hCFMisid[icent][irap], hRightSideCFWeight[icent][irap], Form("hCFMisid_cent%d_rap%d", icent, irap), 1, 1);
		}
	}

	for(int icent = 0; icent < NCent; ++icent) {
		for(int irap = 0; irap < NRap; ++irap) {
			hCFRaw[icent][irap] = (TH1F*)hSameQinv[icent][irap]->Clone();
			hCFRaw[icent][irap]->SetName(Form("hCFRaw_cent%d_rap%d", icent, irap));
			hCFRaw[icent][irap]->Divide(hSameQinv[icent][irap], hMixQinv[icent][irap], 1 / getIntegral(hSameQinv[icent][irap], NormalizeLower, NormalizeUpper), 1 / getIntegral(hMixQinv[icent][irap], NormalizeLower, NormalizeUpper));

			addHist(hCFGenuine[icent][irap], hCFRaw[icent][irap], Form("hCFGenuine_cent%d_rap%d", icent, irap), 1, 1);
			addHist(hCFGenuine[icent][irap], hCFMisid[icent][irap], Form("hCFGenuine_cent%d_rap%d", icent, irap), -1, 1);
			divide(hCFGenuine[icent][irap], hLeftWeight[icent][irap][3], 1);
		}
	}
	//}}}

	//release memory{{{
	for(int icent = 0; icent < NCent9; ++icent) {
		for(int irap = 0; irap < NRap; ++irap) {
			delete hSameQinvCent9[icent][irap];
			delete hMixQinvCent9[icent][irap];
			for(int icase = 0; icase < NCase + 1; ++icase) {
				delete hSameLeftSideKqinvCent9[icent][irap][icase];
				delete hMixLeftSideKqinvCent9[icent][irap][icase];
				delete hMixKqinvLeftWeightCent9[icent][irap][icase];
				delete hSameRightSideKqinvCent9[icent][irap][icase];
				delete hMixRightSideKqinvCent9[icent][irap][icase];
				delete hMixKqinvRightWeightCent9[icent][irap][icase];
			}
		}
	}

	for(int icent = 0; icent < NCent; ++icent) {
		for(int irap = 0; irap < NRap; ++irap) {
			delete hSameQinv[icent][irap];
			delete hMixQinv[icent][irap];
			delete hCFRaw[icent][irap];
			delete hCFGenuine[icent][irap];
			for(int icase = 0; icase < NCase + 1; ++icase) {
				delete hSameLeftSideKqinv[icent][irap][icase];
				delete hMixLeftSideKqinv[icent][irap][icase];
				delete hMixKqinvLeftWeight[icent][irap][icase];
				delete hSameRightSideKqinv[icent][irap][icase];
				delete hMixRightSideKqinv[icent][irap][icase];
				delete hMixKqinvRightWeight[icent][irap][icase];
			}
		}
	}
	//}}}
	return (TH1F***)hCFMisid;
}
