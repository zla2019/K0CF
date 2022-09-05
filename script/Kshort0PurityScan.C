#include <iostream>
#include "utils.h"
#include "func.h"

void Kshort0PurityScan(std::string inputName, std::string cutSetName, std::string nSigma)
{
	//initial{{{
	const int NCent9 = 9;
	const int NRap = 8;
	const int NCent = 3;
	const int NPt = 19;
	const std::string CentName[NCent] = { "0-10%", "10-60%", "0-60%" };
	const float RapEdge[NRap + 1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6};
	//const float PtEdge[NPt + 1] = { 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2 };
	const float PtEdge[NPt + 1] = { 0.10, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0 };
	TFile* ifPlots = new TFile(inputName.c_str());
	//TFile* ifPlots = new TFile("/home/zla/CF/analysis_3p2/3p2WensongAcc.root");
	//TFile* ifPlots = new TFile("/home/zla/CF/data/check/CF5sigma10chi2.root");
	//TFile* ifPlots = new TFile("/home/zla/CF/data/check/production_202206040433.root");
	//TFile* ifPlots = new TFile("/home/zla/CF/data/CFPlots/cutSet1/production_202207050355.root");
	//TFile* ifPlots = new TFile("/home/zla/CF/data/3p2/3p2Raw.root");
	//TFile* ifPlots = new TFile("/home/zla/CF/analysis_3p2/3p2OldAcc.root");
	//TFile* ifPlots = new TFile("/home/zla/CF/data/check/tmp/3sigma11chi2.root");
	//TFile* ifPlots = new TFile("/home/zla/CF/data/check/CFCheck202205142246.root");
	//TFile* ifPlots = new TFile("/home/zla/CF/data/CFPlots/cutSet8/3sigma12chi2.root");
	const float Mean[3][4] = { { 0.4979, 0.4981, 0.4979, 0.4949 }, { 0.4978, 0.4981, 0.4979, 0.4979 }, { 0.4979, 0.4981, 0.4980, 0.4980 } };
	const float Sigma[3][4] = { { 0.0027, 0.0035, 0.0033, 0.0033 }, { 0.0027, 0.0034, 0.0036, 0.0036 }, { 0.0027, 0.0035, 0.0039, 0.0039 } };
	const float NSigma = std::stof(nSigma);
	std::cout << "NSigma window: " << NSigma << std::endl;
	//TFile* ifPlots = new TFile("/home/zla/CF/data/check/cut0.root");
	//TFile* ifPlots = new TFile("/home/zla/CF/analysis/test.root");
	TH3F* hMassRapPt[NCent9];
	TH3F* hTotMassRapPt[NCent];
	TH1F* hMass[NCent][NRap][NPt];
	TH1F* hMassSig[NCent][NRap][NPt];
	TF1* fitFunc[NCent][NRap][NPt];
	TF1* fitFuncBkg[NCent][NRap][NPt];
	TF1* fitFuncSig[NCent][NRap][NPt];
	TH2F* hKAcc[NCent];
	TH2F* hKPurity[NCent];
	float sig[NCent][NRap][NPt];
	float bkg[NCent][NRap][NPt];
	float tot[NCent][NRap][NPt];
	float mean[NCent][NRap][NPt];
	float sigma[NCent][NRap][NPt];
	float sigma_2[NCent][NRap][NPt];
	float quality[NCent][NRap][NPt];
	//}}}

	//get plots{{{
	for(int icent = 0; icent < NCent9; ++icent) {
		hMassRapPt[icent] = (TH3F*)getCopy(ifPlots, Form("hMassCascadeRapidityvsPt_cent%i", icent));
		if(icent == 2) {
			hTotMassRapPt[2] = (TH3F*)hMassRapPt[icent]->Clone();
			hTotMassRapPt[2]->SetName(Form("hMassCascadeRapidityvsPt_2_cent%i", icent));
		} else if(icent > 2) {
			hTotMassRapPt[2]->Add(hMassRapPt[icent]);
		}

		if(icent == 7) {
			hTotMassRapPt[0] = (TH3F*)hMassRapPt[icent]->Clone();
			hTotMassRapPt[0]->SetName(Form("hMassCascadeRapidityvsPt_0_cent%i", icent));
		} else if(icent > 7){
			hTotMassRapPt[0]->Add(hMassRapPt[icent]);
		}

		if(icent == 2) {
			hTotMassRapPt[1] = (TH3F*)hMassRapPt[icent]->Clone();
			hTotMassRapPt[1]->SetName(Form("hMassCascadeRapidityvsPt_1_cent%i", icent));
		} else if(icent < 7 && icent > 2) {
			hTotMassRapPt[1]->Add(hMassRapPt[icent]);
		}
	}
	for(int icent = 0; icent < NCent; ++icent) {
		//hKAcc[icent] = new TH2F(Form("hKAcc_cent%i", icent), Form("K_{s}^{0} Acc. @%s", CentName[icent].c_str()), 7, -1., 0.4, 11, -0.1, 2.1);
		//hKPurity[icent] = new TH2F(Form("hKPurity_cent%i", icent), Form("K_{s}^{0} Purity. @%s", CentName[icent].c_str()), 7, -1., 0.4, 11, -0.1, 2.1);
		hKAcc[icent] = new TH2F(Form("hKAcc_cent%i", icent), Form("K_{s}^{0} Acc. @%s", CentName[icent].c_str()), NRap, RapEdge, NPt, PtEdge);
		hKPurity[icent] = new TH2F(Form("hKPurity_cent%i", icent), Form("K_{s}^{0} Purity. @%s", CentName[icent].c_str()), NRap, RapEdge, NPt, PtEdge);
		for(int irap = 0; irap < NRap; ++irap) {
			for(int ipt = 0; ipt < NPt; ++ipt) {
				hMass[icent][irap][ipt] = projectionX(hTotMassRapPt[icent], RapEdge[irap], RapEdge[irap + 1], PtEdge[ipt], PtEdge[ipt + 1], Form("hMass_cent%i_rap%i_pt%i", icent, irap, ipt));
			}
		}
	}
	//}}}

	//fitting{{{
	for(int icent = 0; icent < NCent; ++icent) {
		for(int irap = 0; irap < NRap; ++irap) {
			//if(irap != 4) {
			//	continue;
			//}
			for(int ipt = 0; ipt < NPt; ++ipt) {
				std::cout << "fitting: icent: " << icent << " irap: " << irap << " ipt: " << ipt << std::endl;
				//hMass[icent][irap][ipt]->Fit(fitFuncTmp, "IMNR", "", 0.45, 0.55);
				TF1* fitFuncTmp = new TF1(Form("fitFuncTmp_rap%i", irap), fitFuncGaus, 0.45, 0.55, /*7*/6);
				//fit tunning{{{
				fitFuncTmp->SetParNames("scale", "mean", "sigma", "p0", "p1", "p2");
				fitFuncTmp->SetParLimits(1, 0.497, 0.499);
				fitFuncTmp->SetParLimits(0, 0., 100000);
				fitFuncTmp->SetParLimits(2, 0.002, 0.012);
				fitFuncTmp->SetParLimits(3, 0., 10000);
				hMass[icent][irap][ipt]->Fit(fitFuncTmp, "INMR", "", 0.45, 0.55);
				//}}}

				fitFunc[icent][irap][ipt] = new TF1(Form("fitFunc_cent%i_rap%i_pt%i", icent, irap, ipt), fitFuncDoubleGaus, 0.42, 0.58, 8);
				//fitFunc[icent][irap][ipt]->SetParNames("scale1", "mean", "sigma1", "p0", "p1", "p2"/*, "p3"*/);
				//fit tunning{{{
				fitFunc[icent][irap][ipt]->SetParNames("scale1", "mean", "sigma1", "scale2", "sigma2", "p0", "p1", "p2"/*, "p3"*/);
				fitFunc[icent][irap][ipt]->FixParameter(1, fitFuncTmp->GetParameter(1));
				fitFunc[icent][irap][ipt]->FixParameter(5, fitFuncTmp->GetParameter(3));
				fitFunc[icent][irap][ipt]->FixParameter(6, fitFuncTmp->GetParameter(4));
				fitFunc[icent][irap][ipt]->FixParameter(7, fitFuncTmp->GetParameter(5));
				fitFunc[icent][irap][ipt]->SetParLimits(2, fitFuncTmp->GetParameter(2) * 0.7, fitFuncTmp->GetParameter(2) * 1.2);
				fitFunc[icent][irap][ipt]->SetParLimits(0, fitFuncTmp->GetParameter(0) * 0.72, fitFuncTmp->GetParameter(0) * 0.8);
				fitFunc[icent][irap][ipt]->SetParLimits(3, 0.33*fitFuncTmp->GetParameter(0), 0.63*fitFuncTmp->GetParameter(0));
				fitFunc[icent][irap][ipt]->SetParLimits(4, fitFuncTmp->GetParameter(2)*1.6, fitFuncTmp->GetParameter(2)*1.9);
				//}}}
				hMass[icent][irap][ipt]->Fit(fitFunc[icent][irap][ipt], "MR", "", 0.45, 0.55);

				fitFuncBkg[icent][irap][ipt] = new TF1(Form("fitFuncBkg_cent%i_rap%i_pt%i", icent, irap, ipt), background, 0.42, 0.58, 3);
				fitFuncBkg[icent][irap][ipt]->SetParameters(&fitFunc[icent][irap][ipt]->GetParameters()[5]);
				fitFuncSig[icent][irap][ipt] = new TF1(Form("fitFuncSig_cent%i_rap%i_pt%i", icent, irap, ipt), doublegaus, 0.42, 0.58, 5);
				fitFuncSig[icent][irap][ipt]->SetParameters(&fitFunc[icent][irap][ipt]->GetParameters()[0]);
				float sigma1 = fitFuncSig[icent][irap][ipt]->GetParameter(2);
				float sigma2 = fitFuncSig[icent][irap][ipt]->GetParameter(4);
				float meanloc = fitFuncSig[icent][irap][ipt]->GetParameter(1);
				quality[icent][irap][ipt] = fitFunc[icent][irap][ipt]->GetChisquare() / fitFunc[icent][irap][ipt]->GetNDF();

				tot[icent][irap][ipt] = getIntegral(hMass[icent][irap][ipt], Mean[2][1] - NSigma * Sigma[2][1], Mean[2][1] + NSigma * Sigma[2][1]);	//nsigma parttern
				//tot[icent][irap][ipt] = getIntegral(hMass[icent][irap][ipt], 0.48, 0.51);
				TH1F* hSig = f2h(fitFuncSig[icent][irap][ipt], hMass[icent][irap][ipt]->GetBinWidth(1), 0.42, 0.58);
				TH1F* hBkg = f2h(fitFuncBkg[icent][irap][ipt], hMass[icent][irap][ipt]->GetBinWidth(1), 0.42, 0.58);
				hMassSig[icent][irap][ipt] = (TH1F*)hMass[icent][irap][ipt]->Clone();
				hMassSig[icent][irap][ipt]->Add(fitFuncBkg[icent][irap][ipt], -1);
				//sig[icent][irap][ipt] = getIntegral(hMassSig[icent][irap][ipt], meanloc - 3 * sigma1, meanloc + 3 * sigma1);
				//bkg[icent][irap][ipt] = getIntegral(hBkg, meanloc - 3 * sigma1, meanloc + 3 * sigma1);
				sig[icent][irap][ipt] = getIntegral(hMassSig[icent][irap][ipt], Mean[2][1] - NSigma * Sigma[2][1], Mean[2][1] + NSigma * Sigma[2][1]);	//nsigma parttern
				//sig[icent][irap][ipt] = getIntegral(hMassSig[icent][irap][ipt], 0.48, 0.51);
				//bkg[icent][irap][ipt] = getIntegral(hBkg, Mean[2][1] - NSigma * Sigma[2][1], Mean[2][1] + NSigma * Sigma[2][1]);
				bkg[icent][irap][ipt] = fitFuncSig[icent][irap][ipt]->Integral(Mean[2][1] - NSigma * Sigma[2][1], Mean[2][1] + NSigma * Sigma[2][1]);		//nsigma parttern
				//bkg[icent][irap][ipt] = fitFuncSig[icent][irap][ipt]->Integral(0.48, 0.51);
				mean[icent][irap][ipt] = meanloc;
				sigma[icent][irap][ipt] = sigma1;
				sigma_2[icent][irap][ipt] = sigma2;
				hKAcc[icent]->SetBinContent(irap + 1, ipt + 1, sig[icent][irap][ipt]);
				float purity = sig[icent][irap][ipt] / tot[icent][irap][ipt] * 100.;
				if(purity < 0 || std::isnan(purity)) {
					purity = 0;
				} else if(purity > 100) {
					purity = 100;
				}
				hKPurity[icent]->SetBinContent(irap + 1, ipt + 1, purity);
			}
		}
	}
	//}}}

	//plotting{{{
	TCanvas* caMass[NCent][NRap];
	for(int icent = 0; icent < NCent; ++icent) {
		for(int irap = 0; irap < NRap; ++irap) {
			//if(irap != 4) {
			//	continue;
			//}
			caMass[icent][irap] = new TCanvas("", "", 1920, 1080);
			caMass[icent][irap]->Divide(5, 4);
			for(int ipt = 0; ipt < NPt; ++ipt) {
				TPad* thisPad = (TPad*)caMass[icent][irap]->cd(ipt + 1);
				thisPad->SetLeftMargin(0.15);
				hMass[icent][irap][ipt]->GetYaxis()->SetRangeUser(0, 1.2 * hMass[icent][irap][ipt]->GetMaximum());
				setStyle(hMass[icent][irap][ipt], 24, 1., kBlack, 1, kBlack);
				hMass[icent][irap][ipt]->SetTitle(Form("K_{S}^{0} Inv. Mass Distribution @cent:%s, rap:%.1f~%.1f, pt:%.1f~%.1f;M_{#pi^{+}#pi^{-}};Counts", CentName[icent].c_str(), RapEdge[irap], RapEdge[irap + 1], PtEdge[ipt], PtEdge[ipt + 1]));
				hMass[icent][irap][ipt]->SetStats(0);
				hMass[icent][irap][ipt]->SetTitle("");
				hMass[icent][irap][ipt]->Draw("ep");
				hMassSig[icent][irap][ipt]->SetLineColor(kBlue);
				hMassSig[icent][irap][ipt]->SetFillColor(kBlue);
				hMassSig[icent][irap][ipt]->SetFillStyle(3352);
				hMassSig[icent][irap][ipt]->GetXaxis()->SetRangeUser(0.45, 0.549);
				hMassSig[icent][irap][ipt]->Draw("same hist");
				fitFuncBkg[icent][irap][ipt]->SetLineColor(kBlack);
				fitFuncBkg[icent][irap][ipt]->Draw("same");
				fitFuncSig[icent][irap][ipt]->SetLineColor(kBlue);
				fitFuncSig[icent][irap][ipt]->SetFillColor(kBlue);
				fitFuncSig[icent][irap][ipt]->SetFillStyle(42);
				fitFuncSig[icent][irap][ipt]->Draw("same");
				drawText(0.2, 0.8, 0.05, Form("S/B = %.2f", sig[icent][irap][ipt] / bkg[icent][irap][ipt]), Form("S/#sqrt{S+B} = %.2f", sig[icent][irap][ipt] / sqrt(sig[icent][irap][ipt] + bkg[icent][irap][ipt])), Form("purity = %.2f%%", sig[icent][irap][ipt] / tot[icent][irap][ipt] * 100.));
				drawText(0.2, 0.55, 0.05, Form("mean: %.4f", mean[icent][irap][ipt]), Form("#sigma_{tight}: %.4f", sigma[icent][irap][ipt]), Form("#sigma_{wide}: %.4f", sigma_2[icent][irap][ipt]), Form("#chi^{2}/NDF = %.3f", quality[icent][irap][ipt]));
				drawText(0.6, 0.8, 0.05, "Au+Au @ #sqrt{s_{NN}} = 3 GeV", Form("y:%.1f~%.1f", RapEdge[irap], RapEdge[irap + 1]), Form("p_{T}:%.1f~%.1f GeV/c", PtEdge[ipt], PtEdge[ipt + 1]), Form("Centrality: %s", CentName[icent].c_str()));
				drawXRegion(Mean[2][1] - NSigma * Sigma[2][1], Mean[2][1] + NSigma * Sigma[2][1], thisPad);
				//drawXRegion(0.48, 0.51, thisPad);
			}
			caMass[icent][irap]->Print(Form("%s/invM_cent%d_rap%d.png", cutSetName.c_str(), icent, irap));
		}
	}
	TCanvas* caAcc;
	caAcc = new TCanvas("", "", 2220, 880);
	caAcc->Divide(3, 1);
	for(int icent = 0; icent < NCent; ++icent) {
		caAcc->cd(icent + 1)->SetLogz();
		setXYTitle(hKAcc[icent], "y", "p_{T}");
		hKAcc[icent]->SetStats(0);
		hKAcc[icent]->Draw("colz text");
	}
	TCanvas* caPurity;
	caPurity = new TCanvas("", "", 2220, 880);
	caPurity->Divide(3, 1);
	for(int icent = 0; icent < NCent; ++icent) {
		caPurity->cd(icent + 1);
		setXYTitle(hKPurity[icent], "y", "p_{T}");
		hKPurity[icent]->SetStats(0);
		hKPurity[icent]->Draw("colz text");
	}
	TFile* ofPurity = new TFile(Form("%s/Purity.root", cutSetName.c_str()), "RECREATE");
	ofPurity->cd();
	for(icent = 0; icent < NCent; ++icent) {
		hKPurity[icent]->Write();
	}
	ofPurity->Close();
	delete ofPurity;
	//}}}
}
