#include <iostream>
#include "utils.h"
#include "func.h"
#include <math.h>

void Kshort0InvM(std::string inputName, std::string cutSetName, std::string nSigma)
{
	//initial{{{
	const int NCent9 = 9;
	const int NRap = 8;
	const int NPt = 20;
	const int NSigma = std::stoi(nSigma);
	const float RapEdge[NRap + 1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6};
	const float PtEdge[NPt + 1] = { 0.0, 0.10, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0 };
	const float Mean[3][4] = { { 0.4979, 0.4981, 0.4979, 0.4949 }, { 0.4978, 0.4981, 0.4979, 0.4979 }, { 0.4979, 0.4981, 0.4980, 0.4980 } };
	const float Sigma[3][4] = { { 0.0027, 0.0035, 0.0033, 0.0033 }, { 0.0027, 0.0034, 0.0036, 0.0036 }, { 0.0027, 0.0035, 0.0039, 0.0039 } };

	TFile* ifPlots = new TFile(inputName.c_str());
	TH3F* hMassRapPt[NCent9];
	TH3F* hRotMassRapPt[NCent9];
	TH1F* hMass[NCent9][NRap][NPt];
	TH1F* hRotMass[NCent9][NRap][NPt];

	TH2F* hPurity[NCent9];
	TH2F* hAcc[NCent9];
	//}}}

	//get plots{{{
	for(int icent = 0; icent < NCent9; ++icent) {
		hMassRapPt[icent] = (TH3F*)getCopy(ifPlots, Form("hMassCascadeRapidityvsPt_cent%i", icent));
		hRotMassRapPt[icent] = (TH3F*)getCopy(ifPlots, Form("hMassRotCascadeRapidityvsPt_cent%i", icent));

		hPurity[icent] = new TH2F(Form("hPurity_cent%d", icent), Form("Purity cent%d", icent), NRap, RapEdge[0], RapEdge[NRap], NPt, PtEdge[0], PtEdge[NPt]);
		hAcc[icent] = new TH2F(Form("hAcc_cent%d", icent), Form("Acc cent%d", icent), NRap, RapEdge[0], RapEdge[NRap], NPt, PtEdge[0], PtEdge[NPt]);
	}
	for(int icent = 0; icent < NCent9; ++icent) {
		for(int irap = 0; irap < NRap; ++irap) {
			for(int ipt = 0; ipt < NPt; ++ipt) {
				hMass[icent][irap][ipt] = projectionX(hMassRapPt[icent], RapEdge[irap], RapEdge[irap + 1], PtEdge[ipt], PtEdge[ipt + 1], Form("hMass_cent%i_rap%i_pt%i", icent, irap, ipt));
				hRotMass[icent][irap][ipt] = projectionX(hRotMassRapPt[icent], RapEdge[irap], RapEdge[irap + 1], PtEdge[ipt], PtEdge[ipt + 1], Form("hRotMass_cent%i_rap%i_pt%i", icent, irap, ipt));
				float mass = getIntegral(hMass[icent][irap][ipt], 0.42, 0.48) + getIntegral(hMass[icent][irap][ipt], 0.52, 0.58);
				float rotMass = getIntegral(hRotMass[icent][irap][ipt], 0.42, 0.48) + getIntegral(hRotMass[icent][irap][ipt], 0.52, 0.58);
				hRotMass[icent][irap][ipt]->Scale(mass / rotMass);
				float tot = getIntegral(hMass[icent][irap][ipt], Mean[2][1] - (NSigma * Sigma[2][1]), Mean[2][1] + (NSigma * Sigma[2][1]));
				float bkg = getIntegral(hRotMass[icent][irap][ipt], Mean[2][1] - (NSigma * Sigma[2][1]), Mean[2][1] + (NSigma * Sigma[2][1]));
				float sig = tot - bkg;
				float purity = sig / tot;
				if(tot == 0 || isnan(purity) || isnan(sig) || isnan(tot) || isnan(bkg) || sig < 0) {
					sig = 0;
					tot = 0;
					purity = 0;
				}
				hAcc[icent]->SetBinContent(irap + 1, ipt + 1, sig);
				hPurity[icent]->SetBinContent(irap + 1, ipt + 1, purity);
			}
		}
	}
	//}}}

	//plotting {{{
	TCanvas* caMass = new TCanvas("caMass", "caMass", 1920, 1080);
	gSystem->Unlink("caMass.gif");
	caMass->Divide(5, 4, 0, 0);
	for(int icent = 0; icent < NCent9; ++icent) {
		for(int irap = 0; irap < NRap; ++irap) {
			for(int ipt = 0; ipt < NPt; ++ipt) {
				caMass->cd(ipt + 1);
				setStyle(hMass[icent][irap][ipt], 20, 0.8, kBlack, 1, kBlack);
				setStyle(hRotMass[icent][irap][ipt], 24, 1., kRed, 1, kRed);
				hMass[icent][irap][ipt]->Draw("ep");
				hRotMass[icent][irap][ipt]->Draw("same ep");
				hMass[icent][irap][ipt]->Draw("same ep");
			}
			caMass->Print("caMass.gif+40");
		}
	}

	TFile* ofPlots = new TFile("Purity.root", "RECREATE");
	ofPlots->cd();
	for(int icent = 0; icent < NCent9; ++icent) {
		hPurity[icent]->Write();
		hAcc[icent]->Write();
	}
	//}}}
}
