#include "Hist.h"
#include <iostream>
#include <ctime>
#include <TMath.h>

void Hist::init()
{
	hSameKPDG = new TH1F("hSameKPDG", "Kaon pdg distribution", 5, 280, 340);
	hVz = new TH1F("hVz", "V_{z} distribution;cm;cnts", 100, 198, 202);
	hVr = new TH2F("hVr", "V_{r} distribution;cm;cm", 100, -2, 2, 100, -4, 0);
	hCent9 = new TH1F("hCent9", "Cent9 dist.", 10, -1, 9);
	hRefMult = new TH1F("hRefMult", "countrefmult", 400, 0, 400);

	hCosTheta = new TH1F("hCosTheta", "cos(#theta) distribution", 50, 0.95, 1.0);
	hDecayLength = new TH1F("hDecayLength", "K^{0}_{s} decay length", 100, 0, 20);
	hDgDCA = new TH1F("hDgDCA", "K_{s}^{0} daughter DCA", 100, 0, 2);
	hDCA = new TH1F("hDCA", "K_{s}^{0} DCA", 100, 0, 2);

	hDedx = new TH2F("hDedx", "dEdx vs p*q;p*q;dEdx", 1000, -5, 5, 500, 1.5, 6.5);
	hMass2 = new TH2F("hMass2", "m^{2} vs p*q;p*q;m^{2}", 1000, -5, 5, 500, -2, 3);

	hAllKRapPt = new TH2F("hAllKRapPt", "K_{s}^{0} Acc.", 200, -1, 1, 300, 0, 3);
	hSameKRapPt = new TH2F("hSameKRapPt", "K_{s}^{0} Acc.", 200, -1, 1, 300, 0, 3);
	hSameKMass = new TH1F("hSameKMass", "K_{s}^{0} mass distribution;M_{inv};cnts", 160, 0.48, 0.52);
	hSameKPhi = new TH1F("hSameKPhi", "K_{s}^{0} #phi distribution;#phi;cnts", 100, -TMath::Pi(), TMath::Pi());
	hSameKPipRapPt = new TH2F("hSameKPipRapPt", "K_{s}^{0} daughter #pi^{+} Acc.", 200, -1, 1, 500, 0, 5);
	hSameKPimRapPt = new TH2F("hSameKPimRapPt", "K_{s}^{0} daughter #pi^{-} Acc.", 200, -1, 1, 500, 0, 5);
	hDaughterPipDCA = new TH1F("hDaughterPipDCA", "K_{s}^{0} daughter #pi^{+} DCA distribution", 50, 0, 10);
	hDaughterPimDCA = new TH1F("hDaughterPimDCA", "K_{s}^{0} daughter #pi^{-} DCA distribution", 50, 0, 10);

	hLeftSideSameKRapPt = new TH2F("hLeftSideSameKRapPt", "K_{s}^{0} Acc.", 200, -1, 1, 300, 0, 3);
	hLeftSideSameKMass = new TH1F("hLeftSideSameKMass", "K_{s}^{0} mass distribution;M_{inv};cnts", 1200, 0.4, 0.70);
	hLeftSideSameKPhi = new TH1F("hLeftSideSameKPhi", "K_{s}^{0} #phi distribution;#phi;cnts", 100, -TMath::Pi(), TMath::Pi());
	hLeftSideSameKPipRapPt = new TH2F("hLeftSideSameKPipRapPt", "K_{s}^{0} daughter #pi^{+} Acc.", 200, -1, 1, 500, 0, 5);
	hLeftSideSameKPimRapPt = new TH2F("hLeftSideSameKPimRapPt", "K_{s}^{0} daughter #pi^{-} Acc.", 200, -1, 1, 500, 0, 5);

	hRightSideSameKRapPt = new TH2F("hRightSideSameKRapPt", "K_{s}^{0} Acc.", 200, -1, 1, 300, 0, 3);
	hRightSideSameKMass = new TH1F("hRightSideSameKMass", "K_{s}^{0} mass distribution;M_{inv};cnts", 1200, 0.4, 0.70);
	hRightSideSameKPhi = new TH1F("hRightSideSameKPhi", "K_{s}^{0} #phi distribution;#phi;cnts", 100, -TMath::Pi(), TMath::Pi());
	hRightSideSameKPipRapPt = new TH2F("hRightSideSameKPipRapPt", "K_{s}^{0} daughter #pi^{+} Acc.", 200, -1, 1, 500, 0, 5);
	hRightSideSameKPimRapPt = new TH2F("hRightSideSameKPimRapPt", "K_{s}^{0} daughter #pi^{-} Acc.", 200, -1, 1, 500, 0, 5);

	hPipNSigma = new TH1F("hPipNSigma", "#pi^{+} n#sigma distribution", 100, -5, 5);
	hPimNSigma = new TH1F("hPimNSigma", "#pi^{-} n#sigma distribution", 100, -5, 5);
	hPipNSigma2 = new TH1F("hPipNSigma2", "#pi^{+} n#sigma distribution", 100, -5, 5);
	hPimNSigma2 = new TH1F("hPimNSigma2", "#pi^{-} n#sigma distribution", 100, -5, 5);

	hDTheta = new TH1F("hDTheta", "K_{s}^{0}-K_{s}^{0} #Delta#theta dist.", 314, 0, TMath::Pi());
	hDPhi = new TH1F("hDPhi", "K_{s}^{0}-K_{s}^{0} #Delta#phi dist.", 314, 0, TMath::Pi());
	hDThetaDPhi = new TH2F("hDThetaDPhi", "K_{s}^{0}-K_{s}^{0} #DeltaTheta vs. #DeltaPhi", 314, 0, TMath::Pi(), 314, 0, TMath::Pi());

	for(int icent = 0; icent < 9; ++icent) {
		hSameKPtRapMass[icent] = new TH3F(Form("hMassCascadeRapidityvsPt_cent%d", icent),Form("hMassCascadeRapidityvsPt_cent%d", icent),160, 0.42, 0.58, 20, -1, 1, 30, 0, 3);
		hSameKqinv[icent] = new TH1F(Form("hSameKqinv_cent%i", icent), Form("K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 500, 0, 1);
		hMixKqinv[icent] = new TH1F(Form("hMixKqinv_cent%i", icent), Form("Mix K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 500, 0, 1);

		hMix0000Qinv[icent] = new TH1F(Form("hMix0000Qinv_cent%d", icent), "", 500, 0, 1);
		hMix0001Qinv[icent] = new TH1F(Form("hMix0001Qinv_cent%d", icent), "", 500, 0, 1);
		hMix0011Qinv[icent] = new TH1F(Form("hMix0011Qinv_cent%d", icent), "", 500, 0, 1);
		hMix0023Qinv[icent] = new TH1F(Form("hMix0023Qinv_cent%d", icent), "", 500, 0, 1);
		hMix0101Qinv[icent] = new TH1F(Form("hMix0101Qinv_cent%d", icent), "", 500, 0, 1);
		hMix0102Qinv[icent] = new TH1F(Form("hMix0102Qinv_cent%d", icent), "", 500, 0, 1);
		hMix0123Qinv[icent] = new TH1F(Form("hMix0123Qinv_cent%d", icent), "", 500, 0, 1);

		for(int icase = 0; icase < 4; ++icase) {
			hMixKqinvWeight[icent][icase] = new TH1F(Form("hMixKqinvWeight_cent%d_case%d", icent, icase), "", 500, 0, 1);
		}
	}

	hMixBkgMass = new TH1F("hMixBkgMass", "", 160, 0.48, 0.52);

	//cut plots
	hCutChi2Topo = new TH1F("hCutChi2Topo", "#chi^{2}_{topo} after cut", 50, 0, 10);
	hCutChi2NDF = new TH1F("hCutChi2NDF", "#chi^{2}_{NDF} after cut", 50, 0, 10);
	hCutChi2PrimA = new TH1F("hCutChi2PrimA", "#chi^{2}_{prim} A after cut", 50, 0, 100);
	hCutChi2PrimB = new TH1F("hCutChi2PrimB", "#chi^{2}_{prim} B after cut", 50, 0, 100);
	hCutNHitsA = new TH1F("hCutNHitsA", "NHitsFit A after cut", 500, 0, 500);
	hCutNHitsB = new TH1F("hCutNHitsB", "NHitsFit B after cut", 500, 0, 500);
	hCutDCAA = new TH1F("hCutDCAA", "DCA A after cut", 50, 0, 10);
	hCutDCAB = new TH1F("hCutDCAB", "DCA B after cut", 50, 0, 10);
	hCutDecayLength = new TH1F("hCutDecayLength", "DecayLength after cut", 50, 0, 10);
	hEtaPt = new TH2F("hEtaPt", "#eta vs. p_{T}", 270, -2.5, 0.2, 300, 0, 3);
}

void Hist::FillAll(MyTree::Particle& p, int cent9)
{
	hDaughterPipDCA->Fill(p.dcaA);
	hDaughterPimDCA->Fill(p.dcaB);
	hDecayLength->Fill(p.decayLength);
	hDgDCA->Fill(p.dgDCA);
	hDCA->Fill(p.dca);
	hSameKPtRapMass[cent9]->Fill(p.mass, p.rap, p.pt);
}

void Hist::Fill(MyTree::Particle& p)
{
	hSameKRapPt->Fill(p.rap, p.pt);
	hSameKPDG->Fill(p.pdg);
	hSameKMass->Fill(p.mass);
	hSameKPipRapPt->Fill(p.rapPip, p.ptPip);
	hSameKPimRapPt->Fill(p.rapPim, p.ptPim);
	hSameKPhi->Fill(p.phi);
	hDedx->Fill(p.pA, p.dEdxA);
	hDedx->Fill(-p.pB, p.dEdxB);
	hMass2->Fill(p.pA, p.m2A);
	hMass2->Fill(-p.pB, p.m2B);
	hPipNSigma2->Fill(p.nSigmaA);
	hPimNSigma2->Fill(p.nSigmaB);
}
void Hist::FillLeft(MyTree::Particle& p)
{
	hLeftSideSameKMass->Fill(p.mass);
	hLeftSideSameKPipRapPt->Fill(p.rapPip, p.ptPip);
	hLeftSideSameKPimRapPt->Fill(p.rapPim, p.ptPim);
	hLeftSideSameKRapPt->Fill(p.rap, p.pt);
	hLeftSideSameKPhi->Fill(p.phi);
}
void Hist::FillRight(MyTree::Particle& p)
{
	hRightSideSameKMass->Fill(p.mass);
	hRightSideSameKPipRapPt->Fill(p.rapPip, p.ptPip);
	hRightSideSameKPimRapPt->Fill(p.rapPim, p.ptPim);
	hRightSideSameKRapPt->Fill(p.rap, p.pt);
	hRightSideSameKPhi->Fill(p.phi);
}

void Hist::FillCut(MyTree::Particle& p)
{
	hCutChi2Topo->Fill(p.chi2Topo);
	hCutChi2NDF->Fill(p.chi2NDF);
	hCutChi2PrimA->Fill(p.chi2PrimPip);
	hCutChi2PrimB->Fill(p.chi2PrimPim);
	hCutNHitsA->Fill(p.nHitsA);
	hCutNHitsB->Fill(p.nHitsB);
	hCutDCAA->Fill(p.dcaA);
	hCutDCAB->Fill(p.dcaB);
	hCutDecayLength->Fill(p.decayLength);
	hEtaPt->Fill(p.etaB, p.ptPim);
}

void Hist::Write(TFile* of)
{
	of->cd();
        hVz->Write();
        hVr->Write();
        hCent9->Write();
	hRefMult->Write();
        hSameKPDG->Write();
        hSameKMass->Write();
        hSameKPhi->Write();
        hAllKRapPt->Write();
        hSameKRapPt->Write();
        hSameKPipRapPt->Write();
        hSameKPimRapPt->Write();

        hLeftSideSameKMass->Write();
        hLeftSideSameKPhi->Write();
        hLeftSideSameKRapPt->Write();
        hLeftSideSameKPipRapPt->Write();
        hLeftSideSameKPimRapPt->Write();
        hRightSideSameKMass->Write();
        hRightSideSameKPhi->Write();
        hRightSideSameKRapPt->Write();
        hRightSideSameKPipRapPt->Write();
        hRightSideSameKPimRapPt->Write();

        hDaughterPipDCA->Write();
        hDaughterPimDCA->Write();
        hCosTheta->Write();
        hDecayLength->Write();
        hDgDCA->Write();
        hDCA->Write();
        hDedx->Write();
        hMass2->Write();
        hPipNSigma->Write();
        hPimNSigma->Write();
        hPipNSigma2->Write();
        hPimNSigma2->Write();

	hMixBkgMass->Write();

	//cut plots
	hCutChi2Topo->Write();
	hCutChi2NDF->Write();
	hCutChi2PrimA->Write();
	hCutChi2PrimB->Write();
	hCutNHitsA->Write();
	hCutNHitsB->Write();
	hCutDCAA->Write();
	hCutDCAB->Write();
	hCutDecayLength->Write();
	hEtaPt->Write();

	//CF plots
        for(int icent = 0; icent < 9; ++icent) {
                hSameKPtRapMass[icent]->Write();
                hSameKqinv[icent]->Write();
                hMixKqinv[icent]->Write();

		hMix0000Qinv[icent]->Write();
		hMix0001Qinv[icent]->Write();
		hMix0011Qinv[icent]->Write();
		hMix0023Qinv[icent]->Write();
		hMix0101Qinv[icent]->Write();
		hMix0102Qinv[icent]->Write();
		hMix0123Qinv[icent]->Write();

		for(int icase = 0; icase < 4; ++icase) {
			hMixKqinvWeight[icent][icase]->Write();
		}
        }
}
