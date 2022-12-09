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

	hPipNSigma = new TH1F("hPipNSigma", "#pi^{+} n#sigma distribution", 100, -5, 5);
	hPimNSigma = new TH1F("hPimNSigma", "#pi^{-} n#sigma distribution", 100, -5, 5);
	hPipNSigma2 = new TH1F("hPipNSigma2", "#pi^{+} n#sigma distribution", 100, -5, 5);
	hPimNSigma2 = new TH1F("hPimNSigma2", "#pi^{-} n#sigma distribution", 100, -5, 5);

	hDTheta = new TH1F("hDTheta", "K_{s}^{0}-K_{s}^{0} #Delta#theta dist.", 314, 0, TMath::Pi());
	hDPhi = new TH1F("hDPhi", "K_{s}^{0}-K_{s}^{0} #Delta#phi dist.", 314, 0, TMath::Pi());
	hDThetaDPhi = new TH2F("hDThetaDPhi", "K_{s}^{0}-K_{s}^{0} #DeltaTheta vs. #DeltaPhi", 314, 0, TMath::Pi(), 314, 0, TMath::Pi());

	for(int icent = 0; icent < 9; ++icent) {
		hSameKPtRapMass[icent] = new TH3F(Form("hMassCascadeRapidityvsPt_cent%d", icent),Form("hMassCascadeRapidityvsPt_cent%d", icent),160, 0.42, 0.58, 20, -1, 1, 30, 0, 3);
		hRotKPtRapMass[icent] = new TH3F(Form("hMassRotCascadeRapidityvsPt_cent%d", icent),Form("hMassRotCascadeRapidityvsPt_cent%d", icent),160, 0.42, 0.58, 20, -1, 1, 30, 0, 3);
		hSameKqinv[icent] = new TH1F(Form("hSameKqinv_cent%i", icent), Form("K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 500, 0, 1);
		hMixKqinv[icent] = new TH1F(Form("hMixKqinv_cent%i", icent), Form("Mix K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 500, 0, 1);
		hSRQinv[icent] = new TH1F(Form("hSRQinv_cent%d", icent), Form("K_{s}^{0} SR q_{inv} @cent%d", icent), 500, 0, 1);
		hMixSRQinv[icent] = new TH1F(Form("hMixSRQinv_cent%d", icent), Form("K_{s}^{0} Mix SR q_{inv} @cent%d", icent), 500, 0, 1);
		hRSQinv[icent] = new TH1F(Form("hRSQinv_cent%d", icent), Form("K_{s}^{0} RS q_{inv} @cent%d", icent), 500, 0, 1);
		hMixRSQinv[icent] = new TH1F(Form("hMixRSQinv_cent%d", icent), Form("K_{s}^{0} Mix RS q_{inv} @cent%d", icent), 500, 0, 1);
		hRRQinv[icent] = new TH1F(Form("hRRQinv_cent%d", icent), Form("K_{s}^{0} RR q_{inv} @cent%d", icent), 500, 0, 1);
		hMixRRQinv[icent] = new TH1F(Form("hMixRRQinv_cent%d", icent), Form("K_{s}^{0} Mix RR q_{inv} @cent%d", icent), 500, 0, 1);

		hSameKqlong[icent] = new TH1F(Form("hSameKqlong_cent%i", icent), Form("K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 2000, -2, 2);
		hMixKqlong[icent] = new TH1F(Form("hMixKqlong_cent%i", icent), Form("Mix K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 2000, -2, 2);
		hSameKqout[icent] = new TH1F(Form("hSameKqout_cent%i", icent), Form("K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 2000, -2, 2);
		hMixKqout[icent] = new TH1F(Form("hMixKqout_cent%i", icent), Form("Mix K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 2000, -2, 2);
		hSameKqside[icent] = new TH1F(Form("hSameKqside_cent%i", icent), Form("K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 2000, -2, 2);
		hMixKqside[icent] = new TH1F(Form("hMixKqside_cent%i", icent), Form("Mix K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 2000, -2, 2);
		for(int icase = 0; icase < 4; ++icase) {
			hMixKqinvWeight[icent][icase] = new TH1F(Form("hMixKqinvWeight_cent%i_case%i", icent, icase), Form("Weighted Mix K_{s}^{0} q_{inv} for SideBand @cent%i, case%i", icent, icase), 500, 0, 1);
		}
	}

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
	//hSameKPtRapMass[cent9]->Fill(p.mass, p.rap, p.pt);
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

void Hist::Fill(float qinv, int cent9)
{
	hSameKqinv[cent9]->Fill(qinv);
}

void Hist::FillMix(float qinv, int cent9, float sideBandWeight[])
{
	hMixKqinv[cent9]->Fill(qinv);
	hMixKqinvWeight[cent9][0]->Fill(qinv, sideBandWeight[0]);
	hMixKqinvWeight[cent9][1]->Fill(qinv, sideBandWeight[1]);
	hMixKqinvWeight[cent9][2]->Fill(qinv, sideBandWeight[2]);
	hMixKqinvWeight[cent9][3]->Fill(qinv, sideBandWeight[3]);
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

        //hDTheta->Write();
        //hDPhi->Write();
        //hDThetaDPhi->Write();

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
                hRotKPtRapMass[icent]->Write();
                hSameKqinv[icent]->Write();
                hMixKqinv[icent]->Write();
		hSRQinv[icent]->Write();
		hMixSRQinv[icent]->Write();
		hRSQinv[icent]->Write();
		hMixRSQinv[icent]->Write();
		hRRQinv[icent]->Write();
		hMixRRQinv[icent]->Write();
		for(int icase = 0; icase < 4; ++icase) {
			hMixKqinvWeight[icent][icase]->Write();
		}
        }
}
