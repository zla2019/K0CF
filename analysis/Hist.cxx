#include "Hist.h"
#include <iostream>
#include <ctime>
#include <TMath.h>

void Hist::init()
{
	hCent9 = new TH1F("hCent9", "Cent9 dist.", 10, -1, 9);
	hRefMult = new TH1F("hRefMult", "countrefmult", 1000, 0, 1000);

	hAllKRapPt = new TH2F("hAllKRapPt", "K_{s}^{0} Acc.", 200, -1, 1, 300, 0, 3);
	hKRapPt = new TH2F("hKRapPt", "K_{s}^{0} Acc.", 200, -1, 1, 300, 0, 3);
	hKPipRapPt = new TH2F("hKPipRapPt", "K_{s}^{0} daughter #pi^{+} Acc.", 200, -1, 1, 500, 0, 5);
	hKPimRapPt = new TH2F("hKPimRapPt", "K_{s}^{0} daughter #pi^{-} Acc.", 200, -1, 1, 500, 0, 5);

	hFrt = new TH1F("hFrt", "Freeze out t", 500, 0, 50);
	hFrx = new TH1F("hFrx", "Freeze out x", 800, -40, 40);
	hFry = new TH1F("hFry", "Freeze out y", 800, -40, 40);
	hFrz = new TH1F("hFrz", "Freeze out z", 800, -40, 40);
	hFrxFry = new TH2F("hFrxFry", "Freeze out x vs. y", 800, -40, 40, 800, -40, 40);

	hQinvCorr = new TH2F("hQinvCorr", "q_{inv} vs. corr(crab)", 200, 0, 1, 200, 0, 2);
	hQdotrCorr = new TH2F("hQdotrCorr", "q_{inv}*r vs. corr(crab)", 200, 0, 400, 200, 0, 2);
	hRCorr = new TH2F("hRCorr", "r vs. corr(crab)", 200, 0, 20, 200, 0, 2);

	hLLWeight = new TH2F("hLLWeight", "q_{inv} vs. LL weight;q_{inv};LL weight", 200, 0, 1, 100, -0.5, 0.5);
	hTotWeight = new TH2F("hTotWeight", "q_{inv} vs. Tot weight;q_{inv};Tot weight", 200, 0, 1, 250, 0.5, 2.0);

	for(int icent = 0; icent < 9; ++icent) {
		hSameKqinv[icent] = new TH1F(Form("hSameKqinv_cent%i", icent), Form("K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 500, 0, 1);
		hSameKqinvSI[icent] = new TH1F(Form("hSameKqinvSI_cent%i", icent), Form("K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 500, 0, 1);
		hSameKqinvQS[icent] = new TH1F(Form("hSameKqinvQS_cent%i", icent), Form("K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 500, 0, 1);
		hSameKqinvWoCrab[icent] = new TH1F(Form("hSameKqinvWoCrab_cent%i", icent), Form("K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 500, 0, 1);
		hMixKqinv[icent] = new TH1F(Form("hMixKqinv_cent%i", icent), Form("Mix K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 500, 0, 1);
	}
}

void Hist::Fill(MyTree::Particle& p)
{
	hKRapPt->Fill(p.y, p.pt);
	hFrt->Fill(p.frt);
	hFrx->Fill(p.frx);
	hFry->Fill(p.fry);
	hFrz->Fill(p.frz);
	hFrxFry->Fill(p.frx, p.fry);
	//hKPipRapPt->Fill(p.rapPip, p.ptPip);
	//hKPimRapPt->Fill(p.rapPim, p.ptPim);
}
void Hist::FillSame(float qinv, int cent9)
{
	hSameKqinv[cent9]->Fill(qinv);
}
void Hist::FillMix(float qinv, int cent9)
{
	hMixKqinv[cent9]->Fill(qinv);
}

void Hist::Write(TFile* of)
{
	of->cd();
        hCent9->Write();
	hRefMult->Write();
        hAllKRapPt->Write();
        hKRapPt->Write();
        hKPipRapPt->Write();
        hKPimRapPt->Write();

	hFrt->Write();
	hFrx->Write();
	hFry->Write();
	hFrz->Write();
	hFrxFry->Write();

	hQinvCorr->Write();
	hQdotrCorr->Write();
	hRCorr->Write();

	hLLWeight->Write();
	hTotWeight->Write();

	//CF plots
        for(int icent = 0; icent < 9; ++icent) {
                hSameKqinv[icent]->Write();
                hSameKqinvSI[icent]->Write();
                hSameKqinvQS[icent]->Write();
                hSameKqinvWoCrab[icent]->Write();
                hMixKqinv[icent]->Write();
        }
}
