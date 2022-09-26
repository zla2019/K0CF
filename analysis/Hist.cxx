#include "Hist.h"
#include <iostream>
#include <ctime>
#include <TMath.h>

void Hist::init()
{
	hCent9 = new TH1F("hCent9", "Cent9 dist.", 10, -1, 9);
	hRefMult = new TH1F("hRefMult", "countrefmult", 400, 0, 400);

	hAllKRapPt = new TH2F("hAllKRapPt", "K_{s}^{0} Acc.", 200, -1, 1, 300, 0, 3);
	hKRapPt = new TH2F("hKRapPt", "K_{s}^{0} Acc.", 200, -1, 1, 300, 0, 3);
	hKPipRapPt = new TH2F("hKPipRapPt", "K_{s}^{0} daughter #pi^{+} Acc.", 200, -1, 1, 500, 0, 5);
	hKPimRapPt = new TH2F("hKPimRapPt", "K_{s}^{0} daughter #pi^{-} Acc.", 200, -1, 1, 500, 0, 5);

	for(int icent = 0; icent < 9; ++icent) {
		hSameKqinv[icent] = new TH1F(Form("hSameKqinv_cent%i", icent), Form("K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 500, 0, 1);
		hMixKqinv[icent] = new TH1F(Form("hMixKqinv_cent%i", icent), Form("Mix K_{s}^{0} q_{inv}(k^{*}) @cent%i", icent), 500, 0, 1);
	}
}

void Hist::Fill(MyTree::Particle& p)
{
	hKRapPt->Fill(p.y, p.pt);
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

	//CF plots
        for(int icent = 0; icent < 9; ++icent) {
                hSameKqinv[icent]->Write();
                hMixKqinv[icent]->Write();
        }
}
