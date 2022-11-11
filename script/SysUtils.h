#ifndef SYSUTILS_H
#define SYSUTILS_H
#include <iostream>
#include "utils.h"
#include "CFutils.h"
#include "CFLLFunc.C"

class SysSource
{
public:
	SysSource() {};
	~SysSource();
	std::vector<TFile*> mSysFile;
	void pushBack(std::string ifName);
	void processCF();
	void processStat(SysSource* sysDefault);
	void processSys(SysSource* sysDefault);
	void processTotSys(SysSource* sysDefault);
	void processSourceSys(SysSource* sysDefault, std::string mode = "Avg");
	void SysAvg(SysSource* sysDefault);
	void SysMax(SysSource* sysDefault);
	void FitAll(float sp);
	void getSysErr(SysSource* sysDefault);
	float getLLPar(int icent, int isys, int ipar);
	float getLLParErr(int icent, int isys, int ipar);
	float getGausPar(int icent, int isys, int ipar);
	float getGausParErr(int icent, int isys, int ipar);
	void close();

	float ssSysErr[3];
	float ldSysErr[3];
	float ssGausSysErr[3];
	float ldGausSysErr[3];

	std::vector<TH1F*> mCF[3];	//3 centrality bin
	std::vector<TF1*> mCFLLFit[3];
	std::vector<TF1*> mCFGausFit[3];
	std::vector<TH1F*> mCFStat[3];
	std::vector<TH1F*> mCFSys[3];
	std::vector<TH1F*> mCFTotSys[3];
	std::vector<TH1F*> mCFRelStat[3];
	std::vector<TH1F*> mCFRelSys[3];
	std::vector<TH1F*> mCFRelTot[3];
	TH1F* hCFRel[3];
	TH1F* hCFSource[3];
private:
};

//function defination
void SysSource::pushBack(std::string ifName)
{
	TFile* iftmp = TFile::Open(ifName.c_str(), "READ");
	mSysFile.push_back(std::move(iftmp));
}

SysSource::~SysSource()
{
	for(int isys = 0; isys < mSysFile.size(); ++isys) {
		mSysFile[isys]->Close();
		delete mSysFile[isys];
		for(int icent = 0; icent < 3; ++icent) {
			delete mCF[icent][isys];
			delete mCFStat[icent][isys];
			delete mCFSys[icent][isys];
			delete mCFTotSys[icent][isys];
			delete mCFRelStat[icent][isys];
			delete mCFRelSys[icent][isys];
			delete mCFRelTot[icent][isys];
			delete mCFLLFit[icent][isys];
			delete mCFGausFit[icent][isys];
		}
	}
	mSysFile.clear();
	for(int icent = 0; icent < 3; ++icent) {
		mCF[icent].clear();
		mCFStat[icent].clear();
		mCFSys[icent].clear();
		mCFTotSys[icent].clear();
		mCFRelStat[icent].clear();
		mCFRelSys[icent].clear();
		mCFRelTot[icent].clear();
		mCFLLFit[icent].clear();
		mCFGausFit[icent].clear();
		delete hCFRel[icent];
		delete hCFSource[icent];
	}
}

void SysSource::close()
{
	for(int isys = 0; isys < mSysFile.size(); ++isys) {
		mSysFile[isys]->Close();
		delete mSysFile[isys];
	}
	mSysFile.clear();
}

void SysSource::processCF()
{
	for(int isys = 0; isys < mSysFile.size(); ++isys) {
		TH1F* hCFRaw[3] = { 0 };
		TH1F* hCFMisid[3] = { 0 };
		TH1F* hPairPurity[3] = { 0 };
		TH1F* hPure[3] = { 0 };
		RawCF(mSysFile[isys], hCFRaw);
		MisidCF(mSysFile[isys], hCFMisid);
		PairPurity(mSysFile[isys], hPairPurity);
		PureCF(hCFRaw, hCFMisid, hPure, hPairPurity);
		mCF[0].push_back(hPure[0]);
		mCF[1].push_back(hPure[1]);
		mCF[2].push_back(hPure[2]);
	}
}

void SysSource::processStat(SysSource* sysDefault)
{
	for(int isys = 0; isys < mSysFile.size(); ++isys) {
		TH1F* hStat[3];
		TH1F* hRelStat[3];
		for(int icent = 0; icent < 3; ++icent) {
			hStat[icent] = (TH1F*)mCF[icent][isys]->Clone(Form("hStat_cent%d", icent));
			hRelStat[icent] = (TH1F*)mCF[icent][isys]->Clone(Form("hRelStat_cent%d", icent));
			for(int ibin = 0; ibin < hStat[icent]->GetNbinsX(); ++ibin) {
				float statErr = hStat[icent]->GetBinError(ibin + 1);
				float statDef = sysDefault->mCF[icent][0]->GetBinError(ibin + 1);
				statErr = sqrt(fabs(statDef*statDef - statErr*statErr));
				hStat[icent]->SetBinError(ibin + 1, statErr);
				hRelStat[icent]->SetBinContent(ibin + 1, statErr / sysDefault->mCF[icent][0]->GetBinContent(ibin + 1));
				hRelStat[icent]->SetBinError(ibin + 1, 0);
			}
			mCFStat[icent].push_back(hStat[icent]);
			mCFRelStat[icent].push_back(hRelStat[icent]);
		}
	}
}

void SysSource::processSys(SysSource* sysDefault)
{
	for(int isys = 0; isys < mSysFile.size(); ++isys) {
		TH1F* hSys[3];
		TH1F* hRelSys[3];
		for(int icent = 0; icent < 3; ++icent) {
			hSys[icent] = (TH1F*)mCF[icent][isys]->Clone(Form("hSys_cent%d", icent));
			hRelSys[icent] = (TH1F*)mCF[icent][isys]->Clone(Form("hRelSys_cent%d", icent));
			for(int ibin = 0; ibin < hSys[icent]->GetNbinsX(); ++ibin) {
				float sys = hSys[icent]->GetBinContent(ibin + 1);
				float sysDef = sysDefault->mCF[icent][0]->GetBinContent(ibin + 1);
				sys = fabs(sys - sysDef);
				hSys[icent]->SetBinError(ibin + 1, sys);
				hRelSys[icent]->SetBinContent(ibin + 1, sys / sysDefault->mCF[icent][0]->GetBinContent(ibin + 1));
				hRelSys[icent]->SetBinError(ibin + 1, 0);
			}
			mCFSys[icent].push_back(hSys[icent]);
			mCFRelSys[icent].push_back(hRelSys[icent]);
		}
	}
}

void SysSource::processTotSys(SysSource* sysDefault)
{
	for(int isys = 0; isys < mSysFile.size(); ++isys) {
		TH1F* hTotSys[3];
		TH1F* hRelTot[3];
		for(int icent =0; icent < 3; ++icent) {
			hTotSys[icent] = (TH1F*)mCF[icent][isys]->Clone(Form("hTotSys_cent%d", icent));
			hRelTot[icent] = (TH1F*)mCF[icent][isys]->Clone(Form("hRelTot_cent%d", icent));
			for(int ibin = 0; ibin < hTotSys[icent]->GetNbinsX(); ++ibin) {
				float totSys = mCFStat[icent][isys]->GetBinError(ibin + 1) > mCFSys[icent][isys]->GetBinError(ibin + 1) ? 0 : mCFSys[icent][isys]->GetBinError(ibin + 1) - mCFStat[icent][isys]->GetBinError(ibin + 1);
				hTotSys[icent]->SetBinError(ibin + 1, totSys);
				//hRelTot[icent]->SetBinContent(ibin + 1, totSys / hRelTot[icent]->GetBinContent(ibin + 1));
				hRelTot[icent]->SetBinContent(ibin + 1, totSys / sysDefault->mCF[icent][0]->GetBinContent(ibin + 1));
				hRelTot[icent]->SetBinError(ibin + 1, 0);
			}
			mCFTotSys[icent].push_back(hTotSys[icent]);
			mCFRelTot[icent].push_back(hRelTot[icent]);
		}
	}
}

void SysSource::processSourceSys(SysSource* sysDefault, std::string mode)
{
	if(mode == "Avg") SysAvg(sysDefault);
	else if(mode == "Max") SysMax(sysDefault);
}

void SysSource::SysAvg(SysSource* sysDefault)
{
	for(int icent = 0; icent < 3; ++icent) {
		hCFSource[icent] = (TH1F*)mCFTotSys[icent][0]->Clone("hSysSource");
		hCFRel[icent] = (TH1F*)mCFTotSys[icent][0]->Clone("hRelSysSource");
		for(int ibin = 0; ibin < hCFSource[icent]->GetNbinsX(); ++ibin) {
			float sysSquare = 0;
			for(int isys = 0; isys < mSysFile.size(); ++isys) {
				float sys = mCFTotSys[icent][isys]->GetBinError(ibin + 1);
				sysSquare += sys*sys;
			}
			hCFSource[icent]->SetBinError(ibin + 1, sqrt(sysSquare / mSysFile.size()));
			hCFRel[icent]->SetBinContent(ibin + 1, sqrt(sysSquare / mSysFile.size()) / sysDefault->mCF[icent][0]->GetBinContent(ibin + 1));
			hCFRel[icent]->SetBinError(ibin + 1, 0);
		}
	}
}

void SysSource::SysMax(SysSource* sysDefault)
{
	for(int icent = 0; icent < 3; ++icent) {
		hCFSource[icent] = (TH1F*)mCFTotSys[icent][0]->Clone("hSysSource");
		hCFRel[icent] = (TH1F*)mCFTotSys[icent][0]->Clone("hRelSysSource");
		for(int ibin = 0; ibin < hCFSource[icent]->GetNbinsX(); ++ibin) {
			float sys = 0;
			for(int isys = 0; isys < mSysFile.size(); ++isys) {
				sys = mCFTotSys[icent][isys]->GetBinError(ibin + 1) > sys ? mCFTotSys[icent][isys]->GetBinError(ibin + 1) : sys;
			}
			hCFSource[icent]->SetBinError(ibin + 1, sys);
			hCFRel[icent]->SetBinContent(ibin + 1, sys / sysDefault->mCF[icent][0]->GetBinContent(ibin + 1));
			hCFRel[icent]->SetBinError(ibin + 1, 0);
		}
	}
}

void SysSource::FitAll(float sp)
{
	for(int icent = 0; icent < 3; ++icent) {
		for(int isys = 0; isys < mSysFile.size(); ++isys) {
			TF1* fLLTmp = new TF1(Form("fLL_%dsys", isys), CFLL, 0, 0.4, 2);
			TF1* fGausTmp = new TF1(Form("fGaus_%dsys", isys), "1 + [0] * exp((-1 * [1]^2 * (x)^2) / 0.038937929230)", 0, 1);
			fLLTmp->SetParameters(0.5, 2.5);
			fLLTmp->SetParLimits(0, 0.2, 1.5);
			fLLTmp->SetParLimits(1, 1.2, 7);
			fGausTmp->SetParameters(0.5, 2.5);
			fGausTmp->SetParLimits(0, 0.2, 1.5);
			fGausTmp->SetParLimits(1, 1.2, 7);
			float start = mCF[icent][isys]->GetBinLowEdge(sp);
			float end = mCF[icent][isys]->GetBinLowEdge(15);
			mCF[icent][isys]->Fit(fLLTmp, "RNQ", "", start, end);
			mCF[icent][isys]->Fit(fGausTmp, "NQ", "", start, 1);
			mCFLLFit[icent].push_back(fLLTmp);
			mCFGausFit[icent].push_back(fGausTmp);
		}
	}
}

void SysSource::getSysErr(SysSource* sysDefault)
{
	for(int icent = 0; icent < 3; ++icent) {
		ssSysErr[icent] = 0;
		ldSysErr[icent] = 0;
		for(int isys = 0; isys < mSysFile.size(); ++isys) {
			float ssStat = sqrt(fabs(getLLParErr(icent, isys, 1)*getLLParErr(icent, isys, 1)
				 - sysDefault->getLLParErr(icent, 0, 1)*sysDefault->getLLParErr(icent, 0, 1)));
			float ssSys = fabs(getLLPar(icent, isys, 1) - sysDefault->getLLPar(icent, 0, 1));
			float ssTotSys = ssStat > ssSys ? 0 : ssSys - ssStat;

			float ldStat = sqrt(fabs(getLLParErr(icent, isys, 0)*getLLParErr(icent, isys, 0)
				 - sysDefault->getLLParErr(icent, 0, 0)*sysDefault->getLLParErr(icent, 0, 0)));
			float ldSys = fabs(getLLPar(icent, 0, 0) - sysDefault->getLLPar(icent, 0, 0));
			float ldTotSys = ldStat > ldSys ? 0 : ldSys - ldStat;
			ssSysErr[icent] += ssTotSys*ssTotSys;
			ldSysErr[icent] += ldTotSys*ldTotSys;

			float ssGausStat = sqrt(fabs(getLLParErr(icent, isys, 1)*getLLParErr(icent, isys, 1)
				- sysDefault->getLLParErr(icent, 0, 1)*sysDefault->getLLParErr(icent, 0, 1)));
			float ssGausSys = fabs(getLLPar(icent, isys, 1) - sysDefault->getLLPar(icent, 0, 1));
			float ssGausTotSys = ssStat > ssSys ? 0 : ssSys - ssStat;

			float ldGausStat = sqrt(fabs(getGausParErr(icent, isys, 0)*getGausParErr(icent, isys, 0)
				- sysDefault->getLLParErr(icent, 0, 0)*sysDefault->getGausParErr(icent, 0, 0)));
			float ldGausSys = fabs(getGausPar(icent, 0, 0) - sysDefault->getGausPar(icent, 0, 0));
			float ldGausTotSys = ldGausStat > ldGausSys ? 0 : ldGausSys - ldGausStat;
			ssGausSysErr[icent] += ssGausTotSys*ssGausTotSys;
			ldGausSysErr[icent] += ldGausTotSys*ldGausTotSys;
		}
		ssSysErr[icent] = sqrt(ssSysErr[icent] / mSysFile.size());
		ldSysErr[icent] = sqrt(ldSysErr[icent] / mSysFile.size());
		ssGausSysErr[icent] = sqrt(ssGausSysErr[icent] / mSysFile.size());
		ldGausSysErr[icent] = sqrt(ldGausSysErr[icent] / mSysFile.size());
	}
}

float SysSource::getLLPar(int icent, int isys, int ipar)
{
	return mCFLLFit[icent][isys]->GetParameter(ipar);
}

float SysSource::getLLParErr(int icent, int isys, int ipar)
{
	return mCFLLFit[icent][isys]->GetParError(ipar);
}

float SysSource::getGausPar(int icent, int isys, int ipar)
{
	return mCFGausFit[icent][isys]->GetParameter(ipar);
}

float SysSource::getGausParErr(int icent, int isys, int ipar)
{
	return mCFGausFit[icent][isys]->GetParError(ipar);
}
#endif
