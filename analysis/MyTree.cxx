#include "MyTree.h"
#include <iostream>
#include <ctime>
#include <TMath.h>

MyTree::MyTree(TTree* tree) : mNEvent(-1), mBufferMass{0}
{
	if(!tree) {
		std::cout << "ERROR: provide tree is a nullptr" << std::endl;
		mTree = nullptr;
	} else {
		mTree = tree;
		setBranchAddress();
	}
	for(int icent = 0; icent < 9; ++icent) {
		mMaxMixEvent[icent] = -1;
	}
}

MyTree::~MyTree()
{
	if(mTree) {
		delete mTree;
	}
}

bool MyTree::setBranchAddress()
{
	if(mTree) {
		mTree->SetBranchAddress("runid", &mBufferRunId);
		mTree->SetBranchAddress("vz", &mBufferVz);
		mTree->SetBranchAddress("vx", &mBufferVx);
		mTree->SetBranchAddress("vy", &mBufferVy);
		mTree->SetBranchAddress("eId", &mBufferEventId);
		mTree->SetBranchAddress("countrefmult", &mBufferRefMult);
		//mTree->SetBranchAddress("ctofmult", &mBufferTofMult);
		//mTree->SetBranchAddress("centid", &mBufferCentId);
		mTree->SetBranchAddress("ntrk", &mBufferNTrack);
		mTree->SetBranchAddress("fCentrality", &mBufferCent9);

		//mTree->SetBranchAddress("pdg", mBufferPDG);
		//mTree->SetBranchAddress("pt", mBufferPt);
		//mTree->SetBranchAddress("phi", mBufferPhi);
		//mTree->SetBranchAddress("eta", mBufferEta);
		//mTree->SetBranchAddress("rap", mBufferRap);
		//mTree->SetBranchAddress("mass", mBufferMass);
		mTree->SetBranchAddress("SIMD_pt", mBufferPt);
		mTree->SetBranchAddress("SIMD_phi", mBufferPhi);
		mTree->SetBranchAddress("SIMD_eta", mBufferEta);
		mTree->SetBranchAddress("SIMD_rapidity", mBufferRap);
		mTree->SetBranchAddress("SIMD_mass", mBufferMass);

		mTree->SetBranchAddress("SIMD_chi2topo", mBufferChi2Topo);
		mTree->SetBranchAddress("SIMD_chi2ndf", mBufferChi2NDF);
		mTree->SetBranchAddress("chi2primaryA", mBufferChi2PrimPip);
		mTree->SetBranchAddress("chi2primaryB", mBufferChi2PrimPim);
		mTree->SetBranchAddress("rapidityA", mBufferRapPip);
		mTree->SetBranchAddress("rapidityB", mBufferRapPim);
		mTree->SetBranchAddress("pt_db", mBufferPtPip);
		mTree->SetBranchAddress("ptB", mBufferPtPim);
		//mTree->SetBranchAddress("bx", mBufferBx);
		//mTree->SetBranchAddress("by", mBufferBy);
		//mTree->SetBranchAddress("bz", mBufferBz);
		mTree->SetBranchAddress("nhitsA", mBufferNHitsA);
		mTree->SetBranchAddress("nhitsB", mBufferNHitsB);
		mTree->SetBranchAddress("m2A", mBufferM2A);
		mTree->SetBranchAddress("m2B", mBufferM2B);
		mTree->SetBranchAddress("eta_db", mBufferEtaA);
		mTree->SetBranchAddress("etaB", mBufferEtaB);
		mTree->SetBranchAddress("dca_db", mBufferDCAA);
		mTree->SetBranchAddress("dcaB", mBufferDCAB);
		//mTree->SetBranchAddress("phDBRF", mBufferTrkIdA);
		//mTree->SetBranchAddress("thDBRF", mBufferTrkIdB);
		//mTree->SetBranchAddress("dgdca", mBufferDgDCA);
		mTree->SetBranchAddress("SIMD_dca", mBufferDCA);
		mTree->SetBranchAddress("SIMD_decaylength", mBufferDecayLength);
		mTree->SetBranchAddress("nsigmaA0", mBufferNSigmaA);
		mTree->SetBranchAddress("nsigmaB0", mBufferNSigmaB);
		//mTree->SetBranchAddress("dedxA", mBufferdEdxA);
		//mTree->SetBranchAddress("dedxB", mBufferdEdxB);
		//mTree->SetBranchAddress("nHitsDedxA", mBufferNHitsDedxA);
		//mTree->SetBranchAddress("nHitsDedxB", mBufferNHitsDedxB);
		mTree->SetBranchAddress("isMC", mBufferIsMC);
		mTree->SetBranchAddress("mcPx", mBufferMCPx);
		mTree->SetBranchAddress("mcPy", mBufferMCPy);
		mTree->SetBranchAddress("mcPz", mBufferMCPz);
		return true;
	} else {
		return false;
	}
}

Long64_t MyTree::getNEvent()
{
	if(mNEvent < 0) {
		mNEvent = mTree->GetEntries();
		return mNEvent;
	} else {
		return mNEvent;
	}
	std::cout << "Get NEvent fail" << std::endl;
	return -1;
}

int MyTree::getEntry(Long64_t entry)
{
	return mTree->GetEntry(entry);
}

void MyTree::copyToBuffer()
{
	if(mBufferCent9 >= 0 && mBufferCent9 <= 8) {
		int replaceEvents = -1;
		srand((int)time(0));
		if(mMaxMixEvent[(int)mBufferCent9] >= 14) {
			replaceEvents = rand()%15;
			if(replaceEvents < 0 || replaceEvents > 14) {
				std::cout << "ERROR: wrong replace events number" << std::endl;
			}
		} else {
			replaceEvents = mMaxMixEvent[(int)mBufferCent9] + 1;
		}
		mMixBuffer[(int)mBufferCent9][replaceEvents].mBufferNTrack = mBufferNTrack;
		mMixBuffer[(int)mBufferCent9][replaceEvents].mBufferEventId = mBufferEventId;
		for(int iK = 0; iK < mBufferNTrack; ++iK) {
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].pdg = mBufferPDG[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].pt = mBufferPt[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].phi = mBufferPhi[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].eta = mBufferEta[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].rap = mBufferRap[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].mass = mBufferMass[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].chi2Topo = mBufferChi2Topo[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].chi2NDF = mBufferChi2NDF[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].chi2PrimPip = mBufferChi2PrimPip[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].chi2PrimPim = mBufferChi2PrimPim[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].rapPip = mBufferRapPip[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].rapPim = mBufferRapPim[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].ptPip = mBufferPtPip[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].ptPim = mBufferPtPim[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].bx = mBufferBx[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].by = mBufferBy[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].bz = mBufferBz[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].nHitsA = mBufferNHitsA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].nHitsB = mBufferNHitsB[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].m2A = mBufferM2A[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].m2B = mBufferM2B[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].etaA = mBufferEtaA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].etaB = mBufferEtaB[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].dcaA = mBufferDCAA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].dcaB = mBufferDCAB[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].trkIdA = mBufferTrkIdA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].trkIdB = mBufferTrkIdB[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].dgDCA = mBufferDgDCA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].dca = mBufferDCA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].decayLength = mBufferDecayLength[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].nSigmaA = mBufferNSigmaA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].nSigmaB = mBufferNSigmaB[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].dEdxA = mBufferdEdxA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].dEdxB = mBufferdEdxB[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].nHitsDedxA = mBufferNHitsDedxA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].nHitsDedxB = mBufferNHitsDedxB[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].isMC = mBufferIsMC[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].mcPx = mBufferMCPx[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].mcPy = mBufferMCPy[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iK].mcPz = mBufferMCPz[iK];
		}
		if(mMaxMixEvent[(int)mBufferCent9] < 13) {
			mMaxMixEvent[(int)mBufferCent9]++;
		} else {
			mMaxMixEvent[(int)mBufferCent9] = 14;
		}
	}
}


void MyTree::copyToBuffer(std::vector<int>& idx)
{
	if(mBufferCent9 >= 0 && mBufferCent9 <= 8) {
		int replaceEvents = -1;
		srand((int)time(0));
		if(mMaxMixEvent[(int)mBufferCent9] >= 24) {
			replaceEvents = rand()%25;
			if(replaceEvents < 0 || replaceEvents > 24) {
				std::cout << "ERROR: wrong replace events number" << std::endl;
			}
		} else {
			replaceEvents = mMaxMixEvent[(int)mBufferCent9] + 1;
		}
		mMixBuffer[(int)mBufferCent9][replaceEvents].mBufferNTrack = idx.size();
		mMixBuffer[(int)mBufferCent9][replaceEvents].mBufferEventId = mBufferEventId;
		for(int iiK = 0; iiK < idx.size(); ++iiK) {
			int iK = idx[iiK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].pdg = mBufferPDG[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].pt = mBufferPt[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].phi = mBufferPhi[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].eta = mBufferEta[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].rap = mBufferRap[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].mass = mBufferMass[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].chi2Topo = mBufferChi2Topo[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].chi2NDF = mBufferChi2NDF[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].chi2PrimPip = mBufferChi2PrimPip[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].chi2PrimPim = mBufferChi2PrimPim[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].rapPip = mBufferRapPip[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].rapPim = mBufferRapPim[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].ptPip = mBufferPtPip[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].ptPim = mBufferPtPim[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].bx = mBufferBx[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].by = mBufferBy[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].bz = mBufferBz[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].nHitsA = mBufferNHitsA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].nHitsB = mBufferNHitsB[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].m2A = mBufferM2A[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].m2B = mBufferM2B[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].etaA = mBufferEtaA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].etaB = mBufferEtaB[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].dcaA = mBufferDCAA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].dcaB = mBufferDCAB[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].trkIdA = mBufferTrkIdA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].trkIdB = mBufferTrkIdB[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].dgDCA = mBufferDgDCA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].dca = mBufferDCA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].decayLength = mBufferDecayLength[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].nSigmaA = mBufferNSigmaA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].nSigmaB = mBufferNSigmaB[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].dEdxA = mBufferdEdxA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].dEdxB = mBufferdEdxB[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].nHitsDedxA = mBufferNHitsDedxA[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].nHitsDedxB = mBufferNHitsDedxB[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].isMC = mBufferIsMC[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].mcPx = mBufferMCPx[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].mcPy = mBufferMCPy[iK];
			mMixBuffer[(int)mBufferCent9][replaceEvents].particle[iiK].mcPz = mBufferMCPz[iK];
		}
		if(mMaxMixEvent[(int)mBufferCent9] < 23) {
			mMaxMixEvent[(int)mBufferCent9]++;
		} else {
			mMaxMixEvent[(int)mBufferCent9] = 24;
		}
	}
}

MyTree::Particle MyTree::getParticle(int iparticle, float beamRapidity)
{
	Particle particle;
	particle.pdg = mBufferPDG[iparticle];
	particle.pt = mBufferPt[iparticle];
	particle.phi = mBufferPhi[iparticle];
	particle.eta = mBufferEta[iparticle];
	particle.rap = -(mBufferRap[iparticle] + beamRapidity);
	particle.mass = mBufferMass[iparticle];
	particle.chi2Topo = mBufferChi2Topo[iparticle];
	particle.chi2NDF = mBufferChi2NDF[iparticle];
	particle.chi2PrimPip = mBufferChi2PrimPip[iparticle];
	particle.chi2PrimPim = mBufferChi2PrimPim[iparticle];
	particle.rapPip = -(mBufferRapPip[iparticle] + beamRapidity);
	particle.rapPim = -(mBufferRapPim[iparticle] + beamRapidity);
	particle.ptPip = mBufferPtPip[iparticle];
	particle.ptPim = mBufferPtPim[iparticle] / 1000.;
	particle.bx = mBufferBx[iparticle];
	particle.by = mBufferBy[iparticle];
	particle.bz = mBufferBz[iparticle];
	particle.nHitsA = mBufferNHitsA[iparticle];
	particle.nHitsB = mBufferNHitsB[iparticle];
	particle.m2A = mBufferM2A[iparticle];
	particle.m2B = mBufferM2B[iparticle];
	particle.etaA = mBufferEtaA[iparticle];
	particle.etaB = mBufferEtaB[iparticle] / 1000.;
	particle.dcaA = mBufferDCAA[iparticle];
	particle.dcaB = mBufferDCAB[iparticle];
	particle.trkIdA = mBufferTrkIdA[iparticle];
	particle.trkIdB = mBufferTrkIdB[iparticle];
	particle.dgDCA = mBufferDgDCA[iparticle];
	particle.dca = mBufferDCA[iparticle];
	particle.decayLength = mBufferDecayLength[iparticle];
	particle.nSigmaA = mBufferNSigmaA[iparticle];
	particle.nSigmaB = mBufferNSigmaB[iparticle];
	particle.dEdxA = mBufferdEdxA[iparticle];
	particle.dEdxB = mBufferdEdxB[iparticle];
	particle.nHitsDedxA = mBufferNHitsDedxA[iparticle];
	particle.nHitsDedxB = mBufferNHitsDedxB[iparticle];
	particle.isMC = mBufferIsMC[iparticle];
	particle.mcPx = mBufferMCPx[iparticle];
	particle.mcPy = mBufferMCPy[iparticle];
	particle.mcPz = mBufferMCPz[iparticle];

	particle.px = particle.pt * TMath::Cos(particle.phi);
	particle.py = particle.pt * TMath::Sin(particle.phi);
	particle.pz = particle.pt * TMath::SinH(particle.eta);
	particle.energy = sqrt(particle.pz*particle.pz + particle.pt*particle.pt + particle.mass*particle.mass);
	particle.pzA = particle.ptPip * TMath::SinH(particle.etaA);
	particle.pzB = particle.ptPim * TMath::SinH(particle.etaB);
	particle.pA = sqrt(particle.pzA*particle.pzA + particle.ptPip*particle.ptPip);
	particle.pB = sqrt(particle.pzB*particle.pzB + particle.ptPim*particle.ptPim);
	return particle;
}

MyTree::Particle MyTree::getMixParticle(int cent, int ievt, int iparticle, float beamRapidity)
{
	MyTree::Particle particle = mMixBuffer[cent][ievt].particle[iparticle];
	particle.rap = -(particle.rap + beamRapidity);
	particle.rapPip = -(particle.rapPip + beamRapidity);
	particle.rapPim = -(particle.rapPim + beamRapidity);
	particle.ptPim /= 1000.;
	particle.etaB /= 1000.;

	particle.px = particle.pt * TMath::Cos(particle.phi);
	particle.py = particle.pt * TMath::Sin(particle.phi);
	particle.pz = particle.pt * TMath::SinH(particle.eta);
	particle.energy = sqrt(particle.pz*particle.pz + particle.pt*particle.pt + particle.mass*particle.mass);
	particle.pzA = particle.ptPip * TMath::SinH(particle.etaA);
	particle.pzB = particle.ptPim * TMath::SinH(particle.etaB);
	particle.pA = sqrt(particle.pzA*particle.pzA + particle.ptPip*particle.ptPip);
	particle.pB = sqrt(particle.pzB*particle.pzB + particle.ptPim*particle.ptPim);
	return particle;
}
