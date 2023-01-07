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

void MyTree::copyToBuffer(std::vector<Particle>& vect, std::vector<Particle>& vectL, std::vector<Particle>& vectR)
{
	int cent9 = (int)mBufferCent9;
	if(cent9 >= 0 && cent9 <= 8) {
		mMixBuffer[cent9].push_front(std::move(vect));
		mMixLBuffer[cent9].push_front(std::move(vectL));
		mMixRBuffer[cent9].push_front(std::move(vectR));

		if(mMixBuffer[cent9].size() > 100) {
			mMixBuffer[cent9].pop_back();
			mMixLBuffer[cent9].pop_back();
			mMixRBuffer[cent9].pop_back();
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
