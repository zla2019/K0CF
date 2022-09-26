#include "MyTree.h"
#include <iostream>
#include <ctime>
#include <TMath.h>
#include <TRandom.h>
#include <TLorentzVector.h>

MyTree::MyTree(TTree* tree) : mNEvent(-1)
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
		mTree->SetBranchAddress("mul", &mMul);
		mTree->SetBranchAddress("b", &mB);
		mTree->SetBranchAddress("Npart", &mNPart);
		mTree->SetBranchAddress("pid", mPid);
		mTree->SetBranchAddress("px", mPx);
		mTree->SetBranchAddress("py", mPy);
		mTree->SetBranchAddress("pz", mPz);
		mTree->SetBranchAddress("frx", mFrx);
		mTree->SetBranchAddress("fry", mFry);
		mTree->SetBranchAddress("frz", mFrz);
		mTree->SetBranchAddress("frt", mFrt);
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
	mPsi = 2 * TMath::Pi() * gRandom->Rndm() - TMath::Pi();
	return mTree->GetEntry(entry);
}

void MyTree::copyToBuffer(std::vector<Particle>& labK, std::vector<Particle>& urqmdK)
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
		mMixBuffer[(int)mBufferCent9][replaceEvents].mBufferNTrack = labK.size();
		mMixBuffer[(int)mBufferCent9][replaceEvents].mBufferEventId = mBufferEventId;
		mMixBuffer[(int)mBufferCent9][replaceEvents].labParticle.swap(labK);
		mMixBuffer[(int)mBufferCent9][replaceEvents].urqmdParticle.swap(urqmdK);
		if(mMaxMixEvent[(int)mBufferCent9] < 23) {
			mMaxMixEvent[(int)mBufferCent9]++;
		} else {
			mMaxMixEvent[(int)mBufferCent9] = 24;
		}
	}
	//std::cout << "mMaxMixEvent in " << (int)mBufferCent9 << " is: " << mMaxMixEvent[(int)mBufferCent9] << std::endl;
}

MyTree::Particle MyTree::getParticle(int iparticle)
{
        Particle particle;
        particle.pid = mPid[iparticle];
        particle.px = mPx[iparticle];
        particle.py = mPy[iparticle];
        particle.pz = mPz[iparticle];
        particle.frx = mFrx[iparticle];
        particle.fry = mFry[iparticle];
        particle.frz = mFrz[iparticle];
        particle.frt = mFrt[iparticle];
        float mass = 0;
        if(fabs(particle.pid) == 321) mass = 0.49368;
        else if(fabs(particle.pid) == 2212) mass = 0.93827;
        else if(fabs(particle.pid) == 211) mass = 0.13957;
        else if(fabs(particle.pid) == 311) mass = 0.49761;
        float energy = sqrt(particle.px*particle.px + particle.py*particle.py + particle.pz*particle.pz + mass*mass);
        float rPhi = TMath::ATan(particle.fry / particle.frx);
        float rt = sqrt(particle.frx*particle.frx + particle.fry*particle.fry);

        TLorentzVector mom_v4(particle.px, particle.py, particle.pz, energy);
        particle.pt = sqrt(particle.px*particle.px + particle.py*particle.py);
        particle.phi = mom_v4.Phi();
        if((particle.phi + mPsi) > TMath::Pi()) {
                particle.phi = particle.phi + mPsi - 2 * TMath::Pi();
        } else if((particle.phi + mPsi) < -TMath::Pi()) {
                particle.phi = particle.phi + mPsi + 2 * TMath::Pi();
        } else {
                particle.phi = particle.phi + mPsi;
        }
        if((rPhi + mPsi) > TMath::Pi()) {
                rPhi = rPhi + mPsi - 2 * TMath::Pi();
        } else if((rPhi + mPsi) < -TMath::Pi()) {
                rPhi = rPhi + mPsi + 2 * TMath::Pi();
        } else {
                rPhi = rPhi + mPsi;
        }

        particle.frx = rt * TMath::Cos(rPhi);
        particle.fry = rt * TMath::Sin(rPhi);
        particle.px = particle.pt * TMath::Cos(particle.phi);
        particle.py = particle.pt * TMath::Sin(particle.phi);
        particle.p = sqrt(particle.pt*particle.pt + particle.pz*particle.pz);
        particle.eta = 0.5 * log((particle.p + particle.pz) / (particle.p - particle.pz));
        particle.energy = sqrt(particle.p*particle.p + mass*mass);
        particle.y = 0.5 * log((particle.energy + particle.pz) / (particle.energy - particle.pz));
        TLorentzVector p_v4(particle.px, particle.py, particle.pz, particle.energy);
	return particle;
}
