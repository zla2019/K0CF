#ifndef MYTREE_H
#define MYTREE_H
#include <TTree.h>

class MyTree
{
public:
	MyTree() : mNEvent(-1) {};
	MyTree(TTree* tree);
	~MyTree();
	bool setBranchAddress();
	Long64_t getNEvent();
	int getEntry(Long64_t entry);
	//private:
	TTree* mTree;
	static const unsigned int nTrackMax = 1000;
	Long64_t mNEvent;
	unsigned int mBufferCent9;
	unsigned int mBufferNTrack;
	unsigned int mBufferEventId;
	int mMaxMixEvent[9];

	int mMul;
	float mB;
	int mNPart;
	float mPsi;
	int mPid[nTrackMax];
	float mPx[nTrackMax];
	float mPy[nTrackMax];
	float mPz[nTrackMax];
	float mFrx[nTrackMax];
	float mFry[nTrackMax];
	float mFrz[nTrackMax];
	float mFrt[nTrackMax];

	struct Particle {
		int pid;
		float px, py, pz;
		float frx, fry, frz, frt;
		float pt, p, eta, y, phi;
		float energy;
	};
        Particle getParticle(int iparticle);

	struct MixBuffer {
		unsigned int mBufferNTrack;
		unsigned int mBufferEventId;
		std::vector<Particle> smearParticle;
		std::vector<Particle> labParticle;
		std::vector<Particle> urqmdParticle;
	};
	MixBuffer mMixBuffer[9][25];
	void copyToBuffer(std::vector<Particle>& labK, std::vector<Particle>& urqmdK);
};
#endif
