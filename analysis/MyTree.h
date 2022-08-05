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

	int mBufferRunId;
	float mBufferVz;
	float mBufferVx;
	float mBufferVy;
	float mBufferCent9;
	int mBufferEventId;
	float mBufferRefMult;
	float mBufferTofMult;
	char mBufferCentId;
	unsigned int mBufferNTrack;
	int mMaxMixEvent[9];

	int mBufferPDG[nTrackMax];
	float mBufferPt[nTrackMax];
	float mBufferPhi[nTrackMax];
	float mBufferEta[nTrackMax];
	float mBufferRap[nTrackMax];
	float mBufferMass[nTrackMax];
	float mBufferChi2Topo[nTrackMax];
	float mBufferChi2NDF[nTrackMax];
	float mBufferChi2PrimPip[nTrackMax];	//pi+
	float mBufferChi2PrimPim[nTrackMax];	//pi-
	float mBufferRapPip[nTrackMax];
	float mBufferRapPim[nTrackMax];
	float mBufferPtPip[nTrackMax];
	short mBufferPtPim[nTrackMax];
	float mBufferBx[nTrackMax];
	float mBufferBy[nTrackMax];
	float mBufferBz[nTrackMax];
	float mBufferNHitsA[nTrackMax];
	float mBufferNHitsB[nTrackMax];
	float mBufferM2A[nTrackMax];
	float mBufferM2B[nTrackMax];
	float mBufferEtaA[nTrackMax];
	short mBufferEtaB[nTrackMax];
	float mBufferDCAA[nTrackMax];
	float mBufferDCAB[nTrackMax];
	float mBufferTrkIdA[nTrackMax];
	float mBufferTrkIdB[nTrackMax];
	float mBufferDgDCA[nTrackMax];
	float mBufferDCA[nTrackMax];
	float mBufferDecayLength[nTrackMax];

	float mBufferNSigmaA[nTrackMax];
	float mBufferNSigmaB[nTrackMax];
	float mBufferdEdxA[nTrackMax];
	float mBufferdEdxB[nTrackMax];
	float mBufferNHitsDedxA[nTrackMax];
	float mBufferNHitsDedxB[nTrackMax];

	struct Particle {
		int pdg;
		float pt;
		float phi;
		float phiA, phiB;
		float eta;
		float px, py, pz;
		float energy;
		float pzA, pzB;
		float pA, pB;
		float rap;
		float mass;
		float chi2Topo;
		float chi2NDF;
		float chi2PrimPip;      //pi+
		float chi2PrimPim;      //pi-
		float rapPip;
		float rapPim;
		float ptPip;
		float ptPim;
		float bx, by, bz;
		float nHitsA;
		float nHitsB;
		float m2A;
		float m2B;
		float etaA;
		float etaB;
		float dcaA;
		float dcaB;
		float trkIdA;
		float trkIdB;
		float dgDCA;
		float dca;
		float decayLength;

		float nSigmaA;
		float nSigmaB;
		float dEdxA;
		float dEdxB;
		float nHitsDedxA;
		float nHitsDedxB;
	};
        Particle getParticle(int iparticle, float beamRapidity);
        Particle getMixParticle(int cent, int ievt, int iparticle, float beamRapidity);

	struct MixBuffer {
		unsigned int mBufferNTrack;
		unsigned int mBufferEventId;
		Particle particle[nTrackMax];
	};
	MixBuffer mMixBuffer[9][25];
	void copyToBuffer();
	void copyToBuffer(std::vector<int>& idx);
};
