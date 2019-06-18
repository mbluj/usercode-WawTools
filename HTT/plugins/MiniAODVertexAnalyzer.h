// M. Bluj
#ifndef WarsawAnalysis_HTT_MiniAODVertexAnalyzer_h
#define WarsawAnalysis_HTT_MiniAODVertexAnalyzer_h

// system include files
#include <memory>
#include <cmath>
//
#include <vector>
#include <string>
#include <algorithm>
#include <set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "WarsawAnalysis/HTTDataFormats/interface/HTTEvent.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

//
// class declaration
//

class MiniAODVertexAnalyzer : public edm::EDAnalyzer {
  public:
    explicit MiniAODVertexAnalyzer(const edm::ParameterSet & iConfig);
    ~MiniAODVertexAnalyzer();

  private:

    ///Fill branches with generator level data.
    bool getGenData(const edm::Event & iEvent, const edm::EventSetup & iSetup);

    ///Fill branches with reconstruction level data.
    bool getRecoData(const edm::Event & iEvent, const edm::EventSetup & iSetup);

    ///Find PV using different discriminators. Return false
    ///if no vertices were found in the event
    ///pfPV  - using PF particles for score calulation
    ///pt2PV - using sum pt^2 os all tracks assigned to vertex.
    ///        as tracks pt<1 do not have errors, we take
    ///        estimate based on TRK-11-001 note    
    bool findPrimaryVertices(const edm::Event & iEvent, const edm::EventSetup & iSetup);

    //Refit PV using track information stored in miniAOD
    ///This has to be run AFTER finding tau candidates
    bool refitPV(const edm::Event & iEvent, const edm::EventSetup & iSetup);

    ///Find reconstructed tau candidates.
    bool findRecoTaus(const edm::Event & iEvent, const edm::EventSetup & iSetup);

    ///Find best mu-tau pair
    std::pair<const reco::Candidate*, const reco::Candidate*> findMTPair(const edm::Handle<std::vector<pat::Muon> > &muColl,
									 const edm::Handle<std::vector<pat::Tau> > &tauColl,
									 const edm::Handle<reco::VertexCollection> &vtxColl);
    ///Find best e-tau pair
    std::pair<const reco::Candidate*, const reco::Candidate*> findETPair(const edm::Handle<std::vector<pat::Electron> > &eColl,
									 const edm::Handle<std::vector<pat::Tau> > &tauColl,
									 const edm::Handle<reco::VertexCollection> &vtxColl);
    ///Find best tau-tau pair
    std::pair<const pat::Tau*, const pat::Tau*> findTTPair(const edm::Handle<std::vector<pat::Tau> > &tauColl);

    bool diMuVeto(const edm::Handle<std::vector<pat::Muon> > &muColl,
		  const edm::Handle<reco::VertexCollection> &vtxColl);
    bool diEVeto(const edm::Handle<std::vector<pat::Electron> > &eColl,
		 const edm::Handle<reco::VertexCollection> &vtxColl);

    size_t recoverStrips(const pat::Tau *aTau,
			 std::vector<reco::Candidate::LorentzVector> &pi0_dau_P4);
    std::set<size_t> buildStrip(const std::vector<reco::Candidate::LorentzVector> &signalGammas,
				std::set<size_t> & usedGammaIds);

    ///Set PCA vector for reco taus. This must be run AFTER
    ///vertex refitting.
    bool setPCAVectors(const edm::Event & iEvent, const edm::EventSetup & iSetup);

    ///Get leading track of tau or track of e/mu
    const reco::Track* getLeadTrack(const edm::Event & iEvent,
				    const reco::Candidate *aTau);

    ///Calculate PCA for given tau candidate, wrt. given vertex.
    TVector3 getPCA(const edm::Event & iEvent, const edm::EventSetup & iSetup,
		    const reco::Track* aTrack,
		    const GlobalPoint & aPoint);

    bool getPVBalance(const edm::Event & iEvent, size_t iPV=0);
    
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    float muRelIso(const pat::Muon *aMu);
    float eRelIso(const pat::Electron *aEle);

    edm::EDGetTokenT<edm::ValueMap<float>> scores_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate> >  cands_, lostCands_, lostCandsEle_;
    edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
    edm::EDGetTokenT<reco::VertexCollection> vertices_;
    edm::EDGetTokenT<reco::BeamSpot> bs_;
    edm::EDGetTokenT<std::vector<pat::Tau> > taus_;
    edm::EDGetTokenT<std::vector<pat::Muon> > mus_;
    edm::EDGetTokenT<std::vector<pat::Electron> > eles_;

    bool useBeamSpot_, useLostCands_, useTauTracks_;
    bool verbose_;

    //std::map<std::string, float> treeVars_;      
    TTree *tree_;

    HTTEvent *myEvent_;

    std::pair<const reco::Candidate*, const reco::Candidate*> thePair_;
    
       
};

#endif

  
