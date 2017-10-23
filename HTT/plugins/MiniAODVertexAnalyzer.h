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
    bool findRecoTau(const edm::Event & iEvent, const edm::EventSetup & iSetup);

    ///Find best tau pair
    std::pair<const pat::Tau*, const pat::Tau*> findTauPair(edm::Handle<std::vector<pat::Tau> > tauColl);

    ///Set PCA vector for reco taus. This must be run AFTER
    ///vertex refitting.
    bool setPCAVectors(const edm::Event & iEvent, const edm::EventSetup & iSetup);

    ///Calculate PCA for given tau candidate, wrt. given vertex.
    TVector3 getPCA(const edm::Event & iEvent, const edm::EventSetup & iSetup,
		    const pat::Tau* aTau,
		    const GlobalPoint & aPoint);

    bool getPVBalance(const edm::Event & iEvent, const edm::EventSetup & iSetup, 
		      size_t iPV=0);
    
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    edm::EDGetTokenT<edm::ValueMap<float>> scores_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate> >  cands_, lostCands_;
    edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
    edm::EDGetTokenT<reco::VertexCollection> vertices_;
    edm::EDGetTokenT<reco::BeamSpot> bs_;
    edm::EDGetTokenT<std::vector<pat::Tau> > taus_;

    bool useBeamSpot_, useLostCands_, useTauTracks_;
    bool verbose_;

    //std::map<std::string, float> treeVars_;      
    TTree *tree_;

    HTTEvent *myEvent_;

    std::pair<const pat::Tau*, const pat::Tau*> thePair_;
    
       
};

#endif

  
