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
    bool refitPV(const edm::Event & iEvent, const edm::EventSetup & iSetup);
   
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    edm::EDGetTokenT<edm::ValueMap<float>> scores_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate> >  cands_;
    edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
    edm::EDGetTokenT<reco::VertexCollection> vertices_;
    edm::EDGetTokenT<reco::BeamSpot> bs_;
    //int threshold_;
    bool useBeamSpot_;
    bool verbose_;

    //std::map<std::string, float> treeVars_;      
    TTree *tree_;

    HTTEvent *myEvent_;
       
};

#endif

  
