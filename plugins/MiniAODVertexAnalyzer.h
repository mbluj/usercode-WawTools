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
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    edm::EDGetTokenT<edm::ValueMap<float>> scores_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate> >  cands_;
    edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
    edm::EDGetTokenT<reco::VertexCollection> vertices_;
    //int threshold_;
    bool verbose_;

    //std::map<std::string, float> treeVars_;      
    TTree *tree_;

    float run_, lumi_, event_;
    int bosonId_, decModeMinus_, decModePlus_;

    TLorentzVector *p4Sum_;
    TLorentzVector *piMinus_, *piPlus_;
    TLorentzVector *tauMinus_, *tauPlus_;
    TLorentzVector *visTauMinus_, *visTauPlus_;

    TVector3 *genPV_, *aodPV_, *pfPV_, *oldPV_;
    TVector3 *svMinus_, *svPlus_;
};

#endif

