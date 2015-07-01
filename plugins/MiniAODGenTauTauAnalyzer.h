// based on a MiniAOD example
// M. Bluj
#ifndef WarsawAnalysis_HTT_MiniAODGenTauTauAnalyzer_h
#define WarsawAnalysis_HTT_MiniAODGenTauTauAnalyzer_h

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

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

//
// class declaration
//

class MiniAODGenTauTauAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MiniAODGenTauTauAnalyzer(const edm::ParameterSet&);
      ~MiniAODGenTauTauAnalyzer();

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      TTree *initTree(edm::Service<TFileService> &fs, std::string name="tree");
      void bookVariable(TTree *t=0, std::string var="foo");      
      void clean();

      //edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
      edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken_;
      edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
      bool verbose_;

      std::map<std::string, float> treeVars_;      
      TTree *tree_; 

      int bosonId_, decModeMinus_, decModePlus_;

      TLorentzVector *p4Sum_, *metNu_, *met_;
      TLorentzVector *piMinus_, *piPlus_;
      TLorentzVector *tauMinus_, *tauPlus_;
      TLorentzVector *visTauMinus_, *visTauPlus_;

      TVector3 *thePV_, *svMinus_, *svPlus_;
      TVector3 *nPiMinus_, *nPiPlus_; 
  
};

#endif

