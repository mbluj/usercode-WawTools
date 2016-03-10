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
#include "DataFormats/PatCandidates/interface/Tau.h"

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

      // auxiliary struct declaration
      struct TVisTau{
	TLorentzVector p4, leadChP4, neutralP4, piZeroP4;
	int nCharged, nNeutral, nPiZero; 
	TVector3 sv;
	int charge, decMode;
	double iso, outerPt; //FIXME: How define it for patTau (iso+1000*outerPt)? Two variables?
      };

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      std::pair<TVisTau,TVisTau> findTauPair(const reco::GenParticleCollection&  particles,
					     const pat::TauCollection& taus);
      bool piZeroQuality(const reco::RecoTauPiZero& piZero, const pat::Tau& tau);
      TTree *initTree(edm::Service<TFileService> &fs, std::string name="tree");
      void bookVariable(TTree *t=0, std::string var="foo");      
      void clean();

      //edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
      edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken_;
      edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
      edm::EDGetTokenT<pat::TauCollection> toyTauToken_;
      bool verbose_;
      bool isSignal_;

      std::map<std::string, float> treeVars_;      
      TTree *tree_; 

      int bosonId_, decModeMinus_, decModePlus_;
      int toyDecModeMinus_, toyDecModePlus_;
      int toyNChargedMinus_, toyNNeutralMinus_, toyNPiZeroMinus_,
	toyNChargedPlus_, toyNNeutralPlus_, toyNPiZeroPlus_;

      TLorentzVector *p4Sum_, *metNu_, *met_;
      TLorentzVector *piMinus_, *piPlus_;
      TLorentzVector *tauMinus_, *tauPlus_;
      TLorentzVector *visTauMinus_, *visTauPlus_;
      TLorentzVector *toyPiMinus_, *toyPiPlus_;
      TLorentzVector *toyNeutralMinus_, *toyNeutralPlus_;
      TLorentzVector *toyPiZeroMinus_, *toyPiZeroPlus_;
      TLorentzVector *toyTauMinus_, *toyTauPlus_;

      TVector3 *thePV_, *svMinus_, *svPlus_, *toySvMinus_, *toySvPlus_;
      TVector3 *nPiMinus_, *nPiPlus_, *toyNPiMinus_, *toyNPiPlus_; 
  
};

#endif

