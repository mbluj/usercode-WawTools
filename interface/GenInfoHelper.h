// M. Bluj

#ifndef WarsawAnalysis_HTT_GenInfoHelper_h
#define WarsawAnalysis_HTT_GenInfoHelper_h

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include <utility>

//Root
#include "TLorentzVector.h"
#include "TVector3.h"

namespace WawGenInfoHelper {
  enum tauDecayModes {kElectron, kMuon, 
		      kOneProng0pi0, kOneProng1pi0, kOneProng2pi0, kOneProng3pi0,
		      kThreeProng0pi0, kThreeProng1pi0,
		      kOther, kUndefined};

  
  enum hZDecayModes {kMuTau,kETau,kTauTau,kMuMu,kEE,kEMu,
		     kEEPrompt, kMMPrompt,
		     kHZOther, kHZUndefined,  kTauTauPrompt};

  typedef reco::GenParticleCollection::const_iterator IG;
  typedef reco::GenParticleRefVector::const_iterator IGR;

  bool isBoson(const reco::GenParticleRef& particle, 
	      bool checkLastCopy=false, bool tauDec=true);
  const reco::GenParticleRef getFinalClone(const reco::GenParticleRef& particle, bool isUnstable=true);
  bool isFinalClone(const reco::GenParticleRef& particle, bool isUnstable=true);
  const reco::GenParticleRef getInitialClone(const reco::GenParticleRef& initialClone);
  bool isInitialClone(const reco::GenParticleRef& particle);
  int getTauDecayMode(const reco::GenParticleRefVector& products);
  int getTauDirDecayMode(const reco::GenParticleRefVector& products);
  int getTausDecays(const reco::GenParticleRef& tau,
		    reco::GenParticleRefVector& products,
		    bool ignoreNus=true,bool direct=false);
  const reco::GenParticleRef getLeadChParticle(const reco::GenParticleRefVector& products);
  TLorentzVector getP4(const reco::GenParticleRef& part){
    return TLorentzVector(part->px(),part->py(),part->pz(),part->energy());
  }
  TLorentzVector getP4(const reco::GenParticle&  part){
    return TLorentzVector(part.px(),part.py(),part.pz(),part.energy());
  }
  void setP4Ptr(const TLorentzVector& p4, TLorentzVector *p4Ptr){
    if(p4Ptr) p4Ptr->SetPxPyPzE(p4.Px(),p4.Py(),p4.Pz(),p4.Energy());
    else p4Ptr = new TLorentzVector(p4);
  }
  void setV3Ptr(const TVector3& v3, TVector3 *v3Ptr){
    if(v3Ptr) v3Ptr->SetXYZ(v3.X(),v3.Y(),v3.Z());
    else v3Ptr = new TVector3(v3);
  }
  TLorentzVector getCombinedP4(const reco::GenParticleRefVector& products){    
    TLorentzVector p4;
    for(IGR id=products.begin(); id!=products.end(); ++id)
      p4 += getP4( *id );
    return p4;
  }

  TLorentzVector getLeadChParticleP4(const reco::GenParticleRefVector& products){
    return getP4( getLeadChParticle(products) );
  }
  TVector3 getVertex(const reco::GenParticleRef& part){
    return TVector3(part->vx(),part->vy(),part->vz());
  }
  void getVertex(const reco::GenParticleRef& part, TVector3 *vtx){
    setV3Ptr(getVertex(part),vtx);
  }

  TVector3 impactParameter(const TVector3& pv, const TVector3& sv, const TLorentzVector& p4){
    TVector3 dir = (p4.Vect()).Unit();
    return (sv-pv) - ((sv-pv)*dir)*dir;
}
  void impactParameter(const TVector3& pv, const TVector3& sv, const TLorentzVector& p4,
		       TVector3 *ip){
    setV3Ptr(impactParameter(pv,sv,p4),ip);
  }
  TLorentzVector getGenMet(const reco::GenParticleRefVector& particles);
  TLorentzVector getGenMet(const reco::GenParticleCollection&  particles);
  TLorentzVector getTauNuMet(const reco::GenParticleRefVector& taus);
  std::pair<float,float> angleBetweenPlanes(const TLorentzVector &prod1, 
					    const TLorentzVector &prod12,
					    const TLorentzVector &prod2, 
					    const TLorentzVector &prod22){
    //Boost to correct frame
    TVector3 boost = (prod1+prod2).BoostVector();
    TLorentzVector prod1Star = prod1; prod1Star.Boost(-boost);
    TLorentzVector prod12Star = prod12; prod12Star.Boost(-boost);
    TLorentzVector prod2Star = prod2; prod2Star.Boost(-boost);
    TLorentzVector prod22Star = prod22; prod22Star.Boost(-boost);
    
    //define common direction and normal vectors to decay planes
    TVector3 direction = prod1Star.Vect().Unit();
    TVector3 n1 = ( direction.Cross( prod12Star.Vect() ) ).Unit(); 
    TVector3 n2 = ( direction.Cross( prod22Star.Vect() ) ).Unit(); 
    
    float phi=TMath::ACos(n1*n2);
    float rho=TMath::ACos( (prod12Star.Vect().Unit() )*(prod22Star.Vect().Unit() ) );

    return std::make_pair(phi,rho);
  }

  //Check recursively if any ancestor of particle is the given one
  bool isAncestor(const reco::Candidate *ancestor, const reco::Candidate *particle);
  //Check recursively if any descendent of particle is the given one
  bool isDescendent(const reco::Candidate *descendent, const reco::Candidate *particle);

  /// find all particles of a given pdgId and status
  /// copy from PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h" which allows ignore status(<=0) or pdgId(==0)
  void findParticles(const reco::GenParticleCollection& sourceParticles,
                     reco::GenParticleRefVector& particleRefs, 
                     int pdgId, int status);
  /// find all descendents of a given status and pdgId (recursive)
  /// copy from PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h" which allows ignore status(<0)
  void findDescendents(const reco::GenParticleRef& base, 
		       reco::GenParticleRefVector& descendents, 
		       int status, int pdgId=0);

  void findAncestors(const reco::GenParticleRef& base, 
		     reco::GenParticleRefVector& ancestors, 
		     int status, int pdgId=0);
}

#endif
