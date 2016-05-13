
#include "WarsawAnalysis/HTT/interface/GenInfoHelper.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include <cmath>

namespace WawGenInfoHelper {
  //////////////
  bool isBoson(const reco::GenParticleRef& particle, 
	       bool checkLastCopy, bool tauDec){

    std::set<int> bosonIds;
    bosonIds.insert(23);//Z
    bosonIds.insert(25);//h
    bosonIds.insert(35);//H
    bosonIds.insert(36);//A

    if(bosonIds.find(std::abs(particle->pdgId()))!=bosonIds.end()){
      bool isFinal =  checkLastCopy ? isFinalClone(particle) : true;
      bool tauDecay = true;
      if(tauDec){
	reco::GenParticleRefVector taus;
	findDescendents( particle, taus, -1, 15);
	tauDecay = (taus.size()>0);
      }
      //std::cout<<"[WawGenInfoHelper::isBoson]: PdgId="<<particle->pdgId()
      //         <<", isFinal: "<<isFinal<<", tauDecay: "<<tauDecay
      //         <<std::endl;
      return (isFinal && tauDecay);
    }
    else
      return false;
  }

  //////////////
  const reco::GenParticleRef getFinalClone(const reco::GenParticleRef& particle, bool isUnstable){
    
    reco::GenParticleRefVector descendents;
    int pdgId = particle->pdgId();
    findDescendents(particle, descendents, -1, pdgId);
    if(descendents.size()==0 || (isUnstable && descendents[0]->numberOfDaughters()==0) ){    
      return particle;
    }
    else{
      return getFinalClone(descendents[0],isUnstable);
    }
  }

  bool isFinalClone(const reco::GenParticleRef& particle, bool isUnstable){
    return particle==getFinalClone(particle,isUnstable);
  }

  //////////////
  const reco::GenParticleRef getInitialClone(const reco::GenParticleRef& particle){
    
    reco::GenParticleRefVector ancestors;
    int pdgId = particle->pdgId();
    findAncestors(particle, ancestors, -1, pdgId);
    if(ancestors.size()==0 ){
      return particle;
    }
    else{
      return getInitialClone(ancestors[0]);
    }
  }

  bool isInitialClone(const reco::GenParticleRef& particle){
    return particle==getInitialClone(particle);
  }

  //////////////
  int getTauDecayMode(const reco::GenParticleRefVector& products){

    int tauDecayMode = tauDecayModes::tauDecayOther;

    int numElectrons      = 0;
    int numMuons          = 0;
    int numChargedPions   = 0;
    int numNeutralPions   = 0;
    int numPhotons        = 0;
    int numNeutrinos      = 0;
    int numOtherParticles = 0;

    for(IGR idr = products.begin(); idr != products.end(); ++idr ) {
      int pdg_id = std::abs((*idr)->pdgId());
      if(pdg_id == 11) numElectrons++;
      else if(pdg_id == 13) numMuons++;
      else if(pdg_id == 211 || pdg_id == 321 ) numChargedPions++; //Count both pi+ and K+
      else if(pdg_id == 111 || pdg_id == 130 || pdg_id == 310 ) numNeutralPions++; //Count both pi0 and K0_L/S
      else if(pdg_id == 12 || 
	      pdg_id == 14 || 
	      pdg_id == 16) {
	numNeutrinos++;
      }
      else if(pdg_id == 22) numPhotons++;
      else {
	numOtherParticles++;
      }
    }
    if(numElectrons>1){//sometimes there are gamma->ee conversions 
      numPhotons += numElectrons/2;
      numElectrons -= 2*(numElectrons/2);
    }

    if( numOtherParticles == 0 ){
      if( numElectrons == 1 ){
	//--- tau decays into electrons
	tauDecayMode = tauDecayModes::tauDecaysElectron;
      } else if( numMuons == 1 ){
	//--- tau decays into muons
	tauDecayMode = tauDecayModes::tauDecayMuon;
      } else {
	//--- hadronic tau decays
	switch ( numChargedPions ){
        case 1 :
          if( numNeutralPions != 0 ){
            tauDecayMode =  tauDecayModes::tauDecayOther;
            break;
	  }
          switch ( numPhotons ){
	  case 0:
	    tauDecayMode = tauDecayModes::tauDecay1ChargedPion0PiZero;
	    break;
	  case 2:
	    tauDecayMode = tauDecayModes::tauDecay1ChargedPion1PiZero;
	    break;
	  case 4:
	    tauDecayMode = tauDecayModes::tauDecay1ChargedPion2PiZero;
	    break;
	  case 6:
	    tauDecayMode = tauDecayModes::tauDecay1ChargedPion3PiZero;
	    break;
	  case 8:
	    tauDecayMode = tauDecayModes::tauDecay1ChargedPion4PiZero;
	    break;
	  default:
	    tauDecayMode = tauDecayModes::tauDecayOther;
	    break;
          }
          break;
        case 3 : 
          if( numNeutralPions != 0 ){
            tauDecayMode = tauDecayModes::tauDecayOther;
            break;
          }
          switch ( numPhotons ){
	  case 0 : 
	    tauDecayMode = tauDecayModes::tauDecay3ChargedPion0PiZero;
	    break;
	  case 2 : 
	    tauDecayMode = tauDecayModes::tauDecay3ChargedPion1PiZero;
	    break;
	  case 4 : 
	    tauDecayMode = tauDecayModes::tauDecay3ChargedPion2PiZero;
	    break;
	  case 6 : 
	    tauDecayMode = tauDecayModes::tauDecay3ChargedPion3PiZero;
	    break;
	  case 8 : 
	    tauDecayMode = tauDecayModes::tauDecay3ChargedPion4PiZero;
	    break;
	  default:
	    tauDecayMode = tauDecayModes::tauDecayOther;
	    break;
          }
          break;
	}
      }
    }
    return tauDecayMode;
  }

  //////////////
  int getTauDirDecayMode(const reco::GenParticleRefVector& products){

    int tauDecayMode =  tauDecayModes::tauDecayOther;

    int numElectrons      = 0;
    int numMuons          = 0;
    int numChargedPions   = 0;
    int numNeutralPions   = 0;
    int numPhotons        = 0;
    int numNeutrinos      = 0;
    int numOtherParticles = 0;
    
    for(IGR idr = products.begin(); idr != products.end(); ++idr ) {
      int pdg_id = std::abs((*idr)->pdgId());
      if(pdg_id == 11) numElectrons++;
      else if(pdg_id == 13) numMuons++;
      else if(pdg_id == 211 || pdg_id == 321 ) numChargedPions++; //Count both pi+ and K+
      else if(pdg_id == 111 || pdg_id == 130 || pdg_id == 310 ) numNeutralPions++; //Count both pi0 and K0_L/S 
      else if(pdg_id == 12 || 
	      pdg_id == 14 || 
	      pdg_id == 16) {
	numNeutrinos++;
      }
      else if(pdg_id == 22) numPhotons++;
      else {
	numOtherParticles++;
      }
    }
    if(numElectrons>1){//sometimes there are gamma->ee conversions 
      numPhotons += numElectrons/2;
      numElectrons -= 2*(numElectrons/2);
    }
    if( numOtherParticles == 0 ){
      if( numElectrons == 1 ){
	//--- tau decays into electrons
	tauDecayMode = tauDecayModes::tauDecaysElectron;
      } else if( numMuons == 1 ){
	//--- tau decays into muons
	tauDecayMode = tauDecayModes::tauDecayMuon;
      } else {
	//--- hadronic tau decays
	switch ( numChargedPions ){
        case 1 :
          if( numOtherParticles != 0 ){
            tauDecayMode = tauDecayModes::tauDecayOther;
            break;
          }
          switch ( numNeutralPions ){
	  case 0:
	    tauDecayMode = tauDecayModes::tauDecay1ChargedPion0PiZero;
	    break;
	  case 1:
	    tauDecayMode = tauDecayModes::tauDecay1ChargedPion1PiZero;
	    break;
	  case 2:
	    tauDecayMode = tauDecayModes::tauDecay1ChargedPion2PiZero;
	    break;
	  case 3:
	    tauDecayMode = tauDecayModes::tauDecay1ChargedPion3PiZero;
	    break;
	  case 4:
	    tauDecayMode = tauDecayModes::tauDecay1ChargedPion4PiZero;
	    break;
	  default:
	    tauDecayMode = tauDecayModes::tauDecayOther;
	    break;
          }
          break;
        case 3 : 
	  if( numOtherParticles != 0 ){
            tauDecayMode = tauDecayModes::tauDecayOther;
            break;
          }
          switch ( numNeutralPions ){
	  case 0 : 
	    tauDecayMode = tauDecayModes::tauDecay3ChargedPion0PiZero;
	    break;
	  case 1 : 
	    tauDecayMode = tauDecayModes::tauDecay3ChargedPion1PiZero;
	    break;
	  case 2 : 
	    tauDecayMode = tauDecayModes::tauDecay3ChargedPion2PiZero;
	    break;
	  case 3 : 
	    tauDecayMode = tauDecayModes::tauDecay3ChargedPion3PiZero;
	    break;
	  case 4 : 
	    tauDecayMode = tauDecayModes::tauDecay3ChargedPion4PiZero;
	    break;
	  default:
	    tauDecayMode = tauDecayModes::tauDecayOther;
	    break;
          }
          break;
	}
      }
    }

    return tauDecayMode;
  }
  //////////////
  int getTausDecays(const reco::GenParticleRef& tau,
		    reco::GenParticleRefVector& products,
		    bool ignoreNus, bool direct){
    
    products.clear();
    if(!direct)
      findDescendents(tau, products, 1, 0);
    else{
      const reco::GenParticleRefVector& daughterRefs = tau->daughterRefVector();
      for(IGR idr = daughterRefs.begin(); idr != daughterRefs.end(); ++idr )
	products.push_back(*idr);
    }
    if(ignoreNus){
      std::set<int> allNus;
      allNus.insert(12);
      allNus.insert(14);
      allNus.insert(16);
      //allNus.insert(18);
      reco::GenParticleRefVector tmp;
      for(IGR idr=products.begin(); idr !=products.end(); ++idr)
	if(allNus.find(std::abs((*idr)->pdgId()))==allNus.end())
	  tmp.push_back((*idr));
      products.swap(tmp);
    }
    
    return direct ? getTauDirDecayMode(products) : getTauDecayMode(products);
  }

  //////////////
  TLorentzVector getGenMet(const reco::GenParticleRefVector& particles){

    float metX=0, metY=0;

    //list of invisible particles for met computation
    std::set<int> invisible;
    //neutrinos
    invisible.insert(12);
    invisible.insert(14);
    invisible.insert(16);
    invisible.insert(18);
    //LQ_ue
    invisible.insert(39);
    //~chi_10
    invisible.insert(1000022);
    //~nu_R (right neutralinos)
    invisible.insert(2000012);
    invisible.insert(2000014);
    invisible.insert(2000016);
    //~gravitino
    invisible.insert(1000039);
    invisible.insert(5000039);
    //?
    invisible.insert(9900012);
    invisible.insert(9900014);
    invisible.insert(9900016);
    //nu*_e0
    invisible.insert(4000012);
    
    for(unsigned int ipp=0; ipp<particles.size(); ++ipp){
      if(particles[ipp]->numberOfDaughters()!=0) continue; //ignore instable
      if( invisible.find( std::abs( particles[ipp]->pdgId() ) )!=invisible.end() ){
	metX += particles[ipp]->px();
	metY += particles[ipp]->py();
      }
    }

    TLorentzVector met;
    met.SetPxPyPzE( metX, metY, 0.0, sqrt(metX*metX + metY*metY) );

    return met;
  }

  //////////////
  TLorentzVector getGenMet(const reco::GenParticleCollection&  particles){

    reco::GenParticleRefVector final;
    findParticles(particles, final, 0, 1);
    return getGenMet(final);
  }

  //////////////
  TLorentzVector getTauNuMet(const reco::GenParticleRefVector& taus){
  
    float metX=0, metY=0;

    std::set<int> allNus;
    allNus.insert(12);
    allNus.insert(14);
    allNus.insert(16);  

    reco::GenParticleRefVector products;
    for(unsigned int itt=0; itt<taus.size(); ++itt)
      findDescendents(taus[itt], products, 1, 0);
  
    for(unsigned int ipp=0; ipp<products.size(); ++ipp){
      if ( allNus.find( std::abs( products[ipp]->pdgId() ) )!=allNus.end() ){
	metX += products[ipp]->px();
	metY += products[ipp]->py();
      }
    }

    TLorentzVector met;
    met.SetPxPyPzE( metX, metY, 0.0, sqrt(metX*metX + metY*metY) );

    return met;
  }

  //////////////
  const reco::GenParticleRef getLeadChParticle(const reco::GenParticleRefVector& products){
    
    float maxPt=0;
    reco::GenParticleRef part;
    /*
    std::set<int> charged;
    charged.insert(211);//pi
    charged.insert(321);//K
    charged.insert(11);//e
    charged.insert(13);//mu
    */
    for(IGR idr = products.begin(); idr != products.end(); ++idr ){
      if( (*idr)->pt() > maxPt &&
	  //charged.find( std::abs( (*idr)->pdgId() ) )!=charged.end() //MB: Logix used in pure Pythia code when charge not defined
	  std::abs( (*idr)->charge() )>0.001 //MB: GenParts have defined charge
	  ){
	maxPt = (*idr)->pt();
	part = (*idr);
      }
    }
    return part;
  }
  //////////////
    TLorentzVector getP4(const reco::GenParticleRef& part){
    return TLorentzVector(part->px(),part->py(),part->pz(),part->energy());
  }
  //////////////
  TLorentzVector getP4(const reco::GenParticle&  part){
    return TLorentzVector(part.px(),part.py(),part.pz(),part.energy());
  }
  //////////////
  void setP4Ptr(const TLorentzVector& p4, TLorentzVector *p4Ptr){
    if(p4Ptr) p4Ptr->SetPxPyPzE(p4.Px(),p4.Py(),p4.Pz(),p4.Energy());
    else p4Ptr = new TLorentzVector(p4);
  }
  //////////////
  void setV3Ptr(const TVector3& v3, TVector3 *v3Ptr){
    if(v3Ptr) v3Ptr->SetXYZ(v3.X(),v3.Y(),v3.Z());
    else v3Ptr = new TVector3(v3);
  }
  //////////////
  TLorentzVector getCombinedP4(const reco::GenParticleRefVector& products){    
    TLorentzVector p4;
    for(IGR id=products.begin(); id!=products.end(); ++id)
      p4 += getP4( *id );
    return p4;
  }
  //////////////
  TLorentzVector getChargedP4(const reco::GenParticleRefVector& products){    
    /*
    std::set<int> charged;
    charged.insert(211);//pi
    charged.insert(321);//K
    charged.insert(11);//e
    charged.insert(13);//mu
    */

    TLorentzVector p4;
    
    for(IGR idr = products.begin(); idr != products.end(); ++idr ){
      if( //charged.find( std::abs( (*idr)->pdgId() ) )!=charged.end() //MB: Logix used in pure Pythia code when charge not defined
	  std::abs( (*idr)->charge() )>0.001 //MB: GenParts have defined charge
	  ){
	p4 += getP4( *idr );
      }
    }
    return p4;
  }
  //////////////
  TLorentzVector getLeadChParticleP4(const reco::GenParticleRefVector& products){
    return getP4( getLeadChParticle(products) );
  }
  //////////////
  TVector3 getVertex(const reco::GenParticleRef& part){
    return TVector3(part->vx(),part->vy(),part->vz());
  }
  //////////////
  TVector3 getVertex(const reco::GenParticle& part){
    return TVector3(part.vx(),part.vy(),part.vz());
  }
  //////////////
  void getVertex(const reco::GenParticleRef& part, TVector3 *vtx){
    setV3Ptr(getVertex(part),vtx);
  }
  //////////////
  TVector3 impactParameter(const TVector3& pv, const TVector3& sv, const TLorentzVector& p4){
    TVector3 dir = (p4.Vect()).Unit();
    return (sv-pv) - ((sv-pv)*dir)*dir;
}
  //////////////
  void impactParameter(const TVector3& pv, const TVector3& sv, const TLorentzVector& p4,
		       TVector3 *ip){
    setV3Ptr(impactParameter(pv,sv,p4),ip);
  }
  //////////////
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
    
    float calO = direction * ( n1.Cross(n2) ); //defines phase
    float phi=TMath::ACos(n1*n2);
    if(calO<0)
      phi = 2.*TMath::Pi()-phi;
    float rho=TMath::ACos( (prod12Star.Vect().Unit() )*(prod22Star.Vect().Unit() ) );

    return std::make_pair(phi,rho);
  }
  //////////////
  //copy of function from PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h" which allows ignore pdgId
  void findParticles(const reco::GenParticleCollection& sourceParticles, 
		     reco::GenParticleRefVector& particleRefs, 
		     int pdgId, int status) {

    //one form status or pdgId has to be specifed!
    if(status<=0 && pdgId==0) return;

    unsigned index = 0;
    for(IG ig = sourceParticles.begin(); 
	ig!= sourceParticles.end(); ++ig, ++index) {
 
      const reco::GenParticle& gen = *ig;
    
      // status has been specified, and this one does not have the correct
      // status
      if(status>0 && gen.status()!=status ) continue;
    
      if(!pdgId || std::abs(gen.pdgId()) == pdgId ) {
	reco::GenParticleRef genref( &sourceParticles, index );
        particleRefs.push_back( genref );
      }
    }
  }

  //////////////
  bool isAncestor(const reco::Candidate *ancestor, 
		  const reco::Candidate *particle){
    //particle is already the ancestor
    if(ancestor == particle) return true;

    //otherwise loop on mothers, if any and return true if the ancestor is found
    for(size_t i=0; i<particle->numberOfMothers(); ++i){
      if( isAncestor(ancestor,particle->mother(i)) ) return true;
    }
    //if we did not return yet, then particle and ancestor are not relatives
    return false;
  }
  //////////////
  bool isDescendent(const reco::Candidate *descendent, 
		    const reco::Candidate *particle){
    //particle is already the descendent
    if(descendent == particle) return true;

    //otherwise loop on daughters, if any and return true if the descendent is found
    for(size_t i=0; i<particle->numberOfDaughters(); ++i){
      if( isDescendent(descendent,particle->daughter(i)) ) return true;
    }
    //if we did not return yet, then particle and ancestor are not relatives
    return false;
  }

  //////////////
  //copy of function from PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h" which allows ignore status
  void findDescendents(const reco::GenParticleRef& base, 
		       reco::GenParticleRefVector& descendents, 
		       int status, int pdgId ) {

    //one form status or pdgId has to be specifed!
    if(status<0 && pdgId==0) return;

    const reco::GenParticleRefVector& daughterRefs = base->daughterRefVector();
  
    for(IGR idr = daughterRefs.begin(); idr != daughterRefs.end(); ++idr ) {
      
      if( (status<0 || (*idr)->status() == status ) && 
          (!pdgId || std::abs((*idr)->pdgId()) == std::abs(pdgId) ) ) {
	
        descendents.push_back(*idr);
      }
      else 
        findDescendents( *idr, descendents, status, pdgId );
    }
  }
  //////////////
  void findAncestors(const reco::GenParticleRef& base, 
		     reco::GenParticleRefVector& ancestors, 
		     int status, int pdgId ) {

    //one form status or pdgId has to be specifed!
    if(status<0 && pdgId==0) return;

    const reco::GenParticleRefVector& motherRefs = base->motherRefVector();
  
    for(IGR idr = motherRefs.begin(); idr != motherRefs.end(); ++idr ) {
      
      if( (status<0 || (*idr)->status() == status ) && 
          (!pdgId || std::abs((*idr)->pdgId()) == std::abs(pdgId) ) ) {
	
        ancestors.push_back(*idr);
      }
      else 
        findAncestors( *idr, ancestors, status, pdgId );
    }
  }
  
  //////////////
  float getGenIso(const TLorentzVector& candidateP4,
		  const reco::GenParticleRefVector& particles,
		  double isoCone, double vetoCone, double minPt){

    float iso=0;

    for(unsigned int ipp=0; ipp<particles.size(); ++ipp){
      if(particles[ipp]->status()!=1 || particles[ipp]->numberOfDaughters()!=0) continue; //ignore unstable
      double dR2=deltaR2(candidateP4.Eta(), candidateP4.Phi(),
			 particles[ipp]->eta(),particles[ipp]->phi() );
      if(dR2<isoCone*isoCone && dR2>vetoCone*vetoCone && 
	 particles[ipp]->pt()>minPt)
	iso += particles[ipp]->pt();
    }

    return iso;
  }

  //////////////
  float getGenIso(const TLorentzVector& candidateP4,
		  const reco::GenParticleCollection&  particles,
		  double isoCone, double vetoCone, double minPt){
    
    reco::GenParticleRefVector final;
    findParticles(particles, final, 0, 1);
    return getGenIso(candidateP4,final,isoCone,vetoCone,minPt);
  }

} //namespace
