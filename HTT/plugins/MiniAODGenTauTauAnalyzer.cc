#include "MiniAODGenTauTauAnalyzer.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h"
#include "WarsawAnalysis/HTT/interface/GenInfoHelper.h"

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGaussQ.h"

MiniAODGenTauTauAnalyzer::MiniAODGenTauTauAnalyzer(const edm::ParameterSet& iConfig):
  //prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  prunedGenToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
  verbose_(iConfig.getUntrackedParameter<bool>("verbose",false)),
  bosonId_(0),
  decModeMinus_(WawGenInfoHelper::tauDecayModes::tauDecayOther),
  decModePlus_(WawGenInfoHelper::tauDecayModes::tauDecayOther)
{
  //edm::Service<edm::RandomNumberGenerator> rng; need to be called here? it is only for check its availability with ng.isAvailable()?

  edm::Service<TFileService> fs;
  tree_ = initTree(fs,"hTTCPTree");
  //Book floats
  bookVariable(tree_,"phi");
  bookVariable(tree_,"rho");
  bookVariable(tree_,"phi2");
  bookVariable(tree_,"yMinus");
  bookVariable(tree_,"yPlus");
  bookVariable(tree_,"yMinus2");
  bookVariable(tree_,"yPlus2");
  bookVariable(tree_,"yMinusLab");
  bookVariable(tree_,"yPlusLab");
  bookVariable(tree_,"yMinusLab2");
  bookVariable(tree_,"yPlusLab2");
  bookVariable(tree_,"phiRho");
  bookVariable(tree_,"mHSmeared");
  //Book integers
  tree_->Branch("bosonId",     &bosonId_,     "bosonId/I");
  tree_->Branch("decModeMinus",&decModeMinus_,"decModeMinus/I");
  tree_->Branch("decModePlus" ,&decModePlus_, "decModePlus/I");
  //Book TLorentzVectors
  p4Sum_ = new TLorentzVector();
  tree_->Branch("p4Sum.","TLorentzVector",p4Sum_);
  metNu_ = new TLorentzVector();
  tree_->Branch("metNu.","TLorentzVector",metNu_);
  met_ = new TLorentzVector();
  tree_->Branch("met.","TLorentzVector",met_);
  piMinus_ = new TLorentzVector();
  tree_->Branch("piMinus.","TLorentzVector",piMinus_);
  piPlus_ = new TLorentzVector();
  tree_->Branch("piPlus.","TLorentzVector",piPlus_);
  tauMinus_ = new TLorentzVector();
  tree_->Branch("tauMinus.","TLorentzVector",tauMinus_);
  tauPlus_ = new TLorentzVector();
  tree_->Branch("tauPlus.","TLorentzVector",tauPlus_);
  visTauMinus_ = new TLorentzVector();
  tree_->Branch("visTauMinus.","TLorentzVector",visTauMinus_);
  visTauPlus_ = new TLorentzVector();
  tree_->Branch("visTauPlus.","TLorentzVector",visTauPlus_);
  //Book TVector3
  thePV_ = new TVector3();
  tree_->Branch("thePV.","TVector3",thePV_);
  svMinus_ = new TVector3();
  tree_->Branch("svMinus.","TVector3",svMinus_);
  svPlus_ = new TVector3();
  tree_->Branch("svPlus.", "TVector3",svPlus_);
  nPiMinus_ = new TVector3();
  tree_->Branch("nPiMinus.","TVector3",nPiMinus_);
  nPiPlus_ = new TVector3();
  tree_->Branch("nPiPlus","TVector3",nPiPlus_);

}


MiniAODGenTauTauAnalyzer::~MiniAODGenTauTauAnalyzer(){

  delete p4Sum_;
  delete metNu_;
  delete met_;
  delete piMinus_;
  delete piPlus_;
  delete tauMinus_;
  delete tauPlus_;
  delete visTauMinus_;
  delete visTauPlus_;

  delete thePV_;
  delete svMinus_;
  delete svPlus_;
  delete nPiMinus_;
  delete nPiPlus_; 

}

void MiniAODGenTauTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  using namespace reco;
  using namespace pat;

  //clean variables
  clean();//FIXME: to be implemented

  // Basic event info
  treeVars_["run"]   = iEvent.id().run();
  treeVars_["lumi"]  = iEvent.id().luminosityBlock();
  treeVars_["event"] = iEvent.id().event();

  // Pruned particles are the one containing "important" stuff
  //edm::Handle<edm::View<reco::GenParticle> > pruned;
  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(prunedGenToken_, pruned);

  // Packed particles are all the status 1, so usable to remake jets
  // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
  edm::Handle<edm::View<pat::PackedGenParticle> > packed;
  iEvent.getByToken(packedGenToken_, packed);

  //Random service
  edm::Service<edm::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine& engine = rng->getEngine(iEvent.streamID());

  //let's try to find all status1 originating directly from a B meson decay 
  reco::GenParticleRef theBoson;
  reco::GenParticleRefVector taus;
  reco::GenParticleRefVector tauProdsPlus, tauProdsMinus;

  for(size_t i=0; i<pruned->size(); ++i){
    if( theBoson.isNonnull() ) break;
    reco::GenParticleRef genref(pruned,i);
    if(theBoson.isNull() && WawGenInfoHelper::isBoson(genref,false,true) ){
      theBoson = WawGenInfoHelper::getFinalClone(genref,true);//.get();
      if(verbose_){
	std::cout<<"Bos_init: "<<std::endl
		 <<"\t pt="<<(*pruned)[i].pt()
		 <<", eta="<<(*pruned)[i].eta()
		 <<", phi="<<(*pruned)[i].phi()
		 <<", pdgId="<<(*pruned)[i].pdgId()<<std::endl;
	std::cout<<"Bos_fin: "<<std::endl
		 <<"\t pt="<<theBoson->pt()
		 <<", eta="<<theBoson->eta()
		 <<", phi="<<theBoson->phi()
		 <<", pdgId="<<theBoson->pdgId()<<std::endl;
      }
    }  
  }
  if(theBoson.isNull() ) return; //FIXME: need to handle cases w/o the boson?
  bosonId_ = theBoson->pdgId();
  WawGenInfoHelper::getVertex(theBoson,thePV_);
  reco::GenParticleRefVector taus_all;
  WawGenInfoHelper::findParticles(*pruned, taus_all, 15, -1);
  for(WawGenInfoHelper::IGR idr = taus_all.begin(); idr != taus_all.end(); ++idr ){
    if(WawGenInfoHelper::isFinalClone((*idr),true) && WawGenInfoHelper::isAncestor(theBoson.get(),(*idr).get()) )
      taus.push_back(*idr);
  }
  if(verbose_){
    std::cout<<"Taus: "<<std::endl;
    for(size_t i=0; i<taus.size(); ++i){
      reco::GenParticleRefVector pp_tmp;
      std::cout<<"\t"<<i<<". "
	       <<"pt="<<taus[i]->pt()
	       <<", eta="<<taus[i]->eta()
	       <<", phi="<<taus[i]->phi()
	       <<", charge="<<taus[i]->charge()
	       <<", pdgId="<<taus[i]->pdgId()
	       <<", decMode w/ direct: "<<WawGenInfoHelper::getTausDecays(taus[i],pp_tmp,true,true)
	       <<", decMode w/ final: "<<WawGenInfoHelper::getTausDecays(taus[i],pp_tmp,true,false)
	       <<std::endl;
    }
  }
  //Less or more than 2 taus from decay of the boson???
  if(taus.size() != 2) return; //FIXME: need to handle cases w/ #taus>2?
  if(taus[0]->charge()*taus[1]->charge()>0) return; //FIXME: need to handle cases w/ same sign
  /* ****************** */
  //order tau pair: nagative charge first
  if(taus[0]->charge()>0){
    reco::GenParticleRefVector taus_tmp;
    taus_tmp.push_back(taus[1]);
    taus_tmp.push_back(taus[0]);
    taus.swap(taus_tmp);
  }
  decModeMinus_ = WawGenInfoHelper::getTausDecays(taus[0],tauProdsMinus,true,false);
  tauMinus_->SetPxPyPzE(taus[0]->px(),taus[0]->py(),taus[0]->pz(),taus[0]->energy());
  reco::GenParticleRef piMinusRef = WawGenInfoHelper::getLeadChParticle(tauProdsMinus);
  piMinus_->SetPxPyPzE(piMinusRef->px(),piMinusRef->py(),piMinusRef->pz(),piMinusRef->energy());
  WawGenInfoHelper::setP4Ptr(WawGenInfoHelper::getCombinedP4(tauProdsMinus),
			     visTauMinus_);
  WawGenInfoHelper::getVertex(WawGenInfoHelper::getLeadChParticle(tauProdsMinus),svMinus_);
  WawGenInfoHelper::impactParameter(*thePV_, *svMinus_, *piMinus_, nPiMinus_);
  decModePlus_ = WawGenInfoHelper::getTausDecays(taus[1],tauProdsPlus,true,false);
  tauPlus_->SetPxPyPzE(taus[1]->px(),taus[1]->py(),taus[1]->pz(),taus[1]->energy());
  reco::GenParticleRef piPlusRef = WawGenInfoHelper::getLeadChParticle(tauProdsPlus);
  piPlus_->SetPxPyPzE(piPlusRef->px(),piPlusRef->py(),piPlusRef->pz(),piPlusRef->energy());
  WawGenInfoHelper::setP4Ptr(WawGenInfoHelper::getCombinedP4(tauProdsPlus),
			     visTauPlus_);
  WawGenInfoHelper::getVertex(WawGenInfoHelper::getLeadChParticle(tauProdsPlus),svPlus_);
  WawGenInfoHelper::impactParameter(*thePV_, *svPlus_, *piPlus_, nPiPlus_);
  WawGenInfoHelper::setP4Ptr((*tauMinus_)+(*tauPlus_), p4Sum_);
  treeVars_["mHSmeared"] = CLHEP::RandGaussQ::shoot(&engine,
						    p4Sum_->M(),
						    p4Sum_->M()*0.15);
  WawGenInfoHelper::setP4Ptr(WawGenInfoHelper::getTauNuMet(taus), metNu_); 
  WawGenInfoHelper::setP4Ptr(WawGenInfoHelper::getGenMet(*pruned), met_);
  //angles for leading charged particle in tau rest frame
  std::pair<float,float> anglesTau = 
    WawGenInfoHelper::angleBetweenPlanes(*tauMinus_, *piMinus_,
					 *tauPlus_, *piPlus_);
  treeVars_["phi"]=anglesTau.first;
  treeVars_["rho"]=anglesTau.second;
  //angles for leading charged particle in charged leading charged rest frame with impact params
  std::pair<float,float> anglesPi = 
    WawGenInfoHelper::angleBetweenPlanes(*piMinus_, TLorentzVector(*nPiMinus_,0),
					 *piPlus_, TLorentzVector(*nPiPlus_,0));
  treeVars_["phi2"]=anglesPi.first;
  //angles for visible tau - leading charge particle system in visible rest frame
  std::pair<float,float> anglesRho = 
    /* Old definition in the rho-rho rest frame
    WawGenInfoHelper::angleBetweenPlanes(*visTauMinus_, *piMinus_,
					 *visTauPlus_, *piPlus_);
    */
    /* New definition in the leading charged partices frame */
    WawGenInfoHelper::angleBetweenPlanes(*piMinus_, (*visTauMinus_)-(*piMinus_),
					 *piPlus_, (*visTauPlus_)-(*piPlus_) );

  treeVars_["phiRho"]=anglesRho.first;
  //energy difference beteen leadig charged and visible tau in different frames
  // tau frames (ideal)
  TVector3 bTauPlus = tauPlus_->BoostVector();
  TVector3 bTauMinus = tauMinus_->BoostVector();
  TLorentzVector piMinusTau, piZero1Tau,
        piPlusTau, piZero2Tau;
  piMinusTau = *piMinus_;
  piMinusTau.Boost(-bTauMinus);
  piPlusTau = *piPlus_;
  piPlusTau.Boost(-bTauPlus);
  piZero1Tau = (*visTauMinus_ - *piMinus_);
  piZero1Tau.Boost(-bTauMinus);
  piZero2Tau = (*visTauPlus_ - *piPlus_);
  piZero2Tau.Boost(-bTauPlus);
  treeVars_["yMinus"] = (piMinusTau.E()-piZero1Tau.E())/(piMinusTau.E()+piZero1Tau.E());
  treeVars_["yPlus"] = (piPlusTau.E()-piZero2Tau.E())/(piPlusTau.E()+piZero2Tau.E());
  // pseudo-tau frames (Was & Worek)
  TVector3 bRhoRho = (*visTauMinus_ + *visTauPlus_).BoostVector();
  TLorentzVector visTauMinusStar, visTauPlusStar;
  visTauMinusStar = *visTauMinus_;
  visTauMinusStar.Boost(-bRhoRho);
  visTauPlusStar = *visTauPlus_;
  visTauPlusStar.Boost(-bRhoRho);
  TLorentzVector pseudoTauMinus, pseudoTauPlus;
  TVector3 pseudoDirMinus = visTauMinusStar.Vect().Unit();
  float mHover2 = p4Sum_->M()/2.;//FIXME: add mass smearing
  //float mHover2 = treeVars_["mHSmeared"]/2;
  float mTau = 1.777;
  pseudoDirMinus.SetMag(sqrt(mHover2*mHover2 - mTau*mTau));//mH/2 and tau mass
  pseudoTauMinus.SetVectM(pseudoDirMinus,mTau);
  TVector3 pseudoDirPlus = visTauPlusStar.Vect().Unit();
  pseudoDirPlus.SetMag(sqrt(mHover2*mHover2 - mTau*mTau));//tau mass
  pseudoTauPlus.SetVectM(pseudoDirPlus,mTau);
  TVector3 bVisTauPlus = pseudoTauPlus.BoostVector();
  TVector3 bVisTauMinus = pseudoTauMinus.BoostVector();
  TLorentzVector piMinusVisTau, piZeroMinusVisTau,
    piPlusVisTau, piZeroPlusVisTau;
  piMinusVisTau = *piMinus_;
  piMinusVisTau.Boost(-bVisTauMinus);
  piPlusVisTau = *piPlus_;
  piPlusVisTau.Boost(-bVisTauPlus);
  piZeroMinusVisTau = (*visTauMinus_ - *piMinus_);
  piZeroMinusVisTau.Boost(-bVisTauMinus);
  piZeroPlusVisTau = (*visTauPlus_ - *piPlus_);
  piZeroPlusVisTau.Boost(-bVisTauPlus);
  treeVars_["yMinus2"]=(piMinusVisTau.E()-piZeroMinusVisTau.E())/(piMinusVisTau.E()+piZeroMinusVisTau.E());
  treeVars_["yPlus2"]=(piPlusVisTau.E()-piZeroPlusVisTau.E())/(piPlusVisTau.E()+piZeroPlusVisTau.E());
  //Lab frame
  treeVars_["yMinusLab"] = (piMinus_->E() - (*visTauMinus_ - *piMinus_).E()) / (piMinus_->E() + (*visTauMinus_ - *piMinus_).E());
  treeVars_["yPlusLab"] = (piPlus_->E() - (*visTauPlus_ - *piPlus_).E()) / (piPlus_->E() + (*visTauPlus_ - *piPlus_).E());
  //Lab frame using Pt
  treeVars_["yMinusLab2"] = 2.* piMinus_->Pt()/visTauMinus_->Pt() - 1.;
  treeVars_["yPlusLab2"] = 2.* piPlus_->Pt()/visTauPlus_->Pt() - 1.;

  /* ************** */

  /*
    if(abs((*pruned)[i].pdgId()) > 500 && abs((*pruned)[i].pdgId()) <600){
      const Candidate * bMeson = &(*pruned)[i];
      std::cout << "PdgID: " << bMeson->pdgId() << " pt " << bMeson->pt() << " eta: " << bMeson->eta() << " phi: " << bMeson->phi() << std::endl;
      std::cout << "  found daugthers: " << std::endl;
      for(size_t j=0; j<packed->size();j++){
	//get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection 
	const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
	if(motherInPrunedCollection != nullptr && isAncestor( bMeson , motherInPrunedCollection)){
	  std::cout << "     PdgID: " << (*packed)[j].pdgId() << " pt " << (*packed)[j].pt() << " eta: " << (*packed)[j].eta() << " phi: " << (*packed)[j].phi() << std::endl;
	}
      }
    }
  }
  */
  //at the end fill tree 
  tree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void MiniAODGenTauTauAnalyzer::beginJob(){}

// ------------ method called once each job just after ending the event loop  ------------
void MiniAODGenTauTauAnalyzer::endJob(){}

//////
TTree * MiniAODGenTauTauAnalyzer::initTree(edm::Service<TFileService> &fs, std::string name)
{
  TTree *aTree = fs->make<TTree>( name.c_str(), name.c_str() );

  treeVars_["run"];
  aTree->Branch("run", &treeVars_["run"], "run/F");
  treeVars_["lumi"];
  aTree->Branch("lumi", &treeVars_["lumi"], "lumi/F");
  treeVars_["event"];
  aTree->Branch("event", &treeVars_["event"], "event/F");

  return aTree;
}

//////
void MiniAODGenTauTauAnalyzer::bookVariable(TTree *t, std::string var)
{
  if(!t) return;
  treeVars_[var];
  t->Branch(var.c_str(), &treeVars_[var], (var + "/F").c_str());

  return;
}
void MiniAODGenTauTauAnalyzer::clean(){
  
  return;
}



