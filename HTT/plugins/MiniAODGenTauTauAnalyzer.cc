#include "MiniAODGenTauTauAnalyzer.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h"
#include "WarsawAnalysis/HTT/interface/GenInfoHelper.h"

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGaussQ.h"

// import LHEEventProduction definition
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

MiniAODGenTauTauAnalyzer::MiniAODGenTauTauAnalyzer(const edm::ParameterSet& iConfig):
  //prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  prunedGenToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
  toyTauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("toyTaus"))),
  verbose_(iConfig.getUntrackedParameter<bool>("verbose",false)),
  isSignal_(iConfig.getUntrackedParameter<bool>("isSignal",true)),
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
  bookVariable(tree_,"yToyMinus");
  bookVariable(tree_,"yToyPlus");
  bookVariable(tree_,"phiRho");
  bookVariable(tree_,"mHSmeared");
  bookVariable(tree_,"isoMinus");
  bookVariable(tree_,"isoPlus");
  bookVariable(tree_,"outerMinus");
  bookVariable(tree_,"outerPlus");
  bookVariable(tree_,"genWeight");
  bookVariable(tree_,"genHt");
  bookVariable(tree_,"genNJet");

  //Book integers
  tree_->Branch("bosonId",     &bosonId_,     "bosonId/I");
  tree_->Branch("decModeMinus",&decModeMinus_,"decModeMinus/I");
  tree_->Branch("decModePlus", &decModePlus_, "decModePlus/I");
  tree_->Branch("toyDecModeMinus",&toyDecModeMinus_,"toyDecModeMinus/I");
  tree_->Branch("toyDecModePlus" ,&toyDecModePlus_, "toyDecModePlus/I");
  tree_->Branch("toyNChargedMinus",&toyNChargedMinus_,"toyNChargedMinus/I");
  tree_->Branch("toyNChargedPlus", &toyNChargedPlus_, "toyNChargedPlus/I");
  tree_->Branch("toyNNeutralMinus",&toyNNeutralMinus_,"toyNNeutralMinus/I");
  tree_->Branch("toyNNeutralPlus", &toyNNeutralPlus_, "toyNNeutralPlus/I");
  tree_->Branch("toyNPiZeroMinus",&toyNPiZeroMinus_,"toyNPiZeroMinus/I");
  tree_->Branch("toyNPiZeroPlus", &toyNPiZeroPlus_, "toyNPiZeroPlus/I");
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
  toyPiMinus_ = new TLorentzVector();
  tree_->Branch("toyPiMinus.","TLorentzVector",toyPiMinus_);
  toyPiPlus_ = new TLorentzVector();
  tree_->Branch("toyPiPlus.","TLorentzVector",toyPiPlus_);
  toyPiZeroMinus_ = new TLorentzVector();
  tree_->Branch("toyPiZeroMinus.","TLorentzVector",toyPiZeroMinus_);
  toyPiZeroPlus_ = new TLorentzVector();
  tree_->Branch("toyPiZeroPlus.","TLorentzVector",toyPiZeroPlus_);
  toyNeutralMinus_ = new TLorentzVector();
  tree_->Branch("toyNeutralMinus.","TLorentzVector",toyNeutralMinus_);
  toyNeutralPlus_ = new TLorentzVector();
  tree_->Branch("toyNeutralPlus.","TLorentzVector",toyNeutralPlus_);
  toyTauMinus_ = new TLorentzVector();
  tree_->Branch("toyTauMinus.","TLorentzVector",toyTauMinus_);
  toyTauPlus_ = new TLorentzVector();
  tree_->Branch("toyTauPlus.","TLorentzVector",toyTauPlus_);

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
  tree_->Branch("nPiPlus.","TVector3",nPiPlus_);
  toySvMinus_ = new TVector3();
  tree_->Branch("toySvMinus.","TVector3",toySvMinus_);
  toySvPlus_ = new TVector3();
  tree_->Branch("toySvPlus.", "TVector3",toySvPlus_);
  toyNPiMinus_ = new TVector3();
  tree_->Branch("toyNPiMinus.","TVector3",toyNPiMinus_);
  toyNPiPlus_ = new TVector3();
  tree_->Branch("toyNPiPlus.","TVector3",toyNPiPlus_);

  if(isSignal_)
    std::cout<<"MiniAODGenTauTauAnalyzer: Running in signal mode"<<std::endl;
  else
    std::cout<<"MiniAODGenTauTauAnalyzer: Running in background mode"<<std::endl;
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
  delete toyPiMinus_;
  delete toyPiPlus_;
  delete toyPiZeroMinus_;
  delete toyPiZeroPlus_;
  delete toyNeutralMinus_;
  delete toyNeutralPlus_;
  delete toyTauMinus_;
  delete toyTauPlus_;

  delete thePV_;
  delete svMinus_;
  delete svPlus_;
  delete nPiMinus_;
  delete nPiPlus_; 
  delete toySvMinus_;
  delete toySvPlus_;
  delete toyNPiMinus_;
  delete toyNPiPlus_; 

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

  //patTaus build from toy PFParticles
  edm::Handle<pat::TauCollection> toyTaus;
  iEvent.getByToken(toyTauToken_, toyTaus);

  //Random service
  edm::Service<edm::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine& engine = rng->getEngine(iEvent.streamID());

  //let's try to find all status1 originating directly from a B meson decay 
  reco::GenParticleRef theBoson;
  reco::GenParticleRefVector taus;
  reco::GenParticleRefVector tauProdsPlus, tauProdsMinus;

  //*** gen information for weighting, stitching, etc.
  // event weight
  treeVars_["genWeight"] = 1;
  edm::Handle<GenEventInfoProduct> genEvt;
  iEvent.getByLabel("generator", genEvt);  
  if( genEvt.isValid() ){
    treeVars_["genWeight"] = genEvt->weight();
  }
  // access generator level HT, npartons
  edm::Handle<LHEEventProduct> lheEventProduct;  
  iEvent.getByLabel("externalLHEProducer", lheEventProduct);
  double lheHt = 0.;
  size_t nOutgoing = 0; // Number of outgoing partons
  if( lheEventProduct.isValid() ){
    const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
    std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
    size_t numParticles = lheParticles.size();
    for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
      int absPdgId = std::abs(lheEvent.IDUP[idxParticle]);
      int status = lheEvent.ISTUP[idxParticle];
      if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) ) { // quarks and gluons
	lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
	++nOutgoing;
      } 
    }
  }
  treeVars_["genHt"] = lheHt;
  treeVars_["genNJet"] = nOutgoing;  
  //***

  if(isSignal_){
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
		 <<", decMode w/ final: "<<WawGenInfoHelper::getTausDecays(taus[i],pp_tmp,true,false);
	std::cout<<", visPt="<<(WawGenInfoHelper::getCombinedP4(pp_tmp)).Pt() //it has to be in new cout to have an access to pp_tmp
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
  }//end of search for signal taus
  else{
    //need to set the PV
    for(size_t i=0; i<pruned->size(); ++i){
      reco::GenParticleRef genref(pruned,i);
      if( genref->isPromptFinalState() || genref->isPromptDecayed() ){
	//if( true ){//FIXME: in older MC statusFlags are not set??
	WawGenInfoHelper::getVertex(genref,thePV_);
	if(verbose_)
	  std::cout<<"Background vertex set"<<std::endl;
	break;
      }
    }
  }
  //
  
  WawGenInfoHelper::setP4Ptr(WawGenInfoHelper::getGenMet(*pruned), met_);

  /* ************** */
  std::pair<MiniAODGenTauTauAnalyzer::TVisTau,MiniAODGenTauTauAnalyzer::TVisTau> tauPair = 
    findTauPair(*pruned, *toyTaus);
  if( (tauPair.first).charge>0 )
    std::cout<<"Charge of negative toyTau>0!"
	     <<" charge="<<(tauPair.first).charge
	     <<" decMode="<<(tauPair.first).decMode
	     <<std::endl;
  if(!isSignal_ && 
     ( (tauPair.first).decMode == WawGenInfoHelper::tauDecayModes::tauDecayOther || 
       (tauPair.second).decMode == WawGenInfoHelper::tauDecayModes::tauDecayOther ) 
     ) return; //do not record non-signal events w/o a candaiate pair
  toyDecModeMinus_ = (tauPair.first).decMode;
  toyNChargedMinus_ = (tauPair.first).nCharged;
  toyNNeutralMinus_ = (tauPair.first).nNeutral;
  toyNPiZeroMinus_ = (tauPair.first).nPiZero;
  treeVars_["isoMinus"]=(tauPair.first).iso;
  treeVars_["outerMinus"]=(tauPair.first).outerPt;
  WawGenInfoHelper::setP4Ptr( (tauPair.first).p4, toyTauMinus_);
  WawGenInfoHelper::setP4Ptr( (tauPair.first).leadChP4, toyPiMinus_);
  WawGenInfoHelper::setP4Ptr( (tauPair.first).neutralP4, toyNeutralMinus_);
  WawGenInfoHelper::setP4Ptr( (tauPair.first).piZeroP4, toyPiZeroMinus_);
  WawGenInfoHelper::setV3Ptr( (tauPair.first).sv, toySvMinus_);
  WawGenInfoHelper::impactParameter(*thePV_, *toySvMinus_, *toyPiMinus_, toyNPiMinus_);
  treeVars_["yToyMinus"] = 2.* toyPiMinus_->Pt()/toyTauMinus_->Pt() - 1.;
  toyDecModePlus_ = (tauPair.second).decMode;
  toyNChargedPlus_ = (tauPair.second).nCharged;
  toyNNeutralPlus_ = (tauPair.second).nNeutral;
  toyNPiZeroPlus_ = (tauPair.second).nPiZero;
  treeVars_["isoPlus"]=(tauPair.second).iso;
  treeVars_["outerPlus"]=(tauPair.second).outerPt;
  WawGenInfoHelper::setP4Ptr( (tauPair.second).p4, toyTauPlus_);
  WawGenInfoHelper::setP4Ptr( (tauPair.second).leadChP4, toyPiPlus_);
  WawGenInfoHelper::setP4Ptr( (tauPair.second).neutralP4, toyNeutralPlus_);
  WawGenInfoHelper::setP4Ptr( (tauPair.second).piZeroP4, toyPiZeroPlus_);
  WawGenInfoHelper::setV3Ptr( (tauPair.second).sv, toySvPlus_);
  WawGenInfoHelper::impactParameter(*thePV_, *toySvPlus_, *toyPiPlus_, toyNPiPlus_);
  treeVars_["yToyPlus"] = 2.* toyPiPlus_->Pt()/toyTauPlus_->Pt() - 1.;

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
//////
void MiniAODGenTauTauAnalyzer::clean(){
  
  return;
}

//////
std::pair<MiniAODGenTauTauAnalyzer::TVisTau,MiniAODGenTauTauAnalyzer::TVisTau> 
MiniAODGenTauTauAnalyzer::findTauPair(const reco::GenParticleCollection&  particles,
				      const pat::TauCollection& taus){

  MiniAODGenTauTauAnalyzer::TVisTau visTau1={}, visTau2={};
  visTau1.decMode=visTau2.decMode=WawGenInfoHelper::tauDecayModes::tauDecayOther;
  //look for good leptons
  reco::GenParticleRefVector lep_all; 
  std::vector<reco::GenParticle> lep_iso;
  WawGenInfoHelper::findParticles(particles, lep_all, 11, 1);//electrons
  WawGenInfoHelper::findParticles(particles, lep_all, 13, 1);//muons
  for(WawGenInfoHelper::IGR idr = lep_all.begin(); idr != lep_all.end(); ++idr ){
    if(WawGenInfoHelper::getGenIso(WawGenInfoHelper::getP4( (*idr) ), particles,0.4,0.001,0.) < 0.5*(*idr)->pt() //loose iso
       && std::abs( (*idr)->eta() )<2.4)//tracker acceptance
      lep_iso.push_back(**idr);
  }
  std::sort(lep_iso.begin(), lep_iso.end(),
	    [](const reco::GenParticle& a, const reco::GenParticle& b){
	      return a.pt() > b.pt();
	    });
  // look for good loose pftaus
  std::vector<pat::Tau> tau_sel;
  for(unsigned int i=0; i<taus.size(); ++i){
    if(taus[i].tauID("decayModeFinding")>0.5 && taus[i].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits")<10){
      //clean pftaus from leptons
      bool overlap=false;
      for(unsigned int j=0; j<lep_iso.size(); ++j){
      	if(deltaR2(taus[i].eta(), taus[i].phi(),
      		   lep_iso[j].eta(), lep_iso[j].phi())<0.1*0.1 &&
      	   std::abs(lep_iso[j].pt()/taus[i].pt()-1.)<0.1 ){
	  if(verbose_)
	    std::cout<<"[MiniAODGenTauTauAnalyzer::findTauPair]: "
		     <<"Reject: lepPt="<<lep_iso[j].pt()<<", id: "<<lep_iso[j].pdgId()
		     <<", tauPt="<<taus[i].pt()<<", decMode: "<<taus[i].decayMode()
		     <<std::endl;
      	  overlap=true;
      	  break;
      	}
      }
      if(!overlap) tau_sel.push_back(taus[i]);
    }
  }
  std::sort(tau_sel.begin(), tau_sel.end(),
	    [](const pat::Tau& a, const pat::Tau& b){
	      return a.pt() > b.pt();
	    });
  /*
  for(unsigned int i=0; i<lep_iso.size(); ++i){
    float iso = WawGenInfoHelper::getGenIso(WawGenInfoHelper::getP4(lep_iso[i]),
                                                                    particles,
                                                                    0.4,0.001,0.);
    std::cout<<i<<". pt="<<lep_iso[i].pt()
	     <<", pt="<<lep_iso[i].eta()
	     <<", pdgId: "<<lep_iso[i].pdgId()
	     <<", charge="<<lep_iso[i].charge()
	     <<", iso="<<iso
	     <<std::endl;
  }
  for(unsigned int i=0; i<tau_sel.size(); ++i)
    std::cout<<i<<". "<<tau_sel[i].pt()
	     <<", decMode: "<<tau_sel[i].decayMode()
	     <<", charge="<<tau_sel[i].charge()
	     <<", iso="<<tau_sel[i].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits")
	     <<", outerPt/Pt="<<tau_sel[i].tauID("photonPtSumOutsideSignalCone")/tau_sel[i].pt()
	     <<std::endl;
  */
  
  //build ll pairs (R>0.3, OS, highest scalarPt)
  float scalarPt=0;
  for(unsigned int i=0; i<lep_iso.size(); ++i){
    for(unsigned int j=i+1; j<lep_iso.size(); ++j){
      if(deltaR2(lep_iso[i].eta(), lep_iso[i].phi(),
  		 lep_iso[j].eta(), lep_iso[j].phi())>0.3*0.3 &&
  	 lep_iso[i].charge()*lep_iso[j].charge()<0 && 
	 lep_iso[i].pt()+lep_iso[j].pt()>scalarPt){
  	scalarPt=lep_iso[i].pt()+lep_iso[j].pt();
	visTau1.decMode=(std::abs(lep_iso[i].pdgId())==11 ? 
			 WawGenInfoHelper::tauDecayModes::tauDecaysElectron :
			 WawGenInfoHelper::tauDecayModes::tauDecayMuon);
	visTau1.p4 = WawGenInfoHelper::getP4(lep_iso[i]);
	visTau1.leadChP4 = visTau1.p4;
	visTau1.sv = WawGenInfoHelper::getVertex(lep_iso[i]);
	visTau1.charge=lep_iso[i].charge();
	visTau1.iso=WawGenInfoHelper::getGenIso(visTau1.p4, particles,0.4,0.001,0.);
	visTau2.decMode=(std::abs(lep_iso[j].pdgId())==11 ? 
			 WawGenInfoHelper::tauDecayModes::tauDecaysElectron :
			 WawGenInfoHelper::tauDecayModes::tauDecayMuon);
	visTau2.p4 = WawGenInfoHelper::getP4(lep_iso[j]);
	visTau2.leadChP4 = visTau2.p4;
	visTau2.sv = WawGenInfoHelper::getVertex(lep_iso[j]);
	visTau2.charge=lep_iso[j].charge();
	visTau2.iso=WawGenInfoHelper::getGenIso(visTau2.p4, particles,0.4,0.001,0.);
      }
    }
  }
  for(unsigned int i=0; i<lep_iso.size(); ++i){
    for(unsigned int j=0; j<tau_sel.size(); ++j){
      if(deltaR2(lep_iso[i].eta(), lep_iso[i].phi(),
  		 tau_sel[j].eta(), tau_sel[j].phi())>0.3*0.3 &&
  	 lep_iso[i].charge()*tau_sel[j].charge()<0 && 
	 lep_iso[i].pt()+tau_sel[j].pt()>scalarPt){
  	scalarPt=lep_iso[i].pt()+tau_sel[j].pt();
	visTau1.decMode=(std::abs(lep_iso[i].pdgId())==11 ? 
			 WawGenInfoHelper::tauDecayModes::tauDecaysElectron :
			 WawGenInfoHelper::tauDecayModes::tauDecayMuon);
	visTau1.p4 = WawGenInfoHelper::getP4(lep_iso[i]);
	visTau1.leadChP4 = visTau1.p4;
	visTau1.sv = WawGenInfoHelper::getVertex(lep_iso[i]);
	visTau1.charge=lep_iso[i].charge();
	visTau1.iso=WawGenInfoHelper::getGenIso(visTau1.p4, particles,0.4,0.001,0.);
	visTau2.decMode=tau_sel[j].decayMode();
	visTau2.p4 = WawGenInfoHelper::getP4(tau_sel[j]);
	visTau2.leadChP4.SetPxPyPzE(tau_sel[j].leadPFChargedHadrCand()->px(),
				    tau_sel[j].leadPFChargedHadrCand()->py(),
				    tau_sel[j].leadPFChargedHadrCand()->pz(),
				    tau_sel[j].leadPFChargedHadrCand()->energy() );
	visTau2.nNeutral=0;
	visTau2.nPiZero=0;
	for(unsigned int iZ=0; iZ<tau_sel[j].signalPiZeroCandidates().size(); ++iZ){
	  
	  TLorentzVector piZeroP4(tau_sel[j].signalPiZeroCandidates()[iZ].px(),
				  tau_sel[j].signalPiZeroCandidates()[iZ].py(),
				  tau_sel[j].signalPiZeroCandidates()[iZ].pz(),
				  tau_sel[j].signalPiZeroCandidates()[iZ].energy());
	  visTau2.neutralP4 += piZeroP4; 
	  visTau2.nNeutral++;	      
	  if( piZeroQuality(tau_sel[j].signalPiZeroCandidates()[iZ], tau_sel[j]) ){
	    visTau2.piZeroP4 += piZeroP4;
	    visTau2.nPiZero++;	    
	  }
	}
	visTau2.sv.SetXYZ( tau_sel[j].leadPFChargedHadrCand()->trackRef()->vx(),
			   tau_sel[j].leadPFChargedHadrCand()->trackRef()->vy(),
			   tau_sel[j].leadPFChargedHadrCand()->trackRef()->vz() );
	visTau2.charge=tau_sel[j].charge();
	visTau2.iso=tau_sel[j].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
	visTau2.outerPt=tau_sel[j].tauID("photonPtSumOutsideSignalCone");
	visTau2.nCharged=tau_sel[j].signalPFChargedHadrCands().size();
      }
    }
  }
  for(unsigned int i=0; i<tau_sel.size(); ++i){
    for(unsigned int j=i+1; j<tau_sel.size(); ++j){
      if(deltaR2(tau_sel[i].eta(), tau_sel[i].phi(),
  		 tau_sel[j].eta(), tau_sel[j].phi())>0.3*0.3 &&
  	 tau_sel[i].charge()*tau_sel[j].charge()<0 && 
	 tau_sel[i].pt()+tau_sel[j].pt()>scalarPt){
  	scalarPt=tau_sel[i].pt()+tau_sel[j].pt();
	visTau1.decMode=tau_sel[i].decayMode();
	visTau1.p4 = WawGenInfoHelper::getP4(tau_sel[i]);
	visTau1.leadChP4.SetPxPyPzE(tau_sel[i].leadPFChargedHadrCand()->px(),
				    tau_sel[i].leadPFChargedHadrCand()->py(),
				    tau_sel[i].leadPFChargedHadrCand()->pz(),
				    tau_sel[i].leadPFChargedHadrCand()->energy() );
	visTau1.nNeutral=0;
	visTau1.nPiZero=0;
	for(unsigned int iZ=0; iZ<tau_sel[i].signalPiZeroCandidates().size(); ++iZ){
	  TLorentzVector piZeroP4(tau_sel[i].signalPiZeroCandidates()[iZ].px(),
				  tau_sel[i].signalPiZeroCandidates()[iZ].py(),
				  tau_sel[i].signalPiZeroCandidates()[iZ].pz(),
				  tau_sel[i].signalPiZeroCandidates()[iZ].energy());
	  visTau1.neutralP4 += piZeroP4; 
	  visTau1.nNeutral++;
	  if( piZeroQuality(tau_sel[i].signalPiZeroCandidates()[iZ], tau_sel[i]) ){
	    visTau1.neutralP4 += piZeroP4;
	    visTau1.nPiZero++;
	  }
	}
	visTau1.sv.SetXYZ( tau_sel[i].leadPFChargedHadrCand()->trackRef()->vx(),
			   tau_sel[i].leadPFChargedHadrCand()->trackRef()->vy(),
			   tau_sel[i].leadPFChargedHadrCand()->trackRef()->vz() );
	visTau1.charge=tau_sel[i].charge();
	visTau1.iso=tau_sel[i].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
	visTau1.outerPt=tau_sel[i].tauID("photonPtSumOutsideSignalCone");
	visTau1.nCharged=tau_sel[i].signalPFChargedHadrCands().size();
	visTau2.decMode=tau_sel[j].decayMode();
	visTau2.p4 = WawGenInfoHelper::getP4(tau_sel[j]);
	visTau2.leadChP4.SetPxPyPzE(tau_sel[j].leadPFChargedHadrCand()->px(),
				    tau_sel[j].leadPFChargedHadrCand()->py(),
				    tau_sel[j].leadPFChargedHadrCand()->pz(),
				    tau_sel[j].leadPFChargedHadrCand()->energy() );
	visTau2.nNeutral=0;
	visTau2.nPiZero=0;
	for(unsigned int iZ=0; iZ<tau_sel[j].signalPiZeroCandidates().size(); ++iZ){
	  TLorentzVector piZeroP4(tau_sel[j].signalPiZeroCandidates()[iZ].px(),
				  tau_sel[j].signalPiZeroCandidates()[iZ].py(),
				  tau_sel[j].signalPiZeroCandidates()[iZ].pz(),
				  tau_sel[j].signalPiZeroCandidates()[iZ].energy());
	  visTau2.neutralP4 += piZeroP4;
	  visTau2.nNeutral++;
	  if( piZeroQuality(tau_sel[j].signalPiZeroCandidates()[iZ], tau_sel[j]) ){
	    visTau2.piZeroP4 += piZeroP4;
	    visTau2.nPiZero++;
	  }
	}
	visTau2.sv.SetXYZ( tau_sel[j].leadPFChargedHadrCand()->trackRef()->vx(),
			   tau_sel[j].leadPFChargedHadrCand()->trackRef()->vy(),
			   tau_sel[j].leadPFChargedHadrCand()->trackRef()->vz() );
	visTau2.charge=tau_sel[j].charge();
	visTau2.iso=tau_sel[j].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
	visTau2.outerPt=tau_sel[j].tauID("photonPtSumOutsideSignalCone");
	visTau2.nCharged=tau_sel[j].signalPFChargedHadrCands().size();
      }
    }
  }
  if(verbose_)
    std::cout<<"[MiniAODGenTauTauAnalyzer::findTauPair]: "
	     <<"ScalarPt="<<scalarPt
	     <<", Pt1="<<visTau1.p4.Pt()
	     <<", decMode1="<<visTau1.decMode<<"("<<visTau1.charge<<")"
	     <<", Pt2="<<visTau2.p4.Pt()
	     <<", decMode2="<<visTau2.decMode<<"("<<visTau2.charge<<")"
	     <<std::endl;
  
  return ( visTau1.charge < 0 ? 
	   std::make_pair(visTau1,visTau2) :
	   std::make_pair(visTau2,visTau1) );
}

bool MiniAODGenTauTauAnalyzer::piZeroQuality(const reco::RecoTauPiZero& piZero, const pat::Tau& tau){

  const double signalConeSize = std::max(std::min(0.1, 3.0/tau.pt() ), 0.05);
  static const double minAbsPhotonSumPt_insideSignalCone = 2.5;
  static const double minRelPhotonSumPt_insideSignalCone = 0.1;
  static const double minAbsPhotonSumPt_outsideSignalCone = 1000000000.0;
  static const double minRelPhotonSumPt_outsideSignalCone = 1000000000.0;

  double photonSumPt_insideSignalCone = 0;
  double photonSumPt_outsideSignalCone = 0;
  int numPhotons = piZero.numberOfDaughters();
  for(int idxPhoton = 0; idxPhoton < numPhotons; ++idxPhoton) {
    const reco::Candidate* photon = piZero.daughter(idxPhoton);
    double dR2 = deltaR2(photon->p4(), tau.p4() );
    if(dR2 < signalConeSize*signalConeSize) {
      photonSumPt_insideSignalCone += photon->pt();
    } else {
      photonSumPt_outsideSignalCone += photon->pt();
    }
  }
  if ( photonSumPt_insideSignalCone  > minAbsPhotonSumPt_insideSignalCone  || 
       photonSumPt_insideSignalCone  > (minRelPhotonSumPt_insideSignalCone*tau.pt())  ||
       photonSumPt_outsideSignalCone > minAbsPhotonSumPt_outsideSignalCone || 
       photonSumPt_outsideSignalCone > (minRelPhotonSumPt_outsideSignalCone*tau.pt()) ) 
    return true;
  else
    return false;
}
