// M. Bluj
#include "MiniAODVertexAnalyzer.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/Provenance/interface/ProductID.h"
//
#include "PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h"
#include "WarsawAnalysis/HTT/interface/GenInfoHelper.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "DataFormats/Math/interface/deltaR.h"

//#include "DataFormats/GeometryVector/interface/LocalPoint.h"


MiniAODVertexAnalyzer::MiniAODVertexAnalyzer(const edm::ParameterSet & iConfig) :
  scores_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("vertexScores"))),
  cands_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("src"))),
  lostCands_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("lostSrc"))),
  lostCandsEle_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("lostEleSrc"))),
  genCands_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genSrc"))),
  vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  bs_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))), 
  taus_(consumes<std::vector<pat::Tau> >(iConfig.getParameter<edm::InputTag>("taus"))),
  mus_(consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"))),
  eles_(consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
  useBeamSpot_(iConfig.getParameter<bool>("useBeamSpot")),
  useLostCands_(iConfig.getParameter<bool>("useLostCands")),
  useTauTracks_(iConfig.getUntrackedParameter<bool>("useTauTracks",false)),
  verbose_(iConfig.getUntrackedParameter<bool>("verbose",false)){

  
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("HTT","HTT");

  myEvent_ = 0; //Important pointer initialisation
  tree_->Branch("HTTEvent", &myEvent_); 

  thePair_ = std::make_pair<const reco::Candidate*, const reco::Candidate*>(0,0); 
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
MiniAODVertexAnalyzer::~MiniAODVertexAnalyzer(){ 

  delete myEvent_;

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void MiniAODVertexAnalyzer::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  myEvent_->clear();
  thePair_ = std::make_pair<const reco::Candidate*, const reco::Candidate*>(0,0);
  
  // Basic event info
  myEvent_->run_ = iEvent.id().run();
  myEvent_->lumi_ = iEvent.id().luminosityBlock();
  myEvent_->event_ = iEvent.id().event();

  bool goodGenData = getGenData(iEvent, iSetup);
  bool goodRecoData = getRecoData(iEvent, iSetup);

  //if(goodGenData && goodRecoData) tree_->Fill();
  if(goodRecoData) tree_->Fill();

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool MiniAODVertexAnalyzer::getGenData(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  edm::Handle<reco::GenParticleCollection> genCands;
  iEvent.getByToken(genCands_, genCands);
  
 //let's try to find initial boson and taus from its decays 
  reco::GenParticleRef theBoson;
  reco::GenParticleRefVector taus;
  reco::GenParticleRefVector tauProdsPlus, tauProdsMinus;

  for(size_t i=0; i<genCands->size(); ++i){
    if( theBoson.isNonnull() ) break;
    reco::GenParticleRef genref(genCands,i);
    if(theBoson.isNull() && WawGenInfoHelper::isBoson(genref,false,true) ){
      theBoson = WawGenInfoHelper::getFinalClone(genref,true);//.get();
      if(verbose_){
        std::cout<<"Bos_init: "<<std::endl
                 <<"\t pt="<<(*genCands)[i].pt()
                 <<", eta="<<(*genCands)[i].eta()
                 <<", phi="<<(*genCands)[i].phi()
                 <<", pdgId="<<(*genCands)[i].pdgId()<<std::endl;
        std::cout<<"Bos_fin: "<<std::endl
                 <<"\t pt="<<theBoson->pt()
                 <<", eta="<<theBoson->eta()
                 <<", phi="<<theBoson->phi()
                 <<", pdgId="<<theBoson->pdgId()<<std::endl;
      }
    }  
  }
  if(theBoson.isNull() ) return false; //FIXME: need to handle cases w/o the boson?
  myEvent_->bosonId_ = theBoson->pdgId();
  WawGenInfoHelper::getVertex(theBoson,&myEvent_->genEvent_.thePV_);

  reco::GenParticleRefVector taus_all;
  WawGenInfoHelper::findParticles(*genCands, taus_all, 15, -1);
  for(WawGenInfoHelper::IGR idr = taus_all.begin(); idr != taus_all.end(); ++idr ){
    if(WawGenInfoHelper::isFinalClone((*idr),true) && WawGenInfoHelper::isAncestor(theBoson.get(),(*idr).get()) )
      //if(WawGenInfoHelper::isFinalClone((*idr),true) )
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
	       <<", energy="<<taus[i]->energy()	
               <<", charge="<<taus[i]->charge()
               <<", pdgId="<<taus[i]->pdgId()
               <<", decMode w/ direct: "<<WawGenInfoHelper::getTausDecays(taus[i],pp_tmp,true,true)
               <<", decMode w/ final: "<<WawGenInfoHelper::getTausDecays(taus[i],pp_tmp,true,false)
               <<std::endl;
    }
  }
  //Less or more than 2 taus from decay of the boson???
  if(taus.size() != 2) return false; //FIXME: need to handle cases w/ #taus>2?
  if(taus[0]->charge()*taus[1]->charge()>0) return false; //FIXME: need to handle cases w/ same sign
   //order tau pair: nagative charge first
  if(taus[0]->charge()>0){
    reco::GenParticleRefVector taus_tmp;
    taus_tmp.push_back(taus[1]);
    taus_tmp.push_back(taus[0]);
    taus.swap(taus_tmp);
  }
  
  //*** Tau minus
  myEvent_->genEvent_.decModeMinus_ = WawGenInfoHelper::getTausDecays(taus[0],tauProdsMinus,true,false); 
  myEvent_->genEvent_.tauMinus_.SetPtEtaPhiE(taus[0]->pt(),taus[0]->eta(),taus[0]->phi(),taus[0]->energy());
  
  reco::GenParticleRef piMinusRef = WawGenInfoHelper::getLeadChParticle(tauProdsMinus);
  myEvent_->genEvent_.piMinus_.SetPtEtaPhiE(piMinusRef->pt(),piMinusRef->eta(),piMinusRef->phi(),piMinusRef->energy());

  WawGenInfoHelper::setP4Ptr(WawGenInfoHelper::getNeutralP4(tauProdsMinus), &myEvent_->genEvent_.pi0Minus_);
  WawGenInfoHelper::setP4Ptr(WawGenInfoHelper::getCombinedP4(tauProdsMinus), &myEvent_->genEvent_.visTauMinus_);
  WawGenInfoHelper::getVertex(WawGenInfoHelper::getLeadChParticle(tauProdsMinus),&myEvent_->genEvent_.svMinus_);
  WawGenInfoHelper::impactParameter(myEvent_->genEvent_.thePV_, myEvent_->genEvent_.svMinus_,
				    myEvent_->genEvent_.piMinus_, &myEvent_->genEvent_.nPiMinus_);
  myEvent_->genEvent_.dzMinus_ = (myEvent_->genEvent_.svMinus_.Z() - myEvent_->genEvent_.thePV_.Z()) - ((myEvent_->genEvent_.svMinus_.X() - myEvent_->genEvent_.thePV_.X()) * myEvent_->genEvent_.piMinus_.Px() + (myEvent_->genEvent_.svMinus_.Y() - myEvent_->genEvent_.thePV_.Y()) * myEvent_->genEvent_.piMinus_.Py()) / myEvent_->genEvent_.piMinus_.Pt() * myEvent_->genEvent_.piMinus_.Pz() / myEvent_->genEvent_.piMinus_.Pt();

  //*** Tau plus
  myEvent_->genEvent_.decModePlus_ = WawGenInfoHelper::getTausDecays(taus[1],tauProdsPlus,true,false);
  myEvent_->genEvent_.tauPlus_.SetPtEtaPhiE(taus[1]->pt(),taus[1]->eta(),taus[1]->phi(),taus[1]->energy());

  reco::GenParticleRef piPlusRef = WawGenInfoHelper::getLeadChParticle(tauProdsPlus);
  myEvent_->genEvent_.piPlus_.SetPtEtaPhiE(piPlusRef->pt(),piPlusRef->eta(),piPlusRef->phi(),piPlusRef->energy());

  WawGenInfoHelper::setP4Ptr(WawGenInfoHelper::getNeutralP4(tauProdsPlus), &myEvent_->genEvent_.pi0Plus_);
  WawGenInfoHelper::setP4Ptr(WawGenInfoHelper::getCombinedP4(tauProdsPlus), &myEvent_->genEvent_.visTauPlus_);
  WawGenInfoHelper::getVertex(WawGenInfoHelper::getLeadChParticle(tauProdsPlus),&myEvent_->genEvent_.svPlus_);
  WawGenInfoHelper::impactParameter(myEvent_->genEvent_.thePV_, myEvent_->genEvent_.svPlus_,
				    myEvent_->genEvent_.piPlus_, &myEvent_->genEvent_.nPiPlus_);
  myEvent_->genEvent_.dzPlus_ = (myEvent_->genEvent_.svPlus_.Z() - myEvent_->genEvent_.thePV_.Z()) - ((myEvent_->genEvent_.svPlus_.X() - myEvent_->genEvent_.thePV_.X()) * myEvent_->genEvent_.piPlus_.Px() + (myEvent_->genEvent_.svPlus_.Y() - myEvent_->genEvent_.thePV_.Y()) * myEvent_->genEvent_.piPlus_.Py()) / myEvent_->genEvent_.piPlus_.Pt() * myEvent_->genEvent_.piPlus_.Pz() / myEvent_->genEvent_.piPlus_.Pt();
  
  WawGenInfoHelper::setP4Ptr(myEvent_->genEvent_.tauMinus_+myEvent_->genEvent_.tauPlus_,
			     &myEvent_->genEvent_.p4Sum_);

  return true;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool  MiniAODVertexAnalyzer::getRecoData(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  bool hasAnyVx = findPrimaryVertices(iEvent, iSetup);
  if(!hasAnyVx) return false;

  bool hasRecoTaus = findRecoTaus(iEvent, iSetup);

  refitPV(iEvent, iSetup);
  setPCAVectors(iEvent, iSetup);
  getPVBalance(iEvent, 0);

  return hasAnyVx && hasRecoTaus;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool  MiniAODVertexAnalyzer::refitPV(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  edm::Handle<edm::View<pat::PackedCandidate> >  cands;
  iEvent.getByToken(cands_, cands);
  edm::Handle<edm::View<pat::PackedCandidate> >  lostCands;
  if(useLostCands_) iEvent.getByToken(lostCands_, lostCands);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertices_, vertices);
  edm::Handle<reco::BeamSpot> beamSpot;
  if(useBeamSpot_) iEvent.getByToken(bs_, beamSpot);

  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);
 
  TransientVertex transVtx, transVtxNoBS;

  //Get tracks associated with pfPV
  //int vtxIdx = myEvent_->recoEvent_.pfPVIndex_; //pfPV
  int vtxIdx = 0; //AOD PV
  reco::TrackCollection pvTracks;
  std::vector<reco::Candidate::LorentzVector> signalTauCandsP4;
  if(!useTauTracks_ && thePair_.first!=nullptr && thePair_.second!=nullptr){
    //t_1
    const pat::Tau* aTau = dynamic_cast<const pat::Tau*>(thePair_.first);
    if(aTau!=nullptr){
      for(size_t j=0 ; j<aTau->signalChargedHadrCands().size(); ++j){
	signalTauCandsP4.push_back(aTau->signalChargedHadrCands()[j]->p4());
      }
    }
    else{
      const pat::Muon* aMu = dynamic_cast<const pat::Muon*>(thePair_.first);
      if(aMu!=nullptr && aMu->originalObjectRef().isNonnull()){
	signalTauCandsP4.push_back(aMu->originalObjectRef()->p4());
      }
      else{
	const pat::Electron* aEle = dynamic_cast<const pat::Electron*>(thePair_.first);
	if(aEle!=nullptr && aEle->associatedPackedPFCandidates().size()>0){
	  reco::Candidate::LorentzVector aEleP4;
	  for(size_t i=0; i<aEle->associatedPackedPFCandidates().size(); ++i){
	    if(aEle->associatedPackedPFCandidates()[i]->pdgId()==aEle->pdgId() &&
	       aEle->associatedPackedPFCandidates()[i]->pt()>aEleP4.pt()){
	      aEleP4 = aEle->associatedPackedPFCandidates()[i]->p4();
	    }
	  }
	  if(aEleP4.pt()>0.5)
	    signalTauCandsP4.push_back(aEleP4);
	}
      }
    }
    //t_2
    aTau = dynamic_cast<const pat::Tau*>(thePair_.second);
    if(aTau!=nullptr){
      for(size_t j=0 ; j<aTau->signalChargedHadrCands().size(); ++j){
	signalTauCandsP4.push_back(aTau->signalChargedHadrCands()[j]->p4());
      }
    }
    else{
      const pat::Muon* aMu = dynamic_cast<const pat::Muon*>(thePair_.second);
      if(aMu!=nullptr && aMu->originalObjectRef().isNonnull()){
	signalTauCandsP4.push_back(aMu->originalObjectRef()->p4());
      }
      else{
	const pat::Electron* aEle = dynamic_cast<const pat::Electron*>(thePair_.second);
	if(aEle!=nullptr && aEle->associatedPackedPFCandidates().size()>0){
	  reco::Candidate::LorentzVector aEleP4;
	  for(size_t i=0; i<aEle->associatedPackedPFCandidates().size(); ++i){
	    if(aEle->associatedPackedPFCandidates()[i]->pdgId()==aEle->pdgId() &&
	       aEle->associatedPackedPFCandidates()[i]->pt()>aEleP4.pt()){
	      aEleP4 = aEle->associatedPackedPFCandidates()[i]->p4();
	    }
	  }
	  if(aEleP4.pt()>0.5)
	    signalTauCandsP4.push_back(aEleP4);
	}
      }
    }
  }
  for(size_t i=0; i<cands->size(); ++i){
    if(!(*cands)[i].hasTrackDetails() || (*cands)[i].vertexRef().isNull()) continue;
    ///Skip tracks comming from tau decay.
    bool skipTrack = false;
    for(size_t j=0 ; j<signalTauCandsP4.size();++j){
      if(deltaR2(signalTauCandsP4[j],(*cands)[i].p4())<0.001*0.001
	 && std::abs(signalTauCandsP4[j].pt()/(*cands)[i].pt()-1.)<0.001
	 ){
	skipTrack = true;
	break;
      }
    }
    if(skipTrack) continue;
    int key = (*cands)[i].vertexRef().key();
    int quality = (*cands)[i].pvAssociationQuality();
    if(key!=vtxIdx ||
       (quality!=pat::PackedCandidate::UsedInFitTight &&
	quality!=pat::PackedCandidate::UsedInFitLoose)) continue;

    pvTracks.push_back(*((*cands)[i].bestTrack()));
  }
  if(useLostCands_){
    for(size_t i=0; i<lostCands->size(); ++i){
      if(!(*lostCands)[i].hasTrackDetails() || (*lostCands)[i].vertexRef().isNull()) continue;
      ///Skip tracks comming from tau decay.
      bool skipTrack = false;
      for(size_t j=0 ; j<signalTauCandsP4.size();++j){
	if(deltaR2(signalTauCandsP4[j],(*lostCands)[i].p4())<0.001*0.001
	   && std::abs(signalTauCandsP4[j].pt()/(*lostCands)[i].pt()-1.)<0.001
	   ){
	  skipTrack = true;
	  break;
	}
      }
      if(skipTrack) continue;
      int key = (*lostCands)[i].vertexRef().key();
      int quality = (*lostCands)[i].pvAssociationQuality();
      if(key!=vtxIdx ||
	 (quality!=pat::PackedCandidate::UsedInFitTight &&
	  quality!=pat::PackedCandidate::UsedInFitLoose)) continue;

      pvTracks.push_back(*((*lostCands)[i].bestTrack()));
    }
  }
  ///Built transient tracks from tracks.
  std::vector<reco::TransientTrack> transTracks;  
  for(auto iter: pvTracks) transTracks.push_back(transTrackBuilder->build(iter));
  myEvent_->recoEvent_.nTracksInRefit_ = transTracks.size();

  bool fitOk = false;  
  if(transTracks.size() >= 2 ) {
    AdaptiveVertexFitter avf;
    avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception    
    try {
      transVtxNoBS = avf.vertex(transTracks);      
      if(useBeamSpot_)
	transVtx = avf.vertex(transTracks, *beamSpot);
      else
	transVtx = transVtxNoBS;
      fitOk = true; 
    } catch (...) {
      fitOk = false; 
      std::cout<<"Vtx fit failed!"<<std::endl;
    }
  }
  
  if(fitOk && transVtx.isValid() && transVtxNoBS.isValid()) { 
    myEvent_->recoEvent_.refitPfPV_.SetXYZ(transVtx.position().x(),transVtx.position().y(),transVtx.position().z());
    myEvent_->recoEvent_.refitPfPVNoBS_.SetXYZ(transVtxNoBS.position().x(),transVtxNoBS.position().y(),transVtxNoBS.position().z());

    //myEvent_->recoEvent_.refitPfPV_.SetXYZ(transVtx.position().x(),transVtx.position().y(),(*vertices)[0].z()); //TEST
    //myEvent_->recoEvent_.refitPfPVNoBS_.SetXYZ(transVtxNoBS.position().x(),transVtxNoBS.position().y(),(*vertices)[0].z());//TEST
     myEvent_->recoEvent_.isRefit_=true;
  }
  else {
     myEvent_->recoEvent_.refitPfPV_.SetXYZ((*vertices)[0].x(),
					    (*vertices)[0].y(),
					    (*vertices)[0].z());
     myEvent_->recoEvent_.refitPfPVNoBS_ = myEvent_->recoEvent_.refitPfPV_;
     myEvent_->recoEvent_.isRefit_=false;
  }
  
  return true;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool MiniAODVertexAnalyzer::findPrimaryVertices(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  edm::Handle<edm::View<pat::PackedCandidate> >  cands;
  iEvent.getByToken(cands_, cands);
  edm::Handle<edm::ValueMap<float> > scores;
  iEvent.getByToken(scores_, scores);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertices_, vertices);

  if(vertices->empty()) return false;   //at least one vertex
  myEvent_->recoEvent_.thePV_.SetXYZ((*vertices)[0].x(),(*vertices)[0].y(),(*vertices)[0].z());

  myEvent_->nPV_ = vertices->size();

  //Find vertex with highest score with PF (miniAOD like)
  //miniAOD uses PF particles instead of tracks
  size_t iPfVtx=0;
  float score=-1;
  for(size_t iVx=0; iVx<vertices->size(); ++iVx){
    reco::VertexRef vtxPrt(vertices,iVx);   
    float aScore = (*scores)[vtxPrt];
    if(iVx < myEvent_->recoEvent_.pfScore_.size()) myEvent_->recoEvent_.pfScore_[iVx] = aScore;
    if(aScore > score){
      score = (*scores)[vtxPrt];
      iPfVtx=iVx;
    }
  }

  myEvent_->recoEvent_.pfPV_.SetXYZ((*vertices)[iPfVtx].x(),(*vertices)[iPfVtx].y(),(*vertices)[iPfVtx].z());
  myEvent_->recoEvent_.pfPVIndex_=iPfVtx;
  ///////////////////////////////
  ///Find PV using sum(pt^2) discriminator.
  ///Calculate weight using tracks with pt>1, and all tracks
  std::map<size_t,float> pt2Scores, pt2ScoresAllTracks;
  for(size_t i=0; i<vertices->size(); ++i){
    pt2Scores[i]=0.;
    pt2ScoresAllTracks[i]=0.;
  }
  
  for(size_t i=0; i<cands->size(); ++i){

    if(verbose_ &&
       (*cands)[i].charge()!=0 &&
       (*cands)[i].bestTrack() &&
       (*cands)[i].bestTrack()->pt()>2000){
      std::cout<<"Here "
	       <<" pt: "<<(*cands)[i].bestTrack()->pt()
	       <<" null: "<<(*cands)[i].vertexRef().isNull()
	//<<" nonull: "<<(*cands)[i].vertexRef().isNonnull()
	       <<" charge: "<<(*cands)[i].charge()
	       <<" pdgId: "<<(*cands)[i].pdgId()
	//<<" key: "<<(*cands)[i].vertexRef().key()
	//<<" fromPV: "<<(*cands)[i].fromPV()
	       <<std::endl;
    }

    if( ( (*cands)[i].charge()==0 && !(*cands)[i].bestTrack() ) || (*cands)[i].vertexRef().isNull()) continue;
    size_t key = (*cands)[i].vertexRef().key();
    int quality = (*cands)[i].pvAssociationQuality();
    if( quality!=pat::PackedCandidate::UsedInFitTight  &&
	quality!=pat::PackedCandidate::UsedInFitLoose) continue;
    if((*cands)[i].bestTrack()){
      float pT = (*cands)[i].bestTrack()->pt();
      float epT=(*cands)[i].bestTrack()->ptError(); 
      pT=pT>epT ? pT-epT : 0;
      pt2Scores[key]+=pT*pT;
      pt2ScoresAllTracks[key]+=pT*pT;	  
    }
    else { //should happen for tracks pT<0.95
      float pT = (*cands)[i].pt();
      if(pT>1 && std::abs((*cands)[i].pdgId())==13) continue;//most probably STA Muons
      //estimate pT error (basing on avarage numbers from TRK-11-001, JINST 9 (2014), P10009
      float epT=0;
      if(std::abs((*cands)[i].eta())<0.9) epT=0.01*pT;
      else if(std::abs((*cands)[i].eta())<1.4) epT=0.02*pT;
      else epT=0.025*pT;
      if(pT>90) epT*=2; //for high-Pt assume bigger error (just to be safe)
      if(pT>1){
	std::cout<<"Something unexpected: Pt="<<pT
		 <<", eta="<<(*cands)[i].eta()
		 <<", charge="<<(*cands)[i].charge()
		 <<", pdgId="<<(*cands)[i].pdgId()
	  //<<", (*cands)[i].bestTrack(): "<<(*cands)[i].bestTrack()
		 <<std::endl;
      }
      pT=pT>epT ? pT-epT : 0;
      pt2ScoresAllTracks[key]+=pT*pT;
    }
  }
  size_t iOldVtx1=0, iOldVtx2=0;
  float pt2Score=-1, pt2ScoreAllTracks=-1;
  for(std::map<size_t,float>::const_iterator it=pt2Scores.begin();it!=pt2Scores.end(); ++it){
    if(it->second > pt2Score){
      pt2Score = it->second;
      iOldVtx1 = it->first;
    }
  }
  for(std::map<size_t,float>::const_iterator it=pt2ScoresAllTracks.begin(); 
       it!=pt2ScoresAllTracks.end(); ++it){
    if(it->second > pt2ScoreAllTracks){
      pt2ScoreAllTracks = it->second;
      iOldVtx2 = it->first;
    }
    if(it->first<myEvent_->recoEvent_.pt2Score_.size())
      myEvent_->recoEvent_.pt2Score_[it->first]=it->second;
  }

  myEvent_->recoEvent_.pt2PV_.SetXYZ((*vertices)[iOldVtx2].x(),(*vertices)[iOldVtx2].y(),(*vertices)[iOldVtx2].z());
  myEvent_->recoEvent_.pt2PVindex_=iOldVtx2;
  iOldVtx1+=0;//hack to avoid compilation errors -Werror=unused-but-set-variable

  //Find vertex with smallest dz distance wrt genPV
  //and store as pfPV in the gen event
  size_t iGenVtx=0;
  float dzMin=9999;
  for(size_t iVx=0; iVx<vertices->size(); ++iVx){
    float dz = std::abs( (*vertices)[iVx].z() - myEvent_->genEvent_.thePV_.Z() );
    if( dz < dzMin){
      dzMin=dz;
      iGenVtx=iVx;
    }
  }
  myEvent_->genEvent_.pfPV_.SetXYZ((*vertices)[iGenVtx].x(),(*vertices)[iGenVtx].y(),(*vertices)[iGenVtx].z());
  myEvent_->genEvent_.pfPVIndex_=iGenVtx;

  //vtx dz
  if(vertices->size()>1) myEvent_->recoEvent_.dzVtx_[0] = (*vertices)[1].z() - (*vertices)[0].z();
  if(vertices->size()>2) myEvent_->recoEvent_.dzVtx_[1] = (*vertices)[2].z() - (*vertices)[0].z();

  return true;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool MiniAODVertexAnalyzer::findRecoTaus(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  edm::Handle<std::vector<pat::Tau> > tauColl;
  iEvent.getByToken(taus_,tauColl);
  edm::Handle<std::vector<pat::Muon> > muColl;
  iEvent.getByToken(mus_,muColl);
  edm::Handle<std::vector<pat::Electron> > eColl;
  iEvent.getByToken(eles_,eColl);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertices_, vertices);

  if( ((muColl->size()<1&&eColl->size()<1) || tauColl->size()<1) && tauColl->size()<2) return false;
  //look for a mu+tau pair firts
  thePair_ = findMTPair(muColl, tauColl, vertices);
  //then for a e+tau pair
  if(thePair_.first==nullptr || thePair_.second==nullptr)
    thePair_ = findETPair(eColl,tauColl,vertices);
  //and then for a di-tau pair
  if(thePair_.first==nullptr || thePair_.second==nullptr)
    thePair_ = findTTPair(tauColl);

  if(thePair_.first==nullptr || thePair_.second==nullptr) return false;

  myEvent_->diMuonVeto_ = diMuVeto(muColl, vertices);
  myEvent_->diEleVeto_ = diEVeto(eColl, vertices);

  //Tau+
  myEvent_->recoEvent_.isoMVAWpPlus_ = 0;
  myEvent_->recoEvent_.antiEWpPlus_ = 0;
  myEvent_->recoEvent_.antiE2WpPlus_ = 0;
  myEvent_->recoEvent_.antiMuWpPlus_ = 0;
  const pat::Tau* aTau = dynamic_cast<const pat::Tau*>(thePair_.first);
  if(aTau!=nullptr){
    if(aTau->tauID("byVVLooseIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5)
      myEvent_->recoEvent_.isoMVAWpPlus_ = 1;
    if(aTau->tauID("byVLooseIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5)
      myEvent_->recoEvent_.isoMVAWpPlus_ = 2;
    if(aTau->tauID("byLooseIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5)
      myEvent_->recoEvent_.isoMVAWpPlus_ = 3;
    if(aTau->tauID("byMediumIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5)
      myEvent_->recoEvent_.isoMVAWpPlus_ = 4;
    if(aTau->tauID("byTightIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5)
      myEvent_->recoEvent_.isoMVAWpPlus_ = 5;
    if(aTau->tauID("byVTightIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5)
      myEvent_->recoEvent_.isoMVAWpPlus_ = 6;
    if(aTau->tauID("byVVTightIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5)
      myEvent_->recoEvent_.isoMVAWpPlus_ = 7;
    myEvent_->recoEvent_.isoPlus_ = aTau->tauID("byIsolationMVArun2v1DBnewDMwLTraw2017v2");
    if(aTau->tauID("againstElectronVLooseMVA6") > 0.5)
      myEvent_->recoEvent_.antiEWpPlus_ = 1;
    if(aTau->tauID("againstElectronLooseMVA6") > 0.5)
      myEvent_->recoEvent_.antiEWpPlus_ = 2;
    if(aTau->tauID("againstElectronMediumMVA6") > 0.5)
      myEvent_->recoEvent_.antiEWpPlus_ = 3;
    if(aTau->tauID("againstElectronTightMVA6") > 0.5)
      myEvent_->recoEvent_.antiEWpPlus_ = 4;
    if(aTau->tauID("againstElectronVTightMVA6") > 0.5)
      myEvent_->recoEvent_.antiEWpPlus_ = 5;
    if(aTau->tauID("againstElectronVLooseMVA62018") > 0.5)
      myEvent_->recoEvent_.antiE2WpPlus_ = 1;
    if(aTau->tauID("againstElectronLooseMVA62018") > 0.5)
      myEvent_->recoEvent_.antiE2WpPlus_ = 2;
    if(aTau->tauID("againstElectronMediumMVA62018") > 0.5)
      myEvent_->recoEvent_.antiE2WpPlus_ = 3;
    if(aTau->tauID("againstElectronTightMVA62018") > 0.5)
      myEvent_->recoEvent_.antiE2WpPlus_ = 4;
    if(aTau->tauID("againstElectronVTightMVA62018") > 0.5)
      myEvent_->recoEvent_.antiE2WpPlus_ = 5;
    if(aTau->tauID("againstMuonLoose3") > 0.5)
      myEvent_->recoEvent_.antiMuWpPlus_ = 1;
    if(aTau->tauID("againstMuonTight3") > 0.5)
      myEvent_->recoEvent_.antiMuWpPlus_ = 2;
    myEvent_->recoEvent_.decModePlus_ = aTau->decayMode();    
    //pi0's
    reco::Candidate::LorentzVector pi0P4;
    myEvent_->recoEvent_.nGammaPlus_ = 0;
    myEvent_->recoEvent_.nGammaInConePlus_ = 0;
    myEvent_->recoEvent_.gammaPtSumInScPlus_ = 0;
    myEvent_->recoEvent_.gammaPtSumOutScPlus_ = 0;
    myEvent_->recoEvent_.dR2midPlus_ = (aTau->signalGammaCands().size()>0 ? 0 : 99);
    float sumPt2 = (aTau->signalGammaCands().size()>0 ? 0 : 1);
    float signalConeR2 = std::clamp(3.0/aTau->pt(),0.05,0.1);
    signalConeR2 *= signalConeR2;
    for(size_t j=0; j<aTau->signalGammaCands().size(); ++j){
      float gammaPt = aTau->signalGammaCands()[j]->pt();
      if(!(gammaPt>0.5)) continue;
      myEvent_->recoEvent_.nGammaPlus_++;
      sumPt2 += gammaPt*gammaPt;
      float dR2 = deltaR2(aTau->p4(),aTau->signalGammaCands()[j]->p4());
      myEvent_->recoEvent_.dR2midPlus_ += dR2*gammaPt*gammaPt;
      pi0P4 += aTau->signalGammaCands()[j]->p4();
      if(dR2 < signalConeR2){
	myEvent_->recoEvent_.nGammaInConePlus_++;
	myEvent_->recoEvent_.gammaPtSumInScPlus_ += gammaPt;
      }
      else {
	myEvent_->recoEvent_.gammaPtSumOutScPlus_ += gammaPt;
      }
    }
    myEvent_->recoEvent_.pi0Plus_.SetXYZT(pi0P4.px(),
					  pi0P4.py(),
					  pi0P4.pz(),
					  pi0P4.e());
    myEvent_->recoEvent_.dR2midPlus_ /= sumPt2;
    const pat::PackedCandidate *leadCharged = dynamic_cast<pat::PackedCandidate const*>(aTau->leadChargedHadrCand().get());
    if(leadCharged!=nullptr){
      myEvent_->recoEvent_.leadIdPlus_ = leadCharged->pdgId();
      myEvent_->recoEvent_.dzPlus_ = leadCharged->dz();
      myEvent_->recoEvent_.piPlus_.SetXYZT(leadCharged->p4().px(),
					   leadCharged->p4().py(),
					   leadCharged->p4().pz(),
					   leadCharged->p4().e());
    }
    //recover mass of 1prong+0pi0 in case of strip in cone
    if(aTau->decayMode()==WawGenInfoHelper::tauDecayModes::tauDecay1ChargedPion0PiZero
       && leadCharged!=nullptr
       && pi0P4.pt()>0.5 && deltaR2(aTau->p4(),pi0P4) < signalConeR2
       ){
      float mass = (leadCharged->p4()+pi0P4).mass();
      myEvent_->recoEvent_.visTauPlus_.SetXYZM(aTau->p4().px(),
					       aTau->p4().py(),
					       aTau->p4().pz(),
					       mass);
    }
    else{
      myEvent_->recoEvent_.visTauPlus_.SetXYZT(aTau->p4().px(),
					       aTau->p4().py(),
					       aTau->p4().pz(),
					       aTau->p4().e());
    }
  }
  else{
    const pat::Muon* aMu = dynamic_cast<const pat::Muon*>(thePair_.first);
    if(aMu!=nullptr){
      myEvent_->recoEvent_.isoPlus_ = muRelIso(aMu);
      myEvent_->recoEvent_.decModePlus_ = WawGenInfoHelper::tauDecayModes::tauDecayMuon;
      myEvent_->recoEvent_.leadIdPlus_ = aMu->pdgId();
      myEvent_->recoEvent_.dzPlus_ = aMu->muonBestTrack()->dz((*vertices)[0].position());
      myEvent_->recoEvent_.piPlus_.SetXYZT(aMu->p4().px(),
					   aMu->p4().py(),
					   aMu->p4().pz(),
					   aMu->p4().e());
    }
    else{
      const pat::Electron* aEle = dynamic_cast<const pat::Electron*>(thePair_.first);
      if(aEle!=nullptr){
	myEvent_->recoEvent_.isoPlus_ = eRelIso(aEle);
	myEvent_->recoEvent_.decModePlus_ = WawGenInfoHelper::tauDecayModes::tauDecaysElectron;
	myEvent_->recoEvent_.leadIdPlus_ = aEle->pdgId();
	myEvent_->recoEvent_.dzPlus_ = aEle->gsfTrack()->dz((*vertices)[0].position());
	myEvent_->recoEvent_.piPlus_.SetXYZT(aEle->p4().px(),
					     aEle->p4().py(),
					     aEle->p4().pz(),
					     aEle->p4().e());
      }
    }
    myEvent_->recoEvent_.visTauPlus_.SetXYZT(thePair_.first->p4().px(),
					     thePair_.first->p4().py(),
					     thePair_.first->p4().pz(),
					     thePair_.first->p4().e());
  }
  myEvent_->recoEvent_.tauPlus_ = myEvent_->recoEvent_.visTauPlus_;
  //gen matching
  // correct matching
  if(deltaR2(myEvent_->recoEvent_.visTauPlus_.Eta(),myEvent_->recoEvent_.visTauPlus_.Phi(),
	     myEvent_->genEvent_.visTauPlus_.Eta(),myEvent_->genEvent_.visTauPlus_.Phi()) < 0.2*0.2){
    if(myEvent_->genEvent_.decModePlus_==WawGenInfoHelper::tauDecayModes::tauDecayMuon)
      myEvent_->recoEvent_.matchedPlus_ = 4;
    else if(myEvent_->genEvent_.decModePlus_==WawGenInfoHelper::tauDecayModes::tauDecaysElectron)
      myEvent_->recoEvent_.matchedPlus_ = 3;
    else if(myEvent_->genEvent_.decModePlus_!=WawGenInfoHelper::tauDecayModes::tauDecayOther || myEvent_->genEvent_.visTauPlus_.Pt()>0)
      myEvent_->recoEvent_.matchedPlus_ = 5;
  }
  // swapped matching
  else if(deltaR2(myEvent_->recoEvent_.visTauPlus_.Eta(),myEvent_->recoEvent_.visTauPlus_.Phi(),
		  myEvent_->genEvent_.visTauMinus_.Eta(),myEvent_->genEvent_.visTauMinus_.Phi()) < 0.2*0.2){
    if(myEvent_->genEvent_.decModeMinus_==WawGenInfoHelper::tauDecayModes::tauDecayMuon)
      myEvent_->recoEvent_.matchedPlus_ = -4;
    else if(myEvent_->genEvent_.decModeMinus_==WawGenInfoHelper::tauDecayModes::tauDecaysElectron)
      myEvent_->recoEvent_.matchedPlus_ = -3;
    else if(myEvent_->genEvent_.decModeMinus_!=WawGenInfoHelper::tauDecayModes::tauDecayOther || myEvent_->genEvent_.visTauMinus_.Pt()>0)
      myEvent_->recoEvent_.matchedPlus_ = -5;
  }
  // matching to prompt e/mu?

  //Tau-
  myEvent_->recoEvent_.isoMVAWpMinus_ = 0;
  myEvent_->recoEvent_.antiEWpMinus_ = 0;
  myEvent_->recoEvent_.antiE2WpMinus_ = 0;
  myEvent_->recoEvent_.antiMuWpMinus_ = 0;
  aTau = dynamic_cast<const pat::Tau*>(thePair_.second);
  if(aTau!=nullptr){
    if(aTau->tauID("byVVLooseIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5)
      myEvent_->recoEvent_.isoMVAWpMinus_ = 1;
    if(aTau->tauID("byVLooseIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5)
      myEvent_->recoEvent_.isoMVAWpMinus_ = 2;
    if(aTau->tauID("byLooseIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5)
      myEvent_->recoEvent_.isoMVAWpMinus_ = 3;
    if(aTau->tauID("byMediumIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5)
      myEvent_->recoEvent_.isoMVAWpMinus_ = 4;
    if(aTau->tauID("byTightIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5)
      myEvent_->recoEvent_.isoMVAWpMinus_ = 5;
    if(aTau->tauID("byVTightIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5)
      myEvent_->recoEvent_.isoMVAWpMinus_ = 6;
    if(aTau->tauID("byVVTightIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5)
      myEvent_->recoEvent_.isoMVAWpMinus_ = 7;
    myEvent_->recoEvent_.isoMinus_ = aTau->tauID("byIsolationMVArun2v1DBnewDMwLTraw2017v2");
    if(aTau->tauID("againstElectronVLooseMVA6") > 0.5)
      myEvent_->recoEvent_.antiEWpMinus_ = 1;
    if(aTau->tauID("againstElectronLooseMVA6") > 0.5)
      myEvent_->recoEvent_.antiEWpMinus_ = 2;
    if(aTau->tauID("againstElectronMediumMVA6") > 0.5)
      myEvent_->recoEvent_.antiEWpMinus_ = 3;
    if(aTau->tauID("againstElectronTightMVA6") > 0.5)
      myEvent_->recoEvent_.antiEWpMinus_ = 4;
    if(aTau->tauID("againstElectronVTightMVA6") > 0.5)
      myEvent_->recoEvent_.antiEWpMinus_ = 5;
    if(aTau->tauID("againstElectronVLooseMVA62018") > 0.5)
      myEvent_->recoEvent_.antiE2WpMinus_ = 1;
    if(aTau->tauID("againstElectronLooseMVA62018") > 0.5)
      myEvent_->recoEvent_.antiE2WpMinus_ = 2;
    if(aTau->tauID("againstElectronMediumMVA62018") > 0.5)
      myEvent_->recoEvent_.antiE2WpMinus_ = 3;
    if(aTau->tauID("againstElectronTightMVA62018") > 0.5)
      myEvent_->recoEvent_.antiE2WpMinus_ = 4;
    if(aTau->tauID("againstElectronVTightMVA62018") > 0.5)
      myEvent_->recoEvent_.antiE2WpMinus_ = 5;
    if(aTau->tauID("againstMuonLoose3") > 0.5)
      myEvent_->recoEvent_.antiMuWpMinus_ = 1;
    if(aTau->tauID("againstMuonTight3") > 0.5)
      myEvent_->recoEvent_.antiMuWpMinus_ = 2;
    myEvent_->recoEvent_.decModeMinus_ = aTau->decayMode();    
    //pi0's
    reco::Candidate::LorentzVector pi0P4;
    myEvent_->recoEvent_.nGammaMinus_ = 0;
    myEvent_->recoEvent_.nGammaInConeMinus_ = 0;
    myEvent_->recoEvent_.gammaPtSumInScMinus_ = 0;
    myEvent_->recoEvent_.gammaPtSumOutScMinus_ = 0;
    myEvent_->recoEvent_.dR2midMinus_ = (aTau->signalGammaCands().size()>0 ? 0 : 99);
    float sumPt2 = (aTau->signalGammaCands().size()>0 ? 0 : 1);
    float signalConeR2 = std::clamp(3.0/aTau->pt(),0.05,0.1);
    signalConeR2 *= signalConeR2;
    for(size_t j=0; j<aTau->signalGammaCands().size(); ++j){
      float gammaPt = aTau->signalGammaCands()[j]->pt();
      if(!(gammaPt>0.5)) continue;
      myEvent_->recoEvent_.nGammaMinus_++;
      sumPt2 += gammaPt*gammaPt;
      float dR2 = deltaR2(aTau->p4(),aTau->signalGammaCands()[j]->p4());
      myEvent_->recoEvent_.dR2midMinus_ += dR2*gammaPt*gammaPt;
      pi0P4 += aTau->signalGammaCands()[j]->p4();
      if(dR2 < signalConeR2){
	myEvent_->recoEvent_.nGammaInConeMinus_++;
	myEvent_->recoEvent_.gammaPtSumInScMinus_ += gammaPt;
      }
      else {
	myEvent_->recoEvent_.gammaPtSumOutScMinus_ += gammaPt;
      }
    }
    myEvent_->recoEvent_.pi0Minus_.SetXYZT(pi0P4.px(),
					   pi0P4.py(),
					   pi0P4.pz(),
					   pi0P4.e());
    myEvent_->recoEvent_.dR2midMinus_ /= sumPt2;
    const pat::PackedCandidate *leadCharged = dynamic_cast<pat::PackedCandidate const*>(aTau->leadChargedHadrCand().get());
    if(leadCharged!=nullptr){
      myEvent_->recoEvent_.leadIdMinus_ = leadCharged->pdgId();
      myEvent_->recoEvent_.dzMinus_ = leadCharged->dz();
      myEvent_->recoEvent_.piMinus_.SetXYZT(leadCharged->p4().px(),
					    leadCharged->p4().py(),
					    leadCharged->p4().pz(),
					    leadCharged->p4().e());
    }
    //recover mass of 1prong+0pi0 in case of strip in cone
    if(aTau->decayMode()==WawGenInfoHelper::tauDecayModes::tauDecay1ChargedPion0PiZero
       && leadCharged!=nullptr
       && pi0P4.pt()>0.5 && deltaR2(aTau->p4(),pi0P4) < signalConeR2
       ){
      float mass = (leadCharged->p4()+pi0P4).mass();
      myEvent_->recoEvent_.visTauMinus_.SetXYZM(aTau->p4().px(),
						aTau->p4().py(),
						aTau->p4().pz(),
						mass);
    }
    else{
      myEvent_->recoEvent_.visTauMinus_.SetXYZT(aTau->p4().px(),
						aTau->p4().py(),
						aTau->p4().pz(),
						aTau->p4().e());
    }
  }
  else{
    const pat::Muon* aMu = dynamic_cast<const pat::Muon*>(thePair_.second);
    if(aMu!=nullptr){
      myEvent_->recoEvent_.isoMinus_ = muRelIso(aMu);
      myEvent_->recoEvent_.decModeMinus_ = WawGenInfoHelper::tauDecayModes::tauDecayMuon;
      myEvent_->recoEvent_.leadIdMinus_ = aMu->pdgId();
      myEvent_->recoEvent_.dzMinus_ = aMu->muonBestTrack()->dz((*vertices)[0].position());
      myEvent_->recoEvent_.piMinus_.SetXYZT(aMu->p4().px(),
					    aMu->p4().py(),
					    aMu->p4().pz(),
					    aMu->p4().e());
    }
    else{
      const pat::Electron* aEle = dynamic_cast<const pat::Electron*>(thePair_.second);
      if(aEle!=nullptr){
	myEvent_->recoEvent_.isoMinus_ = eRelIso(aEle);
	myEvent_->recoEvent_.decModeMinus_ = WawGenInfoHelper::tauDecayModes::tauDecaysElectron;
	myEvent_->recoEvent_.leadIdMinus_ = aEle->pdgId();
	myEvent_->recoEvent_.dzMinus_ = aEle->gsfTrack()->dz((*vertices)[0].position());
	myEvent_->recoEvent_.piMinus_.SetXYZT(aEle->p4().px(),
					      aEle->p4().py(),
					      aEle->p4().pz(),
					      aEle->p4().e());
      }
    }
    myEvent_->recoEvent_.visTauMinus_.SetXYZT(thePair_.second->p4().px(),
					      thePair_.second->p4().py(),
					      thePair_.second->p4().pz(),
					      thePair_.second->p4().e());
  }
  myEvent_->recoEvent_.tauMinus_ = myEvent_->recoEvent_.visTauMinus_;
  //gen matching
  // correct matching
  if(deltaR2(myEvent_->recoEvent_.visTauMinus_.Eta(),myEvent_->recoEvent_.visTauMinus_.Phi(),
	     myEvent_->genEvent_.visTauMinus_.Eta(),myEvent_->genEvent_.visTauMinus_.Phi()) < 0.2*0.2){
    if(myEvent_->genEvent_.decModeMinus_==WawGenInfoHelper::tauDecayModes::tauDecayMuon)
      myEvent_->recoEvent_.matchedMinus_ = 4;
    else if(myEvent_->genEvent_.decModeMinus_==WawGenInfoHelper::tauDecayModes::tauDecaysElectron)
      myEvent_->recoEvent_.matchedMinus_ = 3;
    else if(myEvent_->genEvent_.decModeMinus_!=WawGenInfoHelper::tauDecayModes::tauDecayOther || myEvent_->genEvent_.visTauMinus_.Pt()>0)
      myEvent_->recoEvent_.matchedMinus_ = 5;
  }
  // swapped matching
  else if(deltaR2(myEvent_->recoEvent_.visTauMinus_.Eta(),myEvent_->recoEvent_.visTauMinus_.Phi(),
		  myEvent_->genEvent_.visTauPlus_.Eta(),myEvent_->genEvent_.visTauPlus_.Phi()) < 0.2*0.2){
    if(myEvent_->genEvent_.decModePlus_==WawGenInfoHelper::tauDecayModes::tauDecayMuon)
      myEvent_->recoEvent_.matchedMinus_ = -4;
    else if(myEvent_->genEvent_.decModePlus_==WawGenInfoHelper::tauDecayModes::tauDecaysElectron)
      myEvent_->recoEvent_.matchedMinus_ = -3;
    else if(myEvent_->genEvent_.decModePlus_!=WawGenInfoHelper::tauDecayModes::tauDecayOther || myEvent_->genEvent_.visTauPlus_.Pt()>0)
      myEvent_->recoEvent_.matchedMinus_ = -5;
  }
  // matching to prompt e/mu?

  //the END
  return true;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool MiniAODVertexAnalyzer::setPCAVectors(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  if(thePair_.first==nullptr || thePair_.second==nullptr) return false;

  const reco::Track *leadTrackPlus = getLeadTrack(iEvent, thePair_.first);
  const reco::Track *leadTrackMinus = getLeadTrack(iEvent, thePair_.second);

  GlobalPoint aPoint(myEvent_->recoEvent_.refitPfPV_.X(),
		     myEvent_->recoEvent_.refitPfPV_.Y(),
		     myEvent_->recoEvent_.thePV_.Z());

  myEvent_->recoEvent_.nPiPlus_ = getPCA(iEvent,iSetup, leadTrackPlus, aPoint);
  myEvent_->recoEvent_.nPiMinus_ = getPCA(iEvent,iSetup, leadTrackMinus, aPoint);

   ///Test different PVs.
   ///AOD PV
  aPoint = GlobalPoint(myEvent_->recoEvent_.thePV_.X(),
		       myEvent_->recoEvent_.thePV_.Y(),
		       myEvent_->recoEvent_.thePV_.Z());
  myEvent_->recoEvent_.nPiPlusAODvx_ = getPCA(iEvent,iSetup, leadTrackPlus, aPoint);
  myEvent_->recoEvent_.nPiMinusAODvx_ = getPCA(iEvent,iSetup, leadTrackMinus, aPoint);
  
  ///PV refit (excluding tau tracks?)
  aPoint = GlobalPoint(myEvent_->recoEvent_.refitPfPV_.X(),
		       myEvent_->recoEvent_.refitPfPV_.Y(),
		       myEvent_->recoEvent_.refitPfPV_.Z());
  myEvent_->recoEvent_.nPiPlusRefitvx_ = getPCA(iEvent,iSetup, leadTrackPlus, aPoint);
  myEvent_->recoEvent_.nPiMinusRefitvx_ = getPCA(iEvent,iSetup, leadTrackMinus, aPoint);
  
  ///Generated PV
  aPoint = GlobalPoint(myEvent_->genEvent_.thePV_.X(),
		       myEvent_->genEvent_.thePV_.Y(),
		       myEvent_->genEvent_.thePV_.Z());
  myEvent_->recoEvent_.nPiPlusGenvx_ = getPCA(iEvent,iSetup, leadTrackPlus, aPoint);
  myEvent_->recoEvent_.nPiMinusGenvx_ = getPCA(iEvent,iSetup, leadTrackMinus, aPoint);
  
  return true;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
const reco::Track* MiniAODVertexAnalyzer::getLeadTrack(const edm::Event & iEvent,
						       const reco::Candidate* aTau) {

  const reco::Track *aTrack = nullptr;
  const pat::Tau* aPatTau = dynamic_cast<const pat::Tau*>(aTau);
  if(aPatTau!=nullptr){
    const pat::PackedCandidate *leadCharged = dynamic_cast<pat::PackedCandidate const*>(aPatTau->leadChargedHadrCand().get());
    if(leadCharged!=nullptr && leadCharged->hasTrackDetails()){
      aTrack = leadCharged->bestTrack();
      if(std::abs(leadCharged->pdgId())==11){////look for associated KF track for electrons
	float maxDr2=0.01*0.01, maxDeta=0.01;//MC-based guesstimate
	edm::Handle<edm::View<pat::PackedCandidate> >  eleTracks;
	iEvent.getByToken(lostCandsEle_, eleTracks);
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(vertices_, vertices);
	const reco::Track *aKFTrack = nullptr;
	std::cout<<"leadEle"<<std::endl;
	for(size_t i=0; i<eleTracks->size(); ++i){
	  if(std::abs((*eleTracks)[i].pdgId())!=11 ||//protection againist possible "photon" tracks
	     !(*eleTracks)[i].hasTrackDetails() ||
	     !((*eleTracks)[i].pt()>0.5) ||
	     !(std::abs((*eleTracks)[i].dxy((*vertices)[0].position()))<0.05) || //IP cuts for safety 10% looser than for selection (0.045->0.05, 0.2->0.22) (OK?)
	     !(std::abs((*eleTracks)[i].dz((*vertices)[0].position()))<0.22) ) continue;
	  float dR2=deltaR2(aTrack->eta(),aTrack->phi(),
			    (*eleTracks)[i].eta(),(*eleTracks)[i].phi());
	  float dEta=std::abs(aTrack->eta()-(*eleTracks)[i].eta());
	  if(dR2<maxDr2 && dEta<maxDeta){
	    aKFTrack = (*eleTracks)[i].bestTrack();
	    maxDeta = dEta;
	    std::cout<<"\t\tKF track, dR="<<sqrt(dR2)
		     <<", dEta="<<dEta<<std::endl;
	  }
	}
	if(aKFTrack!=nullptr){
	  std::cout<<"\tKF track, pt="<<aKFTrack->pt()
		   <<", eta="<<aKFTrack->eta()
		   <<", phi="<<aKFTrack->phi()
		   <<std::endl
		   <<"\tGSF track, pt="<<aTrack->pt()
		   <<", eta="<<aTrack->eta()
		   <<", phi="<<aTrack->phi()
		   <<std::endl;
	  aTrack = aKFTrack;
	}
      }
    }
  }
  else{
    const pat::Muon* aPatMu = dynamic_cast<const pat::Muon*>(aTau);
    if(aPatMu!=nullptr){
      aTrack = aPatMu->muonBestTrack().get();
    }
    else{
      const pat::Electron* aPatEle = dynamic_cast<const pat::Electron*>(aTau);
      if(aPatEle!=nullptr){
	aTrack = aPatEle->gsfTrack().get();
      }
    }
  }

  return aTrack;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
TVector3 MiniAODVertexAnalyzer::getPCA(const edm::Event & iEvent, const edm::EventSetup & iSetup,
				       const reco::Track* aTrack,
				       const GlobalPoint & aPoint){

  TVector3 aPCA;
  if(aTrack==nullptr) return aPCA;
  
  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);  
  reco::TransientTrack transTrk=transTrackBuilder->build(aTrack);
     
  //TransverseImpactPointExtrapolator extrapolator(transTrk.field());
  AnalyticalImpactPointExtrapolator extrapolator(transTrk.field());
  GlobalPoint pos  = extrapolator.extrapolate(transTrk.impactPointState(),aPoint).globalPosition();

  aPCA.SetX(pos.x() - aPoint.x());
  aPCA.SetY(pos.y() - aPoint.y());
  aPCA.SetZ(pos.z() - aPoint.z());

  /*
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertices_, vertices);
  GlobalVector direction(aTau->p4().px(), aTau->p4().py(), aTau->p4().pz()); //To compute sign of IP   
  std::pair<bool,Measurement1D> signed_IP2D = IPTools::signedTransverseImpactParameter(transTrk, direction,(*vertices)[0]);

  GlobalVector mom  = extrapolator.extrapolate(transTrk.impactPointState(),aPoint).globalMomentum();
  LocalVector momLoc  = extrapolator.extrapolate(transTrk.impactPointState(),aPoint).localMomentum();
  TVector3 extrapolatedMom(mom.x(), mom.y(), mom.z());  
  GlobalPoint point = extrapolator.extrapolate(transTrk.impactPointState(),aPoint).surface().toGlobal(Local3DPoint(0,0,1));
  TVector3 aTestPoint(point.x() - aPoint.x(),
		      point.y() - aPoint.y(),
		      point.z() - aPoint.z());
  
 
  
   std::cout<<"IP2D: "<<signed_IP2D.second.value()
	   <<" PCA local: "<<extrapolator.extrapolate(transTrk.impactPointState(),aPoint).localPosition()
	    <<" mag local: "<<extrapolator.extrapolate(transTrk.impactPointState(),aPoint).localPosition().mag()
	    <<" PCA global: ("<<aPCA.X()<<" "<<aPCA.Y()<<" "<<aPCA.Z()<<")"
	    <<" mag global: "<<aPCA.Mag()
	    <<std::endl
	    <<"vertex: "<<aPoint
     	    <<" test point: ("<<aTestPoint.X()<<" "<<aTestPoint.Y()<<" "<<aTestPoint.Z()<<")"
	    <<" dot (PCA, extrapol. mom.): "<<aPCA.Unit().Dot(extrapolatedMom.Unit())
	    <<" local mom: "<<momLoc
	    <<" global mom: "<<mom
     //<<" dot (mom, extrapol. mom.): "<<extrapolatedMom.Unit().Dot( myEvent_->recoEvent_.piPlus_.Vect().Unit())
	    <<std::endl;
  */
  return aPCA;
}
////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float MiniAODVertexAnalyzer::muRelIso(const pat::Muon *aMu){
  if(aMu==nullptr) return 999;
  return (aMu->pfIsolationR04().sumChargedHadronPt + 
	  std::max(aMu->pfIsolationR04().sumNeutralHadronEt +
		   aMu->pfIsolationR04().sumPhotonEt - 
		   0.5 * aMu->pfIsolationR04().sumPUPt, 
		   0.0)
	  ) / aMu->pt();
}
////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float MiniAODVertexAnalyzer::eRelIso(const pat::Electron *aEle){
  if(aEle==nullptr) return 999;
  return (aEle->pfIsolationVariables().sumChargedHadronPt +
	  std::max(aEle->pfIsolationVariables().sumNeutralHadronEt +
		   aEle->pfIsolationVariables().sumPhotonEt -
		   0.5 * aEle->pfIsolationVariables().sumPUPt,
		   0.0)
	  ) / aEle->pt();
}
  
////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
std::pair<const reco::Candidate*, const reco::Candidate*>
MiniAODVertexAnalyzer::findMTPair(const edm::Handle<std::vector<pat::Muon> > &muColl,
				  const edm::Handle<std::vector<pat::Tau> > &tauColl,
				  const edm::Handle<reco::VertexCollection> &vtxColl){
  //This section shuld be connected to mu/tau candidate selection from the
  ///analysis code.
  std::pair<const reco::Candidate*, const reco::Candidate*> myPair(0, 0);
  if(muColl->size()<1 || tauColl->size()<1 || vtxColl->empty()) return myPair;

  //select muons
  std::vector<const pat::Muon*> mus;
  for(size_t i=0; i<muColl->size(); ++i){
    if( (*muColl)[i].pt()>20 && std::abs((*muColl)[i].eta())<2.1
	&& std::abs((*muColl)[i].muonBestTrack()->dxy((*vtxColl)[0].position()))<0.045
	&& std::abs((*muColl)[i].muonBestTrack()->dz((*vtxColl)[0].position())) <0.2
	&& (*muColl)[i].isMediumMuon()
	&& muRelIso(&(*muColl)[i])<0.5
	)
      mus.push_back(&(*muColl)[i]);
  }
  if(mus.size()<1) return myPair;
  //sort muons
  std::sort(mus.begin(), mus.end(),
	    [this](const pat::Muon* a, const pat::Muon* b){
	      //return a->pt() > b->pt();
	      float iso_a = muRelIso(a);
	      float iso_b = muRelIso(b);
	      return iso_a < iso_b;
	    });

  //select taus
  std::vector<const pat::Tau*> taus;
  for(size_t i=0; i<tauColl->size(); ++i){
    if( (*tauColl)[i].tauID("decayModeFindingNewDMs")>0.5 &&
	( (*tauColl)[i].decayMode()<WawGenInfoHelper::tauDecayModes::tauDecay2ChargedPion0PiZero ||
	  (*tauColl)[i].decayMode()>WawGenInfoHelper::tauDecayModes::tauDecay2ChargedPion4PiZero )
	&& (*tauColl)[i].pt()>20 && std::abs((*tauColl)[i].eta())<2.3
	&& std::abs(dynamic_cast<pat::PackedCandidate const*>((*tauColl)[i].leadChargedHadrCand().get())->dz()) < 0.2
	&& std::abs((*tauColl)[i].charge()) == 1
	&& (*tauColl)[i].tauID("againstMuonLoose3") > 0.5
	//&& (*tauColl)[i].tauID("photonPtSumOutsideSignalCone") < 0.1 * (*tauColl)[i].pt()
	//&& (*tauColl)[i].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 3.0 //Loose <2.5, Medium<1.5, Tight<0.8
	&& (*tauColl)[i].tauID("byVVLooseIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5
	)
      taus.push_back(&(*tauColl)[i]);
  }
  if(taus.size()<1) return myPair;
  //sort taus
  std::sort(taus.begin(), taus.end(),
	    [](const pat::Tau* a, const pat::Tau* b){
	      //return a->pt() > b->pt();
	      //return a->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") <
	      //       b->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
	      return a->tauID("byIsolationMVArun2v1DBnewDMwLTraw2017v2") >
		     b->tauID("byIsolationMVArun2v1DBnewDMwLTraw2017v2");
	    });

  for(size_t i=0; i<taus.size(); ++i){
    if(deltaR2(mus[0]->p4(),taus[i]->p4())>0.5*0.5 &&
       mus[0]->charge()*taus[i]->charge()==-1){
      if(mus[0]->charge()==1)
	myPair = std::make_pair(mus[0], taus[i]);
      else
	myPair = std::make_pair(taus[i], mus[0]);
      break;
    }
  }

  return myPair;
}
////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
std::pair<const reco::Candidate*, const reco::Candidate*>
MiniAODVertexAnalyzer::findETPair(const edm::Handle<std::vector<pat::Electron> > &eColl,
				  const edm::Handle<std::vector<pat::Tau> > &tauColl,
				  const edm::Handle<reco::VertexCollection> &vtxColl){
  //This section shuld be connected to e/tau candidate selection from the
  ///analysis code.
  std::pair<const reco::Candidate*, const reco::Candidate*> myPair(0, 0);
  if(eColl->size()<1 || tauColl->size()<1 || vtxColl->empty()) return myPair;

  //select electrons
  std::vector<const pat::Electron*> eles;
  for(size_t i=0; i<eColl->size(); ++i){
    if( (*eColl)[i].pt()>20 && std::abs((*eColl)[i].eta())<2.1
	&& std::abs((*eColl)[i].gsfTrack()->dxy((*vtxColl)[0].position()))<0.045
	&& std::abs((*eColl)[i].gsfTrack()->dz((*vtxColl)[0].position())) <0.2
	&& (*eColl)[i].electronID("mvaEleID-Fall17-noIso-V1-wp80")>0.5 //also wp90 and wpLoose
	&& (*eColl)[i].gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) <=1
	&& (*eColl)[i].passConversionVeto()
	&& eRelIso(&(*eColl)[i])<0.5
	)
      eles.push_back(&(*eColl)[i]);
  }
  if(eles.size()<1) return myPair;
  //sort electrons
  std::sort(eles.begin(), eles.end(),
	    [this](const pat::Electron* a, const pat::Electron* b){
	      //return a->pt() > b->pt();
	      float iso_a = eRelIso(a);
	      float iso_b = eRelIso(b);
	      return iso_a < iso_b;
	    });

  //select taus
  std::vector<const pat::Tau*> taus;
  for(size_t i=0; i<tauColl->size(); ++i){
    if( (*tauColl)[i].tauID("decayModeFindingNewDMs")>0.5 &&
	( (*tauColl)[i].decayMode()<WawGenInfoHelper::tauDecayModes::tauDecay2ChargedPion0PiZero ||
	  (*tauColl)[i].decayMode()>WawGenInfoHelper::tauDecayModes::tauDecay2ChargedPion4PiZero )
	&& (*tauColl)[i].pt()>20 && std::abs((*tauColl)[i].eta())<2.3
	&& std::abs(dynamic_cast<pat::PackedCandidate const*>((*tauColl)[i].leadChargedHadrCand().get())->dz()) < 0.2
	&& std::abs((*tauColl)[i].charge()) == 1
	&& (*tauColl)[i].tauID("againstMuonLoose3") > 0.5
	//&& (*tauColl)[i].tauID("photonPtSumOutsideSignalCone") < 0.1 * (*tauColl)[i].pt()
	//&& (*tauColl)[i].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 3.0 //Loose <2.5, Medium<1.5, Tight<0.8
	&& (*tauColl)[i].tauID("byVVLooseIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5
	)
      taus.push_back(&(*tauColl)[i]);
  }
  if(taus.size()<1) return myPair;
  //sort taus
  std::sort(taus.begin(), taus.end(),
	    [](const pat::Tau* a, const pat::Tau* b){
	      //return a->pt() > b->pt();
	      //return a->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") <
	      //       b->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
	      return a->tauID("byIsolationMVArun2v1DBnewDMwLTraw2017v2") >
		     b->tauID("byIsolationMVArun2v1DBnewDMwLTraw2017v2");
	    });

  for(size_t i=0; i<taus.size(); ++i){
    if(deltaR2(eles[0]->p4(),taus[i]->p4())>0.5*0.5 &&
       eles[0]->charge()*taus[i]->charge()==-1){
      if(eles[0]->charge()==1)
	myPair = std::make_pair(eles[0], taus[i]);
      else
	myPair = std::make_pair(taus[i], eles[0]);
      break;
    }
  }

  return myPair;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
std::pair<const pat::Tau*, const pat::Tau*>
MiniAODVertexAnalyzer::findTTPair(const edm::Handle<std::vector<pat::Tau> > &tauColl){
  //This section shuld be connected to tau candidate selection from the
  ///analysis code.
  std::pair<const pat::Tau*, const pat::Tau*> myPair(0, 0);
  if(tauColl->size()<2) return myPair;
  /*
  if(tauColl->size()>=2 && (*tauColl)[0].charge()*(*tauColl)[1].charge()==-1){
    if((*tauColl)[0].charge()==1) myPair = std::make_pair( &(*tauColl)[0], &(*tauColl)[1]);
    else myPair = std::make_pair( &(*tauColl)[1], &(*tauColl)[0]);    
  }
  */
  /*find a most isolated opposite-sign pair */
  std::vector<const pat::Tau*> taus;
  for(size_t i=0; i<tauColl->size(); ++i){
    if( (*tauColl)[i].tauID("decayModeFindingNewDMs")>0.5 &&
	( (*tauColl)[i].decayMode()<WawGenInfoHelper::tauDecayModes::tauDecay2ChargedPion0PiZero ||
	  (*tauColl)[i].decayMode()>WawGenInfoHelper::tauDecayModes::tauDecay2ChargedPion4PiZero )
	&& (*tauColl)[i].pt()>20 && std::abs((*tauColl)[i].eta())<2.1
	&& std::abs(dynamic_cast<pat::PackedCandidate const*>((*tauColl)[i].leadChargedHadrCand().get())->dz()) < 0.2
	&& std::abs((*tauColl)[i].charge()) == 1
	&& (*tauColl)[i].tauID("againstMuonLoose3") > 0.5
	//&& (*tauColl)[i].tauID("photonPtSumOutsideSignalCone") < 0.1 * (*tauColl)[i].pt()
	//&& (*tauColl)[i].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 3.0 //Loose <2.5, Medium<1.5, Tight<0.8
	&& (*tauColl)[i].tauID("byVVLooseIsolationMVArun2v1DBnewDMwLT2017v2") > 0.5
	)
      taus.push_back(&(*tauColl)[i]);
  }
  if(taus.size()<2) return myPair;
  std::sort(taus.begin(), taus.end(),
	    [](const pat::Tau* a, const pat::Tau* b){
	      //return a->pt() > b->pt();
	      //return a->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") <
	      //       b->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
	      return a->tauID("byIsolationMVArun2v1DBnewDMwLTraw2017v2") >
		     b->tauID("byIsolationMVArun2v1DBnewDMwLTraw2017v2");
	    });
  for(size_t i=1; i<taus.size(); ++i){
    if(deltaR2(taus[0]->p4(),taus[i]->p4())>0.5*0.5 &&
       taus[0]->charge()*taus[i]->charge()==-1){
      if(taus[0]->charge()==1)
	myPair = std::make_pair(taus[0], taus[i]);
      else
	myPair = std::make_pair(taus[i], taus[0]);
      break;
    }
  }

  return myPair;
}
////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool MiniAODVertexAnalyzer::diMuVeto(const edm::Handle<std::vector<pat::Muon> > &muColl,
				     const edm::Handle<reco::VertexCollection> &vtxColl){

  if(muColl->size()<2 || vtxColl->empty()) return false;

  //select muons
  std::vector<const pat::Muon*> mus;
  for(size_t i=0; i<muColl->size(); ++i){
    if( (*muColl)[i].pt()>15 && std::abs((*muColl)[i].eta())<2.4
	&& std::abs((*muColl)[i].muonBestTrack()->dxy((*vtxColl)[0].position()))<0.045
	&& std::abs((*muColl)[i].muonBestTrack()->dz((*vtxColl)[0].position())) <0.2
	&& (*muColl)[i].isGlobalMuon()
	&& (*muColl)[i].isTrackerMuon()
	&& (*muColl)[i].isPFMuon()
	&& muRelIso(&(*muColl)[i])<0.3
	)
      mus.push_back(&(*muColl)[i]);
  }
  if(mus.size()<2) return false;
  int nMuPair = 0;
  for(size_t i=0; i<mus.size()-1; ++i){
    for(size_t j=i; j<mus.size(); ++j){
      if(deltaR2(mus[i]->p4(),mus[j]->p4())>0.15*0.15 &&
	 mus[i]->charge()*mus[j]->charge()==-1){
	nMuPair++;
      }
    }
  }
  if(nMuPair<1) return false;

  return true;
}
////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool MiniAODVertexAnalyzer::diEVeto(const edm::Handle<std::vector<pat::Electron> > &eColl,
				    const edm::Handle<reco::VertexCollection> &vtxColl){

  if(eColl->size()<2 || vtxColl->empty()) return false;

  //select electrons
  std::vector<const pat::Electron*> eles;
  for(size_t i=0; i<eColl->size(); ++i){
    if( (*eColl)[i].pt()>15 && std::abs((*eColl)[i].eta())<2.5
	&& std::abs((*eColl)[i].gsfTrack()->dxy((*vtxColl)[0].position()))<0.045
	&& std::abs((*eColl)[i].gsfTrack()->dz((*vtxColl)[0].position())) <0.2
	&& (*eColl)[i].electronID("mvaEleID-Fall17-noIso-V1-wpLoose")>0.5 //by default cut-based is used, the "cutBasedElectronID-Fall17-94X-V1-veto" one?
	&& eRelIso(&(*eColl)[i])<0.3
	)
      eles.push_back(&(*eColl)[i]);
  }
  if(eles.size()<2) return false;
  int nElePair = 0;
  for(size_t i=0; i<eles.size()-1; ++i){
    for(size_t j=i; j<eles.size(); ++j){
      if(deltaR2(eles[i]->p4(),eles[j]->p4())>0.15*0.15 &&
	 eles[i]->charge()*eles[j]->charge()==-1){
	nElePair++;
      }
    }
  }
  if(nElePair<1) return false;

  return true;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool  MiniAODVertexAnalyzer::getPVBalance(const edm::Event & iEvent, size_t iPV){

  edm::Handle<edm::View<pat::PackedCandidate> >  cands;
  iEvent.getByToken(cands_, cands);
  edm::Handle<edm::View<pat::PackedCandidate> >  lostCands;
  if(useLostCands_) iEvent.getByToken(lostCands_, lostCands);
  //edm::Handle<reco::VertexCollection> vertices;
  //iEvent.getByToken(vertices_, vertices);

  //Get tracks associated with iPV
  TVector3 diTauPt;
  if(thePair_.first!=nullptr && thePair_.second!=nullptr){
    float px = thePair_.first->px()+thePair_.second->px();
    float py = thePair_.first->py()+thePair_.second->py();
    diTauPt.SetXYZ(px,py,0);
  }
  float pt2Sum=0, ptSum=0, ptThrust=0;

  std::vector<reco::Candidate::LorentzVector> signalTauCandsP4;
  if(thePair_.first!=nullptr && thePair_.second!=nullptr){
    //t_1
    const pat::Tau* aTau = dynamic_cast<const pat::Tau*>(thePair_.first);
    if(aTau!=nullptr){
      for(size_t j=0 ; j<aTau->signalChargedHadrCands().size(); ++j){
	signalTauCandsP4.push_back(aTau->signalChargedHadrCands()[j]->p4());
      }
    }
    else{
      const pat::Muon* aMu = dynamic_cast<const pat::Muon*>(thePair_.first);
      if(aMu!=nullptr && aMu->originalObjectRef().isNonnull()){
	signalTauCandsP4.push_back(aMu->originalObjectRef()->p4());
      }
      else{
	const pat::Electron* aEle = dynamic_cast<const pat::Electron*>(thePair_.first);
	if(aEle!=nullptr && aEle->associatedPackedPFCandidates().size()>0){
	  reco::Candidate::LorentzVector aEleP4;
	  for(size_t i=0; i<aEle->associatedPackedPFCandidates().size(); ++i){
	    if(aEle->associatedPackedPFCandidates()[i]->pdgId()==aEle->pdgId() &&
	       aEle->associatedPackedPFCandidates()[i]->pt()>aEleP4.pt()){
	      aEleP4 = aEle->associatedPackedPFCandidates()[i]->p4();
	    }
	  }
	  if(aEleP4.pt()>0.5)
	    signalTauCandsP4.push_back(aEleP4);
	}
      }
    }
    //t_2
    aTau = dynamic_cast<const pat::Tau*>(thePair_.second);
    if(aTau!=nullptr){
      for(size_t j=0 ; j<aTau->signalChargedHadrCands().size(); ++j){
	signalTauCandsP4.push_back(aTau->signalChargedHadrCands()[j]->p4());
      }
    }
    else{
      const pat::Muon* aMu = dynamic_cast<const pat::Muon*>(thePair_.second);
      if(aMu!=nullptr && aMu->originalObjectRef().isNonnull()){
	signalTauCandsP4.push_back(aMu->originalObjectRef()->p4());
      }
      else{
	const pat::Electron* aEle = dynamic_cast<const pat::Electron*>(thePair_.second);
	if(aEle!=nullptr && aEle->associatedPackedPFCandidates().size()>0){
	  reco::Candidate::LorentzVector aEleP4;
	  for(size_t i=0; i<aEle->associatedPackedPFCandidates().size(); ++i){
	    if(aEle->associatedPackedPFCandidates()[i]->pdgId()==aEle->pdgId() &&
	       aEle->associatedPackedPFCandidates()[i]->pt()>aEleP4.pt()){
	      aEleP4 = aEle->associatedPackedPFCandidates()[i]->p4();
	    }
	  }
	  if(aEleP4.pt()>0.5)
	    signalTauCandsP4.push_back(aEleP4);
	}
      }
    }
  }
  for(size_t i=0; i<cands->size(); ++i){
    if((*cands)[i].charge()==0 || (*cands)[i].hasTrackDetails() || (*cands)[i].vertexRef().isNull()) continue;
    ///Skip tracks comming from tau decay.
    bool skipTrack = false;
    for(size_t j=0 ; j<signalTauCandsP4.size();++j){
      if(deltaR2(signalTauCandsP4[j],(*cands)[i].p4())<0.001*0.001
	 && std::abs(signalTauCandsP4[j].pt()/(*cands)[i].pt()-1.)<0.001
	 ){
	skipTrack = true;
	break;
      }
    }
    if(skipTrack) continue;
    size_t key = (*cands)[i].vertexRef().key();
    int quality = (*cands)[i].pvAssociationQuality();
    if(key!=iPV ||
       (quality!=pat::PackedCandidate::UsedInFitTight &&
	quality!=pat::PackedCandidate::UsedInFitLoose)) continue;

    pt2Sum += (*cands)[i].pt() * (*cands)[i].pt();
    ptSum += (*cands)[i].pt();
    TVector3 aPt((*cands)[i].px(),(*cands)[i].py(),0);
    ptThrust -= aPt * diTauPt.Unit();
  }
  if(useLostCands_){
    for(size_t i=0; i<lostCands->size(); ++i){
      if( ( (*cands)[i].charge()==0 && !(*lostCands)[i].hasTrackDetails() ) || (*lostCands)[i].vertexRef().isNull()) continue;
      ///Skip tracks comming from tau decay.
      bool skipTrack = false;
      for(size_t j=0 ; j<signalTauCandsP4.size();++j){
	if(deltaR2(signalTauCandsP4[j],(*lostCands)[i].p4())<0.001*0.001
	   && std::abs(signalTauCandsP4[j].pt()/(*lostCands)[i].pt()-1.)<0.001
	   ){
	  skipTrack = true;
	  break;
	}
      }
      if(skipTrack) continue;
      size_t key = (*lostCands)[i].vertexRef().key();
      int quality = (*lostCands)[i].pvAssociationQuality();
      if(key!=iPV ||
	 (quality!=pat::PackedCandidate::UsedInFitTight &&
	  quality!=pat::PackedCandidate::UsedInFitLoose)) continue;

      pt2Sum += (*cands)[i].pt() * (*cands)[i].pt();
      ptSum += (*cands)[i].pt();
      TVector3 aPt((*cands)[i].px(),(*cands)[i].py(),0);
      ptThrust -= aPt * diTauPt.Unit();
    }
  }
  

  myEvent_->recoEvent_.pt2Sum_=pt2Sum;
  myEvent_->recoEvent_.ptThrust_=ptThrust;
  myEvent_->recoEvent_.ptBalance_= (ptSum - diTauPt.Pt())/(ptSum + diTauPt.Pt());
  
  return true;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void MiniAODVertexAnalyzer::beginJob(){}
void MiniAODVertexAnalyzer::endJob(){}

