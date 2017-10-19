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
  genCands_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genSrc"))),
  vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  bs_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))), 
  taus_(consumes<std::vector<pat::Tau> >(iConfig.getParameter<edm::InputTag>("taus"))),
  useBeamSpot_(iConfig.getParameter<bool>("useBeamSpot")),
  useLostCands_(iConfig.getParameter<bool>("useLostCands")),
  useTauTracks_(iConfig.getUntrackedParameter<bool>("useTauTracks",false)),
  verbose_(iConfig.getUntrackedParameter<bool>("verbose",false)){

  
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("HTT","HTT");

  myEvent_ = 0; //Important pointer initialisation
  tree_->Branch("HTTEvent", &myEvent_); 

  thePair_ = std::make_pair<const pat::Tau*, const pat::Tau*>(0,0); 
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
  thePair_ = std::make_pair<const pat::Tau*, const pat::Tau*>(0,0);
  
  // Basic event info
  myEvent_->run_ = iEvent.id().run();
  myEvent_->lumi_ = iEvent.id().luminosityBlock();
  myEvent_->event_ = iEvent.id().event();

  bool goodGenData = getGenData(iEvent, iSetup);
  bool goodRecoData = getRecoData(iEvent, iSetup);

  if(goodGenData && goodRecoData) tree_->Fill();

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
  
  myEvent_->genEvent_.decModeMinus_ = WawGenInfoHelper::getTausDecays(taus[0],tauProdsMinus,true,false); 
  myEvent_->genEvent_.tauMinus_.SetPtEtaPhiE(taus[0]->pt(),taus[0]->eta(),taus[0]->phi(),taus[0]->energy());

  reco::GenParticleRef piMinusRef = WawGenInfoHelper::getLeadChParticle(tauProdsMinus);
  myEvent_->genEvent_.piMinus_.SetPtEtaPhiE(piMinusRef->pt(),piMinusRef->eta(),piMinusRef->phi(),piMinusRef->energy());

  WawGenInfoHelper::setP4Ptr(WawGenInfoHelper::getCombinedP4(tauProdsMinus), &myEvent_->genEvent_.visTauMinus_);
  WawGenInfoHelper::getVertex(WawGenInfoHelper::getLeadChParticle(tauProdsMinus),&myEvent_->genEvent_.svMinus_);
  WawGenInfoHelper::impactParameter(myEvent_->genEvent_.thePV_, myEvent_->genEvent_.svMinus_,
				    myEvent_->genEvent_.piMinus_, &myEvent_->genEvent_.nPiMinus_);
  myEvent_->genEvent_.dzMinus_ = (myEvent_->genEvent_.svMinus_.Z() - myEvent_->genEvent_.thePV_.Z()) - ((myEvent_->genEvent_.svMinus_.X() - myEvent_->genEvent_.thePV_.X()) * myEvent_->genEvent_.piMinus_.Px() + (myEvent_->genEvent_.svMinus_.Y() - myEvent_->genEvent_.thePV_.Y()) * myEvent_->genEvent_.piMinus_.Py()) / myEvent_->genEvent_.piMinus_.Pt() * myEvent_->genEvent_.piMinus_.Pz() / myEvent_->genEvent_.piMinus_.Pt();

  myEvent_->genEvent_.decModePlus_ = WawGenInfoHelper::getTausDecays(taus[1],tauProdsPlus,true,false);
  myEvent_->genEvent_.tauPlus_.SetPtEtaPhiE(taus[1]->pt(),taus[1]->eta(),taus[1]->phi(),taus[1]->energy());

  reco::GenParticleRef piPlusRef = WawGenInfoHelper::getLeadChParticle(tauProdsPlus);
  myEvent_->genEvent_.piPlus_.SetPtEtaPhiE(piPlusRef->pt(),piPlusRef->eta(),piPlusRef->phi(),piPlusRef->energy());

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
  bool hasRecoTau = findRecoTau(iEvent, iSetup);

  refitPV(iEvent, iSetup);
  setPCAVectors(iEvent, iSetup);

  return hasAnyVx && hasRecoTau;
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
  for(size_t i=0; i<cands->size(); ++i){
    if(!(*cands)[i].bestTrack() || (*cands)[i].vertexRef().isNull()) continue;
    ///Skip tracks comming from tau decay.
    bool skipTrack = false;
    if(!useTauTracks_){
      if(thePair_.first && thePair_.second){
	for(size_t j=0 ; j<thePair_.first->signalChargedHadrCands().size(); ++j){
	  if(deltaR2(thePair_.first->signalChargedHadrCands()[j]->p4(),(*cands)[i].p4())<0.001*0.001
	     && std::abs(thePair_.first->signalChargedHadrCands()[j]->pt()/(*cands)[i].pt()-1.)<0.001
	     ){
	    skipTrack = true;
	    break;
	  }
	}
	for(size_t j=0 ; j<thePair_.second->signalChargedHadrCands().size(); ++j){
	  if(deltaR2(thePair_.second->signalChargedHadrCands()[j]->p4(),(*cands)[i].p4())<0.001*0.001
	     && std::abs(thePair_.second->signalChargedHadrCands()[j]->pt()/(*cands)[i].pt()-1.)<0.001
	     ){
	    skipTrack = true;
	    break;
	  }
	}
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
      if(!(*lostCands)[i].bestTrack() || (*lostCands)[i].vertexRef().isNull()) continue;
      ///Skip tracks comming from tau decay.
      bool skipTrack = false;
      if(!useTauTracks_){
	if(thePair_.first && thePair_.second){
	  for(size_t j=0 ; j<thePair_.first->signalChargedHadrCands().size(); ++j){
	    if(deltaR2(thePair_.first->signalChargedHadrCands()[j]->p4(),(*lostCands)[i].p4())<0.001*0.001
	       && std::abs(thePair_.first->signalChargedHadrCands()[j]->pt()/(*lostCands)[i].pt()-1.)<0.001
	       ){
	      skipTrack = true;
	      break;
	    }
	  }
	  for(size_t j=0 ; j<thePair_.second->signalChargedHadrCands().size(); ++j){
	    if(deltaR2(thePair_.second->signalChargedHadrCands()[j]->p4(),(*lostCands)[i].p4())<0.001*0.001
	       && std::abs(thePair_.second->signalChargedHadrCands()[j]->pt()/(*lostCands)[i].pt()-1.)<0.001
	       ){
	      skipTrack = true;
	      break;
	    }
	  }
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

  if(vertices->size()==0) return false;   //at least one vertex
  myEvent_->recoEvent_.thePV_.SetXYZ((*vertices)[0].x(),(*vertices)[0].y(),(*vertices)[0].z());

  //Find vertex with highest score with PF (miniAOD like)
  //miniAOD uses PF particles instead of tracks
  size_t iPfVtx=0;
  float score=-1;
  for(size_t iVx=0; iVx<vertices->size(); ++iVx){
    reco::VertexRef vtxPrt(vertices,iVx);   
    if( (*scores)[vtxPrt] > score){
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

    if((*cands)[i].charge()!=0 &&
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

    if((*cands)[i].charge()==0 || !(*cands)[i].bestTrack() || (*cands)[i].vertexRef().isNull()) continue;
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
		 <<", (*cands)[i].bestTrack(): "<<(*cands)[i].bestTrack()
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

  return true;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool MiniAODVertexAnalyzer::findRecoTau(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  edm::Handle<std::vector<pat::Tau> > tauColl;
  iEvent.getByToken(taus_,tauColl);
  if(tauColl->size()<2) return false;

  thePair_ = findTauPair(tauColl);
  if(!thePair_.first || !thePair_.second) return false;

  myEvent_->recoEvent_.isoMVAWpPlus_ = 0;
  if(thePair_.first->tauID("byVLooseIsolationMVArun2v1DBoldDMwLT") > 0.5)
     myEvent_->recoEvent_.isoMVAWpPlus_ = 1;
  if(thePair_.first->tauID("byLooseIsolationMVArun2v1DBoldDMwLT") > 0.5)
     myEvent_->recoEvent_.isoMVAWpPlus_ = 2;
  if(thePair_.first->tauID("byMediumIsolationMVArun2v1DBoldDMwLT") > 0.5)
     myEvent_->recoEvent_.isoMVAWpPlus_ = 3;
  if(thePair_.first->tauID("byTightIsolationMVArun2v1DBoldDMwLT") > 0.5)
     myEvent_->recoEvent_.isoMVAWpPlus_ = 4;
  if(thePair_.first->tauID("byVTightIsolationMVArun2v1DBoldDMwLT") > 0.5)
     myEvent_->recoEvent_.isoMVAWpPlus_ = 5;
  if(thePair_.first->tauID("byVVTightIsolationMVArun2v1DBoldDMwLT") > 0.5)
     myEvent_->recoEvent_.isoMVAWpPlus_ = 6;
  myEvent_->recoEvent_.isoMVAWpMinus_ = 0;
  if(thePair_.second->tauID("byVLooseIsolationMVArun2v1DBoldDMwLT") > 0.5)
     myEvent_->recoEvent_.isoMVAWpMinus_ = 1;
  if(thePair_.second->tauID("byLooseIsolationMVArun2v1DBoldDMwLT") > 0.5)
     myEvent_->recoEvent_.isoMVAWpMinus_ = 2;
  if(thePair_.second->tauID("byMediumIsolationMVArun2v1DBoldDMwLT") > 0.5)
     myEvent_->recoEvent_.isoMVAWpMinus_ = 3;
  if(thePair_.second->tauID("byTightIsolationMVArun2v1DBoldDMwLT") > 0.5)
     myEvent_->recoEvent_.isoMVAWpMinus_ = 4;
  if(thePair_.second->tauID("byVTightIsolationMVArun2v1DBoldDMwLT") > 0.5)
     myEvent_->recoEvent_.isoMVAWpMinus_ = 5;
  if(thePair_.second->tauID("byVVTightIsolationMVArun2v1DBoldDMwLT") > 0.5)
     myEvent_->recoEvent_.isoMVAWpMinus_ = 6;
  myEvent_->recoEvent_.decModePlus_ = thePair_.first->decayMode();
  myEvent_->recoEvent_.decModeMinus_ = thePair_.second->decayMode();
  
  myEvent_->recoEvent_.dzPlus_ = dynamic_cast<pat::PackedCandidate const*>(thePair_.first->leadChargedHadrCand().get())->dz();
  myEvent_->recoEvent_.dzMinus_ = dynamic_cast<pat::PackedCandidate const*>(thePair_.second->leadChargedHadrCand().get())->dz();
  
  myEvent_->recoEvent_.visTauPlus_ = TLorentzVector(thePair_.first->p4().Px(),
						    thePair_.first->p4().Py(),
						    thePair_.first->p4().Pz(),
						    thePair_.first->p4().E());
  
   myEvent_->recoEvent_.visTauMinus_ = TLorentzVector(thePair_.second->p4().Px(),
						      thePair_.second->p4().Py(),
						      thePair_.second->p4().Pz(),
						      thePair_.second->p4().E());
   
   myEvent_->recoEvent_.tauPlus_ = myEvent_->recoEvent_.visTauPlus_;
   myEvent_->recoEvent_.tauMinus_ = myEvent_->recoEvent_.visTauMinus_;
  
   myEvent_->recoEvent_.piPlus_ = TLorentzVector(thePair_.first->leadChargedHadrCand()->p4().px(),
						 thePair_.first->leadChargedHadrCand()->p4().py(),
						 thePair_.first->leadChargedHadrCand()->p4().pz(),
						 thePair_.first->leadChargedHadrCand()->p4().e());
   
   myEvent_->recoEvent_.piMinus_ = TLorentzVector(thePair_.second->leadChargedHadrCand()->p4().px(),
						  thePair_.second->leadChargedHadrCand()->p4().py(),
						  thePair_.second->leadChargedHadrCand()->p4().pz(),
						  thePair_.second->leadChargedHadrCand()->p4().e());
     
  return true;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool MiniAODVertexAnalyzer::setPCAVectors(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  if(!thePair_.first || !thePair_.second) return false;

  GlobalPoint aPoint(myEvent_->recoEvent_.refitPfPV_.X(),
		     myEvent_->recoEvent_.refitPfPV_.Y(),
		     myEvent_->recoEvent_.refitPfPV_.Z());

  myEvent_->recoEvent_.nPiPlus_ = getPCA(iEvent,iSetup, thePair_.first, aPoint);
  myEvent_->recoEvent_.nPiMinus_ = getPCA(iEvent,iSetup, thePair_.second, aPoint);

   ///Test different PVs.
   ///AOD PV
  aPoint = GlobalPoint(myEvent_->recoEvent_.thePV_.X(),
		       myEvent_->recoEvent_.thePV_.Y(),
		       myEvent_->recoEvent_.thePV_.Z());
  myEvent_->recoEvent_.nPiPlusAODvx_ = getPCA(iEvent,iSetup, thePair_.first, aPoint);
  myEvent_->recoEvent_.nPiMinusAODvx_ = getPCA(iEvent,iSetup, thePair_.second, aPoint);
  
  ///PV refit excluding tau tracks
  aPoint = GlobalPoint(myEvent_->recoEvent_.refitPfPV_.X(),
		       myEvent_->recoEvent_.refitPfPV_.Y(),
		       myEvent_->recoEvent_.refitPfPV_.Z());
  myEvent_->recoEvent_.nPiPlusRefitvx_ = getPCA(iEvent,iSetup, thePair_.first, aPoint);
  myEvent_->recoEvent_.nPiMinusRefitvx_ = getPCA(iEvent,iSetup, thePair_.second, aPoint);
  
  ///Generated PV
   aPoint = GlobalPoint(myEvent_->genEvent_.thePV_.X(),
			myEvent_->genEvent_.thePV_.Y(),
			myEvent_->genEvent_.thePV_.Z());
   myEvent_->recoEvent_.nPiPlusGenvx_ = getPCA(iEvent,iSetup, thePair_.first, aPoint);
   myEvent_->recoEvent_.nPiMinusGenvx_ = getPCA(iEvent,iSetup, thePair_.second, aPoint);
  
   return true;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
TVector3 MiniAODVertexAnalyzer::getPCA(const edm::Event & iEvent, const edm::EventSetup & iSetup,
				       const pat::Tau* aTau,
				       const GlobalPoint & aPoint){

  TVector3 aPCA;
  if(!aTau->leadChargedHadrCand()->bestTrack()) return aPCA;
  
  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);  
  reco::TransientTrack transTrk=transTrackBuilder->build(aTau->leadChargedHadrCand()->bestTrack());
     
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
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
std::pair<const pat::Tau*, const pat::Tau*> MiniAODVertexAnalyzer::findTauPair(edm::Handle<std::vector<pat::Tau> > tauColl){

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
  for(unsigned int i=0; i<tauColl->size(); ++i){
    if( (*tauColl)[i].tauID("decayModeFinding")>0.5
	&& (*tauColl)[i].pt()>20 && std::abs((*tauColl)[i].eta())<2.3
	&& std::abs(dynamic_cast<pat::PackedCandidate const*>((*tauColl)[i].leadChargedHadrCand().get())->dz()) < 0.2
	//&& (*tauColl)[i].tauID("photonPtSumOutsideSignalCone") < 0.1 * (*tauColl)[i].pt()
	//&& (*tauColl)[i].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 3.0 //Loose <2.5, Medium<1.5, Tight<0.8
	&& (*tauColl)[i].tauID("byVLooseIsolationMVArun2v1DBoldDMwLT") > 0.5
	)
      taus.push_back(&(*tauColl)[i]);
  }
  if(taus.size()<2) return myPair;
  std::sort(taus.begin(), taus.end(),
	    [](const pat::Tau* a, const pat::Tau* b){
	      //return a->pt() > b->pt();
	      //return a->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 
	      //       b->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
	      return a->tauID("byIsolationMVArun2v1DBoldDMwLTraw") >
		     b->tauID("byIsolationMVArun2v1DBoldDMwLTraw");
	    });
  for(unsigned int i=1; i<taus.size(); ++i){
    if(deltaR2(taus[0]->p4(),taus[i]->p4())>0.3*0.3 &&
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
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void MiniAODVertexAnalyzer::beginJob(){}
void MiniAODVertexAnalyzer::endJob(){}

