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

MiniAODVertexAnalyzer::MiniAODVertexAnalyzer(const edm::ParameterSet & iConfig) :
  scores_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("vertexScores"))),
  cands_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("src"))),
  genCands_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genSrc"))),
  vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  bs_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  //threshold_(iConfig.getParameter<int>("threshold")),
  useBeamSpot_(iConfig.getParameter<bool>("useBeamSpot")),
  verbose_(iConfig.getUntrackedParameter<bool>("verbose",false)){

  
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("HTT","HTT");

  myEvent_ = 0; //Important pointer initialisation
  tree_->Branch("HTTEvent", &myEvent_); 
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
  WawGenInfoHelper::getVertex(theBoson,&myEvent_->genPV_);

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

  
  
  myEvent_->decModeMinus_ = WawGenInfoHelper::getTausDecays(taus[0],tauProdsMinus,true,false); 
  myEvent_->tauMinus_.SetPtEtaPhiE(taus[0]->pt(),taus[0]->eta(),taus[0]->phi(),taus[0]->energy());

  reco::GenParticleRef piMinusRef = WawGenInfoHelper::getLeadChParticle(tauProdsMinus);
  myEvent_->piMinus_.SetPtEtaPhiE(piMinusRef->pt(),piMinusRef->eta(),piMinusRef->phi(),piMinusRef->energy());

  WawGenInfoHelper::setP4Ptr(WawGenInfoHelper::getCombinedP4(tauProdsMinus), &myEvent_->visTauMinus_);
  WawGenInfoHelper::getVertex(WawGenInfoHelper::getLeadChParticle(tauProdsMinus),&myEvent_->svMinus_);

  myEvent_->decModePlus_ = WawGenInfoHelper::getTausDecays(taus[1],tauProdsPlus,true,false);
  myEvent_->tauPlus_.SetPtEtaPhiE(taus[1]->pt(),taus[1]->eta(),taus[1]->phi(),taus[1]->energy());

  reco::GenParticleRef piPlusRef = WawGenInfoHelper::getLeadChParticle(tauProdsPlus);
  myEvent_->piPlus_.SetPtEtaPhiE(piPlusRef->pt(),piPlusRef->eta(),piPlusRef->phi(),piPlusRef->energy());

  WawGenInfoHelper::setP4Ptr(WawGenInfoHelper::getCombinedP4(tauProdsPlus), &myEvent_->visTauPlus_);
  WawGenInfoHelper::getVertex(WawGenInfoHelper::getLeadChParticle(tauProdsPlus),&myEvent_->svPlus_);
  WawGenInfoHelper::setP4Ptr(myEvent_->tauMinus_+myEvent_->tauPlus_,&myEvent_->p4Sum_);

  return true;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool  MiniAODVertexAnalyzer::getRecoData(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  bool hasAnyVx = findPrimaryVertices(iEvent, iSetup);

  refitPV(iEvent, iSetup);

  return hasAnyVx;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool  MiniAODVertexAnalyzer::refitPV(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  edm::Handle<edm::View<pat::PackedCandidate> >  cands;
  iEvent.getByToken(cands_, cands);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertices_, vertices);
  edm::Handle<reco::BeamSpot> beamSpot;
  if(useBeamSpot_) iEvent.getByToken(bs_, beamSpot);

  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);
 
  TransientVertex transVtx;

  //Get tracks associated wiht pfPV
  reco::TrackCollection pvTracks;
  for(size_t i=0; i<cands->size(); ++i){
    if((*cands)[i].charge()==0 || (*cands)[i].vertexRef().isNull()) continue;      
      int key = (*cands)[i].vertexRef().key();
      int quality = (*cands)[i].pvAssociationQuality();
      if(key!=myEvent_->pfPVIndex_ ||
	 (quality!=pat::PackedCandidate::UsedInFitTight &&
	  quality!=pat::PackedCandidate::UsedInFitLoose)) continue;
      if(!(*cands)[i].bestTrack()) continue;
      pvTracks.push_back(*((*cands)[i].bestTrack()));
  }
  ///Built transient tracks from tracks.
  std::vector<reco::TransientTrack> transTracks;  
  for(auto iter: pvTracks) transTracks.push_back(transTrackBuilder->build(iter));
  myEvent_->nTracksInRefit_ = transTracks.size();

  bool fitOk = false;  
  if(transTracks.size() >= 2 ) {
    AdaptiveVertexFitter avf;
    avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception    
    try {
      if( !useBeamSpot_ ) transVtx = avf.vertex(transTracks);
      else transVtx = avf.vertex(transTracks, *beamSpot);
      fitOk = true; 
    } catch (...) {
      fitOk = false; 
      std::cout<<"Vtx fit failed!"<<std::endl;
    }
  }
  
  if(fitOk && transVtx.isValid() ) { 
     myEvent_->rePfPV_.SetXYZ(transVtx.position().x(),transVtx.position().y(),transVtx.position().z());
     myEvent_->isRefit_=true;
  }
  else {
     myEvent_->rePfPV_.SetXYZ((*vertices)[myEvent_->pfPVIndex_].x(),(*vertices)[myEvent_->pfPVIndex_].y(),(*vertices)[myEvent_->pfPVIndex_].z());
     myEvent_->isRefit_=false;
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
  myEvent_->aodPV_.SetXYZ((*vertices)[0].x(),(*vertices)[0].y(),(*vertices)[0].z());

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

  myEvent_->pfPV_.SetXYZ((*vertices)[iPfVtx].x(),(*vertices)[iPfVtx].y(),(*vertices)[iPfVtx].z());
  myEvent_->pfPVIndex_=iPfVtx;
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
	//<<" key: "<<(*cands)[i].vertexRef().key()
	//<<" fromPV: "<<(*cands)[i].fromPV()
	       <<std::endl;
    }

    if((*cands)[i].charge()==0 || (*cands)[i].vertexRef().isNull()) continue;
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

  myEvent_->pt2PV_.SetXYZ((*vertices)[iOldVtx2].x(),(*vertices)[iOldVtx2].y(),(*vertices)[iOldVtx2].z());
  myEvent_->pt2PVindex_=iOldVtx2;
  iOldVtx1+=0;//hack to avoid compilation errors -Werror=unused-but-set-variable

  return true;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void MiniAODVertexAnalyzer::beginJob(){}
void MiniAODVertexAnalyzer::endJob(){}
