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
  verbose_(iConfig.getUntrackedParameter<bool>("verbose",false)),
  bosonId_(0),
  decModeMinus_(WawGenInfoHelper::kUndefined),
  decModePlus_(WawGenInfoHelper::kUndefined),
  ipfPV_(0), ioldPV_(0), refit_(0)
{
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("hTTVtxTree","hTTVtxTree");
  //Book integers
  tree_->Branch("bosonId", &bosonId_, "bosonId/I");
  tree_->Branch("decModeMinus", &decModeMinus_, "decModeMinus/I");
  tree_->Branch("decModePlus", &decModePlus_, "decModePlus/I");
  tree_->Branch("ipfPV", &ipfPV_, "ipfPV/I");
  tree_->Branch("ioldPV", &ioldPV_, "ioldPV/I");
  tree_->Branch("refit", &refit_, "refit/I");
  //Book floats
  tree_->Branch("run", &run_, "run/F");
  tree_->Branch("lumi", &lumi_, "lumi/F");
  tree_->Branch("event", &event_, "event/F");
  //Book TLorentzVectors
  p4Sum_ = new TLorentzVector();
  tree_->Branch("p4Sum.","TLorentzVector",p4Sum_);
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
  //Book TVector3's
  genPV_ = new TVector3();
  tree_->Branch("genPV.","TVector3",genPV_);
  aodPV_ = new TVector3();
  tree_->Branch("aodPV.","TVector3",aodPV_);
  pfPV_ = new TVector3();
  tree_->Branch("pfPV.","TVector3",pfPV_);
  rePfPV_ = new TVector3();
  tree_->Branch("rePfPV.","TVector3",rePfPV_);
  oldPV_ = new TVector3();
  tree_->Branch("oldPV.","TVector3",oldPV_);
  svMinus_ = new TVector3();
  tree_->Branch("svMinus.","TVector3",svMinus_);
  svPlus_ = new TVector3();
  tree_->Branch("svPlus.", "TVector3",svPlus_);
}

MiniAODVertexAnalyzer::~MiniAODVertexAnalyzer(){ 
  delete genPV_;
  delete aodPV_;
  delete pfPV_;
  delete rePfPV_;
  delete oldPV_;

  delete p4Sum_;
  delete piMinus_;
  delete piPlus_;
  delete tauMinus_;
  delete tauPlus_;
  delete visTauMinus_;
  delete visTauPlus_;

  delete svMinus_;
  delete svPlus_;

}

void MiniAODVertexAnalyzer::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // Basic event info
  run_ = iEvent.id().run();
  lumi_ = iEvent.id().luminosityBlock();
  event_ = iEvent.id().event();

  edm::Handle<edm::View<pat::PackedCandidate> >  cands;
  iEvent.getByToken(cands_, cands);
  edm::Handle<reco::GenParticleCollection> genCands;
  iEvent.getByToken(genCands_, genCands);
  edm::Handle<edm::ValueMap<float> > scores;
  iEvent.getByToken(scores_, scores);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertices_, vertices);
  edm::Handle<reco::BeamSpot> beamSpot;
  if(useBeamSpot_)
    iEvent.getByToken(bs_, beamSpot);

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
  if(theBoson.isNull() ) return; //FIXME: need to handle cases w/o the boson?
  bosonId_ = theBoson->pdgId();
  WawGenInfoHelper::getVertex(theBoson,genPV_);

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
  decModePlus_ = WawGenInfoHelper::getTausDecays(taus[1],tauProdsPlus,true,false);
  tauPlus_->SetPxPyPzE(taus[1]->px(),taus[1]->py(),taus[1]->pz(),taus[1]->energy());
  reco::GenParticleRef piPlusRef = WawGenInfoHelper::getLeadChParticle(tauProdsPlus);
  piPlus_->SetPxPyPzE(piPlusRef->px(),piPlusRef->py(),piPlusRef->pz(),piPlusRef->energy());
  WawGenInfoHelper::setP4Ptr(WawGenInfoHelper::getCombinedP4(tauProdsPlus),
                             visTauPlus_);
  WawGenInfoHelper::getVertex(WawGenInfoHelper::getLeadChParticle(tauProdsPlus),svPlus_);
  WawGenInfoHelper::setP4Ptr((*tauMinus_)+(*tauPlus_), p4Sum_);

  /*************/
  //at least one vertex
  if(vertices->size()==0) return;
  //vertex with higherst score as in AOD
  aodPV_->SetXYZ((*vertices)[0].x(),(*vertices)[0].y(),(*vertices)[0].z());
  //find vertex with highest score with PF (miniAOD like)
  size_t iPfVtx=0;
  float score=-1;
  for(size_t i=0; i<vertices->size(); ++i){
    reco::VertexRef vtxPrt(vertices,i);
    /*
    std::cout<<"No of tracks: "<<vtxPrt->tracksSize()<<std::endl;
    for(std::vector<reco::TrackBaseRef >::const_iterator itt=vtxPrt->tracks_begin();
	itt!=vtxPrt->tracks_end(); ++itt){
      if( (*itt).isNonnull() )
	std::cout<<"\t pt="<<(*itt)->pt()<<"+-"<<(*itt)->ptError()<<std::endl;
      else
	std::cout<<"\t Null\n";
    }
    */
    if( (*scores)[vtxPrt] > score){
      score = (*scores)[vtxPrt];
      iPfVtx=i;
    }
  }
  pfPV_->SetXYZ((*vertices)[iPfVtx].x(),(*vertices)[iPfVtx].y(),(*vertices)[iPfVtx].z());
  ipfPV_=iPfVtx;
  std::map<size_t,float> oldScores1, oldScores2;
  for(size_t i=0; i<vertices->size(); ++i){
    oldScores1[i]=0.;
    oldScores2[i]=0.;
  }
  for(size_t i=0; i<cands->size(); ++i){
    /*
    std::cout<<i<<"."
	     <<" pdgId="<<(*cands)[i].pdgId()
	     <<" q="<<(*cands)[i].charge()
	     <<" pt="<<(*cands)[i].pt()<<std::endl;
    if((*cands)[i].bestTrack())
      std::cout<<"\tTrack:"
	       <<" pt="<<(*cands)[i].bestTrack()->pt()
	       <<"+-"<<(*cands)[i].bestTrack()->ptError()<<std::endl;
    else
      std::cout<<"\tW/o track"<<std::endl;
    */
    if((*cands)[i].charge()!=0 && (*cands)[i].vertexRef().isNonnull() ){
      size_t key = (*cands)[i].vertexRef().key();
      int quality = (*cands)[i].pvAssociationQuality();
      //std::cout<<"asocc vtx: key/Q: "<<key<<"/"<<quality<<std::endl;
      if( quality==pat::PackedCandidate::UsedInFitTight ||
	  quality==pat::PackedCandidate::UsedInFitLoose){
	if((*cands)[i].bestTrack()){
	  float pT = (*cands)[i].bestTrack()->pt();
	  float epT=(*cands)[i].bestTrack()->ptError(); 
	  pT=pT>epT ? pT-epT : 0;
	  oldScores1[key]+=pT*pT;
	  oldScores2[key]+=pT*pT;	  
	}
	else { //should happen for tracks pT<0.95
	  float pT = (*cands)[i].pt();
	  //estimate pT error (basing on avarage numbers from TRK-11-001, JINST 9 (2014), P10009
	  float epT=0;
	  if(std::abs((*cands)[i].eta())<0.9)
	    epT=0.01*pT;
	  else if(std::abs((*cands)[i].eta())<1.4)
	    epT=0.02*pT;
	  else
	    epT=0.025*pT;
	  if(pT>90) epT*=2; //for high-Pt assume bigger error (just to be safe)
	  if(pT>1){
	    std::cout<<"Something unexpected: Pt="<<pT
		     <<", eta="<<(*cands)[i].eta()
		     <<", charge="<<(*cands)[i].charge()
		     <<", pdgId="<<(*cands)[i].pdgId()<<std::endl;
	  }
	  pT=pT>epT ? pT-epT : 0;
	  oldScores2[key]+=pT*pT;
        }
      }
    }
  }
  size_t iOldVtx1=0, iOldVtx2=0;
  float oldScore1=-1, oldScore2=-1;
  for(std::map<size_t,float>::const_iterator it=oldScores1.begin(); 
       it!=oldScores1.end(); ++it){
    //std::cout<<"\toldScores1: "<<it->first<<". "<<it->second<<std::endl;
    if(it->second > oldScore1){
      oldScore1 = it->second;
      iOldVtx1 = it->first;
    }
  }
  for(std::map<size_t,float>::const_iterator it=oldScores2.begin(); 
       it!=oldScores2.end(); ++it){
    //std::cout<<"\toldScores2: "<<it->first<<". "<<it->second<<std::endl;
    if(it->second > oldScore2){
      oldScore2 = it->second;
      iOldVtx2 = it->first;
    }
  }
  if(verbose_ && 
     (iPfVtx!=0 || 
      iOldVtx1!=0 || iOldVtx2!=0 ||
      iOldVtx1!=iPfVtx || iOldVtx2!=iPfVtx) ){
    reco::VertexRef vtxPrt0(vertices,0);
    std::cout<<"Different association of PV with different scoring:\n"
	     <<"\tnew score with PF : "<<score<<" at "<<iPfVtx<<std::endl
	     <<"\tnew score with PF : "<<(*scores)[vtxPrt0]<<" at 0"<<std::endl
	     <<"\told score 1 : "<<oldScore1<<" at "<<iOldVtx1<<std::endl
             <<"\told score 1 : "<<oldScores1[0]<<" at 0"<<std::endl
	     <<"\told score 2 : "<<oldScore2<<" at "<<iOldVtx2<<std::endl
	     <<"\told score 2 : "<<oldScores2[0]<<" at 0"<<std::endl;
  }
  oldPV_->SetXYZ((*vertices)[iOldVtx2].x(),(*vertices)[iOldVtx2].y(),(*vertices)[iOldVtx2].z());
  ioldPV_=iOldVtx2;
  /*************/
  // vertex refit
  TransientVertex transVtx;
  std::vector<reco::TransientTrack> transTracks;
  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);
  //Get tracks from pfPV
  reco::TrackCollection pvTracks;
  for(size_t i=0; i<cands->size(); ++i){
    if((*cands)[i].charge()!=0 && (*cands)[i].vertexRef().isNonnull() ){
      size_t key = (*cands)[i].vertexRef().key();
      int quality = (*cands)[i].pvAssociationQuality();
      if( key==iPfVtx &&
	  (quality==pat::PackedCandidate::UsedInFitTight ||
	   quality==pat::PackedCandidate::UsedInFitLoose) ){
	if((*cands)[i].bestTrack())
	  pvTracks.push_back(*((*cands)[i].bestTrack()));
      }
    }
  }
  for(reco::TrackCollection::iterator iter=pvTracks.begin(); 
      iter!=pvTracks.end(); ++iter){
    transTracks.push_back(transTrackBuilder->build(*iter));
  }

  bool fitOk = true;
  if(transTracks.size() >= 2 ) {
    AdaptiveVertexFitter avf;
    avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception    
    try {
      if( !useBeamSpot_ ){
	transVtx = avf.vertex(transTracks);
      } else {
	transVtx = avf.vertex(transTracks, *beamSpot);
      }
    } catch (...) {
      fitOk = false; 
      std::cout<<"Vtx fit failed!"<<std::endl;
    }
  } else fitOk = false; 
  if( fitOk && transVtx.isValid() ) { 
    rePfPV_->SetXYZ(transVtx.position().x(),transVtx.position().y(),transVtx.position().z());
    //std::cout<<"OK!"<<std::endl;
    refit_=1;
  }
  else {
    rePfPV_->SetXYZ((*vertices)[iPfVtx].x(),(*vertices)[iPfVtx].y(),(*vertices)[iPfVtx].z());
    //std::cout<<"not OK!"<<std::endl;
    refit_=0;
  }

  /*
  if(cands.product()->size()>0) {
    //Find PV with highest score
    int ivtx=0;
    edm::ProductID pid = (*cands.product())[0].vertexRef().id();
    edm::ValueMap<float>::const_iterator itVM = scores->begin();
    for(; itVM!= scores->end();++itVM) {
      if(pid==itVM.id()) {
	float score=-1;
	edm::ValueMap<float>::container::const_iterator itVtx =  itVM.begin();
	for(unsigned int i=0; itVtx!=itVM.end(); ++itVtx,i++){
	  if(*itVtx > score) {ivtx=i; score=*itVtx;}
	  //				    std::cout << " i,imax score,scoremax" << i << " " << ivtx << " " << *itVtx << " "<< score << std::endl;
	}
      }
      //	    std::cout << "Selected " << ivtx << std::endl;
      outPV->push_back((*(*cands.product())[0].vertexRef().product())[ivtx]);
    }
	    //Take only cands with fromPV >= thr relative to the found vertex
	    for (edm::View<pat::PackedCandidate>::const_iterator it = cands->begin(), ed = cands->end(); it != ed; ++it) {
		    if(it->fromPV(ivtx) >= threshold_){
			    out->push_back(cands->ptrAt(it-cands->begin()));
		    }
	    }
    }
  */
  //at the end fill tree 
  tree_->Fill();
}

void MiniAODVertexAnalyzer::beginJob(){}
void MiniAODVertexAnalyzer::endJob(){}

// Define plugin
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MiniAODVertexAnalyzer);
