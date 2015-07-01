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


MiniAODVertexAnalyzer::MiniAODVertexAnalyzer(const edm::ParameterSet & iConfig) :
  scores_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("vertexScores"))),
  cands_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("src"))),
  genCands_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genSrc"))),
  vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  //threshold_(iConfig.getParameter<int>("threshold")),
  verbose_(iConfig.getUntrackedParameter<bool>("verbose",false)),
  bosonId_(0),
  decModeMinus_(WawGenInfoHelper::kUndefined),
  decModePlus_(WawGenInfoHelper::kUndefined)
{
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("hTTVtxTree","hTTVtxTree");
  //Book integers
  tree_->Branch("bosonId", &bosonId_, "bosonId/I");
  tree_->Branch("decModeMinus", &decModeMinus_, "decModeMinus/I");
  tree_->Branch("decModePlus", &decModePlus_, "decModePlus/I");
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
    if( (*scores)[vtxPrt] > score){
      score = (*scores)[vtxPrt];
      iPfVtx=i;
    }
  }
  pfPV_->SetXYZ((*vertices)[iPfVtx].x(),(*vertices)[iPfVtx].y(),(*vertices)[iPfVtx].z());
  if(verbose_ && iPfVtx!=0){
    reco::VertexRef vtxPrt0(vertices,0);
    std::cout<<"PF-sorting gives different result that sorting from AOD: "<<std::endl
	     <<"\tscore,iVtx with PF : "<<score<<", "<<iPfVtx<<std::endl
	     <<"\tscore for iVtx=0 : "<<(*scores)[vtxPrt0]<<std::endl;
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
