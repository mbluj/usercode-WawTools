#ifdef PROJECT_NAME
#include "WarsawAnalysis/HTTDataFormats/interface/HTTEvent.h"
#include "WarsawAnalysis/HTT/interface/GenInfoHelper.h"
#else
#include "HTTEvent.h"
#endif

//////////////////////////////////////////////
//////////////////////////////////////////////
HTTEvent::HTTEvent(){

  clear();
  
}
//////////////////////////////////////////////
//////////////////////////////////////////////
HTTEvent::~HTTEvent(){;}
//////////////////////////////////////////////
//////////////////////////////////////////////
void HTTEvent::clear(){

  run_ = -999;
  lumi_ = -999;
  event_ = -999;
  bosonId_ = 0;

  nPV_ = 0;
  diMuonVeto_ = false;

  genEvent_.clear();
  recoEvent_.clear();

}
//////////////////////////////////////////////
//////////////////////////////////////////////
DiTauData::DiTauData() : 
  pfScore_(3), pt2Score_(2), dzVtx_(2){
  
  clear();
}
//////////////////////////////////////////////
//////////////////////////////////////////////
DiTauData::~DiTauData(){}
//////////////////////////////////////////////
//////////////////////////////////////////////
void DiTauData::clear(){

#ifdef PROJECT_NAME
  decModeMinus_ = WawGenInfoHelper::tauDecayModes::tauDecayOther;
  decModePlus_ = WawGenInfoHelper::tauDecayModes::tauDecayOther;
#else
  decModeMinus_  = 99;
  decModePlus_  = 99;
#endif

  dzPlus_ = 99;
  dzMinus_ = 99;

  thePV_ = TVector3();
  svMinus_ = TVector3();
  svPlus_ = TVector3();
  nPiPlus_ = TVector3();
  nPiMinus_ = TVector3();

  p4Sum_ = TLorentzVector();
  piMinus_ = TLorentzVector();
  piPlus_ = TLorentzVector();
  tauMinus_ = TLorentzVector();
  tauPlus_ = TLorentzVector();
  pi0Minus_ = TLorentzVector();
  pi0Plus_ = TLorentzVector();
  visTauMinus_ = TLorentzVector();
  visTauPlus_ = TLorentzVector();

  pfPV_ = TVector3();
  pt2PV_ = TVector3();
  refitPfPV_ = TVector3();
  refitPfPVNoBS_ = TVector3();

  isRefit_ = false;
  nTracksInRefit_ = 0;
  
  pfPVIndex_  = -1;
  pt2PVindex_ = -1;

  leadIdMinus_ = 0;
  leadIdPlus_ = 0;

  nGammaMinus_ = nGammaPlus_ = 0;
  nGammaInConeMinus_ = nGammaInConePlus_ = 0;

  for(size_t i=0; i<pfScore_.size(); ++i) pfScore_[i]=0; 
  for(size_t i=0; i<pt2Score_.size(); ++i) pt2Score_[i]=0; 
  for(size_t i=0; i<dzVtx_.size(); ++i) dzVtx_[i]=0; 

  isoMinus_ = 999;
  isoPlus_ = 999;

  isoMVAWpMinus_ = 0;
  isoMVAWpPlus_ = 0;

  antiEWpMinus_  = antiEWpPlus_  = 0;
  antiMuWpMinus_ = antiMuWpPlus_ = 0;

  matchedMinus_ = 0;
  matchedPlus_ = 0;

  pt2Sum_ = 0;
  ptThrust_ = 0;
  ptBalance_ = 99;

}
//////////////////////////////////////////////
//////////////////////////////////////////////
