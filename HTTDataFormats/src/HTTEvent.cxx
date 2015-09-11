#include "WarsawAnalysis/HTTDataFormats/interface/HTTEvent.h"

#include "WarsawAnalysis/HTT/interface/GenInfoHelper.h"

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
  decModeMinus_ = WawGenInfoHelper::kUndefined;
  decModePlus_ = WawGenInfoHelper::kUndefined;
  
  pfPVIndex_  = -1;
  pt2PVindex_ = -1;
  
  isRefit_ = false;
  nTracksInRefit_ = 0;
  
  p4Sum_ = TLorentzVector();
  piMinus_ = TLorentzVector();
  piPlus_ = TLorentzVector();
  tauMinus_ = TLorentzVector();
  tauPlus_ = TLorentzVector();
  visTauMinus_ = TLorentzVector();
  visTauPlus_ = TLorentzVector();
  
  genPV_ = TVector3();
  aodPV_ = TVector3();
  pfPV_ = TVector3();
  pt2PV_ = TVector3();
  rePfPV_ = TVector3();
  svMinus_ = TVector3();
  svPlus_ = TVector3();  
}
//////////////////////////////////////////////
//////////////////////////////////////////////
