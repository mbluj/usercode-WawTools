#ifdef PROJECT_NAME
#include "WarsawAnalysis/HTTDataFormats/interface/HTTEvent.h"
#include "WarsawAnalysis/HTT/interface/GenInfoHelper.h"
#else
#include "HTTEvent.h"
#endif

//////////////////////////////////////////////
//////////////////////////////////////////////
HTTEvent::HTTEvent() :
  tauSpinnerWeight_(0) {

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

  for(size_t i=0; i<tauSpinnerWeight_.size(); ++i) tauSpinnerWeight_[i]=1;

  nPV_ = 0;
  diMuonVeto_ = false;
  diEleVeto_ = false;

  genEvent_.clear();
  recoEvent_.clear();

}
//////////////////////////////////////////////
//////////////////////////////////////////////
DiTauData::DiTauData() : 
  pfScore_(3), pt2Score_(2), dzVtx_(3){
  
  clear();
}
//////////////////////////////////////////////
//////////////////////////////////////////////
DiTauData::~DiTauData(){}
//////////////////////////////////////////////
//////////////////////////////////////////////
void DiTauData::clear(){

#ifdef PROJECT_NAME
  decModeMVAMinus_ = decModeMinus_ = WawGenInfoHelper::tauDecayModes::tauDecayOther;
  decModeMVAPlus_ = decModePlus_ = WawGenInfoHelper::tauDecayModes::tauDecayOther;
#else
  decModeMVAMinus_ = decModeMinus_  = 99;
  decModeMVAPlus_ = decModePlus_  = 99;
#endif

  dzPlus_ = 99;
  dzMinus_ = 99;

  thePV_ = TVector3();
  svMinus_ = TVector3();
  svPlus_ = TVector3();
  nPiMinus_ = TVector3();
  nPiPlus_ = TVector3();

  p4Sum_ = TLorentzVector();
  piMinus_ = TLorentzVector();
  piPlus_ = TLorentzVector();
  trkMinus_ = TLorentzVector();
  trkPlus_ = TLorentzVector();
  scMinus_ = TLorentzVector();
  scPlus_ = TLorentzVector();
  neutralMinus_ = TLorentzVector();
  neutralPlus_ = TLorentzVector();
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
  nStripMinus_ = nStripPlus_ = 0;

  gammaPtSumInScMinus_ = gammaPtSumInScPlus_ = 0;
  gammaPtSumOutScMinus_ = gammaPtSumOutScPlus_ = 0;
  maxGammaPtMinus_ = maxGammaPtPlus_ = 0;
  dR2midMinus_ = dR2midPlus_ = 99;


  for(size_t i=0; i<pfScore_.size(); ++i) pfScore_[i]=0; 
  for(size_t i=0; i<pt2Score_.size(); ++i) pt2Score_[i]=0; 
  for(size_t i=0; i<dzVtx_.size(); ++i) dzVtx_[i]=0; 

  isoMinus_ = 999;
  isoPlus_ = 999;

  isoMVAWpMinus_ = 0;
  isoMVAWpPlus_ = 0;

  antiEWpMinus_  = antiEWpPlus_  = 0;
  antiE2WpMinus_  = antiE2WpPlus_  = 0;
  antiMuWpMinus_ = antiMuWpPlus_ = 0;

  matchedMinus_ = 0;
  matchedPlus_ = 0;

  pt2Sum_ = 0;
  ptThrust_ = 0;
  ptBalance_ = 99;

}
//////////////////////////////////////////////
//////////////////////////////////////////////
