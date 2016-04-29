void plotVtx(std::string name, bool noBS=false, bool pfPV=false){

  if(name==""){
    std::cout<<"Provide name of file"<<std::endl;
    return;
  }

  TFile *f = TFile::Open((name+std::string(".root")).c_str());
  if(!f){
    std::cout<<"Invalid file"<<std::endl;
    return;
  }

  TTree *t = (TTree*)f->Get("vtxAna/HTT");

  std::string sel = "";

  TH1F *h[3][2];
  h[0][0] = new TH1F("h00",";PV^{gen}(x) - PV^{reco}(x) [cm]; a.u.",100,-0.025,0.025);
  //h[0][0]->StatOverflows(1);
  h[0][1] = new TH1F("h01",";PV^{gen}(z) - PV^{reco}(z) [cm]; a.u.",100,-0.025,0.025);
  //h[0][1]->StatOverflows(1);
  for(int i=1; i<3; ++i){
    h[i][0] = (TH1F*)h[0][0]->Clone(Form("h%i0",i));
    //h[i][0]->StatOverflows(1);
    h[i][1] = (TH1F*)h[0][1]->Clone(Form("h%i1",i));
    //h[i][1]->StatOverflows(1);
  }
  if(!pfPV){
    t->Project("h00","HTTEvent.genEvent_.thePV_.x()-HTTEvent.recoEvent_.thePV_.x()",sel.c_str());
    t->Project("h01","HTTEvent.genEvent_.thePV_.z()-HTTEvent.recoEvent_.thePV_.z()",sel.c_str());
  }
  else{
    t->Project("h00","HTTEvent.genEvent_.thePV_.x()-HTTEvent.recoEvent_.pfPV_.x()",sel.c_str());
    t->Project("h01","HTTEvent.genEvent_.thePV_.z()-HTTEvent.recoEvent_.pfPV_.z()",sel.c_str());
  }
  t->Project("h10","HTTEvent.genEvent_.thePV_.x()-HTTEvent.recoEvent_.refitPfPV_.x()",sel.c_str());
  t->Project("h11","HTTEvent.genEvent_.thePV_.z()-HTTEvent.recoEvent_.refitPfPV_.z()",sel.c_str());
  t->Project("h20","HTTEvent.genEvent_.thePV_.x()-HTTEvent.recoEvent_.refitPfPVNoBS_.x()",sel.c_str());
  t->Project("h21","HTTEvent.genEvent_.thePV_.z()-HTTEvent.recoEvent_.refitPfPVNoBS_.z()",sel.c_str());
  
  
  int colors[3] = {kRed, kBlue, kGreen+2};
  for(int i=0; i<3; ++i){
    for(int j=0; j<2; ++j){
      h[i][j]->SetStats(0);
      h[i][j]->Scale(1./h[i][j]->Integral(0,h[i][j]->GetNbinsX()+1));
      h[i][j]->SetLineColor(colors[i]);
    }
  }


  TCanvas *c = new TCanvas("c","",600,600);
  c->SetLogy();

  TLegend *leg;
  if(!noBS)
    leg = new TLegend(0.57,0.8,0.9,0.9);
  else
    leg = new TLegend(0.57,0.75,0.9,0.9);

  h[0][0]->Draw();
  if(!pfPV)
    leg->AddEntry(h[0][0],Form("PV^{AOD}: #mu=%.3e, RMS=%f", h[0][0]->GetMean(), h[0][0]->GetStdDev()),"l");
  else
    leg->AddEntry(h[0][0],Form("PV^{PF}: #mu=%.3e, RMS=%f", h[0][0]->GetMean(), h[0][0]->GetStdDev()),"l");
  h[1][0]->Draw("same");
  leg->AddEntry(h[1][0],Form("PV^{refit}: #mu=%.3e, RMS=%f", h[1][0]->GetMean(), h[1][0]->GetStdDev()),"l");
  if(noBS) {
    h[2][0]->Draw("same");
    leg->AddEntry(h[2][0],Form("PV^{refit}_{noBS}: #mu=%.3e, RMS=%f", h[2][0]->GetMean(), h[2][0]->GetStdDev()),"l");
  }
  leg->Draw();
  c->Print((std::string("fig_")+name+std::string("_vtx-x.png")).c_str());
  leg->Clear();

  h[0][1]->Draw();
  if(!pfPV)
    leg->AddEntry(h[0][1],Form("PV^{AOD}: #mu=%.3e, RMS=%f", h[0][1]->GetMean(), h[0][1]->GetStdDev()),"l");
  else
    leg->AddEntry(h[0][1],Form("PV^{PF}: #mu=%.3e, RMS=%f", h[0][1]->GetMean(), h[0][1]->GetStdDev()),"l");
  h[1][1]->Draw("same");
  leg->AddEntry(h[1][1],Form("PV^{refit}: #mu=%.3e, RMS=%f", h[1][1]->GetMean(), h[1][1]->GetStdDev()),"l");
  if(noBS) {
    h[2][1]->Draw("same");
    leg->AddEntry(h[2][1],Form("PV^{refit}_{noBS}: #mu=%.3e, RMS=%f", h[2][1]->GetMean(), h[2][1]->GetStdDev()),"l");
  }
  leg->Draw();
  c->Print((std::string("fig_")+name+std::string("_vtx-z.png")).c_str());
  leg->Clear();

  //
  f->Close();

  return;
}
