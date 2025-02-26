#include "measure_threshold_time.c"

TH1D* getTimeHisto(string filename, double fraction){

  TFile *f = new TFile(filename.c_str());
  TTree *tree = (TTree*)f->Get("tree");
  double ipos[3],fpos[3],ivec[3],fvec[3],time,length;
  int imat, fmat, ftype, nref, npas;
  tree->SetBranchAddress("time", &time);
  tree->SetBranchAddress("ftype", &ftype);

  TH1D* time_histo = new TH1D("time_histo", "photon ToA", 240, 0, 60);

  int nentries = tree->GetEntries();
//  cout<<"File: "<<filename<<"  fraction: "<<fraction<<endl;
  for (int i=0; i<nentries; i++) {
    tree->GetEntry(i);
    if(ftype==4) {
      if(gRandom->Rndm() > fraction) continue;
      time_histo->Fill(time);
    }
  }
  time_histo->SetName(Form("time_histo_%s", filename.c_str()));

  return time_histo;

}



void time_walk_simulation() {


  std::vector<int> distances = {400,600,800,1000,1200,1400,1600,1800};
  std::vector<int> LY = {20,30,40,50,60,70,80,90};

  double threshold = 3.5;
  std::map<pair<int,int>, double> max_SSH;
  std::map<pair<int,int>, double> max_FSH;
  std::map<pair<int,int>, double> Nphotons;
  std::map<pair<int,int>, double> t0;
  std::map<pair<int,int>, double> t_threshold;

  gROOT->SetBatch(kTRUE);

  for(int dist : distances){
    for(int ly : LY){
      string filename = Form("outputs/output%d_%d.root", dist, ly);
      TH1D* time_histo = getTimeHisto(filename, 1);
      TCanvas *c = new TCanvas("c", "c", 800, 600);
      c->cd();
      time_histo->SetFillColor(kBlack);
      time_histo->SetLineColor(kBlack);
      time_histo->Draw("hist");
      double binWidth = time_histo->GetBinWidth(1);
      time_histo->Scale(1./binWidth);
      TH1D* FS = fastShaperOutputLaplace(time_histo);
      TH1D* SS = slowShaperOutput(time_histo);
      FS->SetLineColor(kBlue);
      SS->SetLineColor(kRed);
      FS->SetLineWidth(2);
      SS->SetLineWidth(2);
      FS->Scale(TMath::Exp(1));
      SS->Scale(TMath::Exp(1));
      FS->Draw("same hist");
      SS->Draw("same hist");
      TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
      legend->AddEntry(time_histo, "i(t)", "f");
      legend->AddEntry(FS, "Fast shaper", "l");
      legend->AddEntry(SS, "Slow shaper", "l");
      time_histo->SetTitle(Form("Signal simulation for distance = %d mm, LY = %d", dist, ly));
      time_histo->GetXaxis()->SetTitle("Time [ns]");
      c->SaveAs(Form("time_histo_%s.png", filename.c_str()));

      double max_FS = FS->GetMaximum();
      double max_SS = SS->GetMaximum();
      double this_t0 = computeTime(FS, 0.); // requires fast shaped output as argument
      double this_t_threshold = computeTime(FS, threshold); // requires fast shaped output as argument

      max_SSH[make_pair(dist, ly)] = max_SS;
      max_FSH[make_pair(dist, ly)] = max_FS;
      Nphotons[make_pair(dist, ly)] = time_histo->GetEntries();
      t0[make_pair(dist, ly)] = this_t0;
      t_threshold[make_pair(dist, ly)] = this_t_threshold;
      cout<<"N photons: "<<time_histo->GetEntries()<<" max FS = "<<max_FS<<"  max SS = "<<max_SS<<"  t0 = "<<this_t0<<"  t_threshold = "<<this_t_threshold<<endl;
    }
  }

  gROOT->SetBatch(kFALSE);
//
//  TGraph *g_max_SSH = new TGraph();
//  TGraph *g_max_FSH = new TGraph();
//  TGraph *g_Nphotons = new TGraph();
//  int fix_dist = 400;
//  for(int ly : LY){
//    g_max_SSH->SetPoint(g_max_SSH->GetN(), Nphotons[make_pair(fix_dist, ly)], max_SSH[make_pair(fix_dist, ly)]);
//    g_max_FSH->SetPoint(g_max_FSH->GetN(), Nphotons[make_pair(fix_dist, ly)], max_FSH[make_pair(fix_dist, ly)]);
//    g_Nphotons->SetPoint(g_Nphotons->GetN(), Nphotons[make_pair(fix_dist, ly)], Nphotons[make_pair(fix_dist, ly)]);
//  }
//  TCanvas *c_graphs = new TCanvas("c_graphs", "c_graphs", 800, 600);
//  c_graphs->cd();
//  g_max_SSH->SetMarkerStyle(20);
//  g_max_SSH->SetMarkerSize(2);
//  g_max_SSH->SetMarkerColor(kRed);
//  g_max_SSH->Draw("AP");
//  g_max_SSH->SetTitle("Maximum slow shaper height; LY; Maximum slow shaper height");
//  g_max_FSH->SetMarkerStyle(20);
//  g_max_FSH->SetMarkerSize(1);
//  g_max_FSH->SetMarkerColor(kBlue);
//  g_max_FSH->Draw("P");
//  g_Nphotons->SetMarkerStyle(20);
//  g_Nphotons->SetMarkerSize(1);
//  g_Nphotons->SetMarkerColor(kGreen);
//  g_Nphotons->Draw("P");
//  TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
//  legend->AddEntry(g_max_SSH, "Maximum slow shaper height", "p");
//  legend->AddEntry(g_max_FSH, "Maximum fast shaper height", "p");
//  legend->AddEntry(g_Nphotons, "Number of photons", "p");
//  legend->Draw();

  std::vector<TGraph*> graphs;
  TMultiGraph *mg = new TMultiGraph();
  TGraph *g_speed = new TGraph();
  TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
  for(int dist : distances){
    graphs.push_back(new TGraph());
    graphs.back()->SetName(Form("g_dist%d", dist));
    graphs.back()->SetMarkerStyle(20);
    graphs.back()->SetMarkerSize(1.);
    graphs.back()->SetMarkerColor(graphs.size());
    for(int ly : LY){
      double time_walk = t_threshold[make_pair(dist, ly)] - t0[make_pair(dist, ly)];
      double charge = Nphotons[make_pair(dist, ly)];
      double max_FS = max_FSH[make_pair(dist, ly)];
      double max_SS = max_SSH[make_pair(dist, ly)];
      double Nph = Nphotons[make_pair(dist, ly)];
      graphs.back()->SetPoint(graphs.back()->GetN(), Nph, time_walk);
      if(ly == 90){
        g_speed->SetPoint(g_speed->GetN(), dist, t_threshold[make_pair(dist, ly)]);
      }
    }
    if(dist == 400 || dist == 1800 || dist == 1000)
    mg->Add(graphs.back());
    legend->AddEntry(graphs.back(), Form("Distance = %d mm", 2000-dist), "p");
  }
  mg->SetTitle("; Charge [p.e.]; Time walk [ns]");
  mg->Draw("APL");
  legend->Draw();

  TCanvas *c_speed = new TCanvas("c_speed", "c_speed", 800, 600);
  c_speed->cd();
  g_speed->SetMarkerStyle(20);
  g_speed->SetMarkerSize(1);
  g_speed->SetMarkerColor(kBlack);
  g_speed->SetTitle("speed; Distance [mm]; Time alk [ns]");
  g_speed->Draw("AP");


  int fix_dist = 1000;
  TGraph *g_tw_SSH = new TGraph();
  TGraph *g_tw_FSH = new TGraph();
  for(int ly : LY){
    g_tw_SSH->SetPoint(g_tw_SSH->GetN(), max_SSH[make_pair(fix_dist, ly)], t_threshold[make_pair(fix_dist, ly)] - t0[make_pair(fix_dist, ly)]);
    g_tw_FSH->SetPoint(g_tw_FSH->GetN(), max_FSH[make_pair(fix_dist, ly)], t_threshold[make_pair(fix_dist, ly)] - t0[make_pair(fix_dist, ly)]);
  }

  TCanvas *c_tw = new TCanvas("c_tw", "c_tw", 800, 600);
  c_tw->cd();
  TLegend *legend_tw = new TLegend(0.1, 0.7, 0.3, 0.9);
  g_tw_SSH->SetMarkerStyle(20);
  g_tw_SSH->SetMarkerSize(1);
  g_tw_SSH->SetMarkerColor(kRed);
  g_tw_SSH->SetTitle(Form("Distance = %d mm; peak [p.e.]; Time walk [ns]", 2000 - fix_dist));
  g_tw_SSH->Draw("AP");
  g_tw_FSH->SetMarkerStyle(20);
  g_tw_FSH->SetMarkerSize(1);
  g_tw_FSH->SetMarkerColor(kBlue);
  g_tw_FSH->Draw("P");
  legend_tw->AddEntry(g_tw_SSH, "Slow shaper", "p");
  legend_tw->AddEntry(g_tw_FSH, "Fast shaper", "p");
  legend_tw->Draw();




}


void shapeDeltaFunction(){
  //  // shape a "delta function"
  TH1D* delta = new TH1D("delta", "delta", 400, 0, 100);
  double binWidth = delta->GetBinWidth(1);
  delta->Fill(0);
  delta->Scale(1./binWidth);
  TH1D* FS = fastShaperOutputLaplace(delta);
  TH1D* SS = slowShaperOutput(delta);
  TCanvas *test = new TCanvas("test", "test", 800, 600);
  test->cd();
  delta->SetFillColor(kMagenta);
  delta->Draw();
  FS->SetLineColor(kBlue);
  SS->SetLineColor(kRed);
  FS->SetLineWidth(2);
  SS->SetLineWidth(2);
  FS->Draw("same");
  SS->Draw("same");
  cout<<"FS max = "<<FS->GetMaximum()<<"  SS max = "<<SS->GetMaximum()<<endl;
}
