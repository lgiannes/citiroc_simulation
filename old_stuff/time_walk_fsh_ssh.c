
#include "time_walk_simulation.c"


TGraph* time_walk_graph(string filename, bool SSH = true){

  vector<double> fractions = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.15};

  TGraph* g = new TGraph();

  double threshold = 3.5;

  TGraph *g_charge = new TGraph();

  for (double fraction : fractions){

    TH1D* time_histo = getTimeHisto(filename, 1);

    time_histo->Scale(fraction/time_histo->GetBinWidth(1));

    int NPhotons = time_histo->GetEntries()*fraction;

    TH1D* FS = fastShaperOutputLaplace(time_histo);
    TH1D* SS = slowShaperOutput(time_histo);
    SS->Scale(TMath::Exp(1));
    FS->Scale(TMath::Exp(1));
    double time_walk;
    if ( computeTime(FS, threshold) == -1 ){
      time_walk = 0;
      cout<<" --- Signal does not reach threshold!"<<endl;
    } else {
      time_walk = computeTime(FS, threshold)-computeTime(FS, 0.);
    }

    double charge;
    if(SSH){
      charge = SS->GetMaximum();
      cout<<"SSH: q="<<charge<<"  t_w="<<time_walk<<endl;
    } else {
      charge = FS->GetMaximum();
      cout<<"FSH: q="<<charge<<"  t_w="<<time_walk<<endl;
    }
    if(time_walk > 0) {
      g->SetPoint(g->GetN(), charge, time_walk);
    }
    if(fraction == 0.15) {

      // Draw and save the time histogram
      TCanvas * c = new
      TCanvas("c", "c", 800, 600);
      c->cd();
      time_histo->SetFillColor(kBlack);
      time_histo->SetLineColor(kBlack);
      time_histo->Draw("hist");
      FS->SetLineColor(kBlue);
      SS->SetLineColor(kRed);
      FS->SetLineWidth(3);
      SS->SetLineWidth(2);
      FS->Draw("same hist");
      SS->Draw("same hist");
      TLatex *l = new TLatex();
      l->SetNDC();
      l->SetTextSize(0.04);
      l->DrawLatex(0.1, 0.9, Form("NPhotons = %d", NPhotons));
      l->DrawLatex(0.1, 0.85, Form("Max FS = %.2f", FS->GetMaximum()));
      l->DrawLatex(0.1, 0.8, Form("Max SS = %.2f", SS->GetMaximum()));
      c->SaveAs(Form("time_histo_%s.png", filename.c_str()));

    }
  }


  return g;
}


double ftimewalkCorrection2(double *x, double *p){

  std::vector<double> par = {p[0], p[1], p[2], p[3], p[4]};
//  if(isData){
//    par = {1.003, 0.0042, 1.0136, 114.9, 1.463}; // Data
//  }else{
//    par = {1.37,0.017,1.04,0.0008,-3.0}; // MC
//  }
  return par[0]/log(par[1]*x[0]+par[2]) - par[3]/(x[0] + par[4]);
}

void time_walk_fsh_ssh(int dist = 1800) {

  std::vector<string> filenames = {"outputs/output200.root", "outputs/output400.root", "outputs/output600.root", "outputs/output800.root", "outputs/output1000.root",
                                   "outputs/output1200.root", "outputs/output1400.root", "outputs/output1600.root", "outputs/output1800.root"};


  string filename = Form("outputs/output%d.root", dist);
  TGraph *g_SSH = time_walk_graph(filename, true);
  TGraph *g_FSH = time_walk_graph(filename, false);

  g_SSH->SetMarkerStyle(20);
  g_SSH->SetMarkerSize(1);
  g_SSH->SetMarkerColor(kRed);
  g_SSH->SetTitle(Form("Time walk, D = %d mm; Charge [a.u.]; Time walk [ns]", 2000-dist));
  g_SSH->Draw("AP");
  TF1 *f = new TF1("f", ftimewalkCorrection2, 0, 240, 5);
  f->SetParameters(1.003, 0.0042, 1.0136, 114.9, 1.463);
  g_SSH->Fit("f");
  g_SSH->GetFunction("f")->SetLineColor(kRed);
  g_FSH->SetMarkerStyle(20);
  g_FSH->SetMarkerSize(1);
  g_FSH->SetMarkerColor(kBlue);
  g_FSH->Draw("P");
  TF1 *f2 = new TF1("f2", ftimewalkCorrection2, 0, 240, 5);
  f2->SetParameters(1.003, 0.0042, 1.0136, 114.9, 1.463);
  g_FSH->Fit("f2");
  g_FSH->GetFunction("f2")->SetLineColor(kBlue);

  TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
  legend->AddEntry(g_SSH, "SSH", "p");
  legend->AddEntry(g_FSH, "FSH", "p");
  legend->Draw();

  for(int i=0;i<20;i++){
    double q = 15 + i*5;
    double residual = f->Eval(q) - f2->Eval(q);
    cout<<"q: "<<q<<"  res: "<<residual<<" ns"<<endl;
  }

}

