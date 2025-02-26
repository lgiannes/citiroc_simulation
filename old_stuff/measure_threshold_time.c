
TH1D* getTimeHisto(string filename);
double computeCharge(TH1D* h);
double computeTime(TH1D* h, double threshold);

TH1D* fastShaperOutput(TH1D* h);
int colorArray[] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kPink, kViolet, kRed+1, kBlack};
int I_color_array = 0;
void measure_threshold_time(){

  double threshold = 3.5*0.055; // 3.5 pe

  // Create a canvas
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

  std::vector<string> filenames = {"outputs/output200.root", "outputs/output400.root", "outputs/output600.root", "outputs/output800.root", "outputs/output1000.root",
                                   "outputs/output1200.root", "outputs/output1400.root", "outputs/output1600.root", "outputs/output1800.root"};
  std::vector<double> distances = {1800, 1600, 1400, 1200, 1000, 800, 600, 400, 200};
  std::vector<TH1D*> time_histos;
  TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
  std::vector<double> integrals;
  for (int i=0; i<filenames.size(); i++) {
    time_histos.push_back(getTimeHisto(filenames[i]));
    time_histos.back()->SetLineColor(colorArray[i]);
    time_histos.back()->SetLineWidth(2);
    time_histos.back()->SetName(Form("time_histo_%d", i));
    legend->AddEntry(time_histos.back(), filenames[i].c_str(), "l");
    int bin_at_sht = time_histos.back()->FindBin(12.5);
    integrals.push_back(time_histos.back()->Integral(0, bin_at_sht));
    std::cout<<filenames[i]<<" integral: "<<integrals.back()<<"  total integral: "<<time_histos.back()->Integral()<<std::endl;
  }
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  time_histos[0]->Draw();
  for (int i=1; i<time_histos.size(); i++) {
    time_histos[i]->Draw("same");
  }
  legend->Draw();

  // Define shaping function (e.g., CR-RC shaping)
  double tau = 15;  // Time constant
  int nBins = time_histos[0]->GetNbinsX();
  double binWidth = time_histos[0]->GetBinWidth(1);
  TH1D *hShaper = new TH1D("hShaper", "Shaping Function", nBins, 0, nBins * binWidth);
  for (int i = 1; i <= nBins; ++i) {
    double t = hShaper->GetBinCenter(i);
    double value = (t >= 0) ? (t / tau) * exp(-t / tau) : 0;  // CR-RC example
    hShaper->SetBinContent(i, value);
  }

  std::vector<TH1D*> shapedSignals;
  for(int i=0; i<time_histos.size(); i++) {
    // Perform discrete convolution
    TH1D *hOutput = new TH1D("hOutput", "Convolved Signal", nBins, 0, nBins * binWidth);
    TH1D* hInput = time_histos[i];
    for (int i = 1; i <= nBins; ++i) {
      double sum = 0;
      for (int j = 1; j <= i; ++j) {
        double inputVal = hInput->GetBinContent(j);
        double shaperVal = hShaper->GetBinContent(i - j + 1);
        sum += inputVal * shaperVal;
      }
      hOutput->SetBinContent(i, sum * binWidth);  // Scale by bin width
    }
    hOutput->SetName(Form("shapedSignal_%d", i));
    shapedSignals.push_back(hOutput);
    shapedSignals.back()->SetLineColor(colorArray[i]);
    shapedSignals.back()->SetFillColor(colorArray[i]);
    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    c3->cd();
    shapedSignals.back()->Draw();
    // draw a line at threshold point
    TLine *line = new TLine(0, threshold, 60, threshold);
    line->SetLineColor(kBlack);
    line->Draw();
    c3->SaveAs(Form("time_histo_%s_Shaped.png", filenames[i].c_str()));
  }

  std::vector<double> times;
  TGraph *g_timewalk_vs_peak = new TGraph();
  TGraph *g_timewalk_vs_charge = new TGraph();
  for(int i=0; i<shapedSignals.size(); i++) {
    // find first bin above threshold
    double time_threshold = computeTime(shapedSignals[i], threshold);
    double t0 = computeTime(shapedSignals[i], 0.5);
    // Find peak of shaped signal
    double peak = shapedSignals[i]->GetMaximum();
    std::cout<<filenames[i]<<" - Signal over threshold at "<<time_threshold<<" ns; t_0: "<<t0<<" ns"<<std::endl;
    times.push_back(time_threshold);
    double charge = computeCharge(time_histos[i]);
    std::cout<<filenames[i]<<" - Peak at "<<peak<<" charge = "<<charge<<std::endl;
    g_timewalk_vs_peak->SetPoint(i, time_threshold-t0, peak);
    g_timewalk_vs_charge->SetPoint(i, time_threshold-t0, charge);
  }
  TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
  c4->cd();
  g_timewalk_vs_peak->SetMarkerStyle(20);
  g_timewalk_vs_peak->SetMarkerSize(1.5);
  g_timewalk_vs_peak->SetMarkerColor(kBlack);
  g_timewalk_vs_peak->Draw("APl");
  g_timewalk_vs_charge->SetMarkerStyle(20);
  g_timewalk_vs_charge->SetMarkerSize(1.5);
  g_timewalk_vs_charge->SetMarkerColor(kRed);
  g_timewalk_vs_charge->SetLineColor(kRed);
  g_timewalk_vs_peak->SetTitle("; Time walk [ns];charge [a.u.]");
  g_timewalk_vs_charge->Draw("Pl same");
  TLegend *legend2 = new TLegend(0.1, 0.7, 0.3, 0.9);
  legend2->AddEntry(g_timewalk_vs_peak, "Fast shaper", "p");
  legend2->AddEntry(g_timewalk_vs_charge, "Slow shaper", "p");
  legend2->Draw();


  TGraph *g_distance_time = new TGraph();
  TGraph *g = new TGraph();
  // the following is very hardcoded!
  for(int i=0; i<shapedSignals.size(); i++) {
    g_distance_time->SetPoint(i, distances[i], times[i]);
  }

  TCanvas *c3 = new TCanvas("c3", "c3", 1800, 600);
  c3->Divide(2,1);

  g_distance_time->SetMarkerStyle(20);
  g_distance_time->SetMarkerSize(1.5);
  g_distance_time->SetMarkerColor(kRed);
  g_distance_time->SetTitle("ToA; s [mm]; t [ns]");
  c3->cd(1);
  g_distance_time->Draw("AP");
  g_distance_time->Fit("pol1","","");
  double v = 1/g_distance_time->GetFunction("pol1")->GetParameter(1);
  double t0 = g_distance_time->GetFunction("pol1")->GetParameter(0);
  std::cout<<"Speed of light in the fiber is "<<v<<" mm/ns"<<std::endl;
  for(int i=0; i<shapedSignals.size(); i++) {
    double residual = times[i] - distances[i]/v - t0;
    g->SetPoint(i, distances[i], residual);
  }
  c3->cd(2);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.5);
  g->SetMarkerColor(kBlack);
  g->SetTitle("Residual; s [mm]; t_{residual} [ns]");
  g->Draw("AP");

}


TH1D* getTimeHisto(string filename = "outputs/output200.root"){

  TFile *f = new TFile(filename.c_str());
  TTree *tree = (TTree*)f->Get("tree");
  double ipos[3],fpos[3],ivec[3],fvec[3],time,length;
  int imat, fmat, ftype, nref, npas;
  tree->SetBranchAddress("ipos[3]", ipos);
  tree->SetBranchAddress("fpos[3]", fpos);
  tree->SetBranchAddress("ivec[3]", ivec);
  tree->SetBranchAddress("fvec[3]", fvec);
  tree->SetBranchAddress("time", &time);
  tree->SetBranchAddress("length", &length);
  tree->SetBranchAddress("imat", &imat);
  tree->SetBranchAddress("fmat", &fmat);
  tree->SetBranchAddress("ftype", &ftype);
  tree->SetBranchAddress("nref", &nref);
  tree->SetBranchAddress("npas", &npas);

  TH1D* time_histo = new TH1D("time_histo", "photon ToA", 450, 0, 80);

  int nentries = tree->GetEntries();
  for (int i=0; i<nentries; i++) {
    tree->GetEntry(i);
    if(ftype==4) {
      time_histo->Fill(time);
    }
  }
  time_histo->SetName(Form("time_histo_%s", filename.c_str()));

  // draw and save
   TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
   time_histo->SetFillColor(colorArray[I_color_array]);
    I_color_array++;
   time_histo->Draw();
   time_histo->SetTitle("photon ToA; t [ns]; N_{photons}");
   c1->SaveAs(Form("time_histo_%s.png", filename.c_str()));
   time_histo->SetFillColor(0);

  return time_histo;

}



double computeCharge(TH1D* h) {
  // Define shaping function (e.g., CR-RC shaping)
  double tau = 12.50*4;  // Time constant in SFGD is 4 times larger than the TDC time constant
  int nBins = h->GetNbinsX();
  double binWidth = h->GetBinWidth(1);
  TH1D *hShaper = new TH1D("hShaperCharge", "", nBins, 0, nBins * binWidth);
  for (int i = 1; i <= nBins; ++i) {
    double t = hShaper->GetBinCenter(i);
    double value = (t >= 0) ? (t / tau) * exp(-t / tau) : 0;  // CR-RC example
    hShaper->SetBinContent(i, value);
  }

  // Perform discrete convolution
  TH1D *hOutput = new TH1D("hOutputCharge", "Convolved Signal", nBins, 0, nBins * binWidth);
  TH1D* hInput = h;
  for (int i = 1; i <= nBins; ++i) {
    double sum = 0;
    for (int j = 1; j <= i; ++j) {
      double inputVal = hInput->GetBinContent(j);
      double shaperVal = hShaper->GetBinContent(i - j + 1);
      sum += inputVal * shaperVal;
    }
    hOutput->SetBinContent(i, sum * binWidth);  // Scale by bin width
  }

  double charge = hOutput->GetMaximum();

  delete hShaper;
  delete hOutput;
  return charge;
}

double computeTime(TH1D* h, double threshold){
  // find first bin above threshold
  int bin = 0;
  while(h->GetBinContent(bin) <= threshold) {
    bin++;
    if(bin > h->GetNbinsX()) return -1;
  }
  return h->GetBinCenter(bin);
}

TH1D* fastShaperOutput(TH1D* h){
  // Define shaping function (e.g., CR-RC shaping)
  double tau = 15;  // Time constant
  int nBins = h->GetNbinsX();
  double binWidth = h->GetBinWidth(1);
  TH1D *hShaper = new TH1D("hShaper", "Shaping Function", nBins, 0, nBins * binWidth);
  for (int i = 1; i <= nBins; ++i) {
    double t = hShaper->GetBinCenter(i);
    double value = (t >= 0) ? (t / tau) * exp(-t / tau) : 0;  // CR-RC example
    hShaper->SetBinContent(i, value);
  }

  // Perform discrete convolution
  TH1D *hOutput = new TH1D("hOutput", "Convolved Signal (fast shaper)", nBins, 0, nBins * binWidth);
  TH1D* hInput = h;
  for (int i = 1; i <= nBins; ++i) {
    double sum = 0;
    for (int j = 1; j <= i; ++j) {
      double inputVal = hInput->GetBinContent(j);
      double shaperVal = hShaper->GetBinContent(i - j + 1);
      sum += inputVal * shaperVal;
    }
    hOutput->SetBinContent(i, sum * binWidth);  // Scale by bin width
  }

  return hOutput;
}


TH1D* fastShaperOutputLaplace(TH1D* h){
  // Read shaping function from file
  TFile *f = new TFile("db/FastShaperDeltaResponse.root");
  TH1D* hShaper = (TH1D*)f->Get("deltaResponse");

  int nBins = h->GetNbinsX();
  double binWidth = h->GetBinWidth(1);
  // Perform discrete convolution
  TH1D *hOutput = new TH1D("hOutput", "Convolved Signal (fast shaper)", nBins, 0, nBins * binWidth);
  TH1D* hInput = h;
  for (int i = 1; i <= nBins; ++i) {
    double sum = 0;
    for (int j = 1; j <= i; ++j) {
      double inputVal = hInput->GetBinContent(j);
      int binIndex = hShaper->FindBin( hInput->GetBinCenter(i - j + 1) );
      double shaperVal = hShaper->GetBinContent(binIndex);
      sum += inputVal * shaperVal;
    }
    hOutput->SetBinContent(i, sum * binWidth);  // Scale by bin width
  }

  return hOutput;

}


TH1D* slowShaperOutput(TH1D* h){
  // Define shaping function (e.g., CR-RC shaping)
  double tau = 50;  // Time constant
  int nBins = h->GetNbinsX();
  double binWidth = h->GetBinWidth(1);
  TH1D *hShaper = new TH1D("hShaper", "Shaping Function", nBins, 0, nBins * binWidth);
  for (int i = 1; i <= nBins; ++i) {
    double t = hShaper->GetBinCenter(i);
    double value = (t >= 0) ? (t / tau) * exp(-t / tau) : 0;  // CR-RC example
    hShaper->SetBinContent(i, value);
  }

  // Perform discrete convolution
  TH1D *hOutput = new TH1D("hOutput", "Convolved Signal (slow shaper)", nBins, 0, nBins * binWidth);
  TH1D* hInput = h;
  for (int i = 1; i <= nBins; ++i) {
    double sum = 0;
    for (int j = 1; j <= i; ++j) {
      double inputVal = hInput->GetBinContent(j);
      double shaperVal = hShaper->GetBinContent(i - j + 1);
      sum += inputVal * shaperVal;
    }
    hOutput->SetBinContent(i, sum * binWidth);  // Scale by bin width
  }

  return hOutput;
}
