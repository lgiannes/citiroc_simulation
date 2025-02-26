#include "TransferFunctions.hxx"

void Bode(){

  double range_low = 1;
  double range_high = 1e12;

  // Fast shaper Bode plot and its approximation (all parameters are fixed for the CITIROC fast shaper)
  TF1 *f_actual = new TF1("f_actual", TransferFunctions::TF_FSH_actual, range_low, range_high, 0);
  f_actual->SetNpx(1000);
  f_actual->SetLineColor(kRed);
  f_actual->SetLineWidth(3);
  f_actual->SetTitle("Bode plot of the FSH actual transfer function");
  f_actual->GetXaxis()->SetTitle("Frequency [Hz]");
  f_actual->GetYaxis()->SetTitle("Gain");

  TF1 *f_firstapprox = new TF1("f_firstapprox", TransferFunctions::TF_FSH_firstapprox, range_low, range_high, 0);
  f_firstapprox->SetNpx(1000);
  f_firstapprox->SetLineColor(kBlue);
  f_firstapprox->SetLineWidth(3);
  f_firstapprox->SetTitle("Bode plot of the FSH first approximation transfer function");
  f_firstapprox->GetXaxis()->SetTitle("Frequency [Hz]");
  f_firstapprox->GetYaxis()->SetTitle("Gain");

  TF1 *f_secondapprox = new TF1("f_secondapprox", TransferFunctions::TF_FSH_secondapprox, range_low, range_high, 0);
  f_secondapprox->SetNpx(1000);
  f_secondapprox->SetLineColor(kGreen);
  f_secondapprox->SetLineWidth(3);
  f_secondapprox->SetLineStyle(2);
  f_secondapprox->SetTitle("Bode plot of the FSH second approximation transfer function");
  f_secondapprox->GetXaxis()->SetTitle("Frequency [Hz]");
  f_secondapprox->GetYaxis()->SetTitle("Gain");


  TCanvas *c = new TCanvas("c", "c", 800, 600);
  TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
  c->cd();
  f_actual->Draw();
  legend->AddEntry(f_actual, "Actual transfer function", "l");
//  f_firstapprox->Draw("same");
//  legend->AddEntry(f_firstapprox, "First approximation", "l");
  f_secondapprox->Draw("same");
  legend->AddEntry(f_secondapprox, "Approximation", "l");

  legend->Draw();
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid();

  // HG preamp Bode plot for different DAC values
  vector<double> DACs = {10, 20, 28, 36, 44, 52, 60};
  vector<TF1*> f_HGcls;
  TLegend *legend_HG = new TLegend(0.2, 0.4, 0.4, 0.6);

  TCanvas *c_HG = new TCanvas("c_HG", "c_HG", 800, 600);
  c_HG->cd();
  int color=1;
  for(double DAC : DACs){
    TF1 *f_HGcl = new TF1("f_HGcl", TransferFunctions::HGcl, range_low, range_high, 1);
    f_HGcls.push_back(f_HGcl);
    f_HGcls.back()->SetParameter(0,DAC);
    f_HGcls.back()->SetNpx(1000);
    f_HGcls.back()->SetLineColor(color);
    f_HGcls.back()->SetLineWidth(3);
    f_HGcls.back()->SetTitle(Form("Bode plot of the HG op-amp"));
    f_HGcls.back()->GetXaxis()->SetTitle("Frequency [Hz]");
    f_HGcls.back()->GetYaxis()->SetTitle("Gain");
    if(color==1) {
      f_HGcls.back()->Draw();
    }else {
      f_HGcls.back()->Draw("same");
    }
    legend_HG->AddEntry(f_HGcls.back(), Form("HG DAC = %d", int(DAC)), "l");
    gPad->SetLogx();
    gPad->SetLogy();
    color++;
    if(color==5) color++;
  }
  legend_HG->Draw();

  // repeat, for HG preamp + FSH
  vector<TF1*> f_HGcl_FSHs;
  TLegend *legend_HG_FSH = new TLegend(0.2, 0.4, 0.4, 0.6);
  vector<double> gains;
  vector<double> peak_frequencies;

  TCanvas *c_HG_FSH = new TCanvas("c_HG_FSH", "c_HG_FSH", 800, 600);
  c_HG_FSH->cd();
  color=1;
  for(double DAC : DACs){
    TF1 *f_HGcl_FSH = new TF1("f_HGcl_FSH", TransferFunctions::TF_time_path, range_low, range_high, 1);
    f_HGcl_FSHs.push_back(f_HGcl_FSH);
    f_HGcl_FSHs.back()->SetParameter(0,DAC);
    f_HGcl_FSHs.back()->SetNpx(1000);
    f_HGcl_FSHs.back()->SetLineColor(color);
    f_HGcl_FSHs.back()->SetLineWidth(3);
    f_HGcl_FSHs.back()->SetTitle(Form("Bode plot of the HG op-amp + FSH"));
    f_HGcl_FSHs.back()->GetXaxis()->SetTitle("Frequency [Hz]");
    f_HGcl_FSHs.back()->GetYaxis()->SetTitle("Gain");
    if(color==1) {
      f_HGcl_FSHs.back()->Draw();
    }else {
      f_HGcl_FSHs.back()->Draw("same");
    }
    double peak_frequency = f_HGcl_FSHs.back()->GetMaximumX();
    double gain = f_HGcl_FSHs.back()->Eval(peak_frequency);
    gains.push_back(gain);
    cout<<"Gain for DAC = "<<DAC<<" is "<<gain<<endl;
    peak_frequencies.push_back(peak_frequency);
    legend_HG_FSH->AddEntry(f_HGcl_FSHs.back(), Form("HG DAC = %d", int(DAC)), "l");
    gPad->SetLogx();
    gPad->SetLogy();
    color++;
    if(color==5) color++;
  }
  legend_HG_FSH->Draw();

  TCanvas *c_gains = new TCanvas("c_gains", "c_gains", 800, 600);
  c_gains->cd();
  TMultiGraph *mg = new TMultiGraph();
  TGraph *g_gains = new TGraph();
  TGraph *g_peak_frequencies = new TGraph();
  for(int i=0; i<DACs.size(); i++){
    g_gains->SetPoint(g_gains->GetN(), DACs[i], gains[i]);
    g_peak_frequencies->SetPoint(g_peak_frequencies->GetN(), DACs[i], peak_frequencies[i]*1e-6);
  }
  g_gains->SetLineWidth(4);
  g_gains->SetLineColor(kRed);
  g_gains->SetTitle("Gain vs. DAC");
  g_gains->GetXaxis()->SetTitle("DAC");
  g_gains->GetYaxis()->SetTitle("Gain");
  mg->Add(g_gains);
  g_peak_frequencies->SetLineWidth(4);
  g_peak_frequencies->SetLineColor(kBlue);
  g_peak_frequencies->SetTitle("Peak frequency vs. DAC");
  g_peak_frequencies->GetXaxis()->SetTitle("DAC");
  g_peak_frequencies->GetYaxis()->SetTitle("Peak frequency [MHz]");
  mg->Add(g_peak_frequencies);
  mg->GetXaxis()->SetTitle("DAC");
  mg->Draw("AL");
  TLegend *legend_gains = new TLegend(0.1, 0.7, 0.3, 0.9);
  legend_gains->AddEntry(g_gains, "Gain", "l");
  legend_gains->AddEntry(g_peak_frequencies, "Peak frequency [MHz]", "l");
  legend_gains->Draw();

  TCanvas *c_gains_peaktimes= new TCanvas("c_gains_peaktimes", "c_gains_peaktimes", 800, 600);
  TPad* pad1 = new TPad("pad1", "pad1", 0., 0., 1., 1.);
  pad1->cd();
  TGraph* g_peaktimes = new TGraph();
  std::vector<double> peaktimes = {8.30821643286573e-09, 8.708617234468938e-09, 9.209118236472945e-09, 9.809719438877755e-09, 1.0610521042084168e-08, 1.201192384769539e-08, 1.4914829659318637e-08};
  for(int i=0; i<DACs.size(); i++){
    g_peaktimes->SetPoint(g_peaktimes->GetN(), DACs[i], 1e9*peaktimes[i]);
  }
  g_peaktimes->SetLineWidth(4);
  g_peaktimes->SetLineColor(kBlue);
  g_peaktimes->GetYaxis()->SetAxisColor(kBlue);
  g_peaktimes->GetYaxis()->SetTitleColor(kBlue);
  g_peaktimes->GetYaxis()->SetLabelColor(kBlue);
  g_peaktimes->GetXaxis()->SetTitle("DAC");
  g_peaktimes->GetYaxis()->SetTitle("Peaking time [ns]");
  g_peaktimes->GetYaxis()->SetTitleOffset(1.);
  g_peaktimes->Draw("ALY+");
  g_gains->GetYaxis()->SetAxisColor(kRed);
  g_gains->GetYaxis()->SetTitleColor(kRed);
  g_gains->GetYaxis()->SetLabelColor(kRed);
  pad1->DrawClone();
  // new pad
  TPad *pad = new TPad("pad", "pad", 0., 0., 1., 1.);
  pad->cd();
  g_gains->Draw("AL");
  g_gains->SetTitle("");
  g_gains->SetFillStyle(4100);
  pad->SetFillStyle(4100);
  pad->SetFillColorAlpha(kWhite, 0);
  pad->DrawClone();
  pad->Update();
//  c_gains_peaktimes->Update();
  TLegend *legend_peaktimes = new TLegend(0.1, 0.7, 0.3, 0.9);
  legend_peaktimes->AddEntry(g_peaktimes, "Peaking time", "l");
  legend_peaktimes->AddEntry(g_gains, "Gain", "l");
  legend_peaktimes->Draw();

}
