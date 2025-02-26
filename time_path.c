#include "LT.hxx"

TH1D* time_path(TH1D* input_pulse, int HGDAC){
  // 1. load the transfer function for the desired HG DAC value
  string TF_name = "db/FastShaperDeltaResponse_HG" + to_string(HGDAC) + ".root";
  TFile * TF_file = new TFile(TF_name.c_str());
  TH1D* hShaper = (TH1D *) TF_file->Get("deltaResponse");

  // 2. Perform discrete convolution
  int nBins = input_pulse->GetNbinsX();
  double binWidth = input_pulse->GetBinWidth(1);
  TH1D *hOutput = new TH1D("hOutput", "Convolved Signal (fast shaper)", nBins, 0, nBins * binWidth);
  TH1D *hInput = input_pulse;
  for (int i = 1; i <= nBins; ++i) {
    double sum = 0;
    for (int j = 1; j <= i; ++j) {
      double inputVal = hInput->GetBinContent(j);
      int binIndex = hShaper->FindBin(hInput->GetBinCenter(i - j + 1));
      double shaperVal = hShaper->GetBinContent(binIndex);
      sum += inputVal * shaperVal;
    }
    hOutput->SetBinContent(i, sum * binWidth);  // Scale by bin width
  }

  return hOutput;

}

TH1D* time_path(string input_filename){
  int HGDAC = 44;

  TH1D* input_pulse = LT::getTimeHisto(input_filename);
  TCanvas * c = new TCanvas("c", "c", 1800, 600);
  c->Divide(2,1);
  c->cd(1);
  input_pulse->Draw();
  input_pulse->SetLineColor(kBlack);
  input_pulse->SetLineWidth(3);
  input_pulse->SetTitle("Input pulse;Time [ns];Amplitude [a.u.]");
  TH1D* output_pulse = time_path(input_pulse, HGDAC);
  c->cd(2);
  output_pulse->Draw();
  string title = "Output pulse (HG DAC = " + to_string(HGDAC) + ");Time [ns];Amplitude [a.u.]";
  output_pulse->SetTitle(title.c_str());
  output_pulse->SetStats(0);
  output_pulse->SetLineWidth(3);
  output_pulse->SetLineColor(kBlue);

  return output_pulse;

}