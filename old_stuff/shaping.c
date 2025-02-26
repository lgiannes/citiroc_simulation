#include "time_walk_simulation.c"

void shaping(string filename = "outputs/output200.root"){

  int A = stoi(filename.substr(14,4));

  TH1D* time_histo = getTimeHisto(filename, 1);

  time_histo->Scale(1./time_histo->GetBinWidth(1));
  time_histo->Draw("hist");

  TH1D* FS = fastShaperOutputLaplace(time_histo);

  TCanvas *c = new TCanvas("c", "c", 1600, 450);
  c->Divide(2,1);
  c->cd(1);
  time_histo->SetTitle(Form("Current pulse. Distance: %d mm;Time [ns];i [a.u.]",2000-A));
  time_histo->Draw("hist");
  c->cd(2);
  FS->SetTitle("Fast shaper output;t [ns]; Amplitude [a.u.]");
  FS->Draw("hist");

}
