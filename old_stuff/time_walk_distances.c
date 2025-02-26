
#include "time_walk_fsh_ssh.c"


void time_walk_distances() {

  std::vector<string> filenames = {"outputs/output200.root", "outputs/output400.root", "outputs/output600.root", "outputs/output800.root", "outputs/output1000.root",
                                   "outputs/output1200.root", "outputs/output1400.root", "outputs/output1600.root", "outputs/output1800.root"};


  std::vector<int> distances = {400,600,800,1000,1200,1400,1600,1800};
  TMultiGraph *mg = new TMultiGraph();
  TLegend * legend = new TLegend(0.1, 0.7, 0.3, 0.9);
  int i = 0;
  for(int dist : distances) {
    string filename = Form("outputs/output%d.root", dist);
    TGraph *g_SSH = time_walk_graph(filename, true);

    g_SSH->SetMarkerStyle(20);
    g_SSH->SetMarkerSize(1);
    g_SSH->SetMarkerColor(i+1);
    g_SSH->SetTitle(Form("Time walk, D = %d mm; Charge [a.u.]; Time walk [ns]", 2000 - dist));
    g_SSH->Draw("AP");
    TF1 * f = new
    TF1("f", ftimewalkCorrection2, 0, 240, 5);
    f->SetParameters(1.003, 0.0042, 1.0136, 114.9, 1.463);
    g_SSH->Fit("f");
    g_SSH->GetFunction("f")->SetLineColor(i+1);
    i++;

    legend->AddEntry(g_SSH, Form("Distance = %d mm", 2000-dist), "p");
    mg->Add(g_SSH);
  }

  mg->SetTitle("; Charge [p.e.]; Time walk [ns]");
  mg->Draw("AP");
  legend->Draw();

}

