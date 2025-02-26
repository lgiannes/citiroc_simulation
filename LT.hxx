
namespace LT {


    TH1D* getTimeHisto(string filename){

      TFile *f = new TFile(filename.c_str());
      TTree *tree = (TTree*)f->Get("tree");
      double ipos[3],fpos[3],ivec[3],fvec[3],time,length;
      int imat, fmat, ftype, nref, npas;
      tree->SetBranchAddress("time", &time);
      tree->SetBranchAddress("ftype", &ftype);

      TH1D* time_histo = new TH1D("time_histo", "photon ToA", 240, 0, 60);

      int nentries = tree->GetEntries();
      for (int i=0; i<nentries; i++) {
        tree->GetEntry(i);
        if(ftype==4) {
          time_histo->Fill(time);
        }
      }
      time_histo->SetName(Form("time_histo_%s", filename.c_str()));

      return time_histo;

    }


    TH1D *fastShaperOutputLaplace(TH1D *h) {
      // Read shaping function from file
      TFile *f = new TFile("db/FastShaperDeltaResponse.root");
      TH1D *hShaper = (TH1D *) f->Get("deltaResponse");

      int nBins = h->GetNbinsX();
      double binWidth = h->GetBinWidth(1);
      // Perform discrete convolution
      TH1D *hOutput = new TH1D("hOutput", "Convolved Signal (fast shaper)", nBins, 0, nBins * binWidth);
      TH1D *hInput = h;
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


}

