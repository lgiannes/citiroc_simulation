//
// Created by Lorenzo Giannessi on 02.03.2025.
//

#include "CITIROC_Channel.hxx"
#include "TransferFunctions.hxx"

CITIROC_Channel::CITIROC_Channel(){
  // CITIROC nominal values
  _FSH_tau = 2.5; // ns
  _FSH_Rf = 25000; // Ohm
  _FSH_Rin = 500; // Ohm; beta = Rf/Rin = 50 determines the gain of the fsh closed loop
  _FSH_A0 = 1000;
  _FSH_w0 = 2*TMath::Pi() * 1.2e-3; // rad/ns FSH op-amp cutoff frequency

  _HG_Cin = 15e-3; // nF
  _HG_w0  = 2*ILT::PI*450e-6; // rad/ns
  _HG_A0  = 1778;

  _LG_Cin = 15e-2; // nF
  _LG_w0 = 2*ILT::PI*450e-6; // rad/ns
  _LG_A0 = 1000;

  // Apply configuration (for now, default)
  DefaultConfigure();
  _HG_Cfb = 25e-6*(63-_HG); // nF
  _LG_Cfb = 25e-6*(63-_LG); // nF

}

void CITIROC_Channel::DefaultConfigure(){
  _HG = 44;
  _LG = 44;
  _TimeConstantHG = 4;
  _TimeConstantLG = 4;
  _DAC10b = 220;
  _DAC10b_t = 220;
  _SCAorPeakDetectorLG = false;
  _SCAorPeakDetectorHG = false;
  _FastShaperOnLG = false;
  std::cout << "Default configuration set" << std::endl;
}

ILT::cmplex CITIROC_Channel::TF_HGpreamp(const ILT::cmplex &s) {
  double w0 = _HG_w0;
  double A0 = _HG_A0;
  double Cin = _HG_Cin;
  double Cfb = _HG_Cfb;
  double beta = Cin/Cfb;

  ILT::cmplex Fs_HG = A0 * beta/( 1+beta+A0+s/w0*(1+beta) ); // closed loop formula
  return  Fs_HG;

}

ILT::cmplex CITIROC_Channel::TF_LGpreamp(const ILT::cmplex &s) {
  double w0 = _LG_w0;
  double A0 = _LG_A0;
  double Cin = _LG_Cin;
  double Cfb = _LG_Cfb;
  double beta = Cin/Cfb;

  ILT::cmplex Fs_LG = A0 * beta/( 1+beta+A0+s/w0*(1+beta) ); // closed loop formula
  return  Fs_LG;

}


ILT::cmplex CITIROC_Channel::TF_FSH(const ILT::cmplex &s) {
  double tau = _FSH_tau;
  double beta = _FSH_Rf / _FSH_Rin;
  double A0 = _FSH_A0;
  double w0 = _FSH_w0;

  ILT::cmplex numerator = beta*tau*s ;
  double zero = 1 ;
  double first = 2*tau + beta*tau/A0 + 1/A0/w0 ;
  double second = tau*tau + beta*tau/A0/w0 + 2*tau/A0/w0 ;
  double third = tau*tau/A0/w0 ;
  return numerator / (zero + first*s + second*s*s + third*s*s*s) ;
}

ILT::cmplex CITIROC_Channel::TF_TimePath(const ILT::cmplex &s) {
  if(_FastShaperOnLG){
    return TF_LGpreamp(s) * TF_FSH(s);
  }else{
    return TF_HGpreamp(s) * TF_FSH(s);
  }
}

TH1D* CITIROC_Channel::FastDeltaResponse(double tmin,double tmax, int steps){

  TH1D* delta_response = new TH1D("delta_response","delta_response",steps,tmin,tmax);

  ////////////////////////// debug ////////////////////////////////////
  if(_debug) {
    std::cout << "Components: " << std::endl;
    std::cout << "HG w0: " << _HG_w0 << std::endl;
    std::cout << "HG A0: " << _HG_A0 << std::endl;
    std::cout << "HG Cin: " << _HG_Cin << std::endl;
    std::cout << "HG Cfb: " << _HG_Cfb << std::endl;
    std::cout << "FSH tau: " << _FSH_tau << std::endl;
    std::cout << "FSH Rf: " << _FSH_Rf << std::endl;
    std::cout << "FSH Rin: " << _FSH_Rin << std::endl;
    std::cout << "FSH A0: " << _FSH_A0 << std::endl;
    std::cout << "FSH w0: " << _FSH_w0 << std::endl;
    // Debug: draw transfer function in Laplace space
    TF1 *f_hg = new TF1("f_hg",
                        [this](double *x, double *p) { return abs(TF_HGpreamp(ILT::cmplex(0, x[0]))); },
                        1, 1e12, 0);
    TF1 *f_fsh = new TF1("f_fsh",
                         [this](double *x, double *p) { return abs(TF_FSH(ILT::cmplex(0, x[0]))); },
                         1, 1e12, 0);

    TF1 *f = new TF1("f",
                     [this](double *x, double *p) { return abs(TF_TimePath(ILT::cmplex(0, x[0]))); },
                     1, 1e12, 0);
    TCanvas *c = new TCanvas("c", "c", 1800, 600);
    c->Divide(3, 1);
    c->cd(1);
    f->Draw();
    gPad->SetLogx();
    gPad->SetLogy();
    c->cd(2);
    f_hg->Draw();
    gPad->SetLogx();
    gPad->SetLogy();
    c->cd(3);
    f_fsh->Draw();
    gPad->SetLogx();
    gPad->SetLogy();
    c->SaveAs("TF_TimePath.png");
    c->SaveAs("TF_TimePath.root");
  }
  ////////////////////////// debug ////////////////////////////////////



  // Get inverse laplace transform
  for(int t_step=0;t_step<steps;t_step++){
    double t = tmin + t_step*(tmax-tmin)/steps;

    std::function<ILT::cmplex(const ILT::cmplex&)> Fs =
//            [this](const ILT::cmplex &s) { return 1.0 / (s+1.); }; // test function
              [this](const ILT::cmplex &s) { return TF_TimePath(s); };

    double ft = ILT::zakianInverseLaplace(Fs, t);
    cout << "[Debug] step: "<<t_step<< " t= " << t << " f(t)= " << ft << endl;
    delta_response->SetBinContent(t_step+1,ft);
  }



  return delta_response;
}


