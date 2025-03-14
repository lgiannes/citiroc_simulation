//
// Created by Lorenzo Giannessi on 02.03.2025.
//

#include "CITIROC_Channel.hxx"
#include "TransferFunctions.hxx"

CITIROC_Channel::CITIROC_Channel(){
  // CITIROC nominal values
  _FSH_tau = 2.5e-9; // s
  _FSH_Rf = 25000; // Ohm
  _FSH_Rin = 500; // Ohm; beta = Rf/Rin = 50 determines the gain of the fsh closed loop
  _FSH_A0 = 1000;
  _FSH_w0 = 2*TMath::Pi() * 1.2e6; // rad/s FSH op-amp cutoff frequency

  _HG_Cin = 15e-12; // F
  _HG_w0  = 2*ILT::PI*450e3; // rad/s
  _HG_A0  = 1778;
  _HG_Cfb = 25e-15*(63-_HG); // F

  _LG_Cin = 15e-11; // TODO: check these values
  _LG_w0 = 2*ILT::PI*450e3;
  _LG_A0 = 1000;
  _LG_Cfb = 25e-15*(63-_LG);

  DefaultConfigure();
}

void CITIROC_Channel::DefaultConfigure(){
  _HG = 48;
  _LG = 48;
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

  ILT::cmplex Fs_HG = A0 * beta/( 1+beta+A0+s/w0*(1+beta) ); // HG preamp
  return  Fs_HG;

}

ILT::cmplex CITIROC_Channel::TF_LGpreamp(const ILT::cmplex &s) {
  double w0 = _LG_w0;
  double A0 = _LG_A0;
  double Cin = _LG_Cin;
  double Cfb = _LG_Cfb;
  double beta = Cin/Cfb;

  ILT::cmplex Fs_LG = A0 * beta/( 1+beta+A0+s/w0*(1+beta) ); // HG preamp
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



  // Get inverse laplace transform
  for(int t_step=0;t_step<steps;t_step++){
    double t = tmin + t_step*(tmax-tmin)/steps;

    std::function<ILT::cmplex(const ILT::cmplex&)> Fs =
            [this](const ILT::cmplex &s) { return TF_TimePath(s); };

    double ft = ILT::inverseLaplace(Fs, t, steps);
    cout << "[Debug] step: "<<t_step<< " t= " << t*1e9 << " ns: " << ft << endl;
  }

  return nullptr;
}


