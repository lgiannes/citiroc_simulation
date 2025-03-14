//
// Created by Lorenzo Giannessi on 09.01.2025.
//

#ifndef CITIROC_SIM_TRANSFERFUNCTIONS_HXX
#define CITIROC_SIM_TRANSFERFUNCTIONS_HXX

namespace TransferFunctions{

    double TF_FSH_actual(double *x,double *p){
      //  Transfer function (Laplace space) of FSH op-amp, real model
      double s = x[0];

      double tau = 2.5e-9; // s
      double beta = 50; // R_F/R_in
      double A0 = 1000; // FSH op-amp gain (60 dB)
      double w0 = 2*TMath::Pi()*1.2e6; // rad/s FSH op-amp cutoff frequency

      double numerator = beta*tau*s ;
      double zero = 1 ;
      double first = 2*tau + beta*tau/A0 + 1/A0/w0 ;
      double second = tau*tau + beta*tau/A0/w0 + 2*tau/A0/w0 ;
      double third = tau*tau/A0/w0 ;

      return numerator / (zero + first*s + second*s*s + third*s*s*s) ;
    }

    double TF_FSH_firstapprox(double *x,double *p){
      //  Transfer function (Laplace space) of FSH op-amp, "first" approximation
      double s = x[0];

      double tau = 2.5e-9; // s
      double beta = 50; // R_F/R_in
      double A0 = 1000; // FSH op-amp gain (60 dB)
      double w0 = 2*TMath::Pi()*1.2e6; // rad/s FSH op-amp cutoff frequency

      double numerator = beta*tau*s ;

      return numerator /( (1 + 4*tau*s +  4*tau*tau*s*s)*(1+s/(4*A0*w0)) );
    }

    double TF_FSH_secondapprox(double *x,double *p){
      //  Transfer function (Laplace space) of FSH op-amp, "second" approximation
      double s = x[0];

      double tau = 2.5e-9; // s
      double beta = 50; // R_F/R_in
      double A0 = 1000; // FSH op-amp gain (60 dB)
      double w0 = 2*TMath::Pi()*1.2e6; // rad/s FSH op-amp cutoff frequency

      double numerator = beta*tau*s ;

      return numerator /( (1 + 2*tau*s +  4*tau*tau*s*s)*(1+s/(4*A0*w0)) );
    }

    double HGcl(double *x,double *p){
      // Transfer function (Laplace space) of HG pre-amp closed loop circuit
      double s = x[0];
      double w0 = 2*TMath::Pi()*450e3; // rad/s HG op-amp cutoff frequency
      double A0 = 1778; // HG op-amp gain (65 dB)
      double DAC_HG = p[0];
      double Cin = 15e-12; // 15 pF
      double Cfb = 25e-15*(63-DAC_HG);
      double beta = Cin/Cfb;

      double Fs_HG = A0 * beta/( 1+beta+A0+s/w0*(1+beta) ); // HG preamp
      return  Fs_HG;
    }

    double TF_time_path(double *x,double *p){
      // Transfer function (Laplace space) of full time-stamping circuit: HG preamp + Fast Shaper
      double s = x[0];

      return HGcl(x,p) * TF_FSH_actual(x,p);
    }

};

#endif //CITIROC_SIM_TRANSFERFUNCTIONS_HXX
