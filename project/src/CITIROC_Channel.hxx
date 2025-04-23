//
// Created by Lorenzo Giannessi on 02.03.2025.
//

#ifndef CITIROC_SIMULATION_CITIROC_CHANNEL_HXX
#define CITIROC_SIMULATION_CITIROC_CHANNEL_HXX


// ROOT includes
#include <TH1D.h>
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>

// others
#include "ILT.hxx"



class CITIROC_Channel {
private:
    bool _debug{true};
    // Configuration
    bool _isConfigured{false};
    int _HG, _LG; // Gain (6-bit: 0-63)
    int _TimeConstantHG; // Shaping time constant for HG shaper (3-bit: 0-7)
    int _TimeConstantLG; // Shaping time constant for LG shaper (3-bit: 0-7)
    int _DAC10b, _DAC10b_t; // 10-bit: 0-1023
    bool _SCAorPeakDetectorLG; // Switch between Peak Detecor and SCA for charge measurement (0: PD, 1: SCA)
    bool _SCAorPeakDetectorHG; // Switch between Peak Detecor and SCA for charge measurement (0: PD, 1: SCA)
    bool _FastShaperOnLG; // Switch preamplifiers to FSH circuit (0: HG, 1: LG)

    // Components
    // HG Preamplifier
    double _HG_w0;
    double _HG_A0;
    double _HG_Cin;
    double _HG_Cfb; // F
    // LG Preamplifier
    double _LG_w0;
    double _LG_A0;
    double _LG_Cin;
    double _LG_Cfb; // F


    // Fast Shaper
    double _FSH_tau; // seconds
    double _FSH_Rf; // Ohm
    double _FSH_Rin; // Ohm; beta = Rf/Rin = 50 determines the gain of the fsh closed loop
    double _FSH_A0; // FSH op-amp gain (60 dB)
    double _FSH_w0; // rad/s FSH op-amp cutoff frequency
    // Slow shaper

public:

    CITIROC_Channel();
//    ~CITIROC_Channel();

    void DefaultConfigure();
    void Configure(std::string configFile){DefaultConfigure();} // To be implemented

    // Transfer functions
    ILT::cmplex TF_FSH(const ILT::cmplex &s);
    ILT::cmplex TF_HGpreamp(const ILT::cmplex &s);
    ILT::cmplex TF_LGpreamp(const ILT::cmplex &s);
    ILT::cmplex TF_TimePath(const ILT::cmplex &s);

    TH1D* FastDeltaResponse(double tmin,double tmax, int steps);
    TH1D* SlowDeltaResponse(double tmin,double tmax, int steps);
    TH1D* FastStepResponse(double tmin,double tmax, int steps);
    TH1D* SlowStepResponse(double tmin,double tmax, int steps);
    TH1D* FastPulseResponse(TH1D* inputPulse); // Response of the time-measuring (fast shaper) part of the citiroc circuit
    TH1D* SlowPulseResponse(TH1D* inputPulse); // Response of the charge-measuring (slow shaper) part of the citiroc circuit



};


#endif //CITIROC_SIMULATION_CITIROC_CHANNEL_HXX
