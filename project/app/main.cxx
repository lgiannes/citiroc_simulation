//
// Created by Lorenzo Giannessi on 02.03.2025.
//

#include <iostream>
#include "CITIROC_Channel.hxx"


int main(int argc, char** argv){

  CITIROC_Channel channel1;

  std::cout << "Channel 1 configured (default config)" << std::endl;

  TH1D* delta_response = channel1.FastDeltaResponse(0.1,100,1000); // in ns

  if(!delta_response){
    std::cerr << "Error in FastDeltaResponse" << std::endl;
    return 1;
  }

  TCanvas* c1 = new TCanvas("c1","c1",800,600);
  delta_response->Draw();
//  gPad->SetLogy();
  TF1* f = new TF1("f","[0]*x*exp(-[1]*x)",0.1,100);
  f->SetParameter(0,1);
  f->SetParameter(1,0.2);
  f->Draw("same");
  c1->SaveAs("delta_response_timepath.png");
  c1->SaveAs("delta_response_timepath.root");

  for(int iBin = 1; iBin <= delta_response->GetNbinsX(); ++iBin){
    std::cout << "Bin " << iBin << " "<<delta_response->GetBinContent(iBin)-f->Eval(delta_response->GetBinCenter(iBin)) << std::endl;
  }


  return 0;
}