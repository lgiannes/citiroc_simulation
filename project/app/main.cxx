//
// Created by Lorenzo Giannessi on 02.03.2025.
//

#include <iostream>
#include "CITIROC_Channel.hxx"


int main(int argc, char** argv){

  CITIROC_Channel channel1;

  std::cout << "Channel 1 configured (default config)" << std::endl;

  channel1.FastDeltaResponse(0,10e-9,10);


  return 0;
}