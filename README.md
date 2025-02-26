# citiroc_simulation
Analytical simulation of signal processing in the CITIROC read-out chip by Omega. For physicists

## THIS IS A WORK IN PROGRESS!!!


The `Bode.c` macro produces Bode diagrams (gain as a function of the input frequency) for different HG DAC values, for the Fast shaper part of the CITIROC.

### Usage

the main function so far is `time_path.c`. It's a simple root macro that takes as input a discrete waveform in the form of a 
root TH1D (current as a function of time) and produces the time-domain response of the system, for the Fast shaper part. 
It needs the HG_DAC value as argument as well. 

run with:
`root "time_path.c(<waveform_filename> <HG_DAC value>)"`

