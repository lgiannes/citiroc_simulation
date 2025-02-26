
## inverse laplace transform with Python

This small subpackage is used to perform inverse numerical Laplace transform of the Bode gain functions. 
It produces a TH1D (easily processed with ROOT libraries) containing the time domain delta response of the system.

This uses the `functions.py` module , that contains the response functions in teh Laplace space.  
For a certain HG_DAC value, the response in the Laplace space is computed and the inverse Laplace transform is performed.

### Usage

The main function is `make_time_domain_response.py` which takes the HG_DAC value as an argument and produces the predicted
delta-response of the FAST SHAPER circuit of the CITIROC.

`python3 make_time_domain_response.py <HG_DAC value>`

### Work in progress

So far, only the fast shaper is implemented. Slow shaper coming soon.