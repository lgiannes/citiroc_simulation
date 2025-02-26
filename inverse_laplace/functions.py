
import numpy as np
def make_response_function(HGDAC):
    # HG preamp
    A_HG = pow(10,3.25)
    w_HG = 2*np.pi*450e3
    beta = 15e3 /( 25*(63-HGDAC) )
    # fast shaper
    A_FSH = pow(10,3)
    w_FSH = 2*np.pi*1.2e6
    ratio = 50 #R_F/R_in in the FSH loop
    tau = 2.5e-9 #R_F*C_F=R_in*C_in in the FSH loop

    Fs_HG = lambda s: A_HG * beta / ( A_HG  + (1+s/w_HG)*(1+beta) )
    Fs_FSH = lambda s: A_FSH*ratio*tau*s/( A_FSH*(tau*s+1)**2 + (1 + (ratio+2)*tau*s + tau*tau*s*s)*(1+s/w_FSH) )
    #multiply lambdas
    Fs_tot = lambda s: Fs_HG(s)*Fs_FSH(s)
    return Fs_tot
