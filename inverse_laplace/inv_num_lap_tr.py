
import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt
from functions import make_response_function

mp.dps = 15
mp.pretty = True



# tau_CRRC = 2.5e-9
# CRRC = lambda s: tau_CRRC*s/pow((1 + tau_CRRC*s),3)


# for a range of DAC_HG values, plot the step responses in the same plot
Fs_tot_list = []
HGDAC_list = [10,20,28,36,44,52,60]
for HGDAC in HGDAC_list:
    Fs_tot = make_response_function(HGDAC)
    Fs_tot_list.append(Fs_tot)
    # draw_function(Fs_tot)

# draw all the step responses in the same plot
x = np.linspace(0.1e-9, 100e-9, 999)
plt.figure(figsize=(8, 6))
for Fs_tot,HGDAC in zip(Fs_tot_list,HGDAC_list):
    step_response = lambda s: Fs_tot(s)*1/s
    y = [(mp.invertlaplace(step_response,val,method='talbot')) for val in x]
    plt.plot(x, y, label='HGDAC = ' + str(HGDAC))
    max_value = max(y)
    max_time = x[y.index(max_value)]
    print("Maximum value of step response: ", max_value)
    print("Time at which maximum value occurs: ", max_time)

plt.xlabel('Time (s)')
plt.ylabel('V_out')
plt.title('Step response')
plt.grid()
plt.legend()
plt.show()

