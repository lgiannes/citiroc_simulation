
import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt
import ROOT
from functions import make_response_function
import sys

HGDAC = int(sys.argv[1])
if(len(sys.argv) != 2):
    print("Please provide HGDAC value as argument")
if(HGDAC < 0 or HGDAC > 62):
    print("HGDAC value must be between 0 and 63")
    sys.exit()

Fs = make_response_function(HGDAC)


x = np.linspace(0.1e-9, 100e-9, 999)
y_deltaResponse = [(mp.invertlaplace(Fs,val,method='talbot')) for val in x]
# make TH1D (root) from the delta response
h_deltaResponse = ROOT.TH1D("deltaResponse", "deltaResponse", len(x), 1e9*x[0], 1e9*x[-1])
h_deltaResponse.GetXaxis().SetTitle("Time (ns)")
h_deltaResponse.GetYaxis().SetTitle("V_out")
h_deltaResponse.Title = "Response of the HG+FSH circuit to a Dirac delta pulse"
for i in range(len(x)):
    h_deltaResponse.SetBinContent(i+1, y_deltaResponse[i])

h_deltaResponse.Draw()
h_deltaResponse.SaveAs("FastShaperDeltaResponse_HG"+str(HGDAC)+".root")











def draw_function(Fs):
    # Draw function in laplace space
    # log-log plot
    x = np.logspace(1, 12, 10000)  # From 10^-2 to 10^2
    y = [(Fs(val)) for val in x]
    # Convert y to float for matplotlib plotting
    y = np.array(y, dtype=float)
    # Plot with matplotlib on log-log scale
    plt.figure(figsize=(8, 6))
    plt.loglog(x, y, color="blue")
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('V_out/V_in')
    plt.title('Bode plot')
    plt.grid()
    # Show the plot
    plt.show()

    step_response = lambda s: Fs(s)*1/s

    x = np.linspace(0.1e-9, 100e-9, 999)
    y = [(mp.invertlaplace(step_response,val,method='talbot')) for val in x]
    y_deltaResponse = [(mp.invertlaplace(Fs,val,method='talbot')) for val in x]
    # make TH1D (root) from the delta response
    h_deltaResponse = ROOT.TH1D("deltaResponse", "deltaResponse", len(x), 1e9*x[0], 1e9*x[-1])
    h_deltaResponse.GetXaxis().SetTitle("Time (ns)")
    h_deltaResponse.GetYaxis().SetTitle("V_out")
    h_deltaResponse.Title = "Response of the HG+FSH circuit to a Dirac delta pulse"
    for i in range(len(x)):
        h_deltaResponse.SetBinContent(i+1, y_deltaResponse[i])

    h_deltaResponse.SaveAs("FastShaperDeltaResponse.root")

    plt.figure(figsize=(8, 6))
    plt.plot(x, y, color="blue")
    plt.xlabel('Time (s)')
    plt.ylabel('V_out')
    plt.title('Step response')
    plt.grid()
    plt.show()
    # find maximum of the step response
    max_value = max(y)
    max_time = x[y.index(max_value)]
    print("Maximum value of step response: ", max_value)
    print("Time at which maximum value occurs: ", max_time)