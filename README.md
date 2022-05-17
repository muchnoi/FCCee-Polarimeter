# FCCee-Polarimeter
MC simulations for the project of FCC-ee Compton beam polarimeter.

The method is discussed in the publication https://arxiv.org/abs/1803.09595.

Tested with python 3.8.10 [GCC 9.4.0] and ROOT 6.22/03 (https://root.cern.ch/) on Linux Mint kernel 5.4.0-109-generic.
The packages numpy and scipy are required (recent versions).

1) Simulation parameters may be altered in the <b>hardware.py</b>. 
2) Use <b>mc-compton.py</b> for Monte-Carlo generation of photon and electron distributions.
3) Simulation parameters and resulting histograms will be saved to the file with a name, similar to <b>Tue_May_10_130538_2022.MC</b>.
4) To look at the data and fit functions with initial parameters use command <pre>allfits.py filename</pre>.
5) For fitting use command <pre>allfits.py --fit filename</pre>.
5) Bear in mind that the fitting procedure is quite long at this stage.



