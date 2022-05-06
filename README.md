# FCCee-Polarimeter
MC simulations for the project of FCC-ee Compton beam polarimeter.

The method is discussed in the publication https://arxiv.org/abs/1803.09595.

Tested with python 3.85 [GCC 9.3.0] and ROOT 6.22/03 (https://root.cern.ch/) on Linux Mint kernel 5.4.0-104-generic

1) Simulation parameters may be altered in the <b>hardware.py</b>. 
2) Use <b>mc-compton.py</b> for Monte-Carlo generation of photon and electron distributions.
3) The result will be saved to a root file with a name, similar to <b>Sun Dec 13 10:49:17 2020.root</b>.
4) To look at the data and functions with initial parameters use command <pre>allfits.py filename</pre>.
5) For fitting use command <pre>allfits.py --fit filename</pre>.
5) Bear in mind that the fitting procedure is quite long at this stage.



