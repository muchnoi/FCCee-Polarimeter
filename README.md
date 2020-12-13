# FCCee-Polarimeter
MC simulations for the project of FCC-ee Compton beam polarimeter.

The method is discussed in the publication https://arxiv.org/abs/1803.09595.

Tested with python 3.85 [GCC 9.3.0] and ROOT 6.22/03 on Linux Mint kernel 5.4.0-58-generic

1) Simulation parameters may be <pre>altered</pre> in the <b>hardware.py</b>. 
2) Use <b>mc-compton.py</b> for Monte-Carlo generation of photon and electron distributions.
3) The result will be saved to a root file with a name, similar to <b>Sun Dec 13 10:49:17 2020.root</b>.
4) For fitting, alter this filename in <b>eFit.py</b> (somewhere around 140 line) and run the script.
5) Bear in mind that the fitting procedure is quite long at this stage. On my Intel(R) Core(TM) i3-6100U CPU @ 2.30GHz laptop the fit (with default parameter set) converges after iterations, which takes about real time. If this method is ever used, the fitting procedure should definitely be improved.


