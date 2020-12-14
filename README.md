# FCCee-Polarimeter
MC simulations for the project of FCC-ee Compton beam polarimeter.

The method is discussed in the publication https://arxiv.org/abs/1803.09595.

Tested with python 3.85 [GCC 9.3.0] and ROOT 6.22/03 (https://root.cern.ch/) on Linux Mint kernel 5.4.0-58-generic

1) Simulation parameters may be altered in the <b>hardware.py</b>. 
2) Use <b>mc-compton.py</b> for Monte-Carlo generation of photon and electron distributions.
3) The result will be saved to a root file with a name, similar to <b>Sun Dec 13 10:49:17 2020.root</b>.
4) For fitting, alter the filename in <b>eFit.py</b> (somewhere around 140 line) and run the script.
5) Bear in mind that the fitting procedure is quite long at this stage.<br> On my Intel(R) Core(TM) i3-6100U CPU @ 2.30GHz laptop the fit (with default parameter set) converges after 197 iterations, which takes 
<pre>Real Time = 4115.06 seconds Cpu Time = 4122.41 seconds.</pre>

If this method is ever used, the fitting procedure should definitely be improved.

The fit result with default settings is:
<pre>
 FCN=6447.49 FROM MIGRAD    STATUS=CONVERGED     208 CALLS         209 TOTAL
                     EDM=2.08144e-09    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   0.8 per cent
  EXT PARAMETER                  PARABOLIC         MINOS ERRORS
  NO.   NAME      VALUE            ERROR      NEGATIVE      POSITIVE
   1  X_{0}        1.15679e-02   1.00231e-02
   2  X_{1}        3.47635e+02   3.26152e-03
   3  Y_{0}       -1.06824e+00   4.15379e-05
   4  Y_{1}        1.06820e+00   4.07565e-05
   5  #kappa       1.62800e+00     fixed
   6  P_{#parallel}   0.00000e+00     fixed
   7  P_{#perp}    1.02652e-01   1.65003e-03
   8  #sigma_{x}   1.96462e-01   4.54456e-03
   9  #sigma_{y}   2.36739e-02   2.45517e-05
  10  norm         5.22868e+07   1.41892e+04
Chi2: 6447.487180
NDF: 6119
Probability: 0.001729
A = 347.6238 +/- 0.0105 (0.00003)
1.627839909428406 0.00012432482133096817
E = 45.59708 +/- 0.00348
</pre>

