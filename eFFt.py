#!/usr/bin/env python3
import sys, ROOT
import numpy as np; π = np.pi
from scipy             import fft                 as F_F_T
from scipy.interpolate import RectBivariateSpline as R_B_S
from hardware          import Constants, Pixels, Laser

class ConvolutionFFT:                                       # Integration with convolution theorem
  Rx, Ry  = 1.1, 1.1                                        # ratio of calculation range to unit radius
  Nx, Ny  = 4096, 4096                                      # number of x,y divisions
  Dx, Dy  = 2*Rx/Nx, 2*Ry/Ny                                # distance between neighbour knots
  x    = np.linspace(Dx/2-Rx, Rx-Dx/2, num=Nx)              # knots positions in x
  y    = np.linspace(Dy/2-Ry, Ry-Dy/2, num=Ny)              # knots positions in y
  u    = F_F_T.fftshift(F_F_T.fftfreq(Nx, d=Dx))            # FFT frequencies in x 
  v    = F_F_T.fftshift(F_F_T.fftfreq(Ny, d=Dy))            # FFT frequencies in y
  u,v  = np.meshgrid(u, v)                                  # FFT frequencies in x and y
  fim  = np.sinc(2*(u**2 + v**2)**0.5)                      # Fourier image of 1/π/sqrt(1-r**2)

  def __init__(self, σx, σy):                               # Calculates the Convolution
    J = self.fim * np.exp(-2*π*π*((σx*self.v)**2 + (σy*self.u)**2))
    R = np.abs(F_F_T.ifftshift(F_F_T.irfft2(J, s=[self.Ny, self.Nx], workers=2, norm='forward')))
    limits = [None, None, None, None]
    self.SP =  R_B_S(self.y, self.x, R, bbox=limits, kx=3, ky=3 )              # this is the RectBivariateSpline
    print('fft: {:8.6f} {:8.6f} '.format(σx, σy), end = '')

  def __call__(self, xa, xb, ya, yb):                       # Call with integration range [xa,xb] [ya,yb] 
    return self.SP.integral(xa, xb, ya, yb)                 # Integration by Interpolation


class XSection(Pixels):
  def __init__(self):
    self.iteration = 0
    self.xmid_n    = np.zeros((self.X_npix             ), dtype=np.float64)
    self.ymid_n    = np.zeros((self.Y_npix             ), dtype=np.float64)
    self.Compton   = np.zeros((self.X_npix, self.Y_npix), dtype=np.float64)
    self.ellipse, self.compton, self.spreads = [], [], []

  def xSection(self):
    σx, σy         = self.spreads
    k, zx, zy, zz  = self.compton
    Xmax, Ymax     = 1 + σx/2, 1 + σy/2
    Rmax           = (Xmax**2 + Ymax**2)**0.5
    for i in range(self.X_npix):
      X   = self.xmid_n[i];  X2  = X*X
      if X > Rmax:
        for j in range(self.Y_npix): self.Compton[i][j] = 0.0
      else:
        uok = (1 + X)/2;                 u   = k*uok;             v   = 1 + u
        CSU = 1 + v*v - 4*uok*v*(1-uok); CSL = u*(u+2)*(1-2*uok); iCS = 1/(k*v*v*v)
        for j in range(self.Y_npix):
          Y = self.ymid_n[j]
          R = (X2+Y**2)**0.5
          if R > Rmax:  self.Compton[i][j] = 0.0
          else:         self.Compton[i][j] = (CSU + zy*u*Y + zz*CSL)*iCS


  def __call__(self,x,p):
    ellipse = [p[i] for i in (0,1,2,3)]
    compton = [p[i] for i in (4,5,6,7)]
    spreads = σx, σy = [p[i] for i in (8,9)]
    parameters_has_changed = False
    if self.spreads != spreads:          # if emittance parameters changed
      self.spreads   = spreads
      self.BLUR = ConvolutionFFT(σx, σy) # σx, σy are in arbitrary units
      if self.ellipse == ellipse:  self.SPLINE = self.BLUR.SP(self.xmid_n, self.ymid_n)
      parameters_has_changed = True
    if self.ellipse != ellipse: # if geometry parameters changed
      self.ellipse = X0, X1, Y0, Y1 = ellipse
      # 1) center position [mm]        2) unit radius [1/mm]
      CX = self.X_beam - (X1 + X0)/2;  RX = 2/(X1-X0)
      CY = self.Y_beam - (Y1 + Y0)/2;  RY = 2/(Y1-Y0)
      for xpix in range(self.X_npix):  self.xmid_n[xpix] = RX*(CX + (xpix+0.5)*self.X_pix)
      for ypix in range(self.Y_npix):  self.ymid_n[ypix] = RY*(CY + (ypix+0.5)*self.Y_pix)
      self.SPLINE = self.BLUR.SP(self.xmid_n, self.ymid_n)
      if self.compton == compton: self.xSection()
      parameters_has_changed = True
    if self.compton != compton: # if Compton parameters changed
      self.compton = compton
      self.xSection()
      parameters_has_changed = True
    if parameters_has_changed:
      self.Results = self.SPLINE*self.Compton
      self.iteration += 1;   print('iteration %d ' % (self.iteration))
    xpix = int((x[0]-self.X_beam) / self.X_pix)
    ypix = int((x[1]-self.Y_beam) / self.Y_pix)
    return p[10]*self.Results[xpix][ypix]

def cpuinfo():
  with open('/proc/cpuinfo', 'r') as f: info = f.readlines()
  model = [el for el in info if 'model name' in el]
  return model[0].strip().split(': ')[1]
  

def main(argv):
  FIT = (len(argv)==1)
  ROOT.gStyle.SetOptFit(1111)
  ROOT.gStyle.SetOptStat('ne')
  ROOT.gStyle.SetFitFormat("8.3g")
  ROOT.gROOT.SetStyle("Plain")
  ROOT.gROOT.ForceStyle()
#  ROOT.gStyle.SetPalette(56)
  XYD   = ROOT.TList()
  pure  = {
          'unpol'              :'Thu Dec  2 12:58:32 2021.root',
          'x=0.5'              :'Thu Dec  2 13:08:59 2021.root',
          'x=0.5, y=0.5'       :'Thu Dec  2 13:13:28 2021.root',
          'x=0.5, y=0.5, z=0.5':'Thu Dec  2 13:18:11 2021.root',
          }
  beam  = {
          'unpol-4x'           :'Thu Jan 27 19:28:31 2022.root',
          'unpol'              :'Thu Dec  2 14:16:20 2021.root',
          'x=0.5'              :'Thu Dec  2 14:13:04 2021.root',
          'x=0.5, y=0.5'       :'Thu Dec  2 14:06:51 2021.root',
          'x=0.5, y=0.5, z=0.5':'Thu Dec  2 13:53:13 2021.root',
          }
  hfile = ROOT.TFile(beam['unpol'])
  name  = hfile.GetListOfKeys()[0].GetName()
  lo, Eo, k, gf, to, sx, sy = [float(elem) for elem in name.split()]
  XYD = hfile.Get(name)
  hfile.Close()
  HDe, HDp = XYD[1], XYD[0]
  HDe.SetStats(0)
  X_gamma, dX_gamma = HDp.GetMean(1), HDp.GetMeanError(1)
#  X_gamma, dX_gamma = -213.538, 0.001
  HDd = HDe.Clone()
  HDd.SetStats(0)

  X0 = Pixels.X_beam; X1 = X0 + Pixels.X_size
  Y0 = Pixels.Y_beam; Y1 = Y0 + Pixels.Y_size
  
  fitf = XSection()
  FXY = ROOT.TF2('FXY', fitf , X0, X1, Y0, Y1, 11)
  FXY.SetParName(0, 'X_{0}');         FXY.SetParameter(0,   -0.004)  # beam position x, mm
  FXY.SetParName(1, 'X_{1}');         FXY.SetParameter(1,  347.66)  # edge position x, mm
  FXY.SetParName(2, 'Y_{0}');         FXY.SetParameter(2,   -1.068)  # y_min, mm
  FXY.SetParName(3, 'Y_{1}');         FXY.SetParameter(3,    1.068)  # y_max, mm
  FXY.SetParName(4, '#kappa');        FXY.SetParameter(4,       k )  # κ (kappa)
  FXY.SetParName(5, '#zeta_{x}');     FXY.SetParameter(5,      0.0); FXY.SetParLimits(5, -1.0, 1.0)
  FXY.SetParName(6, '#zeta_{y}');     FXY.SetParameter(6,      0.0); FXY.SetParLimits(6, -1.0, 1.0)
  FXY.SetParName(7, '#zeta_{z}');     FXY.SetParameter(7,      0.0); FXY.SetParLimits(7, -1.0, 1.0)
  FXY.SetParName(8, '#sigma_{x}');    FXY.SetParameter(8,    0.001); FXY.SetParLimits(8, 0.0001, 0.1) # σ_x [a.u.]
  FXY.SetParName(9, '#sigma_{y}');    FXY.SetParameter(9,    0.025); FXY.SetParLimits(9, 0.0001, 0.1) # σ_y [a.u.]
  FXY.SetParName(10,'norm');          FXY.SetParameter(10,     3e+2)  # amplitude

  cv = ROOT.TCanvas('cv','cv',10,10,1200,1000)
  cv.Divide(2,2)
  cv.cd(1); cv.GetPad(1).SetGrid()
  HDe.Draw('COLZ')
  cv.Modified()
  cv.Update()

  success = False
  ROOT.gBenchmark.Start( '2DFit' )
  if FIT:
    fixed_parameters = [4,5]
    for p in fixed_parameters: FXY.FixParameter(p, FXY.GetParameter(p))
    Result = HDe.Fit(FXY, 'SVNP')
    success = not Result.Status()
    for p in fixed_parameters: FXY.ReleaseParameter(p)

  if success:
    chi2 = Result.Chi2(); print('Chi2: %f'        % (chi2))
    NDF  = Result.Ndf();  print('NDF: %d'         % (NDF))
    prob = Result.Prob(); print('Probability: %f' % (prob))

    X0 = X_gamma;             dX0 = dX_gamma
    X1 = FXY.GetParameter(0); dX1 = FXY.GetParError(0) 
    X2 = FXY.GetParameter(1); dX2 = FXY.GetParError(1)
    Y1 = FXY.GetParameter(2); dY1 = FXY.GetParError(2)
    Y2 = FXY.GetParameter(3); dY2 = FXY.GetParError(3)
    A  = X2-X1;               dA  = (dX1**2+dX2**2)**0.5
    B  = Y2-Y1;               dB  = (dY1**2+dY2**2)**0.5
    print('A = %.4f +/- %.4f (%.5f)' % (A, dA, dA/A))

    R   = (X2-X1)/(X1-X0)
    R   = R/(1.0 + 0.5/to**2) # Correction !!!
    dR  = R * ( (dX0/(X1-X0))**2 + (dX1*(X0-X2)/(X1-X0)/(X2-X1))**2 + (dX2/(X2-X1))**2 )**0.5
    Eo  = 0.25*Constants.me**2/Laser.wo * R 
    dEo = 0.25*Constants.me**2/Laser.wo * dR

    print(R,  dR)
    print('E = %.5f +/- %.5f' % (Eo*1.e-9, dEo*1.e-9))
    Pt  = FXY.GetParameter(6); dPt = FXY.GetParError(6)

  cv.cd(3); cv.GetPad(3).SetGrid()
  FXY.SetNpx(Pixels.X_npix)
  FXY.SetNpy(Pixels.Y_npix)
  FXY.SetTitle('F(x,y)')
  FXY.Draw('COLZ1');  FXY.GetXaxis().SetTitle('X, mm')

  cv.cd(4); cv.GetPad(4).SetGrid()
  NZ = 0
  for binx in range(1, Pixels.X_npix):
    for biny in range(1, Pixels.Y_npix):
      H = HDd.GetBinContent(binx, biny)
      if H:
        F = FXY.Eval(Pixels.X_beam + (binx-0.5)*Pixels.X_pix, Pixels.Y_beam + (biny-0.5)*Pixels.Y_pix)
        HDd.SetBinContent(binx, biny, (F - H)/(1+abs(F))**0.5)
      else:
        NZ += 1
  HDd.SetTitle('(F(x,y) - H(x,y)) / (1+F(x,y))^{1/2}')
  HDd.Draw('COLZ1');  HDd.GetXaxis().SetTitle('X, mm')
  ROOT.gBenchmark.Show( '2DFit' )
  bench = 'Real time {:.0f} s, CPU time {:.0f} s'.format(ROOT.gBenchmark.GetRealTime('2DFit'), ROOT.gBenchmark.GetCpuTime('2DFit'))

  if success:
    NDF -= NZ
    prob = ROOT.TMath.Prob(chi2,NDF)
    cv.cd(2)
    pt = ROOT.TPaveText(.05, .05, .95, .96)
    pt.AddText(cpuinfo())
    pt.AddText(bench)
    pt.AddText('#chi^{2}/NDF = %.1f/%d | Prob = %.4f' % (chi2,NDF,prob))
    pt.AddText('X_{0} = %08.3f #pm %5.3f mm' % (X0, dX0))
    pt.AddText('X_{1} = %08.3f #pm %5.3f mm' % (X1, dX1))
    pt.AddText('X_{2} = %08.3f #pm %5.3f mm' % (X2, dX2))
#    pt.AddText('#zeta_{x} = %05.3f #pm %5.3f' % (FXY.GetParameter(5), FXY.GetParError(5)))
    pt.AddText('#zeta_{y} = %05.3f #pm %5.3f' % (FXY.GetParameter(6), FXY.GetParError(6)))
    pt.AddText('#zeta_{z} = %05.3f #pm %5.3f' % (FXY.GetParameter(7), FXY.GetParError(7)))
    pt.AddText('#sigma_{x} = %5.1f #pm %4.1f #mum' % (500*FXY.GetParameter(8)*A, 500*FXY.GetParError(8)*A))
    pt.AddText('#sigma_{y} = %5.2f #pm %4.2f #mum' % (500*FXY.GetParameter(9)*B, 500*FXY.GetParError(9)*B))
    pt.AddText('E_{beam} = %7.4f #pm %6.4f GeV. ' % (Eo*1.e-9, dEo*1.e-9))
    pt.Draw()

  cv.Modified()
  cv.Update()

  input()

  exit()

if __name__ == "__main__": main(sys.argv)


