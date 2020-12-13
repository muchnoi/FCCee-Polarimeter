#!/usr/bin/env python3
import math, sys, ROOT
import numpy as np

from hardware import Constants, Pixels, Laser

class Convolution_Kernel_2D:

  def Intexy(self, x0, x1, y0, y1):
    q = 2.**-0.5
    return 0.25*(math.erf(x1*q)-math.erf(x0*q))*(math.erf(y1*q)-math.erf(y0*q))

  def Meanxy(self, vmin, vmax):
#    return 0.5 * (vmin+vmax)
    q = 2.**-0.5
    return 2.*(2.*math.pi)**-0.5*(math.exp(-0.5*vmin**2) - math.exp(-0.5*vmax**2))/(math.erf(vmax*q) - math.erf(vmin*q))

  def Prepare(self, nX = 10, nY = 10, ScaleX = 1.0, ScaleY = 1.0, Threshold = 1.e-4):
    convolution_kernel, S = [], 0.0
    i = 100.; RgX, RgY = [-i, i], [-i, i] # (this is infinity)
    for i in range(-nX//2, nX//2+1): RgX.insert(-1, i*ScaleX)
    for j in range(-nY//2, nY//2+1): RgY.insert(-1, j*ScaleY)
    for i in range(len(RgX) - 1):
      print()
      for j in range(len(RgY) - 1):
        xmean = self.Meanxy( RgX[i], RgX[i+1])
        ymean = self.Meanxy( RgY[j], RgY[j+1] )
        value = self.Intexy( RgX[i], RgX[i+1], RgY[j], RgY[j+1] )
        if value >Threshold:
          convolution_kernel.append([xmean, ymean, value])
          S += value
          print(' %7.5f ' % value, end=' ')
        else: print('- - - - -', end=' ')
    print()
    print('The integral is: %.3f' % S)
    for el in convolution_kernel: el[2] = el[2] / S
    input()
    return convolution_kernel

class Integrator(Pixels):
  def __init__(self, convolve, Rebinx):
    self.i2pi = 0.5/math.pi
    if convolve:
       self.convolve = True
       CK = Convolution_Kernel_2D()
#       self.simple_convolution = CK.Prepare(nX = 14, nY = 14, ScaleX = 0.5, ScaleY=0.5, Threshold = 1.e-4)
       self.simple_convolution = CK.Prepare(nX = 16, nY = 16, ScaleX = 0.4, ScaleY=0.4, Threshold = 5.e-5)
    else: self.convolve = False
    self.X_pix  = self.X_pix  * Rebinx
    self.X_npix = self.X_npix // Rebinx
    self.iteration = 0
    self.result = np.zeros((self.X_npix, self.Y_npix), dtype=np.float64)
    self.xmid_n = np.zeros((self.X_npix), dtype=np.float64)
    self.ymid_n = np.zeros((self.Y_npix), dtype=np.float64)
    self.parset    = []


  def Convolute(self, nR, X0, X1, Y0, Y1):
    S = 0.0
    if   nR > 1.2:  return 0.0
    for p in self.simple_convolution:
      dx, dy = p[0]*self.sx_n, p[1]*self.sy_n
      S += p[2]*self.Integrate(X0-dx, X1-dx, Y0-dy, Y1-dy)
    return S

  def Integrate(self, X0, X1, Y0, Y1):
    I = 0.0
    X02, X12, IY02, IY12  = X0*X0, X1*X1, 1.0-Y0*Y0, 1.0-Y1*Y1
    D00, D10 = math.sqrt(max(0.0, IY02 - X02)), math.sqrt(max(0.0, IY02 - X12))
    D01, D11 = math.sqrt(max(0.0, IY12 - X02)), math.sqrt(max(0.0, IY12 - X12))
    if not (D00 + D10 + D01 + D11): return 0.0
    Y0Y1   = Y0*Y1;   X0X1   = X0*X1
    D11D10 = D11*D10; D00D01 = D00*D01

    NUM, DEN = Y1*D10 - Y0*D11, Y0Y1 + D11D10
    if NUM: I = I + X1 * math.atan2(NUM, DEN)
    NUM, DEN = Y1*D00 - Y0*D01, Y0Y1 + D00D01
    if NUM: I = I - X0 * math.atan2(NUM, DEN)
    NUM, DEN = X1*D01 - X0*D11, X0X1 + D11*D01
    if NUM: I = I + Y1 * math.atan2(NUM, DEN)
    NUM, DEN = X1*D00 - X0*D10, X0X1 + D00*D10
    if NUM: I = I - Y0 * math.atan2(NUM, DEN)
    NUM, DEN = X1*(Y0*D11 - Y1*D10), X12*Y0Y1 + D10*D11
    if NUM: I = I +      math.atan2(NUM, DEN)
    NUM, DEN = X0*(Y1*D00 - Y0*D01), X02*Y0Y1 + D00*D01
    if NUM: I = I +      math.atan2(NUM, DEN)

    I *= self.i2pi
    if   I >  0.0: pass
    elif I < -0.4: I +=  0.50
    else: print('I = %f' % I); input()
    return I

  def __call__(self,x,p):
    X, Y = x[0], x[1]
    self.curpar = [p[i] for i in range(9)]
    if self.parset != self.curpar:
      X0, X1, Y0, Y1, self.k, self.Pl, self.Pt, sx, sy = self.curpar
      self.ik = 1./self.k
      CX, RX = self.X_beam-0.5*(X1+X0), 2.0/(X1-X0); PX = self.X_pix*RX; self.sx_n = sx*RX
      CY, RY = self.Y_beam-0.5*(Y1+Y0), 2.0/(Y1-Y0); PY = self.Y_pix*RY; self.sy_n = sy*RY
      self.iteration += 1;   print('iteration %d ' % (self.iteration))
      self.parset = self.curpar
      for xpix in range(self.X_npix): self.xmid_n[xpix] = (CX + (xpix+0.5)*self.X_pix)*RX
      for ypix in range(self.Y_npix): self.ymid_n[ypix] = (CY + (ypix+0.5)*self.Y_pix)*RY
      PX, PY = PX*0.499, PY*0.499
      for i in range(self.X_npix):
        xx = self.xmid_n[i]*self.xmid_n[i]
        for j in range(self.Y_npix):
          yy = self.ymid_n[j]*self.ymid_n[j]
          nR = math.sqrt(xx+yy)
          if self.convolve:
            self.result[i][j] = self.Convolute(nR, self.xmid_n[i]-PX, self.xmid_n[i]+PX, self.ymid_n[j]-PY, self.ymid_n[j]+PY)
          else:
            self.result[i][j] = self.Integrate(    self.xmid_n[i]-PX, self.xmid_n[i]+PX, self.ymid_n[j]-PY, self.ymid_n[j]+PY)
    xpix = int((X-self.X_beam) / self.X_pix)
    ypix = int((Y-self.Y_beam) / self.Y_pix)
#    return p[9]*self.result[xpix][ypix]
    u = 0.5*self.k*(1.+self.xmid_n[xpix])
    v = u + 1.
    iCS = self.k*v*v*v
    CSU = 1.+ v*v - 4.*u*self.ik*v*(1.-u*self.ik)
    CST = u*self.ymid_n[ypix]
#   CSL = u*(u+2.)*(1.-2.*u*self.ik)
    return p[9]*self.result[xpix][ypix]*(CSU + self.Pt*CST)/iCS


def main(argv):

  ROOT.gStyle.SetOptFit(1111)
  ROOT.gStyle.SetOptStat('ne')
  ROOT.gStyle.SetFitFormat("8.3g")
#  ROOT.gROOT.SetStyle("Plain")
#  ROOT.gROOT.ForceStyle()
  ROOT.gStyle.SetPalette(56)
  XYD   = ROOT.TList()
#  hfile = ROOT.TFile('zero-emittance.root') # 1600 x 80
#  hfile = ROOT.TFile('x4-true-emittance.root') # 1600 x 80
  hfile = ROOT.TFile('Sun Dec 13 10:49:17 2020.root') # 1600 x 80
  name  = hfile.GetListOfKeys()[0].GetName()
  lo, Eo, k, gf, to, sx, sy = [float(elem) for elem in name.split()]
  XYD = hfile.Get(name)
  hfile.Close()
  HDe, HDp = XYD[1], XYD[0]
  X_gamma, dX_gamma = HDp.GetMean(1), HDp.GetMeanError(1)
  rebin_x = 8; HDe.RebinX(rebin_x) # use to decrease the number of x-pixels
  HDd = HDe.Clone()

  X0 = Pixels.X_beam; X1 = X0 + Pixels.X_size
  Y0 = Pixels.Y_beam; Y1 = Y0 + Pixels.Y_size
  
  function = Integrator(convolve=True, Rebinx = rebin_x)
  FXY = ROOT.TF2('FXY', function , X0, X1, Y0, Y1, 10)
  FXY.SetParName(0, 'X_{0}');         FXY.SetParameter(0,    0.000)  # beam position x, mm
  FXY.SetParName(1, 'X_{1}');         FXY.SetParameter(1,  347.640)  # edge position x, mm
  FXY.SetParName(2, 'Y_{0}');         FXY.SetParameter(2,   -1.068)  # y_min, mm
  FXY.SetParName(3, 'Y_{1}');         FXY.SetParameter(3,    1.068)  # y_max, mm
  FXY.SetParName(4, '#kappa');        FXY.SetParameter(4,    1.628)  # κ (kappa)
  FXY.SetParName(5, 'P_{#parallel}'); FXY.SetParameter(5,      0.0)  # longitudinal polarization degree, a. u.
  FXY.SetParName(6, 'P_{#perp}');     FXY.SetParameter(6,      0.1); FXY.SetParLimits(6, 0.0, 1.00)
  FXY.SetParName(7, '#sigma_{x}');    FXY.SetParameter(7,    0.204); FXY.SetParLimits(7, 0.0, 1.00) # σ_x,   mm
  FXY.SetParName(8, '#sigma_{y}');    FXY.SetParameter(8,    0.024); FXY.SetParLimits(8, 0.0, 1.00) # σ_y,   mm
  FXY.SetParName(9, 'norm');          FXY.SetParameter(9, 4*1.305e+7)  # amplitude

  cv = ROOT.TCanvas('cv','cv',10,10,1200,1000)
  cv.Divide(2,2)
  cv.cd(1); cv.GetPad(1).SetGrid()
  HDe.Draw('COLZ')
  cv.Modified()
  cv.Update()

  success = False
  ROOT.gBenchmark.Start( '2DFit' )
  fixed_parameters = [4,5]
  for p in fixed_parameters: FXY.FixParameter(p, FXY.GetParameter(p))
  Result = HDe.Fit(FXY, 'SVN')
  success = not Result.Status()
  for p in fixed_parameters: FXY.ReleaseParameter(p)

  if success:
    chi2 = Result.Chi2(); print('Chi2: %f'        % (chi2))
    NDF  = Result.Ndf();  print('NDF: %d'         % (NDF))
    prob = Result.Prob(); print('Probability: %f' % (prob))

    X0 = X_gamma;             dX0 = dX_gamma
    X1 = FXY.GetParameter(0); dX1 = FXY.GetParError(0)
    X2 = FXY.GetParameter(1); dX2 = FXY.GetParError(1)
    A  = X2-X1;               dA  = (dX1**2+dX2**2)**0.5
    print('A = %.4f +/- %.4f (%.5f)' % (A, dA, dA/A))

    R   = (X2-X1)/(X1-X0)
    dR  = R * ( (dX0/(X1-X0))**2 + (dX1*(X0-X2)/(X1-X0)/(X2-X1))**2 + (dX2/(X2-X1))**2 )**0.5
    Eo  = 0.25*Constants.me**2/Laser.wo * R
    dEo = 0.25*Constants.me**2/Laser.wo * dR

    print(R,  dR)
    print('E = %.5f +/- %.5f' % (Eo*1.e-9, dEo*1.e-9))
    Pt  = FXY.GetParameter(6); dPt = FXY.GetParError(6)

    cv.cd(2)
    pt = ROOT.TPaveText(.05, .05, .95, .96)
    pt.AddText('')
    pt.AddText('#chi^{2}/NDF = %.1f/%d | Prob = %.4f' % (chi2,NDF,prob))
    pt.AddText('X_{0} = %08.3f #pm %5.3f mm' % (X0, dX0))
    pt.AddText('X_{1} = %08.3f #pm %5.3f mm' % (X1, dX1))
    pt.AddText('X_{2} = %08.3f #pm %5.3f mm' % (X2, dX2))
    pt.AddText('#sigma_{x} = %5.1f #pm %4.1f #mum' % (1000*FXY.GetParameter(7), 1000*FXY.GetParError(7)))
    pt.AddText('#sigma_{y} = %5.2f #pm %4.2f #mum' % (1000*FXY.GetParameter(8), 1000*FXY.GetParError(8)))
    pt.AddText('E_{beam} = %7.4f #pm %6.4f GeV. ' % (Eo*1.e-9, dEo*1.e-9))
    pt.AddText('P_{#perp}  = %6.4f #pm %6.4f' % (Pt, dPt))
    pt.AddText('')
    pt.Draw()

  cv.cd(3); cv.GetPad(3).SetGrid()
  FXY.SetNpx(Pixels.X_npix)
  FXY.SetNpy(Pixels.Y_npix)
  FXY.SetTitle('F(x,y)')
  FXY.Draw('COLZ');  FXY.GetXaxis().SetTitle('X, mm')


  cv.cd(4); cv.GetPad(4).SetGrid()
  for binx in range(1, Pixels.X_npix//rebin_x):
    for biny in range(1, Pixels.Y_npix):
      H = HDd.GetBinContent(binx, biny)
      if H:
        F = FXY.Eval(Pixels.X_beam + (binx-0.5)*rebin_x*Pixels.X_pix, Pixels.Y_beam + (biny-0.5)*Pixels.Y_pix)
        HDd.SetBinContent(binx, biny, (F - H)/H**0.5)
  HDd.SetTitle('(F(x,y) - H(x,y)) / H(x,y)^{1/2}')
  HDd.Draw('COLZ');  HDd.GetXaxis().SetTitle('X, mm')


  cv.Modified()
  cv.Update()
  ROOT.gBenchmark.Show( '2DFit' )

  input()

  exit()

if __name__ == "__main__": main(sys.argv)


