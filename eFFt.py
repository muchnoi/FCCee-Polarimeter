#!/usr/bin/env python3
import sys, ROOT
import numpy as np; π = np.pi
from scipy             import fft                 as F_F_T
from scipy.interpolate import RectBivariateSpline as R_B_S
from hardware          import Constants, Pixels, Laser

class xSection:
  Rx, Ry     = 1.5, 1.5                                     # ratio of calculation range to unit radius
  Nx, Ny     = 2048, 2048                                   # number of x,y divisions
  Dx, Dy     = 2*Rx/Nx, 2*Ry/Ny                             # distance between neighbour knots
  x          = np.linspace(Dx/2-Rx, Rx-Dx/2, num=Nx)        # knots positions in x
  y          = np.linspace(Dy/2-Ry, Ry-Dy/2, num=Ny)        # knots positions in y
  Y, X       = np.meshgrid(y, x)                            # knots grid in x and y
  νx         = F_F_T.fftshift(F_F_T.fftfreq(Nx, d=Dx))      # FFT frequencies in x 
  νy         = F_F_T.fftshift(F_F_T.fftfreq(Ny, d=Dy))      # FFT frequencies in y
  νy, νx     = np.meshgrid(νy, νx)                          # FFT grid in νx and νy
  img        = np.sinc(2*(νx**2 + νy**2)**0.5)              # Fourier image of 1/π/sqrt(1-r**2)
  emittance  = [None, None]
  stockspar  = [None, None, None]

  def __init__(self, κ):
    uoκ      = (1 + self.X)/2                               # u over kappa
    u        = κ*uoκ                                        # u
    idn      = 1/(κ*(1+u)**3)                               # inverse denominator
    self.xso = idn*(1 + (1+u)**2 - 4*uoκ*(1+u)*(1-uoκ))     # unpolarized xSection
    self.xsy = idn*u*self.Y                                 # transverse vertical electron polarization
    self.xsz = idn*u*(u + 2)*(1 - 2*uoκ)                    # longitudinal electron polarization

  def Prepare(self, parameters):
    ζx, ζy, ζz, σx, σy, norm = parameters
    if [σx, σy]     != self.emittance:
      self.emittance = [σx, σy]
      fg             = self.img * np.exp(-2*π*π*((σx*self.νx)**2 + (σy*self.νy)**2))             # Fourier image by convolution theorem
      self.cvl       = np.abs(F_F_T.ifftshift(F_F_T.irfft2(fg, s=[self.Ny, self.Nx], workers=-1, norm='forward'))) # convolution result
      print('emittance: σx={:8.6f} σy={:8.6f} '.format(σx, σy))#, end = '')
    if [ζx, ζy, ζz] != self.stockspar:
      self.stockspar = [ζx, ζy, ζz]
      self.xst       = self.xso + ζy*self.xsy + ζz*self.xsz
      print('stockspar: ζx={:8.6f} ζy={:8.6f} ζz={:8.6f} '.format(ζx, ζy, ζz))#, end = '')
    self.BLUR      =  R_B_S(self.y, self.x, norm*self.cvl*self.xst, kx=3, ky=3 )  # this is the RectBivariateSpline


class Electrons2D(Pixels):
  def __init__(self, k):
    self.iteration, self.ellipse, self.compton = 0, [], []
    self.xmid_n      = np.zeros((self.X_npix), dtype=np.double)
    self.ymid_n      = np.zeros((self.Y_npix), dtype=np.double)
    self.XS          = xSection(k)
    self.Ax, self.Bx = 1/self.X_pix, - self.X_beam/self.X_pix
    self.Ay, self.By = 1/self.Y_pix, - self.Y_beam/self.Y_pix

  def __call__(self,x,p):
    ellipse, compton, change = [p[i] for i in (0,1,2,3)], [p[i] for i in (4,5,6,7,8,9)], False
    if self.compton != compton: # if Compton parameters changed
      self.compton, change = compton, True
      self.XS.Prepare(compton)
    if self.ellipse != ellipse: # if geometry parameters changed
      self.ellipse = X0, X1, Y0, Y1 = ellipse; change = True
      # 1) center position [mm]        2) unit radius [1/mm]
      CX = self.X_beam - (X1 + X0)/2;  RX = 2/(X1-X0)
      CY = self.Y_beam - (Y1 + Y0)/2;  RY = 2/(Y1-Y0)
      for xpix in range(self.X_npix):  self.xmid_n[xpix] = RX*(CX + (xpix+0.5)*self.X_pix)
      for ypix in range(self.Y_npix):  self.ymid_n[ypix] = RY*(CY + (ypix+0.5)*self.Y_pix)
    if change:
      self.SPLINE = self.XS.BLUR(self.xmid_n, self.ymid_n)
      self.iteration += 1;   print('iteration: %d ' % (self.iteration))
    return self.SPLINE[int(self.Ax*x[0] + self.Bx)][int(self.Ay*x[1] + self.By)]


def cpuinfo():
  with open('/proc/cpuinfo', 'r') as f: info = f.readlines()
  model = [el for el in info if 'model name' in el]
  return model[0].strip().split(': ')[1]


class POLARIMETER(Pixels):
  DataMC = {}

  def __init__(self):
    beam  = {
            'unpol-4x'           :'Thu Jan 27 19:28:31 2022.root',
            'unpol'              :'Thu Dec  2 14:16:20 2021.root',
            'x=0.5'              :'Thu Dec  2 14:13:04 2021.root',
            'x=0.5, y=0.5'       :'Thu Dec  2 14:06:51 2021.root',
            'x=0.5, y=0.5, z=0.5':'Thu Dec  2 13:53:13 2021.root',
            }
    XYD   = ROOT.TList()
    hfile = ROOT.TFile(beam['unpol'])
    name  = hfile.GetListOfKeys()[0].GetName()
    XYD   = hfile.Get(name)
    hfile.Close()
    for i in range(7): self.DataMC[['λo','Eo','κ','γ','θo','σ`x','σ`y'][i]]=[float(p) for p in name.split()][i]
    self.HDe, self.HDp, self.HDd = XYD[1], XYD[0], XYD[1].Clone()
    self.HDe.SetStats(0);     self.HDd.SetStats(0)

  def SetFunction(self):
    self.fitf = Electrons2D(self.DataMC['κ'])
    Xfrom = Pixels.X_beam; Xto = Xfrom + Pixels.X_size
    Yfrom = Pixels.Y_beam; Yto = Yfrom + Pixels.Y_size
    self.FXY = ROOT.TF2('FXY', self.fitf , Xfrom, Xto, Yfrom, Yto, 10)
    self.FXY.SetParName(0, 'X_{1}');       self.FXY.SetParameter(0,  -0.01)  # beam position x, mm
    self.FXY.SetParName(1, 'X_{2}');       self.FXY.SetParameter(1,  347.5)  # edge position x, mm
    self.FXY.SetParName(2, 'Y_{0}');       self.FXY.SetParameter(2, -1.068)  # y_min, mm
    self.FXY.SetParName(3, 'Y_{1}');       self.FXY.SetParameter(3,  1.068)  # y_max, mm
    self.FXY.SetParName(4, '#zeta_{x}');   self.FXY.SetParameter(4,    0.0); self.FXY.SetParLimits(4, -1.0, 1.0)
    self.FXY.SetParName(5, '#zeta_{y}');   self.FXY.SetParameter(5,    0.0); self.FXY.SetParLimits(5, -1.0, 1.0)
    self.FXY.SetParName(6, '#zeta_{z}');   self.FXY.SetParameter(6,    0.0); self.FXY.SetParLimits(6, -1.0, 1.0)
    self.FXY.SetParName(7, '#sigma_{x}');  self.FXY.SetParameter(7,  0.001); self.FXY.SetParLimits(7, 0.0001, 0.1) # σ_x [a.u.]
    self.FXY.SetParName(8, '#sigma_{y}');  self.FXY.SetParameter(8,  0.026); self.FXY.SetParLimits(8, 0.0001, 0.1) # σ_y [a.u.]
    self.FXY.SetParName(9, 'norm');        self.FXY.SetParameter(9,   2e+2)  # amplitude
    self.FXY.SetNpx(Pixels.X_npix);        self.FXY.SetNpy(Pixels.Y_npix)
    self.FXY.SetTitle('F(x,y)');           self.FXY.GetXaxis().SetTitle('X, mm')

  def FitElectrons2D(self):
    fixed_parameters = [4]
    for p in fixed_parameters: self.FXY.FixParameter(p, self.FXY.GetParameter(p))
    self.FitElectronsResult = self.HDe.Fit(self.FXY, 'SVNP')
    for p in fixed_parameters: self.FXY.ReleaseParameter(p)
    return not self.FitElectronsResult.Status()

  def Residuals(self):
    self.NZeros = 0
    for binx   in range(1, Pixels.X_npix):
      for biny in range(1, Pixels.Y_npix):
        H = self.HDd.GetBinContent(binx, biny)
        if H:
          F = self.FXY.Eval(Pixels.X_beam + (binx-0.5)*Pixels.X_pix, Pixels.Y_beam + (biny-0.5)*Pixels.Y_pix)
          self.HDd.SetBinContent(binx, biny, (F - H)/(1+abs(F))**0.5)
        else:
          self.NZeros += 1
    self.HDd.SetTitle('(F(x,y) - H(x,y)) / (1+F(x,y))^{1/2}')


  def ElectronsResults(self):
    chi2 = self.FitElectronsResult.Chi2(); print('Chi2: %f'        % (chi2))
    NDF  = self.FitElectronsResult.Ndf();  print('NDF: %d'         % (NDF) )
    prob = self.FitElectronsResult.Prob(); print('Probability: %f' % (prob))
    NDF -= self.NZeros;                    prob = ROOT.TMath.Prob(chi2, NDF)

    X0 = self.HDp.GetMean(1);      dX0 = self.HDp.GetMeanError(1)
    X1 = self.FXY.GetParameter(0); dX1 = self.FXY.GetParError(0) 
    X2 = self.FXY.GetParameter(1); dX2 = self.FXY.GetParError(1)
    Y1 = self.FXY.GetParameter(2); dY1 = self.FXY.GetParError(2)
    Y2 = self.FXY.GetParameter(3); dY2 = self.FXY.GetParError(3)
    A  = X2-X1;                    dA  = (dX1**2+dX2**2)**0.5
    B  = Y2-Y1;                    dB  = (dY1**2+dY2**2)**0.5
    print('A = %.4f +/- %.4f (%.5f)' % (A, dA, dA/A))
    R   = (X2-X1)/(X1-X0)
    R  /= 1.0 + 0.5/self.DataMC['θo']**2 # Correction !!!
    dR  = R * ( (dX0/(X1-X0))**2 + (dX1*(X0-X2)/(X1-X0)/(X2-X1))**2 + (dX2/(X2-X1))**2 )**0.5
    Eo  = 0.25*Constants.me**2/Laser.wo * R 
    dEo = 0.25*Constants.me**2/Laser.wo * dR
    print(R,  dR)
    print('E = %.5f +/- %.5f' % (Eo*1.e-9, dEo*1.e-9))
    self.ResultsTable = ROOT.TPaveText(.05, .05, .95, .96)
    self.ResultsTable.AddText(cpuinfo())
    self.ResultsTable.AddText('Real time {:.0f} s, CPU time {:.0f} s'.format(ROOT.gBenchmark.GetRealTime('2DFit'), ROOT.gBenchmark.GetCpuTime('2DFit')))
    self.ResultsTable.AddText('#chi^{2}/NDF = %.1f/%d | Prob = %.4f' % (chi2,NDF,prob))
    self.ResultsTable.AddText('X_{0} = %08.3f #pm %5.3f mm' % (X0, dX0))
    self.ResultsTable.AddText('X_{1} = %08.3f #pm %5.3f mm' % (X1, dX1))
    self.ResultsTable.AddText('X_{2} = %08.3f #pm %5.3f mm' % (X2, dX2))
    self.ResultsTable.AddText('#zeta_{y} = %05.3f #pm %5.3f' % (self.FXY.GetParameter(5), self.FXY.GetParError(5)))
    self.ResultsTable.AddText('#zeta_{z} = %05.3f #pm %5.3f' % (self.FXY.GetParameter(6), self.FXY.GetParError(6)))
    self.ResultsTable.AddText('#sigma_{x} = %5.1f #pm %4.1f #mum' % (500*self.FXY.GetParameter(7)*A, 500*self.FXY.GetParError(7)*A))
    self.ResultsTable.AddText('#sigma_{y} = %5.2f #pm %4.2f #mum' % (500*self.FXY.GetParameter(8)*B, 500*self.FXY.GetParError(8)*B))
    self.ResultsTable.AddText('E_{beam} = %7.4f #pm %6.4f GeV. ' % (Eo*1.e-9, dEo*1.e-9))


class DISPLAY:
  def __init__(self):
    ROOT.gStyle.SetOptFit(1111);    ROOT.gStyle.SetOptStat('ne');    ROOT.gStyle.SetFitFormat("8.3g")
    ROOT.gROOT.SetStyle("Plain");   ROOT.gROOT.ForceStyle()      #   ROOT.gStyle.SetPalette(56)
    self.cv = ROOT.TCanvas('cv','cv',10,10,1200,1000);               self.cv.Divide(2,2)

  def ShowOnPad(self, nPad, entity, grid=False, goption=''):
    self.cv.cd(nPad) 
    if grid: self.cv.GetPad(nPad).SetGrid()
    entity.Draw(goption)
    self.cv.Modified();    self.cv.Update()


def main(argv):
  DATA = POLARIMETER()
  LOOK = DISPLAY()
  LOOK.ShowOnPad(nPad=1, entity = DATA.HDe, grid = True, goption='COLZ')
  DATA.SetFunction()
  ROOT.gBenchmark.Start('2DFit')
  if len(argv)==1:      success = DATA.FitElectrons2D()
  else:                 success = False
  LOOK.ShowOnPad(nPad=3, entity = DATA.FXY, grid = True, goption='COLZ1')
  DATA.Residuals()
  LOOK.ShowOnPad(nPad=4, entity = DATA.HDd, grid = True, goption='COLZ1')
  ROOT.gBenchmark.Show('2DFit')
  if success: 
    DATA.ElectronsResults()
    LOOK.ShowOnPad(nPad=2, entity = DATA.ResultsTable)
  input()
  exit()

if __name__ == "__main__": main(sys.argv)

