#!/usr/bin/env python3
import sys, ROOT, argparse, pickle
from hardware import Constants
from xSection import pPIXELS, ePIXELS


def cpuinfo():
  with open('/proc/cpuinfo', 'r') as f: info = f.readlines()
  model = [el for el in info if 'model name' in el]
  return model[0].strip().split(': ')[1]


class POLARIMETER:
  def __init__(self, filename=None, grebin=1):
    if filename:
      with open(filename, 'rb') as fp: self.MC = pickle.load(fp)
    else:
      print('Please specify filename!'); exit()
    self.grebin = grebin
    self.HDp = self.MC['XYp']
    self.HDp.RebinX(grebin)
    self.HDp.RebinY(grebin)
    self.DDp = self.HDp.Clone()
    self.HDe = self.MC['XYe']
    self.DDe = self.HDe.Clone()
    self.HDe.SetStats(0); self.HDe.SetTitle('Electrons: MC')
    self.HDp.SetStats(0); self.HDp.SetTitle('Photons: MC')
    self.DDe.SetStats(0); self.DDe.SetTitle('Electrons: residuald')
    self.DDp.SetStats(0); self.DDp.SetTitle('Photons: residuals')
    
  def ParametersMC(self):
    self.ParametersTable = ROOT.TPaveText(.05, .05, .95, .96)
    self.ParametersTable.AddText('Monte-Carlo Parameters:')
    self.ParametersTable.AddText('Laser #lambda_{0} = %5.3f um' %       (1e+4*self.MC['Laser'].λo))
    self.ParametersTable.AddText('Electron E_{0} = %6.3f GeV' %         (1e-9*self.MC['Spectrometer'].Eo))
    self.ParametersTable.AddText('Electron #gamma = %5.3f#times10^{3}'% (1e-3*self.MC['Spectrometer'].γ))
    self.ParametersTable.AddText('Compton #kappa = %5.3f'      %        (     self.MC['Spectrometer'].κ))
    self.ParametersTable.AddText('Bend: #gamma#theta_{0} = %5.3f'  %    (     self.MC['Spectrometer'].γ*self.MC['Spectrometer'].θo))
    outstr = '(#xi_{1}, #xi_{2}, #xi_{3}) = (%6.3f, %6.3f, %6.3f)'
    self.ParametersTable.AddText(outstr % (self.MC['Laser'].ξ1, self.MC['Laser'].ξ2, self.MC['Laser'].ξ3))
    outstr = '(#zeta_{x}, #zeta_{y}, #zeta_{z}) = (%6.3f, %6.3f, %6.3f)'
    self.ParametersTable.AddText(outstr % (self.MC['Spectrometer'].ζx, self.MC['Spectrometer'].ζy, self.MC['Spectrometer'].ζz))

  def ElectronsSetFunction(self):
    self.efit = ePIXELS(EPD = self.MC['EPD'], setup = self.MC['Spectrometer'])
    X0 = self.MC['EPD'].X_beam; X1 = X0 + self.MC['EPD'].X_size
    Y0 = self.MC['EPD'].Y_beam; Y1 = Y0 + self.MC['EPD'].Y_size
    ξ1 = self.MC['Laser'].ξ1
    ξ2 = self.MC['Laser'].ξ2
    ζx = self.MC['Laser'].ξ3 * self.MC['Spectrometer'].ζx
    ζy = self.MC['Laser'].ξ3 * self.MC['Spectrometer'].ζy
    ζz = self.MC['Laser'].ξ3 * self.MC['Spectrometer'].ζz
    self.EXY = ROOT.TF2('EXY', self.efit , X0, X1, Y0, Y1, 12)
    self.EXY.SetParName(0,  'X1');      self.EXY.SetParameter(0,   -0.15)  # beam position x, mm
    self.EXY.SetParName(1,  'X2');      self.EXY.SetParameter(1,   347.75)  # edge position x, mm
    self.EXY.SetParName(2,  'Y0');      self.EXY.SetParameter(2,  -1.068)  # y_min, mm
    self.EXY.SetParName(3,  'Y1');      self.EXY.SetParameter(3,   1.068)  # y_max, mm
    self.EXY.SetParName(4,  'ξ1');      self.EXY.SetParameter(4,      ξ1); #self.EXY.SetParLimits(4,  -1.0, 1.0)
    self.EXY.SetParName(5,  'ξ2');      self.EXY.SetParameter(5,      ξ2); #self.EXY.SetParLimits(5,  -1.0, 1.0)
    self.EXY.SetParName(6,  'ξ3*ζx');   self.EXY.SetParameter(6,      ζx); #self.EXY.SetParLimits(6,  -1.0, 1.0)
    self.EXY.SetParName(7,  'ξ3*ζy');   self.EXY.SetParameter(7,      ζy); #self.EXY.SetParLimits(7,  -1.0, 1.0)
    self.EXY.SetParName(8,  'ξ3*ζz');   self.EXY.SetParameter(8,      ζz); #self.EXY.SetParLimits(8,  -1.0, 1.0)
    self.EXY.SetParName(9,  'σx');      self.EXY.SetParameter(9,  1.2e-3); self.EXY.SetParLimits( 9, 0.0001, 0.1) # σ_x [a.u.]
    self.EXY.SetParName(10, 'σy');      self.EXY.SetParameter(10, 2.6e-2); self.EXY.SetParLimits(10, 0.0001, 0.1) # σ_y [a.u.]
    self.EXY.SetParName(11, 'norm');    self.EXY.SetParameter(11, 100.0)  # amplitude
    self.EXY.SetNpx(self.MC['EPD'].X_npix)
    self.EXY.SetNpy(self.MC['EPD'].Y_npix)
    self.EXY.SetTitle('Electrons: Fit'); self.EXY.GetXaxis().SetTitle('X, mm')

  def FitElectrons(self):
    fixed_parameters = [4,5,6,7,8]
    for p in fixed_parameters: self.EXY.FixParameter(p, self.EXY.GetParameter(p))
    self.FitElectronsResult =  self.HDe.Fit(self.EXY, 'SVNP') # Use Pearsons chi-square method
    for p in fixed_parameters: self.EXY.ReleaseParameter(p)
    return not self.FitElectronsResult.Status()

  def PhotonsSetFunction(self):
    self.pfit = pPIXELS(PPD = self.MC['PPD'], setup = self.MC['Spectrometer'], rebin = self.grebin)
    X0 = self.MC['PPD'].X_beam + self.MC['PPD'].X_pix; X1 = X0 + self.MC['PPD'].X_size - self.MC['PPD'].X_pix
    Y0 = self.MC['PPD'].Y_beam + self.MC['PPD'].Y_pix; Y1 = Y0 + self.MC['PPD'].Y_size - self.MC['PPD'].Y_pix
    ξ1 = self.MC['Laser'].ξ1
    ξ2 = self.MC['Laser'].ξ2
    ζx = self.MC['Laser'].ξ3 * self.MC['Spectrometer'].ζx
    ζy = self.MC['Laser'].ξ3 * self.MC['Spectrometer'].ζy
    ζz = self.MC['Laser'].ξ3 * self.MC['Spectrometer'].ζz
    self.PXY = ROOT.TF2('PXY', self.pfit , X0, X1, Y0, Y1, 10)
    self.PXY.SetParName(0, 'X0');      self.PXY.SetParameter(0,    self.HDp.GetMean(1))
    self.PXY.SetParName(1, 'Y0');      self.PXY.SetParameter(1,    self.HDp.GetMean(2))
    self.PXY.SetParName(2, 'ξ1');      self.PXY.SetParameter(2,     ξ1); self.PXY.SetParLimits(2,  -1.1, 1.1)
    self.PXY.SetParName(3, 'ξ2');      self.PXY.SetParameter(3,     ξ2); self.PXY.SetParLimits(3,  -1.1, 1.1)
    self.PXY.SetParName(4, 'ξ3*ζx');   self.PXY.SetParameter(4,     ζx); self.PXY.SetParLimits(4,  -1.1, 1.1)
    self.PXY.SetParName(5, 'ξ3*ζy');   self.PXY.SetParameter(5,     ζy); self.PXY.SetParLimits(5,  -1.1, 1.1)
    self.PXY.SetParName(6, 'ξ3*ζz');   self.PXY.SetParameter(6,     ζz); self.PXY.SetParLimits(6,  -1.1, 1.1)
    self.PXY.SetParName(7, 'σx');      self.PXY.SetParameter(7,    0.2); self.PXY.SetParLimits(7, 0.001, 1.0) # σ_x [mm]
    self.PXY.SetParName(8, 'σy');      self.PXY.SetParameter(8,   0.05); self.PXY.SetParLimits(8, 0.001, 1.0) # σ_y [mm]
    self.PXY.SetParName(9,' norm');    self.PXY.SetParameter(9,   9e+3)  # amplitude
    self.PXY.SetNpx(self.MC['PPD'].X_npix)
    self.PXY.SetNpy(self.MC['PPD'].Y_npix)
    self.PXY.SetTitle('Photons: Fit'); self.PXY.GetXaxis().SetTitle('X, mm')

  def FitPhotons(self):
    fixed_parameters = []
    for p in fixed_parameters: self.PXY.FixParameter(p, self.PXY.GetParameter(p))
    if self.grebin < 4: self.FitPhotonsResult =  self.HDp.Fit(self.PXY, 'SVNP') 
    else:               self.FitPhotonsResult =  self.HDp.Fit(self.PXY, 'SVNI') 
    for p in fixed_parameters: self.PXY.ReleaseParameter(p)
    return not self.FitPhotonsResult.Status()

  def ElectronsResiduals(self):
    self.eZeros = 0
    for binx   in range(1, self.MC['EPD'].X_npix+1):
      for biny in range(1, self.MC['EPD'].Y_npix+1):
        H = self.DDe.GetBinContent(binx, biny)
        if H:
          F = self.EXY.Eval(self.MC['EPD'].X_beam + (binx-0.5)*self.MC['EPD'].X_pix, self.MC['EPD'].Y_beam + (biny-0.5)*self.MC['EPD'].Y_pix)
          self.DDe.SetBinContent(binx, biny, (F - H)/(1+abs(F))**0.5)
        else: self.eZeros += 1
    print ('NeZero=', self.eZeros)
    self.DDe.SetTitle('Electrons: (Fit - MC)/(1+Fit)^{1/2}')

  def PhotonsResiduals(self):
    self.pZeros = 0
    for binx   in range(1, self.MC['PPD'].X_npix//self.grebin+1):
      for biny in range(1, self.MC['PPD'].Y_npix//self.grebin+1):
        H = self.DDp.GetBinContent(binx, biny)
        if H:
          F = self.PXY.Eval(self.MC['PPD'].X_beam + (binx-0.5)*self.MC['PPD'].X_pix*self.grebin, self.MC['PPD'].Y_beam + (biny-0.5)*self.MC['PPD'].Y_pix*self.grebin)
          self.DDp.SetBinContent(binx, biny, (F - H)/(1+abs(F))**0.5)
        else: self.pZeros += 1
    self.DDp.SetTitle('Photons: (Fit - MC)/(1+Fit)^{1/2}')

  def ElectronsResults(self):
    X0 = 1000*self.MC['Spectrometer'].θo*self.MC['Spectrometer'].spec_L
    X2 = X0*self.MC['Spectrometer'].κ
    print ('expectation: X0 = {:7.5f} mm'.format(-X0))
    print ('expectation: X2 = {:7.5f} mm'.format(X2))

    chi2 = self.FitElectronsResult.Chi2(); print('Chi2: %f'        % (chi2))
    NDF  = self.FitElectronsResult.Ndf();  print('NDF: %d'         % (NDF) )
    prob = self.FitElectronsResult.Prob(); print('Probability: %f' % (prob))
    NDF -= self.eZeros;                    prob = ROOT.TMath.Prob(chi2, NDF)
    X0   = self.PXY.GetParameter(0);       dX0  = self.PXY.GetParError(0) 
#    X0   = self.HDp.GetMean(1);            dX0 = self.HDp.GetMeanError(1)
    X1   = self.EXY.GetParameter(0);       dX1 = self.EXY.GetParError(0) 
    X2   = self.EXY.GetParameter(1);       dX2 = self.EXY.GetParError(1)
    Y1   = self.EXY.GetParameter(2);       dY1 = self.EXY.GetParError(2)
    Y2   = self.EXY.GetParameter(3);       dY2 = self.EXY.GetParError(3)
    A    = X2-X1;                          dA  = (dX1**2+dX2**2)**0.5
    B    = Y2-Y1;                          dB  = (dY1**2+dY2**2)**0.5
    print('A = %.4f +/- %.4f (%.5f)' % (A, dA, dA/A))
    Δ    = 1e+3*self.MC['Laser'].ωo/self.MC['Spectrometer'].θo/self.MC['Spectrometer'].Eo*self.MC['Spectrometer'].leip_L # mm
    print ('Δ = {:7.5f} mm'.format(Δ))
    Δ    = 250.*self.MC['Spectrometer'].κ/self.MC['Spectrometer'].θo/self.MC['Spectrometer'].γ**2*self.MC['Spectrometer'].leip_L # mm
    print ('Δ = {:7.5f} mm'.format(Δ))
    print ('Beam X = {:7.5f} ± {:7.5f} mm'.format(X1+Δ, dX1))
    R    = (X2-X1-2*Δ)/(X1-X0+Δ)
    dR   = R * ( (dX0/(X1-X0))**2 + (dX1*(X0-X2)/(X1-X0)/(X2-X1))**2 + (dX2/(X2-X1))**2 )**0.5
    Eo   = 0.25*Constants.me**2/self.MC['Laser'].ωo * R
    dEo  = 0.25*Constants.me**2/self.MC['Laser'].ωo * dR
    print(R,  dR)
    print('E = %.5f +/- %.5f GeV' % (Eo*1.e-9, dEo*1.e-9))
    print('Y1 = %.5f +/- %.5f mm' % (Y1, dY1))
    print('Y2 = %.5f +/- %.5f mm' % (Y2, dY2))
    print(0.25e-3*(Y2-Y1)/self.MC['Spectrometer'].leip_L*Constants.me/self.MC['Laser'].ωo)
    self.ElectronsTable = ROOT.TPaveText(.05, .05, .95, .96)
    self.ElectronsTable.AddText(cpuinfo())
    self.ElectronsTable.AddText('Electrons fit: t = {:.0f} s (CPU {:.0f} s)'.format(ROOT.gBenchmark.GetRealTime('EFit'), ROOT.gBenchmark.GetCpuTime('EFit')))
    self.ElectronsTable.AddText('#chi^{2}/NDF = %.1f/%d | Prob = %.4f' % (chi2,NDF,prob))
    self.ElectronsTable.AddText('X_{1} = %08.3f #pm %5.3f mm' % (X1, dX1))
    self.ElectronsTable.AddText('X_{2} = %08.3f #pm %5.3f mm' % (X2, dX2))
    self.ElectronsTable.AddText('#xi_{1} = %05.3f #pm %5.3f'  % (self.EXY.GetParameter(4), self.EXY.GetParError(4)))
    self.ElectronsTable.AddText('#xi_{2} = %05.3f #pm %5.3f'  % (self.EXY.GetParameter(5), self.EXY.GetParError(5)))
    self.ElectronsTable.AddText('#xi_{3}#zeta_{x} = %05.3f #pm %5.3f' % (self.EXY.GetParameter(6), self.EXY.GetParError(6)))
    self.ElectronsTable.AddText('#xi_{3}#zeta_{y} = %05.3f #pm %5.3f' % (self.EXY.GetParameter(7), self.EXY.GetParError(7)))
    self.ElectronsTable.AddText('#xi_{3}#zeta_{z} = %05.3f #pm %5.3f' % (self.EXY.GetParameter(8), self.EXY.GetParError(8)))
    self.ElectronsTable.AddText('#sigma_{x} = %5.1f #pm %4.1f #mum' % (500*self.EXY.GetParameter(9)*A, 500*self.EXY.GetParError(9)*A))
    self.ElectronsTable.AddText('#sigma_{y} = %5.2f #pm %4.2f #mum' % (500*self.EXY.GetParameter(10)*B, 500*self.EXY.GetParError(10)*B))
    self.ElectronsTable.AddText('E_{beam} = %7.4f #pm %6.4f GeV. ' % (Eo*1.e-9, dEo*1.e-9))

  def PhotonsResults(self):
    chi2 = self.FitPhotonsResult.Chi2(); print('Chi2: %f'        % (chi2))
    NDF  = self.FitPhotonsResult.Ndf();  print('NDF: %d'         % (NDF) )
    prob = self.FitPhotonsResult.Prob(); print('Probability: %f' % (prob))
    NDF -= self.pZeros;                  prob = ROOT.TMath.Prob(chi2, NDF)
    X0   = self.PXY.GetParameter(0);     dX0  = self.PXY.GetParError(0) 
    Y0   = self.PXY.GetParameter(1);     dY0  = self.PXY.GetParError(1)
    self.PhotonsTable = ROOT.TPaveText(.05, .05, .95, .96)
    self.PhotonsTable.AddText(cpuinfo())
    self.PhotonsTable.AddText('Photons fit: t = {:.0f} s (CPU {:.0f} s)'.format(ROOT.gBenchmark.GetRealTime('PFit'), ROOT.gBenchmark.GetCpuTime('PFit')))
    self.PhotonsTable.AddText('#chi^{2}/NDF = %.1f/%d | Prob = %.4f' % (chi2,NDF,prob))
    self.PhotonsTable.AddText('X_{0} = %08.3f #pm %5.3f mm' % (X0, dX0))
    self.PhotonsTable.AddText('#xi_{1} = %05.3f #pm %5.3f'  % (self.PXY.GetParameter(2), self.PXY.GetParError(2)))
    self.PhotonsTable.AddText('#xi_{2} = %05.3f #pm %5.3f'  % (self.PXY.GetParameter(3), self.PXY.GetParError(3)))
    self.PhotonsTable.AddText('#xi_{3}#zeta_{x} = %05.3f #pm %5.3f'  % (self.PXY.GetParameter(4), self.PXY.GetParError(4)))
    self.PhotonsTable.AddText('#xi_{3}#zeta_{y} = %05.3f #pm %5.3f'  % (self.PXY.GetParameter(5), self.PXY.GetParError(5)))
    self.PhotonsTable.AddText('#xi_{3}#zeta_{z} = %05.3f #pm %5.3f'  % (self.PXY.GetParameter(6), self.PXY.GetParError(6)))
    self.PhotonsTable.AddText('#sigma_{x} = %5.1f #pm %4.1f um'      % (1000*self.PXY.GetParameter(7), 1000*self.PXY.GetParError(7)))
    self.PhotonsTable.AddText('#sigma_{y} = %5.2f #pm %4.2f um'      % (1000*self.PXY.GetParameter(8), 1000*self.PXY.GetParError(8)))


class DISPLAY:
  def __init__(self):
    ROOT.gStyle.SetOptFit(1111);    ROOT.gStyle.SetOptStat('ne');    ROOT.gStyle.SetFitFormat("8.3g")
    ROOT.gROOT.SetStyle("Plain");   ROOT.gROOT.ForceStyle()      #   ROOT.gStyle.SetPalette(56)
    self.cv = ROOT.TCanvas('cv','cv',0,0,1600,1200);               self.cv.Divide(3,3)

  def ShowOnPad(self, nPad, entity, grid=False, goption=''):
    self.cv.cd(nPad) 
    if grid: self.cv.GetPad(nPad).SetGrid()
    entity.Draw(goption)
    self.cv.Modified();    self.cv.Update()


def main(argv):
  parser = argparse.ArgumentParser(description='Process cli arguments.')
  parser.add_argument('--fit', action='store_const', const=True, default=False, help = 'Fit the data. Default is NO.')
  parser.add_argument('filename', type=str, help='Name of the datafile.')
  args = parser.parse_args()

  DATA = POLARIMETER(filename=args.filename, grebin=1)
  LOOK = DISPLAY()
  DATA.ParametersMC()
  LOOK.ShowOnPad(nPad=1, entity = DATA.HDp, grid = True, goption='COLZ')
  LOOK.ShowOnPad(nPad=2, entity = DATA.HDe, grid = True, goption='COLZ')
  LOOK.ShowOnPad(nPad=3, entity = DATA.ParametersTable)

  ROOT.gBenchmark.Start('PFit')
  DATA.PhotonsSetFunction()
  if args.fit:  psuccess = DATA.FitPhotons()
  else:         psuccess = False
  DATA.PhotonsResiduals()
  ROOT.gBenchmark.Stop('PFit')
  LOOK.ShowOnPad(nPad=4, entity = DATA.PXY, grid = True, goption='COLZ1')
  LOOK.ShowOnPad(nPad=5, entity = DATA.DDp, grid = True, goption='COLZ')
  if psuccess: 
    DATA.PhotonsResults()
    LOOK.ShowOnPad(nPad=6, entity = DATA.PhotonsTable)
  input('Fit electrons ?')
  ROOT.gBenchmark.Start('EFit')
  DATA.ElectronsSetFunction()
  if args.fit:  esuccess = DATA.FitElectrons()
  else:         esuccess = False
  DATA.ElectronsResiduals()
  ROOT.gBenchmark.Stop('EFit')

  LOOK.ShowOnPad(nPad=7, entity = DATA.EXY, grid = True, goption='COLZ')
  LOOK.ShowOnPad(nPad=8, entity = DATA.DDe, grid = True, goption='COLZ1')
  if psuccess and esuccess: 
    DATA.ElectronsResults()
    LOOK.ShowOnPad(nPad=9, entity = DATA.ElectronsTable)

  input()
  exit()

if __name__ == "__main__": main(sys.argv)

