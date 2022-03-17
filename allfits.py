#!/usr/bin/env python3
import sys, ROOT, argparse
from hardware import Constants, Laser, EPD, PPD
from xSection import pPIXELS, ePIXELS


def cpuinfo():
  with open('/proc/cpuinfo', 'r') as f: info = f.readlines()
  model = [el for el in info if 'model name' in el]
  return model[0].strip().split(': ')[1]


class POLARIMETER:
  DataMC = {}
  def __init__(self, filename=None):
    beam  = {
            'test'               :'Tue Feb 22 16:45:20 2022.root',
            'mpol'               :'Wed Mar 16 13:55:28 2022.root',
            'polm'               :'Wed Mar 16 14:20:16 2022.root',
            }
    if filename: hfile = ROOT.TFile(filename)
    else:        hfile = ROOT.TFile(beam['test'])
    name  = hfile.GetListOfKeys()[0].GetName()
    XYD   = hfile.Get(name) #    XYD   = ROOT.TList()
    hfile.Close()
    for i in range(11): self.DataMC[['λo','Eo','κ','γ','θo','ξ1','ξ2','ξ3','ζx','ζy','ζz'][i]] = [float(p) for p in name.split()][i]
    self.HDp, self.DDp = XYD[0], XYD[0].Clone()
    self.HDe, self.DDe = XYD[1], XYD[1].Clone()
    self.HDe.SetStats(0); self.HDe.SetTitle('Electrons: MC')
    self.HDp.SetStats(0); self.HDp.SetTitle('Photons: MC')
    self.DDe.SetStats(0); self.DDe.SetTitle('Electrons: residuald')
    self.DDp.SetStats(0); self.DDp.SetTitle('Photons: residuals')

  def ParametersMC(self):
#  'λo','Eo','κ','γ','θo','ξ1','ξ2','ξ3','ζx','ζy','ζz'
    self.ParametersTable = ROOT.TPaveText(.05, .05, .95, .96)
    self.ParametersTable.AddText('Monte-Carlo Parameters:')
    self.ParametersTable.AddText('Electron E_{0} = %6.3f GeV' %         (1e-9*self.DataMC['Eo']))
    self.ParametersTable.AddText('Laser #lambda_{0} = %5.3f um' % (1e+4*self.DataMC['λo']))
    self.ParametersTable.AddText('Electron #gamma = %5.3f#times10^{3}'% (1e-3*self.DataMC['γ']))
    self.ParametersTable.AddText('Compton #kappa = %5.3f'      %       (     self.DataMC['κ']))
    self.ParametersTable.AddText('Bend: #gamma#theta_{0} = %5.3f'  % (     self.DataMC['θo']))
    self.ParametersTable.AddText('(#xi_{1}, #xi_{2}, #xi_{3}) = (%6.3f, %6.3f, %6.3f)'%(self.DataMC['ξ1'],self.DataMC['ξ2'],self.DataMC['ξ3']))
    self.ParametersTable.AddText('(#zeta_{x}, #zeta_{y}, #zeta_{z}) = (%6.3f, %6.3f, %6.3f)'%(self.DataMC['ζx'],self.DataMC['ζy'],self.DataMC['ζz']))

  def ElectronsSetFunction(self):
    self.efit = ePIXELS()
    X0 = EPD.X_beam; X1 = X0 + EPD.X_size
    Y0 = EPD.Y_beam; Y1 = Y0 + EPD.Y_size
    self.EXY = ROOT.TF2('EXY', self.efit , X0, X1, Y0, Y1, 10)
    self.EXY.SetParName(0, 'X1');        self.EXY.SetParameter(0,  -0.01)  # beam position x, mm
    self.EXY.SetParName(1, 'X2');        self.EXY.SetParameter(1,  347.5)  # edge position x, mm
    self.EXY.SetParName(2, 'Y0');        self.EXY.SetParameter(2, -1.068)  # y_min, mm
    self.EXY.SetParName(3, 'Y1');        self.EXY.SetParameter(3,  1.068)  # y_max, mm
    self.EXY.SetParName(4, 'ξ1');        self.EXY.SetParameter(4,    0.0); self.EXY.SetParLimits(4,  -1.0, 1.0)
    self.EXY.SetParName(5, 'ξ3*ζy');     self.EXY.SetParameter(5,    0.0); self.EXY.SetParLimits(5,  -1.0, 1.0)
    self.EXY.SetParName(6, 'ξ3*ζz');     self.EXY.SetParameter(6,    0.0); self.EXY.SetParLimits(6,  -1.0, 1.0)
    self.EXY.SetParName(7, 'σx');        self.EXY.SetParameter(7, 1.2e-3); self.EXY.SetParLimits(7, 0.0001, 0.1) # σ_x [a.u.]
    self.EXY.SetParName(8, 'σy');        self.EXY.SetParameter(8, 2.6e-2); self.EXY.SetParLimits(8, 0.0001, 0.1) # σ_y [a.u.]
    self.EXY.SetParName(9,' norm');      self.EXY.SetParameter(9, 1.2e+2)  # amplitude
    self.EXY.SetNpx(EPD.X_npix);         self.EXY.SetNpy(EPD.Y_npix)
    self.EXY.SetTitle('Electrons: Fit'); self.EXY.GetXaxis().SetTitle('X, mm')

  def FitElectrons(self):
    fixed_parameters = []
    for p in fixed_parameters: self.EXY.FixParameter(p, self.EXY.GetParameter(p))
    self.FitElectronsResult =  self.HDe.Fit(self.EXY, 'SVNP') # Use Pearsons chi-square method
    for p in fixed_parameters: self.EXY.ReleaseParameter(p)
    return not self.FitElectronsResult.Status()

  def PhotonsSetFunction(self):
    self.pfit = pPIXELS()
    X0 = PPD.X_beam; X1 = X0 + PPD.X_size 
    Y0 = PPD.Y_beam; Y1 = Y0 + PPD.Y_size
    self.PXY = ROOT.TF2('PXY', self.pfit , X0, X1, Y0, Y1, 10)
    self.PXY.SetParName(0, 'X0');      self.PXY.SetParameter(0,    self.HDp.GetMean(1))
    self.PXY.SetParName(1, 'Y0');      self.PXY.SetParameter(1,    self.HDp.GetMean(2))
    self.PXY.SetParName(2, 'ξ1');      self.PXY.SetParameter(2,    0.0); self.PXY.SetParLimits(2,  -1.0, 1.0)
    self.PXY.SetParName(3, 'ξ2');      self.PXY.SetParameter(3,    0.0); self.PXY.SetParLimits(3,  -1.0, 1.0)
    self.PXY.SetParName(4, 'ξ3*ζx');   self.PXY.SetParameter(4,    0.0); self.PXY.SetParLimits(4,  -1.0, 1.0)
    self.PXY.SetParName(5, 'ξ3*ζy');   self.PXY.SetParameter(5,    0.0); self.PXY.SetParLimits(5,  -1.0, 1.0)
    self.PXY.SetParName(6, 'ξ3*ζz');   self.PXY.SetParameter(6,    0.0); self.PXY.SetParLimits(6,  -1.0, 1.0)
    self.PXY.SetParName(7, 'σx');      self.PXY.SetParameter(7,    0.2); self.PXY.SetParLimits(7, 0.001, 1.0) # σ_x [mm]
    self.PXY.SetParName(8, 'σy');      self.PXY.SetParameter(8,   0.04); self.PXY.SetParLimits(8, 0.001, 1.0) # σ_y [mm]
    self.PXY.SetParName(9,' norm');    self.PXY.SetParameter(9,   9e+3)  # amplitude
    self.PXY.SetNpx(PPD.X_npix);       self.PXY.SetNpy(PPD.Y_npix)
    self.PXY.SetTitle('Photons: Fit'); self.PXY.GetXaxis().SetTitle('X, mm')

  def FitPhotons(self):
    fixed_parameters = []
    for p in fixed_parameters: self.PXY.FixParameter(p, self.PXY.GetParameter(p))
    self.FitPhotonsResult =  self.HDp.Fit(self.PXY, 'SVNP') # Use Pearsons chi-square method
    for p in fixed_parameters: self.PXY.ReleaseParameter(p)
    return not self.FitPhotonsResult.Status()

  def ElectronsResiduals(self):
    self.eZeros = 0
    for binx   in range(1, EPD.X_npix+1):
      for biny in range(1, EPD.Y_npix+1):
        F = self.EXY.Eval(EPD.X_beam + (binx-0.5)*EPD.X_pix, EPD.Y_beam + (biny-0.5)*EPD.Y_pix)
        if F>0.01:
          H = self.DDe.GetBinContent(binx, biny)
          self.DDe.SetBinContent(binx, biny, (F - H)/(1+abs(F))**0.5)
        else: self.eZeros += 1
    print ('NeZero=', self.eZeros)
    self.DDe.SetTitle('Electrons: (Fit - MC)/(1+Fit)^{1/2}')

  def PhotonsResiduals(self):
    self.pZeros = 0
    for binx   in range(1, PPD.X_npix+1):
      for biny in range(1, PPD.Y_npix+1):
        H = self.DDp.GetBinContent(binx, biny)
        if H:
          F = self.PXY.Eval(PPD.X_beam + (binx-0.5)*PPD.X_pix, PPD.Y_beam + (biny-0.5)*PPD.Y_pix)
          self.DDp.SetBinContent(binx, biny, (F - H)/(1+abs(F))**0.5)
        else: self.pZeros += 1
    self.DDp.SetTitle('Photons: (Fit - MC)/(1+Fit)^{1/2}')

  def ElectronsResults(self):
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
    R   = (X2-X1)/(X1-X0)
    R  /= 1.0 + 0.5/self.DataMC['θo']**2 # Correction !!!
    dR  = R * ( (dX0/(X1-X0))**2 + (dX1*(X0-X2)/(X1-X0)/(X2-X1))**2 + (dX2/(X2-X1))**2 )**0.5
    Eo  = 0.25*Constants.me**2/Laser.ωo * R 
    dEo = 0.25*Constants.me**2/Laser.ωo * dR
    print(R,  dR)
    print('E = %.5f +/- %.5f' % (Eo*1.e-9, dEo*1.e-9))
    self.ElectronsTable = ROOT.TPaveText(.05, .05, .95, .96)
    self.ElectronsTable.AddText(cpuinfo())
    self.ElectronsTable.AddText('Electrons fit: t = {:.0f} s (CPU {:.0f} s)'.format(ROOT.gBenchmark.GetRealTime('EFit'), ROOT.gBenchmark.GetCpuTime('EFit')))
    self.ElectronsTable.AddText('#chi^{2}/NDF = %.1f/%d | Prob = %.4f' % (chi2,NDF,prob))
    self.ElectronsTable.AddText('X_{1} = %08.3f #pm %5.3f mm' % (X1, dX1))
    self.ElectronsTable.AddText('X_{2} = %08.3f #pm %5.3f mm' % (X2, dX2))
    self.ElectronsTable.AddText('#xi_{1} = %05.3f #pm %5.3f'  % (self.EXY.GetParameter(4), self.EXY.GetParError(4)))
    self.ElectronsTable.AddText('#xi_{3}#zeta_{y} = %05.3f #pm %5.3f' % (self.EXY.GetParameter(5), self.EXY.GetParError(5)))
    self.ElectronsTable.AddText('#xi_{3}#zeta_{z} = %05.3f #pm %5.3f' % (self.EXY.GetParameter(6), self.EXY.GetParError(6)))
    self.ElectronsTable.AddText('#sigma_{x} = %5.1f #pm %4.1f #mum' % (500*self.EXY.GetParameter(7)*A, 500*self.EXY.GetParError(7)*A))
    self.ElectronsTable.AddText('#sigma_{y} = %5.2f #pm %4.2f #mum' % (500*self.EXY.GetParameter(8)*B, 500*self.EXY.GetParError(8)*B))
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

  DATA = POLARIMETER(filename=args.filename)
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
  LOOK.ShowOnPad(nPad=5, entity = DATA.DDp, grid = True, goption='COLZ1')
  if psuccess: 
    DATA.PhotonsResults()
    LOOK.ShowOnPad(nPad=6, entity = DATA.PhotonsTable)

  ROOT.gBenchmark.Start('EFit')
  DATA.ElectronsSetFunction()
  if args.fit:  esuccess = DATA.FitElectrons()
  else:         esuccess = False
  DATA.ElectronsResiduals()
  ROOT.gBenchmark.Stop('EFit')

  LOOK.ShowOnPad(nPad=7, entity = DATA.EXY, grid = True, goption='COLZ1')
  LOOK.ShowOnPad(nPad=8, entity = DATA.DDe, grid = True, goption='COLZ1')
  if psuccess and esuccess: 
    DATA.ElectronsResults()
    LOOK.ShowOnPad(nPad=9, entity = DATA.ElectronsTable)

  input()
  exit()

if __name__ == "__main__": main(sys.argv)

