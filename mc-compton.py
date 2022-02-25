#!/usr/bin/env python3
# Α Β Γ Δ Ε Ζ Η Θ Ι Κ Λ Μ Ν Ξ Ο Π Ρ Σ Τ Υ Φ Χ Ψ Ω
# α β γ δ ε ζ η θ ι κ λ μ ν ξ ο π ρ σ τ υ φ χ ψ ω

import ROOT, sys
from ctypes import c_double
from time import asctime
from hardware import EPD, PPD, Laser, Spectrometer
from xSection import XSMC

class MonteCarlo(XSMC):
  def __init__(self):
    cDouble = c_double*1
    ROOT.gRandom.SetSeed()
    XSMC.__init__(self)
    self.rannor    = ROOT.gRandom.Rannor
    self.u, self.φ = cDouble(), cDouble()
    self.x, self.y = cDouble(), cDouble()

  def game(self):
    self.U.GetRandom2(self.u, self.φ)        # Get u,φ
    self.rannor(      self.x, self.y)        # Get x,y from normal distribution
    θex = self.x[0] * Spectrometer.nsx       # electron betatron x-angle [1/γ]
    θey = self.y[0] * Spectrometer.nsy       # electron betatron y-angle [1/γ]
    sin = ROOT.TMath.Sin(self.φ[0])
    cos = ROOT.TMath.Cos(self.φ[0])
    θp  = (Spectrometer.κ/self.u[0]-1)**0.5  # photon   scattering angle [1/γ]
    θe  = θp * self.u[0]                     # electron scattering angle [1/γ]
    θpx = θex + θp*cos                       # photon            x-angle [1/γ]
    θpy = θey + θp*sin                       # photon            y-angle [1/γ]
    θex = θex - θe*cos                       # electron          x-angle [1/γ]
    θey = θey - θe*sin                       # electron          y-angle [1/γ]
    return self.u[0], θpx, θpy, θex, θey


def main(argv):
  S      = Spectrometer()
  print('Laser ωo:                %.3f eV'   % (Laser.ωo))
  print('Comptony κ-parameter:    %.3f'      % (S.κ))
  print('beam γ-factor:           %.3f'      % (S.γ))
  print('beam bending angle:      %.3f mrad' % (1000*S.bend))
  print('beam bending angle:      %.3f 1/γ'  % (S.no))
  print('γ-to-beam distance :     %.3f mm'   % (1000   * S.bend*S.spec_L))
  print('emin-to-beam distance :  %.3f mm'   % (1000*S.κ*S.bend*S.spec_L))
  print('horizontal spread:       %.3f mm'   % (1000      *S.sx*S.leip_L))
  print('vertical spread:         %.3f mm'   % (1000      *S.sy*S.leip_L))
  input('Continue, OK?')

  Xemin, Xemax = EPD.X_beam, EPD.X_beam + EPD.X_size
  Yemin, Yemax = EPD.Y_beam, EPD.Y_beam + EPD.Y_size
  XYe = ROOT.TH2I('XYe', 'electrons XY', EPD.X_npix, Xemin, Xemax, EPD.Y_npix, Yemin, Yemax)
  Ye  = ROOT.TH1I('Ye',  'electrons  Y', EPD.Y_npix, Yemin, Yemax)
  XYe.GetXaxis().SetTitle('X, mm'); XYe.GetZaxis().SetTitle('e^{#pm}/pixel');
  XYe.GetYaxis().SetTitle('Y, mm'); Ye.GetXaxis().SetTitle('Y, mm')

  Xpmin, Xpmax = PPD.X_beam, PPD.X_beam + PPD.X_size
  Ypmin, Ypmax = PPD.Y_beam, PPD.Y_beam + PPD.Y_size
  XYp = ROOT.TH2I('XYp', 'photons XY', PPD.X_npix, Xpmin, Xpmax, PPD.Y_npix, Ypmin, Ypmax)
  Yp  = ROOT.TH1I('Ye',  'photons  Y', PPD.Y_npix, Ypmin, Ypmax)
  XYp.GetXaxis().SetTitle('X, mm');   XYp.GetYaxis().SetTitle('Y, mm'); Yp.GetXaxis().SetTitle('Y, mm')

  G = MonteCarlo()

  cv = ROOT.TCanvas('cv','cv',0,0,1600,1200)
  cv.Divide(2,2)
  for p in range(4): cv.GetPad(p+1).SetGrid()
  netomm = 1000*S.leip_L/S.γ # scattering angle in 1/γ to mm
  notomm = 1000*S.spec_L/S.γ #    bending angle in 1/γ to mm
  for i in range(100):
    for j in range(100000):
      u, txp, typ, txe, tye = G.game()
      xe, ye = netomm*txe + u*notomm*S.no, netomm*tye
      xp, yp = netomm*txp -   notomm*S.no, netomm*typ
      if xe > EPD.X_beam:
        XYe.Fill(xe,ye); Ye.Fill(ye)
      if Xpmin < xp < Xpmax and Ypmin < yp < Ypmax:
        XYp.Fill(xp,yp); Yp.Fill(yp)
    print('\rProgress:  %3d%%' % (i+1), end=' ')
    cv.cd(1);   XYp.Draw('COLZ'); cv.cd(2);   XYe.Draw('COLZ')
    cv.cd(3);    Yp.Draw();       cv.cd(4);    Ye.Draw()
    cv.Modified(); cv.Update()
    sys.stdout.flush()
  print()

  mean, rms, N = Yp.GetMean(), Yp.GetRMS(), Yp.Integral()
  print(mean, rms, N)
  print('Accuracy: %f percent' % (100*rms/mean/N**0.5))

  q = input('Save data? (Yes/no)')
  if not 'n' in q:
    T = '%.3e %.3e %.6f %.3e %.3f' % (Laser.λo, S.Eo, S.κ, S.γ, S.no)
    T+= ' %6.3f %6.3f %6.3f' % (Laser.ξ1, Laser.ξ2, Laser.ξ3)
    T+= ' %6.3f %6.3f %6.3f' % (    S.ζx,     S.ζy,     S.ζz)
    XYD = ROOT.TList(); XYD.Add(XYp); XYD.Add(XYe)
    fname = asctime()+'.root'
    f = ROOT.TFile(fname,'new'); f.WriteObject(XYD, T); f.Close()
    print ('Data saved to ' + fname)

  exit(0)


if __name__ == "__main__": main(sys.argv)
