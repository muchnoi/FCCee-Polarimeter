#!/usr/bin/env python3
# Α Β Γ Δ Ε Ζ Η Θ Ι Κ Λ Μ Ν Ξ Ο Π Ρ Σ Τ Υ Φ Χ Ψ Ω
# α β γ δ ε ζ η θ ι κ λ μ ν ξ ο π ρ σ τ υ φ χ ψ ω

import ROOT, sys, pickle
from ctypes import c_double
from time import asctime
from hardware import EPD, PPD, Laser, Spectrometer as Sptr
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
    self.U.GetRandom2(self.u, self.φ)         # get u,φ
    Δ    = ROOT.gRandom.Gaus(0.0, Sptr.σE)    # electron energy deviation [a.u.]
    x    = Δ*Sptr.Dx                          # synchrotron x-coordinate [mm]
    self.rannor(      self.x, self.y)         # get x,y from normal distribution
    x   += self.x[0] * Sptr.σx                # electron betatron x-coordinate [mm]
    y    = self.y[0] * Sptr.σy                # electron betatron y-coordinate [mm]
    self.rannor(      self.x, self.y)         # Get x,y from normal distribution
    θex  = self.x[0] * Sptr.ηx                # electron betatron x-angle [rad]
    θey  = self.y[0] * Sptr.ηy                # electron betatron y-angle [rad]
    sin  = ROOT.TMath.Sin(self.φ[0])
    cos  = ROOT.TMath.Cos(self.φ[0])
    θp   = (Sptr.κ/self.u[0]-1)**0.5/Sptr.γ   # photon   scattering   angle [rad]
    θpx  = θex + θp*cos/(1+Δ)                 # photon   scattering x-angle [rad]
    θpy  = θey + θp*sin/(1+Δ)                 # photon   scattering y-angle [rad]
    θe   = θp * self.u[0]                     # electron scattering   angle [rad]
    θex  = θex - θe*cos                       # electron scattering x-angle [rad]
    θey  = θey - θe*sin                       # electron scattering y-angle [rad]
    uΔ   = self.u[0] - Δ/(1+Δ)                # electron bending angle / θo 
    return x, y, θpx, θpy, θex, θey, uΔ


def main(argv):
  S      = Sptr()
  print('Laser ωo:                %.3f eV'   % (Laser.ωo))
  print('Comptony κ-parameter:    %.4f'      % (S.κ))
  print('beam γ-factor:           %.3f'      % (S.γ))
  print('beam bending angle:      %.3f mrad' % (1000*S.θo))
  print('beam bending angle:      %.3f 1/γ'  % (S.γ*S.θo))
  print('γ-to-beam distance :     %.3f mm'   % (1000    *S.θo*S.spec_L))
  print('emin-to-beam distance :  %.3f mm'   % (1000*S.κ*S.θo*S.spec_L))
  print('horizontal size:         %.3f mm'   % (S.σx))
  print('vertical size:           %.3f mm'   % (S.σy))
  print('horizontal spread add:   %.3f mm'   % (1000*S.ηx*S.leip_L))
  print('vertical spread add:     %.3f mm'   % (1000*S.ηy*S.leip_L))
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
  L1     = 1000*S.leip_L                         # scattering base [mm]
  L2     = 1000*S.spec_L*ROOT.TMath.Tan(S.θo)    # tan(θo) * bending base [mm]
  icos   = 1./ROOT.TMath.Cos(S.θo)
  for i in range(100):
    for j in range(100000):
      x, y, θpx, θpy, θex, θey, uΔ = G.game()
      xp = icos*(x + L1*ROOT.TMath.Tan(θpx)) - L2
      yp =       y + L1*ROOT.TMath.Tan(θpy)
      xe = icos*(x + L1*ROOT.TMath.Tan(θex)) + L2*uΔ
      ye =       y + L1*ROOT.TMath.Tan(θey)
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

  if not ('n' or 'N') in input('Save data? (Yes/no)'):
    SAVE = {};         SAVE['Laser'] = Laser;       SAVE['Spectrometer'] = Sptr
    SAVE['EPD'] = EPD; SAVE['PPD'] = PPD; SAVE['XYe'] = XYe; SAVE['XYp'] = XYp
    fname = asctime().replace(' ','_').replace(':','') + '.MC'
    with open(fname, 'wb') as fp:  pickle.dump(SAVE, fp, -1)
    print('Data saved to %s' % fname)
  exit(0)


if __name__ == "__main__": main(sys.argv)
