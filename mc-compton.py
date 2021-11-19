#!/usr/bin/env python3
# Α Β Γ Δ Ε Ζ Η Θ Ι Κ Λ Μ Ν Ξ Ο Π Ρ Σ Τ Υ Φ Χ Ψ Ω
# α β γ δ ε ζ η θ ι κ λ μ ν ξ ο π ρ σ τ υ φ χ ψ ω

import ROOT, sys
from ctypes import c_double
from math import pi, cos, sin
from time import asctime
from hardware import Pixels, Laser, Spectrometer

class MonteCarlo(Spectrometer):
  def __init__(self):
    ROOT.gRandom.SetSeed()
    self.rannor  = ROOT.gRandom.Rannor
    CCS = '([0]*(1+(1+x)**2) - 4*x/[0]*(1+x)*([0]-x)*(1-[1]*cos(2*(y-[2]))) + \
            [3]*x*(x+2)*([0]-2*x) - \
            [4]*2*x*sqrt(x*([0]-x))*sin(y)) / [0]**2/(1+x)**3'
    self.U = ROOT.TF2('U', CCS)
    self.U.SetRange( 0, 0, Spectrometer.k, 2*pi)
    self.U.SetParameter(0, Spectrometer.k)
    self.U.SetParameter(1, Spectrometer.pp1)
    self.U.SetParameter(2, Spectrometer.pp2)
    self.U.SetParameter(3, Spectrometer.pp3)
    self.U.SetParameter(4, Spectrometer.pp4)
    self.U.SetNpx(1000)
    self.U.SetNpy(1000)
    cDouble = c_double*1
    self.u, self.phi = cDouble(), cDouble()
    self.bx, self.by = cDouble(), cDouble()

  def game(self):
    self.U.GetRandom2(self.u,  self.phi)
    self.rannor(      self.bx, self.by)
    ex = self.bx[0] * Spectrometer.nsx
    ey = self.by[0] * Spectrometer.nsy
    sinus, cosinus = sin(self.phi[0]), cos(self.phi[0])
    np  =  (Spectrometer.k/self.u[0]-1.)**0.5
    ne  =  np * self.u[0]
    npx =  np*cosinus + ex
    npy =  np*sinus   + ey
    nex = -ne*cosinus + ex
    ney = -ne*sinus   + ey
    return self.u[0], npx, npy, nex, ney


def main(argv):
  S      = Spectrometer()
  print('Comptony κ-parameter:    %.3f'      % (S.k))
  print('beam γ-factor:           %.3f'      % (S.gf))
  print('beam bending angle:      %.3f mrad' % (1000*S.bend))
  print('beam bending angle:      %.3f 1/γ'  % (S.no))
  print('γ-to-beam distance :     %.3f mm'   % (1000   * S.bend*S.spec_L))
  print('emin-to-beam distance :  %.3f mm'   % (1000*S.k*S.bend*S.spec_L))
  print('horizontal spread:       %.3f mm'   % (1000      *S.sx*S.leip_L))
  print('vertical spread:         %.3f mm'   % (1000      *S.sy*S.leip_L))
  input('Continue, OK?')

  NameTitle = 'E_{0} = %.1f GeV, #lambda_{0} = %.1f nm, #kappa = %5.3f, P_{#perp}  =%.2f' % (S.Eo*1e-9, Laser.lo*1.e+7, S.k, S.pp4)

  D = Pixels()
  Xemin, Xemax = D.X_beam, D.X_beam + D.X_size; nbinsx = Pixels.X_npix
  Yemin, Yemax = D.Y_beam, D.Y_beam + D.Y_size; nbinsy = Pixels.Y_npix
  XYe = ROOT.TH2I('XYe', 'electrons ' + NameTitle, nbinsx, Xemin, Xemax, nbinsy, Yemin, Yemax)
  Ye  = ROOT.TH1I('Ye',  'electrons Y', nbinsy, Yemin, Yemax)
  XYe.GetXaxis().SetTitle('X, mm'); XYe.GetZaxis().SetTitle('e^{#pm}/pixel');
  XYe.GetYaxis().SetTitle('Y, mm'); Ye.GetXaxis().SetTitle('Y, mm')

  Xpcen = - 1000*S.bend*S.spec_L
  Xpmin, Xpmax = Xpcen - 5.0, Xpcen + 5.0
  Ypmin, Ypmax = -5.0, 5.0
  XYp = ROOT.TH2I('XYp', ' photons ' + NameTitle, nbinsx, Xpmin, Xpmax, nbinsy, Ypmin, Ypmax)
  Yp  = ROOT.TH1I('Ye',  ' photons Y', nbinsy, Ypmin, Ypmax)
  XYp.GetXaxis().SetTitle('X, mm');   XYp.GetYaxis().SetTitle('Y, mm'); Yp.GetXaxis().SetTitle('Y, mm')

  G = MonteCarlo()

  cv = ROOT.TCanvas('cv','cv',0,0,1600,1200)
  cv.Divide(2,2)
  for p in range(4): cv.GetPad(p+1).SetGrid()
  netomm = 1000*S.leip_L/S.gf # scattering angle in 1/gf to mm
  notomm = 1000*S.spec_L/S.gf #    bending angle in 1/gf to mm
  for i in range(100):
    for j in range(200000):
      u, txp, typ, txe, tye = G.game()
      xe, ye = netomm*txe + u*notomm*S.no, netomm*tye
      xp, yp = netomm*txp -   notomm*S.no, netomm*typ
      if xe> 40.0:
        XYe.Fill(xe,ye); Ye.Fill(ye)
      if Xpmin<xp<Xpmax and Ypmin<yp<Ypmax:
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

  T = '%.3e %.3e %.6f %.1e %.3f %.3f %.3f' % (Laser.lo, S.Eo, S.k, S.gf, S.no, S.sx, S.sy)
  XYD = ROOT.TList(); XYD.Add(XYp); XYD.Add(XYe)
  f = ROOT.TFile(asctime()+'.root','new'); f.WriteObject(XYD, T); f.Close()

  input()

  exit(0)


if __name__ == "__main__": main(sys.argv)
