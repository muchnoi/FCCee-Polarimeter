# -*- coding: utf-8 -*-
import math

class Constants:
  c      = 2.99792458e+10   # speed of light, cm/s
  me     = 0.510998928e+6   # electron rest energy, eV/c^2
  hbar   = 6.58211928e-16   # Plank constant, eV*s
  re     = 2.8179403267e-13 # classical electron radius, cm
  barn   = 1.0e-24          # cm

class Pixels:
  X_beam = 40.0   # beam-detector horizontal space, mm
  Y_beam = -2.0   # beam-detector  vertical  space, mm
  X_size = 400.0  # detector X-size (horizontal)  , mm
  Y_size = 4.0    # detector Y-size   (vertical)  , mm
  X_npix = 1600            # number of pixels in X
  Y_npix = 80              # number of pixels in Y
  X_pix  = X_size/float(X_npix)   # one pixel in X, mm
  Y_pix  = Y_size/float(Y_npix)   # one pixel in Y, mm

class Laser:
  lo     = 0.532e-4                                # λo, cm
  wo     = 2*math.pi*Constants.hbar*Constants.c/lo # ωo, eV

class Spectrometer:
  Eo     = 45.6e+9                     # electron (beam) energy, eV
  gf     = Eo/Constants.me             # electron Lorentz factor
  k      = 4*gf*Laser.wo/Constants.me  # 2 * rest frame momentum of laser photon
  bend   = .0021341                    # bending angle, rad
  bend_l = 24.12                       # dipole length, m
  spec_L = 0.5*bend_l + 88.0           # spectrometer arc length, m
  leip_L = 0.5*bend_l + spec_L + 5.0   # laser electron i.p. distance, m
  no     = gf*bend                     # bending angle in units 1/gf
  ex     = 0.27e-9                     # emittance x, m*rad
  ey     = 1.0e-12                     # emittance y, m*rad
  Bx     = 100.0                       # βx, m
  By     = 25.00                       # βy, m
  sx     = (ex/Bx)**0.5                # betatron energy spread in x, rad
  sy     = (ey/By)**0.5                # betatron energy spread in y, rad
#  sx     = 0.001*sx
#  sy     = 0.01*sy
  nsx    = gf*sx                       # betatron energy spread in x, 1/gf
  nsy    = gf*sy                       # betatron energy spread in y, 1/gf
  # polarization parameters:
  pp1    = 0.0 # №1: laser linear polarization degree ξ_1 [0:1]
  pp2    = 0.0 # №2: laser linear polarization angle φl   [0:2pi]
  pp3    = 0.0 # №3: e ζ_parallel      * laser ξ_3        [0:1]
  pp4    = 0.1 # №4: e ζ_perpendicular * laser ξ_3        [0:1]

