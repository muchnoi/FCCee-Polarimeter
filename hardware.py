class Constants:
  c      = 2.99792458e+10              # speed of light, cm/s
  me     = 0.510998928e+6              # electron rest energy, eV/c^2
  h      = 6.62607015e-34              # Plank constant, J*s
  e      = 1.602176634e-19             # electron charge magnitude, C
  re     = 2.8179403267e-13            # classical electron radius, cm
  barn   = 1.0e-24                     # one barn, m^2
  hc     = h*c/e                       # hc, eV*cm

class Laser:
  λo     = 0.532e-4                    # λo, cm
  ωo     = Constants.hc/λo             # ωo, eV
  ξ1     = 0.0                        # laser linear polarization horizontal (ξ1=1) or vertical ξ1=-1
  ξ2     = 0.0                        # laser linear polarization (φ=π/4 ξ2=1) or vertical (φ=-π/4 ξ2=-1)
  ξ3     = (1 - ξ1**2 - ξ2**2)**0.5    # laser circular polarization degree ξ3 [-1:1]

class Spectrometer:
  Eo     = 45.6e+9                     # electron (beam) energy, eV
  σE     = 1.e-3                       # beam energy spread, a.u.
  Dx     = 25.0                        # horizontal dispersion at i.p., mm
  γ      = Eo/Constants.me             # electron Lorentz factor
  κ      = 4*γ*Laser.ωo/Constants.me   # 2 * rest frame momentum of laser photon
  θo     = .0021341                    # bending angle, rad
  bend_l = 24.12                       # dipole length, m
  spec_L = 0.5*bend_l + 88.0           # spectrometer arc length, m
  leip_L = 0.5*bend_l + spec_L + 5.0   # laser electron i.p. distance, m
  εx     = 0.27e-9                     # emittance x, m*rad
  εy     = 1.0e-12                     # emittance y, m*rad
  βx     = 100.0                       # βx, m
  βy     = 20.00                       # βy, m
  σx     = 1000*(εx*βx)**0.5           # betatron size in x, mm
  σy     = 1000*(εy*βy)**0.5           # betatron size in y, mm
  ηx     = (εx/βx)**0.5                # betatron angular spread in x, rad
  ηy     = (εy/βy)**0.5                # betatron angular spread in y, rad
  ζx     = 0.0                         # transverse electron spin polariztion along x
  ζy     = 0.0                         # transverse electron spin polariztion along y
  ζz     = 0.0                         # longitudinal electron spin polariztion along z

class EPD:                             # Electron Pixel Detector
  X_beam = 15.0                        # beam-detector horizontal space, mm
  Y_beam = -2.0                        # beam-detector  vertical  space, mm
  X_size = 400.0                       # detector X-size (horizontal)  , mm
  Y_size = 4.0                         # detector Y-size   (vertical)  , mm
  X_npix = 1600                        # number of pixels in X
  Y_npix = 80                          # number of pixels in Y
  X_pix  = X_size/float(X_npix)        # one pixel in X, mm
  Y_pix  = Y_size/float(Y_npix)        # one pixel in Y, mm

class PPD:                             # Photon Pixel Detector
  X_beam = -218.5                      # beam-detector horizontal space, mm
  Y_beam = -5.0                        # beam-detector  vertical  space, mm
  X_size = 10.                         # detector X-size (horizontal)  , mm
  Y_size = 10.                         # detector Y-size (horizontal)  , mm
  X_npix = 100                         # number of pixels in X
  Y_npix = 100                         # number of pixels in Y
  X_pix  = X_size/float(X_npix)        # one pixel in X, mm
  Y_pix  = Y_size/float(Y_npix)        # one pixel in Y, mm

