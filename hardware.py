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
  ξ1     = 0.0                         # laser linear polarization horizontal (ξ1=1) or vertical ξ1=-1
  ξ2     = 0.0                         # laser linear polarization (φ=π/4 ξ2=1) or vertical (φ=-π/4 ξ2=-1)
  ξ3     = (1 - ξ1**2 - ξ2**2)**0.5    # laser circular polarization degree ξ3 [-1:1]

class Spectrometer:
  Eo     = 45.6e+9                     # electron (beam) energy, eV
  γ      = Eo/Constants.me             # electron Lorentz factor
  κ      = 4*γ*Laser.ωo/Constants.me   # 2 * rest frame momentum of laser photon
  bend   = .0021341                    # bending angle, rad
  bend_l = 24.12                       # dipole length, m
  spec_L = 0.5*bend_l + 88.0           # spectrometer arc length, m
  leip_L = 0.5*bend_l + spec_L + 5.0   # laser electron i.p. distance, m
  no     = γ*bend                      # bending angle in units 1/γ
  ex     = 0.27e-9                     # emittance x, m*rad
  ey     = 1.0e-12                     # emittance y, m*rad
  Bx     = 100.0                       # βx, m
  By     = 25.00                       # βy, m
  sx     = (ex/Bx)**0.5                # betatron angular spread in x, rad
  sy     = (ey/By)**0.5                # betatron angular spread in y, rad
  nsx    = γ*sx                        # betatron angular spread in x, 1/γ
  nsy    = γ*sy                        # betatron angular spread in y, 1/γ
  ζx     = 0.0                         # transverse electron spin polariztion along x
  ζy     = 0.0                         # transverse electron spin polariztion along y
  ζz     = 0.0                         # longitudinal electron spin polariztion along z

class EPD:                             # Electron Pixel Detector
  X_beam = 40.0                        # beam-detector horizontal space, mm
  Y_beam = -2.0                        # beam-detector  vertical  space, mm
  X_size = 400.0                       # detector X-size (horizontal)  , mm
  Y_size = 4.0                         # detector Y-size   (vertical)  , mm
  X_npix = 1600                        # number of pixels in X
  Y_npix = 80                          # number of pixels in Y
  X_pix  = X_size/float(X_npix)        # one pixel in X, mm
  Y_pix  = Y_size/float(Y_npix)        # one pixel in Y, mm

class PPD:                             # Photon Pixel Detector
  X_beam = -219.0                      # beam-detector horizontal space, mm
  Y_beam = -5.0                        # beam-detector  vertical  space, mm
  X_size = 10.0                        # detector X-size (horizontal)  , mm
  Y_size = 10.0                        # detector Y-size (horizontal)  , mm
  X_npix = 128                         # number of pixels in X
  Y_npix = 128                         # number of pixels in Y
  X_pix  = X_size/float(X_npix)        # one pixel in X, mm
  Y_pix  = Y_size/float(Y_npix)        # one pixel in Y, mm

