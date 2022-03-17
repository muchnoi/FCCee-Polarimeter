import ROOT
from hardware import Laser, Spectrometer, EPD, PPD
import numpy as np; π = np.pi
from scipy             import fft                 as F_F_T
from scipy.interpolate import RectBivariateSpline as R_B_S
from scipy.ndimage     import gaussian_filter     as GAUSS


class XSMC: # DIFFERENTIAL CROSS SECTIONS FOR MONTE-CARLO (dσ / du dφ) +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
  def __init__(self): # x is u, y is φ
    A, B, C, D = '(1+x)', 'x/[0]', '(1-x/[0])', 'x/[0]*(1-x/[0])'
    T1 = '(1 + {}**2 - 4*{}*{}*(1 - [1]*cos(2*y) - [2]*sin(2*y)))/[0]'.format(A,A,D)
    T2 = ' - 2*[3]*{}*sqrt({})*([4]*cos(y) + [5]*sin(y))'.format(B,D)
    T3 = ' + [3]*[6]*{}*(x+2)*(1-2*{})'.format(B,B)
    XS = '({})/{}**3'.format(T1+T2+T3, A)
    self.U = ROOT.TF2('U', XS)
    self.U.SetRange( 0, 0, Spectrometer.κ, 2*π)
    self.U.SetParameter(0, Spectrometer.κ)
    self.U.SetParameter(1, Laser.ξ1)
    self.U.SetParameter(2, Laser.ξ2)
    self.U.SetParameter(3, Laser.ξ3)
    self.U.SetParameter(4, Spectrometer.ζx)
    self.U.SetParameter(5, Spectrometer.ζy)
    self.U.SetParameter(6, Spectrometer.ζz)
    self.U.SetNpx(1000); self.U.SetNpy(720)
#    print(XS)

class Photons: # THIS IS THE CROSS SECTION FOR SCATTERED PHOTONS (dσ / dηx dηy) +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

  def __init__(self):
    self.xs    = [None for i in range(6)]
    self.ηx    =  None
    self.ηy    =  None

  def Components(self):
    uoκ        = 1/(1 + self.ηx**2 + self.ηy**2)         # u over kappa
    uoκ2       = uoκ*uoκ                                 # u over kappa ^2
    uoκ3       = uoκ*uoκ2                                # u over kappa ^3
    uoκ4       = uoκ*uoκ3                                # u over kappa ^4
    u          = uoκ*Spectrometer.κ                      # u
    u1         = 1 + u                                   # 1 + u
    idn2       = 1/u1**2                                 # inverse denominator
    idn3       = idn2/u1                                 # inverse denominator
    self.xs[0] = 2 *idn3*uoκ2*(1+u1**2-4*u1*uoκ*(1-uoκ)) # unpolarized xSection
    self.xs[1] = 8 *idn2*uoκ4 *(self.ηx**2 - self.ηy**2) # vert/hor linear laser polarization ξ1
    self.xs[2] = 16*idn2*uoκ4 * self.ηx    * self.ηy     # diagonal linear laser polarization ξ2
    self.xs[3] = -4*idn3*uoκ3 * u * self.ηx              # transverse horizontal electron polarization
    self.xs[4] = -4*idn3*uoκ3 * u * self.ηy              # transverse vertical electron polarization
    self.xs[5] =  2*idn3*uoκ3*(u+2)*(Spectrometer.κ-2*u) # longitudinal electron polarization

  def Total(self, parameters):
    ξ1, ξ2, ζx, ζy, ζz = parameters
    self.xst   = (self.xs[0] + ξ1*self.xs[1] + ξ2*self.xs[2] + ζx*self.xs[3] + ζy*self.xs[4] + ζz*self.xs[5])
#    print('stocksp: ξ1={:8.6f} ξ2={:8.6f}'.format(ξ1, ξ2))#, end = '')
#    print('beampol: ζx={:8.6f} ζy={:8.6f} ζz={:8.6f}'.format(ζx, ζy, ζz))#, end = '')


class pPIXELS: # THIS IS THE PIXEL DETECTOR FOR SCATTERED PHOTONS  +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
  def __init__(self):
    self.iteration, self.position, self.compton, self.emittance = 0, [], [], []
    self.xmid_n = np.zeros(PPD.X_npix, dtype=np.double)
    self.ymid_n = np.zeros(PPD.Y_npix, dtype=np.double)
    self.ibase  = 1.e-3*Spectrometer.γ/Spectrometer.leip_L # from mm to 1/γ
    self.XS     = Photons()
    self.Ax, self.Bx = 1/PPD.X_pix, - PPD.X_beam/PPD.X_pix
    self.Ay, self.By = 1/PPD.Y_pix, - PPD.Y_beam/PPD.Y_pix

  def __call__(self, x, p):
    position, compton, change = [p[i] for i in (0,1)], [p[i] for i in (2,3,4,5,6)], False
    emittance = [p[7]/PPD.X_pix, p[8]/PPD.Y_pix]

    if self.position != position: # if spot center parameters changed
      self.position = X0, Y0 = position; change = True
      for xpix in range(PPD.X_npix): self.xmid_n[xpix] = self.ibase*(PPD.X_beam + (xpix+0.5)*PPD.X_pix - X0)
      for ypix in range(PPD.Y_npix): self.ymid_n[ypix] = self.ibase*(PPD.Y_beam + (ypix+0.5)*PPD.Y_pix - Y0)
      self.XS.ηy, self.XS.ηx = np.meshgrid(self.ymid_n, self.xmid_n, sparse=True)
      self.XS.Components()
    if (self.compton != compton) or change: # if Compton parameters changed
      self.compton = compton;            change = True
      self.XS.Total(compton)
    if (self.emittance != emittance) or change:  # if emittance parameters changed
      self.emittance = emittance
      self.convolution         = GAUSS(self.XS.xst, sigma = emittance, mode='reflect', truncate=5.0)
#      print('emittance: σx={:8.6f} σy={:8.6f} '.format(emittance[0],emittance[1]))#, end = '')
      self.iteration += 1;   print('iteration: %d ' % (self.iteration))
    return p[9]*self.convolution[int(self.Ax*x[0] + self.Bx)][int(self.Ay*x[1] + self.By)]


class Electrons: # THIS IS THE CROSS SECTION FOR SCATTERED ELECTRONS  (dσ / dx dy) +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
  R          = 1.5                                          # calculation range in unit radius
  N          = 4096                                         # number of x,y divisions
  D          = 2*R/N                                        # distance between knots
  x          = np.linspace(D/2-R, R-D/2, num=N)             # knots positions in x
  y          = x.copy()                                     # knots positions in y
  Y, X       = np.meshgrid(y, x, sparse=True)               # knots grid in x and y
  νx         = F_F_T.fftshift(F_F_T.fftfreq(N, d=D))        # FFT frequencies in x 
  νy         = νx.copy()                                    # FFT frequencies in y
  νy, νx     = np.meshgrid(νy, νx, sparse=True)             # FFT grid in νx and νy
  ρ          = 2*π*(νx**2 + νy**2)**0.5 + np.finfo(float).tiny**0.5
  Hankel     = np.sin(ρ)/ρ                                  # numpy function sinc(x) = sin(πx)/(πx)
  emittance  = [None, None]                                 # σx, σy
  stockspar  = [None, None, None]                           # ξ1, ζy, ζz

  def __init__(self):
    uoκ      = (1 + self.X)/2                               # u over kappa
    u        = uoκ*Spectrometer.κ                           # u
    u1       = 1+u                                          # 1 + u
    idn      = 1/(1+u)**3                                   # inverse denominator
    self.xso = idn*(1 + u1**2 - u1*(1-self.X**2))           # unpolarized xSection
    self.xs1 = idn*(1 - self.X**2 - 2*self.Y**2)*u1         # vert/hor linear laser polarization ξ1
    self.xsy = idn*u*self.Y                                 # transverse vertical electron polarization
    self.xsz = -idn*u*(u + 2)*self.X                        # longitudinal electron polarization

  def Prepare(self, parameters):
    ξ1, ζy, ζz, σx, σy = parameters
    if [σx, σy]     != self.emittance:
      self.emittance = [σx, σy]
      HA             = self.Hankel * np.exp(-2*π*π*((σx*self.νx)**2 + (σy*self.νy)**2))        # Fourier image by convolution theorem
      self.cvl       = np.abs(F_F_T.ifftshift(F_F_T.irfft2(HA, s=[self.N, self.N], workers=-1, norm='forward'))) # convolution result
      print('emittance: σx={:8.6f} σy={:8.6f} '.format(σx, σy))#, end = '')
    if [ξ1, ζy, ζz] != self.stockspar:
      self.stockspar = [ξ1, ζy, ζz]
      self.xst        = self.xso + ξ1*self.xs1 + ζy*self.xsy + ζz*self.xsz
      print('stockspar: ξ={:8.6f} ζy={:8.6f} ζz={:8.6f} '.format(ξ1, ζy, ζz))#, end = '')
    self.SPLINE      =  R_B_S(self.y, self.x, self.xst*self.cvl, kx=3, ky=3, s = 0 )  # this is the RectBivariateSpline


class ePIXELS: # THIS IS THE PIXEL DETECTOR FOR SCATTERED ELECTRONS  +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
  def __init__(self):
    self.iteration, self.ellipse, self.compton = 0, [], []
    self.xmid_n      = np.zeros(EPD.X_npix, dtype = np.double)
    self.ymid_n      = np.zeros(EPD.Y_npix, dtype = np.double)
    self.XS          = Electrons()
    self.Ax, self.Bx = 1/EPD.X_pix, - EPD.X_beam/EPD.X_pix
    self.Ay, self.By = 1/EPD.Y_pix, - EPD.Y_beam/EPD.Y_pix

  def __call__(self,x,p):
    ellipse, compton, change = [p[i] for i in (0,1,2,3)], [p[i] for i in (4,5,6,7,8)], False
    if self.compton != compton: # if Compton parameters changed
      self.compton = compton;                  change = True
      self.XS.Prepare(compton)
    if self.ellipse != ellipse: # if geometry parameters changed
      self.ellipse = X0, X1, Y0, Y1 = ellipse; change = True
      # 1) center position [mm]        2) unit radius [1/mm]
      CX = EPD.X_beam - (X1 + X0)/2;  RX = 2/(X1-X0)
      CY = EPD.Y_beam - (Y1 + Y0)/2;  RY = 2/(Y1-Y0)
      for xpix in range(EPD.X_npix):  self.xmid_n[xpix] = RX*(CX + (xpix+0.5)*EPD.X_pix)
      for ypix in range(EPD.Y_npix):  self.ymid_n[ypix] = RY*(CY + (ypix+0.5)*EPD.Y_pix)
    if change:
      self.BLUR = self.XS.SPLINE(self.xmid_n, self.ymid_n)
      self.iteration += 1;   print('iteration: %d ' % (self.iteration))
    return p[9]*self.BLUR[int(self.Ax*x[0] + self.Bx)][int(self.Ay*x[1] + self.By)]


