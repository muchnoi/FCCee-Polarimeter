import ROOT
from hardware import Laser, Spectrometer
import numpy as np; π = np.pi
from scipy.ndimage     import gaussian_filter     as GAUSS
#from scipy             import fft                 as F_F_T
#from scipy.interpolate import RectBivariateSpline as R_B_S
#from scipy.ndimage     import fourier_gaussian    as GAUSS


class XSMC: # DIFFERENTIAL CROSS SECTIONS FOR MONTE-CARLO (dσ / du dφ) +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
  def __init__(self): # x is u, y is φ
    A, B, C, D = '(1+x)', 'x/[0]', '(1-x/[0])', 'x/[0]*(1-x/[0])'
    T1 = '(1 + {}**2 - 4*{}*{}*(1 - [1]*cos(2*y) - [2]*sin(2*y)))/[0]'.format(A,A,D)
    T2 = ' - 2*[3]*{}*sqrt({})*([4]*cos(y) + [5]*sin(y))'.format(B,D)
    T3 = ' + [3]*[6]*{}*(x+2)*(1-2*{})'.format(B,B)
    XS = '({})/{}**3'.format(T1+T2+T3, A); #    print(XS)
    self.U = ROOT.TF2('U', XS)
    self.U.SetRange( 0, 0, Spectrometer.κ, 2*π)
    self.U.SetParameter(0, Spectrometer.κ)
    self.U.SetParameter(1, Laser.ξ1)
    self.U.SetParameter(2, Laser.ξ2)
    self.U.SetParameter(3, Laser.ξ3)
    self.U.SetParameter(4, Spectrometer.ζx)
    self.U.SetParameter(5, Spectrometer.ζy)
    self.U.SetParameter(6, Spectrometer.ζz)
    self.U.SetNpx(1000); self.U.SetNpy(3600)


class Photons: # THIS IS THE CROSS SECTION FOR SCATTERED PHOTONS (dσ / dηx dηy) +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

  def __init__(self, setup):
    self.κ     = setup.κ
    self.xs    = [None for i in range(6)]
    self.ηx    =  None
    self.ηy    =  None

  def Components(self):
    uoκ        = 1/(1 + self.ηx**2 + self.ηy**2)         # u over kappa
    uoκ2       = uoκ*uoκ                                 # u over kappa ^2
    uoκ3       = uoκ*uoκ2                                # u over kappa ^3
    uoκ4       = uoκ*uoκ3                                # u over kappa ^4
    u          = uoκ*self.κ                              # u
    u1         = 1 + u                                   # 1 + u
    idn2       = 1/u1**2                                 # inverse denominator
    idn3       = idn2/u1                                 # inverse denominator
    self.xs[0] =  2*idn3*uoκ2*(1+u1**2-4*u1*uoκ*(1-uoκ)) # unpolarized xSection
    self.xs[1] =  8*idn2*uoκ4 *(self.ηx**2 - self.ηy**2) # vert/hor linear laser polarization ξ1
    self.xs[2] = 16*idn2*uoκ4 * self.ηx    * self.ηy     # diagonal linear laser polarization ξ2
    self.xs[3] = -4*idn3*uoκ3 * u * self.ηx              # transverse horizontal electron polarization
    self.xs[4] = -4*idn3*uoκ3 * u * self.ηy              # transverse vertical electron polarization
    self.xs[5] =  2*idn3*uoκ3*(u+2)*(self.κ-2*u)         # longitudinal electron polarization

  def Total(self, parameters):
    ξ1, ξ2, ζx, ζy, ζz = parameters
    self.xst   = (self.xs[0] + ξ1*self.xs[1] + ξ2*self.xs[2] + ζx*self.xs[3] + ζy*self.xs[4] + ζz*self.xs[5])
#    print('stocksp: ξ1={:8.6f} ξ2={:8.6f}'.format(ξ1, ξ2))#, end = '')
#    print('beampol: ζx={:8.6f} ζy={:8.6f} ζz={:8.6f}'.format(ζx, ζy, ζz))#, end = '')


class pPIXELS: # THIS IS THE PIXEL DETECTOR FOR SCATTERED PHOTONS  +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
  def __init__(self, PPD, setup):
    self.Xm = 2 # number of each pixel subdivisions in X (for integration)
    self.Ym = 4 # number of each pixel subdivisions in Y (for integration)

    self.iteration, self.position, self.compton, self.emittance = 0, [], [], []
    self.X_npix, self.X_pix, self.X_beam = PPD.X_npix*self.Xm, PPD.X_pix/self.Xm, PPD.X_beam
    self.Y_npix, self.Y_pix, self.Y_beam = PPD.Y_npix*self.Ym, PPD.Y_pix/self.Ym, PPD.Y_beam
    self.xmid_n = np.zeros(self.X_npix, dtype=np.double)
    self.ymid_n = np.zeros(self.Y_npix, dtype=np.double)
    self.ibase  = 1.e-3*setup.γ/setup.leip_L # from mm to 1/γ
    self.XS     = Photons(setup)
    self.Ax, self.Bx = 1/PPD.X_pix, - PPD.X_beam/PPD.X_pix
    self.Ay, self.By = 1/PPD.Y_pix, - PPD.Y_beam/PPD.Y_pix

  def __call__(self, x, p):
    position, compton, change = [p[i] for i in (0,1)], [p[i] for i in (2,3,4,5,6)], False
    emittance = [p[7]/self.X_pix, p[8]/self.Y_pix]

    if self.position != position: # if spot center parameters changed
      self.position = X0, Y0 = position; change = True
      for xpix in range(self.X_npix): self.xmid_n[xpix] = self.ibase*(self.X_beam + (xpix+0.5)*self.X_pix - X0)
      for ypix in range(self.Y_npix): self.ymid_n[ypix] = self.ibase*(self.Y_beam + (ypix+0.5)*self.Y_pix - Y0)
      self.XS.ηy, self.XS.ηx = np.meshgrid(self.ymid_n, self.xmid_n, sparse=True)
      self.XS.Components()
    if (self.compton != compton) or change: # if Compton parameters changed
      self.compton = compton;            change = True
      self.XS.Total(compton)
    if (self.emittance != emittance) or change:  # if emittance parameters changed
      self.emittance = emittance
      self.convolution         = GAUSS(self.XS.xst, sigma = emittance, truncate=5.0, mode='nearest')
#      print('emittance: σx={:8.6f} σy={:8.6f} '.format(emittance[0],emittance[1]))#, end = '')
      self.iteration += 1;
      if np.random.rand()<0.05:  print('iteration: %d ' % (self.iteration))
    res = 0.0
    for i in range(self.Xm):
      xbin = int(self.Ax*x[0] + self.Bx)*self.Xm + i
      for j in range(self.Ym):
        res += self.convolution[xbin][int(self.Ay*x[1] + self.By)*self.Ym + j]
    return p[9]*res
#    return p[9]*self.convolution[int(self.Ax*x[0] + self.Bx)][int(self.Ay*x[1] + self.By)]


class Electrons: # THIS IS THE CROSS SECTION FOR SCATTERED ELECTRONS  (dσ / dx dy) +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

  def __init__(self, setup):
    self.κ     = setup.κ
    θo         = setup.θo*setup.γ
    self.A     = 1/(1 + θo**2)**0.5
    self.B     = θo*self.A
    self.xs    = [None for i in range(6)]
    self.x     =  None
    self.y     =  None

  def Components(self):
    R          = np.real(np.emath.sqrt(1.0 - self.x**2 - self.y**2))             # 'Root' - square root
    Δ          = {'+': self.B*self.x + self.A*R,  '-': self.B*self.x - self.A*R} # Δ±
    Δ2         = {'+': Δ['+']**2,                 '-': Δ['-']**2               } # (Δ±)^2
    δ          = {'+':-self.A*self.x + self.B*R,  '-':-self.A*self.x - self.B*R} # δ±
    u          = {'+':  0.5*self.κ*(1 + Δ['+']),  '-':  0.5*self.κ*(1 + Δ['-'])} # u±
    u1         = {'+': 1 + u['+'],                '-': 1 + u['-']              } # (1 + u±)
    u12        = {'+': u1['+']**2,                '-': u1['-']**2              } # (1 + u±)^2
    u13        = {'+': u1['+']*u12['+'],          '-': u1['-']*u12['-']        } # (1 + u±)^3
    self.xs[0] = ((1+u12['+']-u1['+']*(1-Δ2['+']))/u13['+'] + (1+u12['-']-u1['-']*(1-Δ2['-']))/u13['-'])*0.5 # unpolarized part
    self.xs[1] = ((δ['+']**2 - self.y**2)         /u12['+'] + (δ['-']**2 - self.y**2)         /u12['-'])*0.5 # hor/vert linear laser polarization ξ1
    self.xs[2] = -(δ['+']    * self.y             /u12['+'] +  δ['-']    * self.y             /u12['-'])     # diagonal linear laser polarization ξ2
    self.xs[3] = -(δ['+']    * u['+']             /u13['+'] +  δ['-']    * u['-']             /u13['-'])*0.5 # x electron polarization
    self.xs[4] =  (u['+']    * self.y             /u13['+'] +  u['-']    * self.y             /u13['-'])*0.5 # y electron polarization
    self.xs[5] = -(u['+'] * (2 + u['+']) * Δ['+'] /u13['+'] +  u['-'] * (2 + u['-']) * Δ['-'] /u13['-'])*0.5 # z electron polarization

  def Total(self, parameters):
    ξ1, ξ2, ζx, ζy, ζz = parameters
    self.xst  = 1e+3*(self.xs[0] + ξ1*self.xs[1] + ξ2*self.xs[2] + ζx*self.xs[3] + ζy*self.xs[4] + ζz*self.xs[5])


class ePIXELS: # THIS IS THE PIXEL DETECTOR FOR SCATTERED ELECTRONS  +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
  def __init__(self, EPD, setup):
    self.Xm = 2 # number of each pixel subdivisions in X (for integration)
    self.Ym = 5 # number of each pixel subdivisions in Y (for integration)
    self.iteration, self.ellipse, self.compton, self.emittance = 0, [], [], []
    self.X_npix, self.X_pix, self.X_beam = EPD.X_npix*self.Xm, EPD.X_pix/self.Xm, EPD.X_beam
    self.Y_npix, self.Y_pix, self.Y_beam = EPD.Y_npix*self.Ym, EPD.Y_pix/self.Ym, EPD.Y_beam
    self.gridx       = np.ones(  self.X_npix + 1,             dtype = np.double) # x grid knots
    self.gridy       = np.ones(  self.Y_npix + 1,             dtype = np.double) # y grid knots
    self.pixcx       = np.zeros( self.X_npix,                 dtype = np.double) # x pixel centers
    self.pixcy       = np.zeros( self.Y_npix,                 dtype = np.double) # y pixel centers
    self.Ax, self.Bx = 1/EPD.X_pix, - self.X_beam/EPD.X_pix
    self.Ay, self.By = 1/EPD.Y_pix, - self.Y_beam/EPD.Y_pix
    self.XS          = Electrons(setup)

  def __call__(self,x,p):
    ellipse,  compton = [p[i] for i in (0,1,2,3)], [p[i] for i in (4,5,6,7,8)]
    emittance, change = [p[9]/self.X_pix, p[10]/self.Y_pix], False

    if self.ellipse != ellipse: # if geometry parameters changed
      self.ellipse = X0, X1, Y0, Y1 = ellipse; change = True
      print ('X0={:8.3f} X1={:8.3f} Y0={:6.3f} Y1={:6.3f}'.format(X0, X1, Y0, Y1))
      # 1) center position [mm]        2) unit radius [1/mm]
      CX = self.X_beam - (X1 + X0)/2;  RX = 2/(X1-X0)
      CY = self.Y_beam - (Y1 + Y0)/2;  RY = 2/(Y1-Y0)
      for xpix in range(self.X_npix): 
        g = RX*(CX + xpix*self.X_pix)
        self.gridx[xpix] = -1.*(g<-1) + g*(-1<g<1) + 1.*(g>1)
        self.pixcx[xpix] = g + 0.5*RX*self.X_pix
      for ypix in range(self.Y_npix):
        g = RY*(CY + ypix*self.Y_pix) 
        self.gridy[ypix] = -1.*(g<-1) + g*(-1<g<1) + 1.*(g>1)
        self.pixcy[ypix] = g + 0.5*RY*self.Y_pix
      X, Y  = np.meshgrid(self.gridx, self.gridy, indexing='ij')
      R = np.real(np.emath.sqrt(1.0 - X**2 - Y**2)) 
      I = X*np.arctan2(Y,R) + Y*np.arctan2(X,R) - np.arctan2(X*Y,R)
      self.dxdy = I[:-1,:-1] + I[1:,1:] - I[1:,:-1] - I[:-1,1:]
      self.dxdy = np.where(self.dxdy<1e-10, 0, self.dxdy)
      self.XS.x, self.XS.y = np.meshgrid(self.pixcx, self.pixcy, indexing='ij')
      if self.iteration == 0: self.XS.Components()
    if self.compton != compton or change: # if Compton parameters changed
      self.compton = ξ1, ξ2, ζx, ζy, ζz = compton; change = True
      self.XS.Total(compton)
      self.cs = self.XS.xst*self.dxdy
      print ('ξ1={:6.3f} ξ2={:6.3f} ζx={:6.3f} ζy={:6.3f} ζz={:6.3f}'.format(ξ1, ξ2, ζx, ζy, ζz))
    if self.emittance != emittance or change:  # if emittance parameters changed
      self.emittance = σx, σy = emittance; change = True
      print ('σx = {:8.5f} pix, σy = {:8.5f} pix'.format(σx, σy))
      self.convolution         = GAUSS(self.cs, sigma = emittance, truncate=6.0, mode='nearest')
    if change:
      self.iteration += 1;   print('iteration: %d ' % (self.iteration))
    res = 0.0
    for i in range(self.Xm):
      xbin = int(self.Ax*x[0] + self.Bx)*self.Xm + i
      for j in range(self.Ym):
        res += self.convolution[xbin][int(self.Ay*x[1] + self.By)*self.Ym + j]
    return p[11]*res


