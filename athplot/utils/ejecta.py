# AthPlot: manipulate and plot AthenaK binary data
# Copyright (C) 2026, Jacob Fields <fields@ias.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from .eos import EquationOfState
from .sph_integrate import *

from .units import *

class Ejecta:
  '''
    A class representing an ejecta snapshot. Parts of the functionality have been adapted
    from Ruocheng Zhai's ejecta scripts.
    Inputs:
      fname_mhd: the path to the MHD data
      fname_adm: the path to the ADM data
      eos: the EquationOfState object for this ejecta profile
      fname_gauge: the path to the gauge variables. If None (default), the gauge is
                   read in from the ADM data. If a single string is provided, it is
                   assumed to be Z4c data. Alternatively, one can also provide a tuple in
                   the form (alpha, betax, betay, betaz) to provide separate (Z4c) paths
                   for each variable, such as if they are output separately from the
                   complete Z4c output.
  '''
  def __init__(self, fname_mhd, fname_adm, eos, fname_gauge=None):
    mhd = SphericalData(fname_mhd)
    adm = SphericalData(fname_adm)
    self.r = mhd.r
    self.theta = mhd.theta
    self.phi = mhd.phi
    self.time = mhd.time
    self.eos = eos
    
    self.dens  = mhd.fields['dens']
    self.press = mhd.fields['press']
    self.yq    = mhd.fields['s_00']
    self.temperature = mhd.fields['temperature']

    self.velx = mhd.fields['velx']
    self.vely = mhd.fields['vely']
    self.velz = mhd.fields['velz']

    self.bcc1 = mhd.fields['bcc1']
    self.bcc2 = mhd.fields['bcc2']
    self.bcc3 = mhd.fields['bcc3']

    self.gxx = adm.fields['adm_gxx']
    self.gxy = adm.fields['adm_gxy']
    self.gxz = adm.fields['adm_gxz']
    self.gyy = adm.fields['adm_gyy']
    self.gyz = adm.fields['adm_gyz']
    self.gzz = adm.fields['adm_gzz']

    # Note that there are three different distinct ways the gauge can be provided:
    # 1. In prescribed spacetimes, the gauge variables will be stored inside the ADM
    #    variables.
    # 2. In dynamical spacetimes, the gauge variables are part of the Z4c variables and
    #    can be provided via the Z4c output.
    # 3. The gauge variables can also be output separately from the other Z4c variables,
    #    which is commonly done if none of the auxiliary spacetime information is wanted.
    if fname_gauge is None:
      self.alpha = adm.fields['adm_alpha']
      self.betax= adm.fields['adm_betax']
      self.betay= adm.fields['adm_betay']
      self.betaz= adm.fields['adm_betaz']
    elif type(fname_gauge) is tuple:
      (fname_alpha, fname_betax, fname_betay, fname_betaz) = fname_gauge
      data_alpha = SphericalData(fname_alpha)
      data_betax = SphericalData(fname_betax)
      data_betay = SphericalData(fname_betay)
      data_betaz = SphericalData(fname_betaz)

      self.alpha = data_alpha.fields['z4c_alpha']
      self.betax = data_betax.fields['z4c_betax']
      self.betay = data_betay.fields['z4c_betay']
      self.betaz = data_betaz.fields['z4c_betaz']
    else:
      z4c = SphericalData(fname_gauge)

      self.alpha = z4c.fields['z4c_alpha']
      self.betax = z4c.fields['z4c_betax']
      self.betay = z4c.fields['z4c_betay']
      self.betaz = z4c.fields['z4c_betaz']

    # There was a bug in AthenaK at one point that led to double-counting in the data for
    # cells near mesh-block boundaries. This is a kluge that takes care of the
    # double-counting by recognizing that the lapse should be close to 1 but not more than
    # 1.
    z4c_alpha_raw = self.alpha
    phi_zero    = (self.phi == self.phi.min())
    theta_north = (self.theta == self.theta.min())
    alpha_phi0  = z4c_alpha_raw[phi_zero & ~theta_north]
    alpha_north = z4c_alpha_raw[theta_north]

    # -- Build duplication correction mask (only if broken data detected) --
    dup_corr = np.ones(z4c_alpha_raw.shape)
    if alpha_phi0.max() > 1.0:
        dup_corr[phi_zero] = 0.5          # phi=0: 2x duplication
    if alpha_north.max() > 1.0:
        dup_corr[theta_north] = 0.25      # north pole: 4x, overrides phi=0 there

    self.dens *= dup_corr
    self.press *= dup_corr
    self.yq *= dup_corr
    self.temperature *= dup_corr

    self.velx *= dup_corr
    self.vely *= dup_corr
    self.velz *= dup_corr

    self.bcc1 *= dup_corr
    self.bcc2 *= dup_corr
    self.bcc3 *= dup_corr

    self.gxx *= dup_corr
    self.gxy *= dup_corr
    self.gxz *= dup_corr
    self.gyy *= dup_corr
    self.gyz *= dup_corr
    self.gxz *= dup_corr

    self.alpha *= dup_corr
    self.betax *= dup_corr
    self.betay *= dup_corr
    self.betaz *= dup_corr

    # Further corrections in case the interpolation does not restrict the physical
    # constraints of the data.
    self.alpha = np.clip(self.alpha, 0.0, 1.0)

    # Note that the EOS is provided rather than loaded separately. This is done because
    # copies of the EOS can be rather large objects and are inconvenient to carry around.
    # This also means that the Ejecta class doesn't need to know anything about how the
    # EOS is stored.
    self.eos = eos

  '''
    Compute the Poynting flux (T_em)^r_t across a spherical shell.
    Outputs: flux
  '''
  def calc_Poynting_flux(self):
    # Compute the four-velocity
    W = np.sqrt(1.0 + (self.gxx*self.velx**2 + 2.0*self.gxy*self.velx*self.vely +
                       2.0*self.gxz*self.velx*self.velz + self.gyy*self.vely**2 +
                       2.0*self.gyz*self.vely*self.velz + self.gzz*self.velz**2))

    ut_u = W/self.alpha
    ux_u = self.velx - self.betax*ut_u
    uy_u = self.vely - self.betay*ut_u
    uz_u = self.velz - self.betaz*ut_u

    # Lower the time-like component
    ut_d = -self.alpha*W \
           + self.gxx*self.betax*self.velx + self.gyy*self.betay*self.vely \
           + self.gzz*self.betaz*self.velz \
           + self.gxy*(self.betax*self.vely + self.betay*self.velx) \
           + self.gxz*(self.betax*self.velz + self.betaz*self.velx) \
           + self.gyz*(self.betay*self.velz + self.betaz*self.vely)

    # Undensitize the magnetic field
    sqrt_det = sqrt_det_metric_adm(self.gxx, self.gxy, self.gxz, self.gyy,
                                   self.gyz, self.gzz, 1.0)
    Bx = self.bcc1 / sqrt_det
    By = self.bcc2 / sqrt_det
    Bz = self.bcc3 / sqrt_det

    # Compute the magnetic field in the fluid frame
    Bsq = self.gxx*Bx*Bx + 2.0*self.gxy*Bx*By + 2.0*self.gxz*Bx*Bz + \
          self.gyy*By*By + 2.0*self.gyz*By*Bz + self.gzz*Bz*Bz

    BWv = self.gxx*Bx*self.velx + self.gyy*By*self.vely + self.gzz*Bz*self.velz + \
          self.gxy*(Bx*self.vely + self.velx*By) + \
          self.gxz*(Bx*self.velz + self.velx*Bz) + \
          self.gyz*(By*self.velz + self.vely*Bz)

    bsq = (Bsq + BWv*BWv)/(W*W)

    bt_u = BWv/self.alpha
    bx_u = (Bx + BWv*ux_u)/W
    by_u = (By + BWv*uy_u)/W
    bz_u = (Bz + BWv*uz_u)/W

    # Lower the time-like component
    betax_d = self.gxx*self.betax + self.gxy*self.betay + self.gxz*self.betaz
    betay_d = self.gxy*self.betax + self.gyy*self.betay + self.gyz*self.betaz
    betaz_d = self.gxz*self.betax + self.gyz*self.betay + self.gzz*self.betaz

    betasq = self.gxx*self.betax**2 + 2.0*self.gxy*self.betax*self.betay + \
             2.0*self.gxz*self.betax*self.betaz + self.gyy*self.betay**2 + \
             2.0*self.gyz*self.betay*self.betaz + self.gzz*self.betaz**2

    bt_d = (-self.alpha*self.alpha + betasq)*bt_u + betax_d*bx_u + betay_d*by_u + \
            betaz_d*bz_u

    # Compute the radial components of the velocity and magnetic field
    ur_u = calc_ur(ux_u, uy_u, uz_u, self.theta, self.phi)
    br_u = calc_ur(bx_u, by_u, bz_u, self.theta, self.phi)

    # Compute the Ponyting flux
    integrand = (bsq*ur_u*ut_d - br_u*bt_d)*self.alpha*sqrt_det*self.r**2*np.sin(self.theta)

    flux = integrate_over_sphere(integrand, self.theta[0,:], self.phi[:,0])

    return flux

  '''
    Compute the mass ejection rate, including estimates of the unbound ejecta using
    both the geodesic and Bernoulli criteria.
    Outputs: Mej_rate, Mej_rate_Ber, Mej_rate_geo
      Mej_rate: the total outgoing mass flux across the spherical surface
      Mej_rate_Ber: the estimated unbound mass flux using the Bernoulli criterion
      Mej_rate_geo: the estimated unbound mass flux using the geodesic criterion
  '''
  def calc_Mej_rate(self):
    # Compute the number density (in fm^-3) from the rest-mass density (in Msun^-2)
    mn = self.eos.mn
    nb = conv_dens(cactus, cgs, self.dens)*1e-39/(mn*cgs.eV*1e6/cgs.light_speed**2)

    # Note that the data has been interpolated, so the temperature is *not* consistent
    # with the pressure inside the output. We need to recompute the pressure using that
    # temperature to get a consistent enthalpy, so it's just easier to compute the
    # enthalpy using the EOS, even though it would technically be cheaper to compute
    # the energy density and add it directly to the pressure.
    h = self.eos.calc_enthalpy(nb, self.yq, self.temperature)

    W = np.sqrt(1.0 + (self.gxx*self.velx**2 + 2.0*self.gxy*self.velx*self.vely +
                       2.0*self.gxz*self.velx*self.velz + self.gyy*self.vely**2 +
                       2.0*self.gyz*self.vely*self.velz + self.gzz*self.velz**2))

    ut_u = W/self.alpha
    ux_u = self.velx - self.betax*ut_u
    uy_u = self.vely - self.betay*ut_u
    uz_u = self.velz - self.betaz*ut_u

    ut_d = -self.alpha*W \
           + self.gxx*self.betax*self.velx + self.gyy*self.betay*self.vely \
           + self.gzz*self.betaz*self.velz \
           + self.gxy*(self.betax*self.vely + self.betay*self.velx) \
           + self.gxz*(self.betax*self.velz + self.betaz*self.velx) \
           + self.gyz*(self.betay*self.velz + self.betaz*self.vely)

    ur_u = calc_ur(ux_u, uy_u, uz_u, self.theta, self.phi)

    mask = ur_u > 0
    mask_Ber = ((h*ut_d) < -self.eos.hmin) & (ur_u > 0)
    mask_geo = (ut_d < -1) & (ur_u > 0)
    dens = np.where(~mask, 0.0, self.dens)
    dens_Ber = np.where(~mask_Ber, 0.0, dens)
    dens_geo = np.where(~mask_geo, 0.0, dens)

    sqrt_det = sqrt_det_metric_adm(self.gxx, self.gxy, self.gxz, self.gyy,
                                   self.gyz, self.gzz, self.alpha)

    integrand = dens*ur_u*sqrt_det*self.r**2*np.sin(self.theta)
    integrand_Ber = dens_Ber*ur_u*sqrt_det*self.r**2*np.sin(self.theta)
    integrand_geo = dens_geo*ur_u*sqrt_det*self.r**2*np.sin(self.theta)

    Mej_rate = integrate_over_sphere(integrand, self.theta[0,:], self.phi[:,0])
    Mej_rate_Ber = integrate_over_sphere(integrand_Ber, self.theta[0,:], self.phi[:,0])
    Mej_rate_geo = integrate_over_sphere(integrand_geo, self.theta[0,:], self.phi[:,0])

    return Mej_rate, Mej_rate_Ber, Mej_rate_geo
