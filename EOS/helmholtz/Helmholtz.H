#ifndef EOS_HELMHOLTZ_H
#define EOS_HELMHOLTZ_H

#include "Eos.H"
#include "Array1D.H"
#include "Array2D.H"

/*
  Frank Timmes Helmholtz based Equation of State
  http://cococubed.asu.edu/

  given a temperature temp [K], density den [g/cm**3], and a composition
  characterized by abar and zbar, this routine returns most of the other
  thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
  specific thermal energy [erg/gr], the entropy [erg/g/K], along with
  their derivatives with respect to temperature, density, abar, and zbar.
  other quantites such the normalized chemical potential eta (plus its
  derivatives), number density of electrons and positron pair (along
  with their derivatives), adiabatic indices, specific heats, and
  relativistically correct sound speed are also returned.
    
  this routine assumes planckian photons, an ideal gas of ions,
  and an electron-positron gas with an arbitrary degree of relativity
  and degeneracy. interpolation in a table of the helmholtz free energy
  is used to return the electron-positron thermodynamic quantities.
  all other derivatives are analytic.
    
  references: cox & giuli chapter 24 ; timmes & swesty apj 1999
*/

class Helmholtz : public Eos, public UnifiedMemoryClass {

  // Runtime parameters  
  bool do_coulomb;
  bool input_is_constant;

  // for the tables, in general
  static const int imax = 541;
  static const int jmax = 201;

  int itmax, jtmax;

  static Array1D<double, imax> d;
  static Array1D<double, jmax> t;

  double tlo, thi, tstp, tstpi;
  double dlo, dhi, dstp, dstpi;
  double ttol, dtol;

  // ..for the helmholtz free energy tables (2D)
  static Array2D<double, imax, jmax> f;
  static Array2D<double, imax, jmax> fd;
  static Array2D<double, imax, jmax> ft;
  static Array2D<double, imax, jmax> fdd;
  static Array2D<double, imax, jmax> ftt;
  static Array2D<double, imax, jmax> fdt;
  static Array2D<double, imax, jmax> fddt;
  static Array2D<double, imax, jmax> fdtt;
  static Array2D<double, imax, jmax> fddtt;

  // ..for the pressure derivative with density ables
  static Array2D<double, imax, jmax> dpdf;
  static Array2D<double, imax, jmax> dpdfd;
  static Array2D<double, imax, jmax> dpdft;
  static Array2D<double, imax, jmax> dpdfdt;

  // ..for chemical potential tables
  static Array2D<double, imax, jmax> ef;
  static Array2D<double, imax, jmax> efd;
  static Array2D<double, imax, jmax> eft;
  static Array2D<double, imax, jmax> efdt;

  // ..for the number density tables
  static Array2D<double, imax, jmax> xf;
  static Array2D<double, imax, jmax> xfd;
  static Array2D<double, imax, jmax> xft;
  static Array2D<double, imax, jmax> xfdt;

  // ..for storing the differences
  static Array1D<double, jmax> dt_sav;
  static Array1D<double, jmax> dt2_sav;
  static Array1D<double, jmax> dti_sav;
  static Array1D<double, jmax> dt2i_sav;
  
  static Array1D<double, imax> dd_sav;
  static Array1D<double, imax> dd2_sav;
  static Array1D<double, imax> ddi_sav;
  static Array1D<double, imax> dd2i_sav;

  static const int max_newton = 100;

  // 2006 CODATA physical constants

  // Math constants
  static const double pi       = 3.1415926535897932384e0;

  // Physical constants
  static const double h       = 6.6260689633e-27;
  static const double hbar    = 0.5e0 * h/pi;
  static const double qe      = 4.8032042712e-10;
  static const double avo_eos = 6.0221417930e23;
  static const double clight  = 2.99792458e10;
  static const double kerg    = 1.380650424e-16;
  static const double ev2erg_eos  = 1.60217648740e-12;
  static const double kev     = kerg/ev2erg_eos;
  static const double amu     = 1.66053878283e-24;
  static const double me_eos  = 9.1093821545e-28;
  static const double rbohr   = hbar*hbar/(me_eos * qe * qe);
  static const double fine    = qe*qe/(hbar*clight);

#ifdef RADIATION
  static const double ssol    = 0.0e0;
#else
  static const double ssol    = 5.67051e-5;
#endif
  static const double asol    = 4.0d0 * ssol / clight;
  static const double weinlam = h*clight/(kerg * 4.965114232e0);
  static const double weinfre = 2.821439372e0*kerg/h;

  // Astronomical constants
  static const double ly      = 9.460528e17;
  static const double pc      = 3.261633e0 * ly;

  // Some other useful combinations of the constants
  static const double sioncon = (2.0e0 * pi * amu * kerg)/(h*h);
  static const double forth   = 4.0e0/3.0e0;
  static const double forpi   = 4.0e0 * pi;
  static const double kergavo = kerg * avo_eos;
  static const double ikavo   = 1.0e0/kergavo;
  static const double asoli3  = asol/3.0e0;
  static const double light2  = clight * clight;

  // Constants used for the Coulomb corrections
  static const double a1    = -0.898004e0;
  static const double b1    =  0.96786e0;
  static const double c1    =  0.220703e0;
  static const double d1    = -0.86097e0;
  static const double e1    =  2.5269e0;
  static const double a2    =  0.29561e0;
  static const double b2    =  1.9885e0;
  static const double c2    =  0.288675e0;
  static const double onethird = 1.0e0/3.0e0;
  static const double esqu = qe * qe;

public:

  void Helmholtz()
  {
    eos_name = "helmholtz";
  }

  void ~Helmholtz() {};

  void eos(const EosInput input, EosState& state)
  {

    const double ZERO = 0.0;
    const double HALF = 0.5;
    const double TWO = 2.0;

    double temp_row, den_row, abar_row, zbar_row;
    double ye_row, etot_row, ptot_row;
    double cv_row, cp_row;
    double xne_row, xnp_row;
    double etaele_row, pele_row, ppos_row;
    double dpd_row, dpt_row, dpa_row, dpz_row;
    double ded_row, det_row, dea_row, dez_row;
    double stot_row, dsd_row, dst_row, htot_row;
    double dhd_row, dht_row, dpe_row, dpdr_e_row;
    double gam1_row, cs_row;

    bool single_iter, double_iter, converged;
    int var, dvar, var1, var2, iter;
    double v_want;
    double v1_want, v2_want;
    double xnew, xtol, dvdx, smallx, error, v;
    double v1, v2, dv1dt, dv1dr, dv2dt,dv2dr, delr;
    double error1, error2, told, rold, tnew, rnew, v1i, v2i;

    double x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida;
    double dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt;
    double dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt;
    double deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt;
    double dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion;
    double sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd;
    double dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp;
    double gam1,gam2,gam3,chit,chid,nabad,sound,etaele;
    double detadt,detadd,xnefer,dxnedt,dxnedd,s;
    double temp,den,abar,zbar,ytot1,ye;

#ifdef EXTRA_THERMO
    // ..for the abar derivatives
    double dpradda,deradda,dsradda;
    double dpionda,deionda,dsionda;
    double dpepda,deepda,dsepda;
    double dpresda,denerda,dentrda;
    double detada,dxneda;


    // ..for the zbar derivatives
    double dpraddz,deraddz,dsraddz;
    double dpiondz,deiondz,dsiondz;
    double dpepdz,deepdz,dsepdz;
    double dpresdz,denerdz,dentrdz;
    double detadz,dxnedz;
#endif


    // ..for the interpolations
    int iat,jat;
    double free,df_d,df_t,df_tt,df_dt;
    double xt,xd,mxt,mxd;
    double si0t,si1t,si2t,si0mt,si1mt,si2mt;
    double si0d,si1d,si2d,si0md,si1md,si2md;
    double dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt;
    double dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md;
    double ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt;
    double z,din;
    Array1D<double, 36> fi;

    // ..for the coulomb corrections
    double dsdd,dsda,lami,inv_lami,lamida,lamidd;
    double plasg,plasgdd,plasgdt,plasgda,plasgdz;
    double ecoul,decouldd,decouldt,decoulda,decouldz;
    double pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz;
    double scoul,dscouldd,dscouldt,dscoulda,dscouldz;

    double p_temp, e_temp;

    temp_row = state.T;
    den_row  = state.rho;
    abar_row = state.abar;
    zbar_row = state.zbar;
    ye_row   = state.y_e;

    // Initial setup for iterations

    single_iter = false;
    double_iter = false;

    if (input == EosInput::rt) {

      // Nothing to do here.

    } else if (input == EosInput::rh) {

      single_iter = true;
      v_want = state.h;
      var  = ienth;
      dvar = itemp;

    }  else if (input == EosInput::tp) {

      single_iter = true;
      v_want = state.p;
      var  = ipres;
      dvar = idens;

    } else if (input == EosInput::rp) {

      single_iter = true;
      v_want = state.p;
      var  = ipres;
      dvar = itemp;

    } else if (input == EosInput::re) {

      single_iter = true;
      v_want = state.e;
      var  = iener;
      dvar = itemp;

    } else if (input == EosInput::ps) {

      double_iter = true;
      v1_want = state.p;
      v2_want = state.s;
      var1 = ipres;
      var2 = ientr;

    } else if (input == EosInput::ph) {

      double_iter = true;
      v1_want = state.p;
      v2_want = state.h;
      var1 = ipres;
      var2 = ienth;

    } else if (input == EosInput::th) {

      single_iter = true;
      v_want = state.h;
      var  = ienth;
      dvar = idens;

    }

    ptot_row = ZERO;
    dpt_row = ZERO;
    dpd_row = ZERO;
    dpa_row = ZERO;
    dpz_row = ZERO;
    dpe_row = ZERO;
    dpdr_e_row = ZERO;

    etot_row = ZERO;
    det_row = ZERO;
    ded_row = ZERO;
    dea_row = ZERO;
    dez_row = ZERO;

    stot_row = ZERO;
    dst_row = ZERO;
    dsd_row = ZERO;

    htot_row = ZERO;
    dhd_row = ZERO;
    dht_row = ZERO;

    pele_row = ZERO;
    ppos_row = ZERO;

    xne_row = ZERO;
    xnp_row = ZERO;

    etaele_row = ZERO;

    cv_row = ZERO;
    cp_row = ZERO;
    cs_row = ZERO;
    gam1_row = ZERO;

    converged = false;

    if (input == EosInput::rt) converged = true;

    for (int iter = 1; iter <= max_newton; iter++) {

      temp  = temp_row;
      den   =  den_row;
      abar  = abar_row;
      zbar  = zbar_row;

      ytot1 = 1.0e0 / abar;
      ye    = ye_row;
      din   = ye * den;

      // ..initialize
      deni    = 1.0e0/den;
      tempi   = 1.0e0/temp;
      kt      = kerg * temp;
      ktinv   = 1.0e0/kt;

      // ..radiation section:
      prad    = asoli3 * temp * temp * temp * temp;
      dpraddd = ZERO;
      dpraddt = 4.0e0 * prad*tempi;
#ifdef EXTRA_THERMO
      dpradda = ZERO;
      dpraddz = ZERO;
#endif

      erad    = 3.0e0 * prad*deni;
      deraddd = -erad*deni;
      deraddt = 3.0e0 * dpraddt*deni;
#ifdef EXTRA_THERMO
      deradda = ZERO;
      deraddz = ZERO;
#endif

      srad    = (prad*deni + erad)*tempi;
      dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi;
      dsraddt = (dpraddt*deni + deraddt - srad)*tempi;
#ifdef EXTRA_THERMO
      dsradda = ZERO;
      dsraddz = ZERO;
#endif

      // ..ion section:
      xni     = avo_eos * ytot1 * den;
      dxnidd  = avo_eos * ytot1;
      dxnida  = -xni * ytot1;

      pion    = xni * kt;
      dpiondd = dxnidd * kt;
      dpiondt = xni * kerg;
#ifdef EXTRA_THERMO
      dpionda = dxnida * kt;
      dpiondz = ZERO;
#endif

      eion    = 1.5e0 * pion*deni;
      deiondd = (1.5e0 * dpiondd - eion)*deni;
      deiondt = 1.5e0 * dpiondt*deni;
#ifdef EXTRA_THERMO
      deionda = 1.5e0 * dpionda*deni;
      deiondz = ZERO;
#endif

      x       = abar*abar*sqrt(abar) * deni/avo_eos;
      s       = sioncon * temp;
      z       = x * s * sqrt(s);
      y       = log(z);
      sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y;
      dsiondd = ((dpiondd*deni - pion*deni*deni + deiondd)*tempi
                 - kergavo * deni * ytot1);
      dsiondt = ((dpiondt*deni + deiondt)*tempi -  
                 (pion*deni + eion) * tempi*tempi  
                 + 1.5e0 * kergavo * tempi*ytot1);
      x       = avo_eos*kerg/abar;
#ifdef EXTRA_THERMO
      dsionda = ((dpionda*deni + deionda)*tempi  
                 + kergavo*ytot1*ytot1* (2.5e0 - y));
      dsiondz = ZERO;
#endif

      // ..electron-positron section:
      // ..assume complete ionization
      xnem    = xni * zbar;

      // ..enter the table with ye*den
      din = ye*den;

      // ..hash locate this temperature and density
      jat = int((log10(temp) - tlo)*tstpi) + 1;
      jat = max(1,min(jat,jtmax-1));
      iat = int((log10(din) - dlo)*dstpi) + 1;
      iat = max(1,min(iat,itmax-1));

      // ..access the table locations only once
      fi(1)  = f(iat,jat);
      fi(2)  = f(iat+1,jat);
      fi(3)  = f(iat,jat+1);
      fi(4)  = f(iat+1,jat+1);
      fi(5)  = ft(iat,jat);
      fi(6)  = ft(iat+1,jat);
      fi(7)  = ft(iat,jat+1);
      fi(8)  = ft(iat+1,jat+1);
      fi(9)  = ftt(iat,jat);
      fi(10) = ftt(iat+1,jat);
      fi(11) = ftt(iat,jat+1);
      fi(12) = ftt(iat+1,jat+1);
      fi(13) = fd(iat,jat);
      fi(14) = fd(iat+1,jat);
      fi(15) = fd(iat,jat+1);
      fi(16) = fd(iat+1,jat+1);
      fi(17) = fdd(iat,jat);
      fi(18) = fdd(iat+1,jat);
      fi(19) = fdd(iat,jat+1);
      fi(20) = fdd(iat+1,jat+1);
      fi(21) = fdt(iat,jat);
      fi(22) = fdt(iat+1,jat);
      fi(23) = fdt(iat,jat+1);
      fi(24) = fdt(iat+1,jat+1);
      fi(25) = fddt(iat,jat);
      fi(26) = fddt(iat+1,jat);
      fi(27) = fddt(iat,jat+1);
      fi(28) = fddt(iat+1,jat+1);
      fi(29) = fdtt(iat,jat);
      fi(30) = fdtt(iat+1,jat);
      fi(31) = fdtt(iat,jat+1);
      fi(32) = fdtt(iat+1,jat+1);
      fi(33) = fddtt(iat,jat);
      fi(34) = fddtt(iat+1,jat);
      fi(35) = fddtt(iat,jat+1);
      fi(36) = fddtt(iat+1,jat+1);

      // ..various differences
      xt  = max( (temp - t(jat))*dti_sav(jat), ZERO);
      xd  = max( (din - d(iat))*ddi_sav(iat), ZERO);
      mxt = 1.0e0 - xt;
      mxd = 1.0e0 - xd;

      // ..the six density and six temperature basis functions
      si0t =   psi0(xt);
      si1t =   psi1(xt)*dt_sav(jat);
      si2t =   psi2(xt)*dt2_sav(jat);

      si0mt =  psi0(mxt);
      si1mt = -psi1(mxt)*dt_sav(jat);
      si2mt =  psi2(mxt)*dt2_sav(jat);

      si0d =   psi0(xd);
      si1d =   psi1(xd)*dd_sav(iat);
      si2d =   psi2(xd)*dd2_sav(iat);

      si0md =  psi0(mxd);
      si1md = -psi1(mxd)*dd_sav(iat);
      si2md =  psi2(mxd)*dd2_sav(iat);

      // ..derivatives of the weight functions
      dsi0t =   dpsi0(xt)*dti_sav(jat);
      dsi1t =   dpsi1(xt);
      dsi2t =   dpsi2(xt)*dt_sav(jat);

      dsi0mt = -dpsi0(mxt)*dti_sav(jat);
      dsi1mt =  dpsi1(mxt);
      dsi2mt = -dpsi2(mxt)*dt_sav(jat);

      dsi0d =   dpsi0(xd)*ddi_sav(iat);
      dsi1d =   dpsi1(xd);
      dsi2d =   dpsi2(xd)*dd_sav(iat);

      dsi0md = -dpsi0(mxd)*ddi_sav(iat);
      dsi1md =  dpsi1(mxd);
      dsi2md = -dpsi2(mxd)*dd_sav(iat);

      // ..second derivatives of the weight functions
      ddsi0t =   ddpsi0(xt)*dt2i_sav(jat);
      ddsi1t =   ddpsi1(xt)*dti_sav(jat);
      ddsi2t =   ddpsi2(xt);

      ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat);
      ddsi1mt = -ddpsi1(mxt)*dti_sav(jat);
      ddsi2mt =  ddpsi2(mxt);

      // comment begins
      ddsi0d =   ddpsi0(xd)*dd2i_sav(iat);
      ddsi1d =   ddpsi1(xd)*ddi_sav(iat);
      ddsi2d =   ddpsi2(xd);

      ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat);
      ddsi1md = -ddpsi1(mxd)*ddi_sav(iat);
      ddsi2md =  ddpsi2(mxd);


      // ..the free energy
      free  = h5(fi,
                 si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
                 si0d,   si1d,   si2d,   si0md,   si1md,   si2md);

      // ..derivative with respect to density
      df_d  = h5(fi,
                 si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
                 dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md);

      // ..derivative with respect to temperature
      df_t = h5(fi,
                dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md);

      // ..derivative with respect to density**2
      df_dd = h5(si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
                 ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md);

      // ..derivative with respect to temperature**2
      df_tt = h5(fi,
                 ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt,
                 si0d,   si1d,   si2d,   si0md,   si1md,   si2md);

      // ..derivative with respect to temperature and density
      df_dt = h5(fi,
                 dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
                 dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md);

      // ..now get the pressure derivative with density, chemical potential, and
      // ..electron positron number densities
      // ..get the interpolation weight functions
      si0t   =  xpsi0(xt);
      si1t   =  xpsi1(xt)*dt_sav(jat);

      si0mt  =  xpsi0(mxt);
      si1mt  =  -xpsi1(mxt)*dt_sav(jat);

      si0d   =  xpsi0(xd);
      si1d   =  xpsi1(xd)*dd_sav(iat);

      si0md  =  xpsi0(mxd);
      si1md  =  -xpsi1(mxd)*dd_sav(iat);

      // ..derivatives of weight functions
      dsi0t  = xdpsi0(xt)*dti_sav(jat);
      dsi1t  = xdpsi1(xt);

      dsi0mt = -xdpsi0(mxt)*dti_sav(jat);
      dsi1mt = xdpsi1(mxt);

      dsi0d  = xdpsi0(xd)*ddi_sav(iat);
      dsi1d  = xdpsi1(xd);

      dsi0md = -xdpsi0(mxd)*ddi_sav(iat);
      dsi1md = xdpsi1(mxd);

      // ..look in the pressure derivative only once
      fi(1)  = dpdf(iat,jat);
      fi(2)  = dpdf(iat+1,jat);
      fi(3)  = dpdf(iat,jat+1);
      fi(4)  = dpdf(iat+1,jat+1);
      fi(5)  = dpdft(iat,jat);
      fi(6)  = dpdft(iat+1,jat);
      fi(7)  = dpdft(iat,jat+1);
      fi(8)  = dpdft(iat+1,jat+1);
      fi(9)  = dpdfd(iat,jat);
      fi(10) = dpdfd(iat+1,jat);
      fi(11) = dpdfd(iat,jat+1);
      fi(12) = dpdfd(iat+1,jat+1);
      fi(13) = dpdfdt(iat,jat);
      fi(14) = dpdfdt(iat+1,jat);
      fi(15) = dpdfdt(iat,jat+1);
      fi(16) = dpdfdt(iat+1,jat+1);

      // ..pressure derivative with density
      dpepdd  = h3(fi,
                   si0t,   si1t,   si0mt,   si1mt,
                   si0d,   si1d,   si0md,   si1md);
      dpepdd  = max(ye * dpepdd,ZERO);

      // ..look in the electron chemical potential table only once
      fi(1)  = ef(iat,jat);
      fi(2)  = ef(iat+1,jat);
      fi(3)  = ef(iat,jat+1);
      fi(4)  = ef(iat+1,jat+1);
      fi(5)  = eft(iat,jat);
      fi(6)  = eft(iat+1,jat);
      fi(7)  = eft(iat,jat+1);
      fi(8)  = eft(iat+1,jat+1);
      fi(9)  = efd(iat,jat);
      fi(10) = efd(iat+1,jat);
      fi(11) = efd(iat,jat+1);
      fi(12) = efd(iat+1,jat+1);
      fi(13) = efdt(iat,jat);
      fi(14) = efdt(iat+1,jat);
      fi(15) = efdt(iat,jat+1);
      fi(16) = efdt(iat+1,jat+1);

      // ..electron chemical potential etaele
      etaele  = h3(fi,
                   si0t,   si1t,   si0mt,   si1mt,
                   si0d,   si1d,   si0md,   si1md);

      // ..derivative with respect to density
      x       = h3(fi,
                   si0t,   si1t,   si0mt,   si1mt,
                   dsi0d,  dsi1d,  dsi0md,  dsi1md);
      detadd  = ye * x;

      // ..derivative with respect to temperature
      detadt  = h3(fi,
                   dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
                   si0d,   si1d,   si0md,   si1md);

#ifdef EXTRA_THERMO
      // ..derivative with respect to abar and zbar
      detada = -x * din * ytot1;
      detadz =  x * den * ytot1;
#endif

      // ..look in the number density table only once
      fi(1)  = xf(iat,jat);
      fi(2)  = xf(iat+1,jat);
      fi(3)  = xf(iat,jat+1);
      fi(4)  = xf(iat+1,jat+1);
      fi(5)  = xft(iat,jat);
      fi(6)  = xft(iat+1,jat);
      fi(7)  = xft(iat,jat+1);
      fi(8)  = xft(iat+1,jat+1);
      fi(9)  = xfd(iat,jat);
      fi(10) = xfd(iat+1,jat);
      fi(11) = xfd(iat,jat+1);
      fi(12) = xfd(iat+1,jat+1);
      fi(13) = xfdt(iat,jat);
      fi(14) = xfdt(iat+1,jat);
      fi(15) = xfdt(iat,jat+1);
      fi(16) = xfdt(iat+1,jat+1);

      // ..electron + positron number densities
      xnefer   = h3(fi,
                    si0t,   si1t,   si0mt,   si1mt,
                    si0d,   si1d,   si0md,   si1md);

      // ..derivative with respect to density
      x        = h3(fi,
                    si0t,   si1t,   si0mt,   si1mt,
                    dsi0d,  dsi1d,  dsi0md,  dsi1md);
      x = max(x,ZERO);
      dxnedd   = ye * x;

      // ..derivative with respect to temperature
      dxnedt   = h3(fi,
                    dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
                    si0d,   si1d,   si0md,   si1md);

#ifdef EXTRA_THERMO
      // ..derivative with respect to abar and zbar
      dxneda = -x * din * ytot1;
      dxnedz =  x  * den * ytot1;
#endif

      // ..the desired electron-positron thermodynamic quantities

      // ..dpepdd at high temperatures and low densities is below the
      // ..floating point limit of the subtraction of two large terms.
      // ..since dpresdd doesn't enter the maxwell relations at all, use the
      // ..bicubic interpolation done above instead of this one
      x       = din * din;
      pele    = x * df_d;
      dpepdt  = x * df_dt;
      //     dpepdd  = ye * (x * df_dd + 2.0e0 * din * df_d);
      s       = dpepdd/ye - 2.0e0 * din * df_d;
#ifdef EXTRA_THERMO
      dpepda  = -ytot1 * (2.0e0 * pele + s * din);
      dpepdz  = den*ytot1*(2.0e0 * din * df_d  +  s);
#endif

      x       = ye * ye;
      sele    = -df_t * ye;
      dsepdt  = -df_tt * ye;
      dsepdd  = -df_dt * x;
#ifdef EXTRA_THERMO
      dsepda  = ytot1 * (ye * df_dt * din - sele);
      dsepdz  = -ytot1 * (ye * df_dt * den  + df_t);
#endif

      eele    = ye*free + temp * sele;
      deepdt  = temp * dsepdt;
      deepdd  = x * df_d + temp * dsepdd;
#ifdef EXTRA_THERMO
      deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda;
      deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz;
#endif

      // ..coulomb section:
      // ..initialize
      pcoul    = ZERO;
      dpcouldd = ZERO;
      dpcouldt = ZERO;
      dpcoulda = ZERO;
      dpcouldz = ZERO;
      ecoul    = ZERO;
      decouldd = ZERO;
      decouldt = ZERO;
      decoulda = ZERO;
      decouldz = ZERO;
      scoul    = ZERO;
      dscouldd = ZERO;
      dscouldt = ZERO;
      dscoulda = ZERO;
      dscouldz = ZERO;

      // ..uniform background corrections only
      // ..from yakovlev & shalybkov 1989
      // ..lami is the average ion seperation
      // ..plasg is the plasma coupling parameter
      z        = forth * pi;
      s        = z * xni;
      dsdd     = z * dxnidd;
      dsda     = z * dxnida;

      lami     = 1.0e0/s**onethird;
      inv_lami = 1.0e0/lami;
      z        = -onethird * lami;
      lamidd   = z * dsdd/s;
      lamida   = z * dsda/s;

      plasg    = zbar*zbar*esqu*ktinv*inv_lami;
      z        = -plasg * inv_lami;
      plasgdd  = z * lamidd;
      plasgda  = z * lamida;
      plasgdt  = -plasg*ktinv * kerg;
      plasgdz  = 2.0e0 * plasg/zbar;

      //     TURN ON/OFF COULOMB
      if ( do_coulomb ) {
        // ...yakovlev & shalybkov 1989 equations 82, 85, 86, 87
        if (plasg >= 1.0e0) {
          x        = plasg**(0.25e0);
          y        = avo_eos * ytot1 * kerg;
          ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1);
          pcoul    = onethird * den * ecoul;
          scoul    = (-y * (3.0e0*b1*x - 5.0e0*c1/x
                            + d1 * (log(plasg) - 1.0e0) - e1));

          y        = avo_eos*ytot1*kt*(a1 + 0.25e0/plasg*(b1*x - c1/x));
          decouldd = y * plasgdd;
          decouldt = y * plasgdt + ecoul/temp;
          decoulda = y * plasgda - ecoul/abar;
          decouldz = y * plasgdz;

          y        = onethird * den;
          dpcouldd = onethird * ecoul + y*decouldd;
          dpcouldt = y * decouldt;
          dpcoulda = y * decoulda;
          dpcouldz = y * decouldz;

          y        = (-avo_eos*kerg/(abar*plasg)* 
                      (0.75d0*b1*x+1.25d0*c1/x+d1));
          dscouldd = y * plasgdd;
          dscouldt = y * plasgdt;
          dscoulda = y * plasgda - scoul/abar;
          dscouldz = y * plasgdz;

          // ...yakovlev & shalybkov 1989 equations 102, 103, 104
        } else if (plasg < 1.0e0) {
          x        = plasg*sqrt(plasg);
          y        = plasg**b2;
          z        = c2 * x - onethird * a2 * y;
          pcoul    = -pion * z;
          ecoul    = 3.0e0 * pcoul/den;
          scoul    = -avo_eos/abar*kerg*(c2*x -a2*(b2-1.0e0)/b2*y);

          s        = 1.5e0*c2*x/plasg - onethird*a2*b2*y/plasg;
          dpcouldd = -dpiondd*z - pion*s*plasgdd;
          dpcouldt = -dpiondt*z - pion*s*plasgdt;
#ifdef EXTRA_THERMO
          dpcoulda = -dpionda*z - pion*s*plasgda;
          dpcouldz = -dpiondz*z - pion*s*plasgdz;
#endif

          s        = 3.0e0/den;
          decouldd = s * dpcouldd - ecoul/den;
          decouldt = s * dpcouldt;
          decoulda = s * dpcoulda;
          decouldz = s * dpcouldz;

          s        = (-avo_eos*kerg/(abar*plasg)*
                      (1.5e0*c2*x-a2*(b2-1.0e0)*y));
          dscouldd = s * plasgdd;
          dscouldt = s * plasgdt;
          dscoulda = s * plasgda - scoul/abar;
          dscouldz = s * plasgdz;
        }

        // Disable Coulomb corrections if they cause
        // the energy or pressure to go negative.

        p_temp = prad + pion + pele + pcoul;
        e_temp = erad + eion + eele + ecoul;

        if (p_temp <= ZERO .or. e_temp <= ZERO) {

                 pcoul    = ZERO;
                 dpcouldd = ZERO;
                 dpcouldt = ZERO;
                 dpcoulda = ZERO;
                 dpcouldz = ZERO;
                 ecoul    = ZERO;
                 decouldd = ZERO;
                 decouldt = ZERO;
                 decoulda = ZERO;
                 decouldz = ZERO;
                 scoul    = ZERO;
                 dscouldd = ZERO;
                 dscouldt = ZERO;
                 dscoulda = ZERO;
                 dscouldz = ZERO;

        }
      }

      // ..sum all the components
      pres    = prad + pion + pele + pcoul;
      ener    = erad + eion + eele + ecoul;
      entr    = srad + sion + sele + scoul;

      dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd;
      dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt;
#ifdef EXTRA_THERMO
      dpresda = dpradda + dpionda + dpepda + dpcoulda;
      dpresdz = dpraddz + dpiondz + dpepdz + dpcouldz;
#endif
      denerdd = deraddd + deiondd + deepdd + decouldd;
      denerdt = deraddt + deiondt + deepdt + decouldt;
#ifdef EXTRA_THERMO
      denerda = deradda + deionda + deepda + decoulda;
      denerdz = deraddz + deiondz + deepdz + decouldz;
#endif

      dentrdd = dsraddd + dsiondd + dsepdd + dscouldd;
      dentrdt = dsraddt + dsiondt + dsepdt + dscouldt;
#ifdef EXTRA_THERMO
      dentrda = dsradda + dsionda + dsepda + dscoulda;
      dentrdz = dsraddz + dsiondz + dsepdz + dscouldz;
#endif

      // ..the temperature and density exponents (c&g 9.81 9.82)
      // ..the specific heat at constant volume (c&g 9.92)
      // ..the third adiabatic exponent (c&g 9.93)
      // ..the first adiabatic exponent (c&g 9.97)
      // ..the second adiabatic exponent (c&g 9.105)
      // ..the specific heat at constant pressure (c&g 9.98)
      // ..and relativistic formula for the sound speed (c&g 14.29)
      zz    = pres*deni;
      zzi   = den/pres;
      chit  = temp/pres * dpresdt;
      chid  = dpresdd*zzi;
      cv    = denerdt;
      x     = zz * chit/(temp * cv);
      gam3  = x + 1.0e0;
      gam1  = chit*x + chid;
      nabad = x/gam1;
      gam2  = 1.0e0/(1.0e0 - nabad);
      cp    = cv * gam1/chid;
      z     = 1.0e0 + (ener + light2)*zzi;
      sound = clight * sqrt(gam1/z);

      // ..maxwell relations; each is zero if the consistency is perfect
      x   = den * den;
      dse = temp*dentrdt/denerdt - 1.0e0;
      dpe = (denerdd*x + temp*dpresdt)/pres - 1.0e0;
      dsp = -dentrdd*x/dpresdt - 1.0e0;

      ptot_row = pres;
      dpt_row = dpresdt;
      dpd_row = dpresdd;
#ifdef EXTRA_THERMO
      dpa_row = dpresda;
      dpz_row = dpresdz;
#endif
      dpe_row = dpresdt / denerdt;
      dpdr_e_row = dpresdd - dpresdt * denerdd / denerdt;

      etot_row = ener;
      det_row = denerdt;
      ded_row = denerdd;
#ifdef EXTRA_THERMO
      dea_row = denerda;
      dez_row = denerdz;
#endif

      stot_row = entr;
      dst_row = dentrdt;
      dsd_row = dentrdd;

      htot_row = ener + pres / den;
      dhd_row = denerdd + dpresdd / den - pres / den**2;
      dht_row = denerdt + dpresdt / den;

      pele_row = pele;
      ppos_row = ZERO;

      xne_row = xnefer;
      xnp_row = ZERO;

      etaele_row = etaele;

      cv_row = cv;
      cp_row = cp;
      cs_row = sound;
      gam1_row = gam1;

      if (converged)

        break;

      else if (single_iter) {

        if (dvar == itemp) {

          x = temp_row;
          smallx = smallt;
          xtol = ttol;

          if (var == ipres) {
            v    = ptot_row;
            dvdx = dpt_row;
          } else if (var == iener) {
            v    = etot_row;
            dvdx = det_row;
          } else if (var == ientr) {
            v    = stot_row;
            dvdx = dst_row;
          } else if (var == ienth) {
            v    = htot_row;
            dvdx = dht_row;
          } else {
            break;
          }

        } else { // dvar == density

          x = den_row;
          smallx = smalld;
          xtol = dtol;

          if (var == ipres) {
            v    = ptot_row;
            dvdx = dpd_row;
          } else if (var == iener) {
            v    = etot_row;
            dvdx = ded_row;
          } else if (var == ientr) {
            v    = stot_row;
            dvdx = dsd_row;
          } else if (var == ienth) {
            v    = htot_row;
            dvdx = dhd_row;
          } else {
            break;
          }

        }

        // Now do the calculation for the next guess for T/rho

        xnew = x - (v - v_want) / dvdx;

        // Don't let the temperature/density change by more than a factor of two
        xnew = max(0.5 * x, min(xnew, 2.0 * x));

        // Don't let us freeze/evacuate
        xnew = max(smallx, xnew);

        // Store the new temperature/density

        if (dvar == itemp) {
          temp_row = xnew;
        } else {
          den_row  = xnew;
        }

        // Compute the error from the last iteration

        error = abs( (xnew - x) / x );

        if (error < xtol) converged = true;
        else if (double_iter) {

          // Figure out which variables we're using

          told = temp_row;
          rold = den_row;

          if (var1 == ipres) {
            v1    = ptot_row;
            dv1dt = dpt_row;
            dv1dr = dpd_row;
          } else if (var1 == iener) {
            v1    = etot_row;
            dv1dt = det_row;
            dv1dr = ded_row;
          } else if (var1 == ientr) {
            v1    = stot_row;
            dv1dt = dst_row;
            dv1dr = dsd_row;
          } else if (var1 == ienth) {
            v1    = htot_row;
            dv1dt = dht_row;
            dv1dr = dhd_row;
          } else {
            break;
          }

          if (var2 == ipres) {
            v2    = ptot_row;
            dv2dt = dpt_row;
            dv2dr = dpd_row;
          } else if (var2 == iener) {
            v2    = etot_row;
            dv2dt = det_row;
            dv2dr = ded_row;
          } else if (var2 == ientr) {
            v2    = stot_row;
            dv2dt = dst_row;
            dv2dr = dsd_row;
          } else if (var2 == ienth) {
            v2    = htot_row;
            dv2dt = dht_row;
            dv2dr = dhd_row;
          } else {
            break;
          }

          // Two functions, f and g, to iterate over
          v1i = v1_want - v1;
          v2i = v2_want - v2;

          //
          // 0 = f + dfdr * delr + dfdt * delt
          // 0 = g + dgdr * delr + dgdt * delt
          //

          // note that dfi/dT = - df/dT
          delr = (-v1i*dv2dt + v2i*dv1dt) / (dv2dr*dv1dt - dv2dt*dv1dr);

          rnew = rold + delr;

          tnew = told + (v1i - dv1dr*delr) / dv1dt;

          // Don't let the temperature or density change by more
          // than a factor of two
          tnew = max(HALF * told, min(tnew, TWO * told));
          rnew = max(HALF * rold, min(rnew, TWO * rold));

          // Don't let us freeze or evacuate
          tnew = max(smallt, tnew);
          rnew = max(smalld, rnew);

          // Store the new temperature and density
          den_row  = rnew;
          temp_row = tnew;

          // Compute the errors
          error1 = abs( (rnew - rold) / rold );
          error2 = abs( (tnew - told) / told );

          if (error1 < dtol && error2 < ttol) converged = true;

        }
      }
    }

    state.T    = temp_row;
    state.rho  = den_row;

    state.p    = ptot_row;
    state.dpdT = dpt_row;
    state.dpdr = dpd_row;

#ifdef EXTRA_THERMO
    state.dpdA = dpa_row;
    state.dpdZ = dpz_row;
#endif

    state.dpde = dpe_row;
    state.dpdr_e = dpdr_e_row;

    state.e    = etot_row;
    state.dedT = det_row;
    state.dedr = ded_row;

#ifdef EXTRA_THERMO
    state.dedA = dea_row;
    state.dedZ = dez_row;
#endif

    state.s    = stot_row;
    state.dsdT = dst_row;
    state.dsdr = dsd_row;

    state.h    = htot_row;
    state.dhdR = dhd_row;
    state.dhdT = dht_row;

    state.pele = pele_row;
    state.ppos = ppos_row;

    state.xne = xne_row;
    state.xnp = xnp_row;

    state.eta = etaele_row;

    state.cv   = cv_row;
    state.cp   = cp_row;
    state.gam1 = gam1_row;
    // state.cs   = cs_row;

    // Take care of final housekeeping.

    // Count the positron contribution in the electron quantities.

    state.xne  = state.xne  + state.xnp;
    state.pele = state.pele + state.ppos;

    // Use the non-relativistic version of the sound speed, cs = sqrt(gam_1 * P / rho).
    // This replaces the relativistic version that comes out of helmeos.

    state.cs = sqrt(state.gam1 * state.p / state.rho);

    if (input_is_constant) {

      if (input == EosInput::rh) {

        state.h = v_want;

      } else if (input == EosInput::tp) {

        state.p = v_want;

      } else if (input == EosInput::rp) {

        state.p = v_want;

      } else if (input == EosInput::re) {

        state.e = v_want;

      } else if (input == EosInput::ps) {

        state.p = v1_want;
        state.s = v2_want;

      } else if (input == EosInput::ph) {

        state.p = v1_want;
        state.h = v2_want;

      } else if (input == EosInput::th) {

        state.h = v_want;

      }
    }
  }

  void actual_eos_init() {

    //  use extern_probin_module, only: eos_input_is_constant, use_eos_coulomb, eos_ttol, eos_dtol

    double dth, dt2, dti, dt2i;
    double dd, dd2, ddi, dd2i;
    double tsav, dsav;
    int i, j;
    int status;

    // Read in the runtime parameters (TODO)
    /*
        input_is_constant = eos_input_is_constant
        do_coulomb = use_eos_coulomb
        ttol = eos_ttol
        dtol = eos_dtol
    */

    if (ParallelDescriptor::IOProcessor()) {
      if (do_coulomb) {
        Print() << std::endl << "Initializing Helmholtz EOS and using Coulomb corrections." << std::endl;
      } else {
        Print() << std::endl << "Initializing Helmholtz EOS without using Coulomb corrections." << std::endl;
      }
    }

    // ..   read the helmholtz free energy table
    itmax = imax;
    jtmax = jmax;
    tlo   = 3.0e0;
    thi   = 13.0e0;
    tstp  = (thi - tlo)/float(jmax-1);
    tstpi = 1.0e0/tstp;
    dlo   = -12.0e0;
    dhi   = 15.0e0;
    dstp  = (dhi - dlo)/float(imax-1);
    dstpi = 1.0e0/dstp;

    for (int j = 1; j <= jmax; j++) {
      tsav = tlo + (j-1)*tstp;
      t(j) = pow(10.0e0, tsav);
      for (int i = 1; i <= imax; i++) {
        dsav = dlo + (i-1)*dstp;
        d(i) = pow(10.0e0, dsav);
      }
    }

    
    if (ParallelDescriptor::IOProcessor()) {

      // ..   open the table
      /*
        open(unit=2,file='helm_table.dat',status='old',iostat=status,action='read')
        if (status > 0) then

        call amrex_error('actual_eos_init: Failed to open helm_table.dat')

        endif

        !...  read in the free energy table
        do j=1,jmax
        do i=1,imax
        read(2,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
        fddt(i,j),fdtt(i,j),fddtt(i,j)
        end do
        end do

        !..   read the pressure derivative with density table
        do j = 1, jmax
        do i = 1, imax
        read(2,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
        end do
        end do

        !..   read the electron chemical potential table
        do j = 1, jmax
        do i = 1, imax
        read(2,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
        end do
        end do

        !..   read the number density table
        do j = 1, jmax
        do i = 1, imax
        read(2,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
        end do
        end do

        close(unit=2)

      */

    }

    auto BcastArray = [&](auto array){
                        ParallelDescriptor::Bcast(array.dataPtr(), array.size(),
                                                  ParallelDescriptor::IOProcessorNumber());
                      };

    BcastArray(f);
    BcastArray(fd);
    BcastArray(ft);
    BcastArray(fdd);
    BcastArray(ftt);
    BcastArray(fdt);
    BcastArray(fddt);
    BcastArray(fdtt);
    BcastArray(fddtt);
    BcastArray(dpdf);
    BcastArray(dpdfd);
    BcastArray(dpdft);
    BcastArray(dpdfdt);
    BcastArray(ef);
    BcastArray(efd);
    BcastArray(eft);
    BcastArray(efdt);
    BcastArray(xf);
    BcastArray(xfd);
    BcastArray(xft);
    BcastArray(xfdt);

    // ..   construct the temperature and density deltas and their inverses
    for (int j = 1; j < jmax; j++) {
      dth         = t(j+1) - t(j);
      dt2         = dth * dth;
      dti         = 1.0e0/dth;
      dt2i        = 1.0e0/dt2;
      dt_sav(j)   = dth;
      dt2_sav(j)  = dt2;
      dti_sav(j)  = dti;
      dt2i_sav(j) = dt2i;
    }

    for (int i = 1; i < imax; i++) {
      dd          = d(i+1) - d(i);
      dd2         = dd * dd;
      ddi         = 1.0e0/dd;
      dd2i        = 1.0e0/dd2;
      dd_sav(i)   = dd;
      dd2_sav(i)  = dd2;
      ddi_sav(i)  = ddi;
      dd2i_sav(i) = dd2i;
    }

    // Set up the minimum and maximum possible densities.

    mintemp = pow(10.e0, tlo);
    maxtemp = pow(10.e0, thi);
    mindens = pow(10.e0, dlo);
    maxdens = pow(10.e0, dhi);
  }

  // quintic hermite polynomial functions
  // psi0 and its derivatives
  AMREX_GPU_HOST_DEVICE AMREX_INLINE
  double psi0(double z)
  {
    double psi0r = pow(z, 3) * ( z * (-6.0e0*z + 15.0e0) -10.0e0) + 1.0e0;
    return psi0r;
  }

  AMREX_GPU_HOST_DEVICE AMREX_INLINE  
  double dpsi0(double z)
  {
    double dpsi0r = pow(z, 2) * ( z * (-30.0e0*z + 60.0e0) - 30.0e0);
    return dpsi0r;
  }

  AMREX_GPU_HOST_DEVICE AMREX_INLINE  
  double ddpsi0(z)
  {
    double ddpsi0r = z* ( z*( -120.0e0*z + 180.0e0) -60.0e0);
    return ddpsi0r;
  }

  // psi1 and its derivatives
  AMREX_GPU_HOST_DEVICE AMREX_INLINE    
  double psi1(z)
  {
    double psi1r = z* ( pow(z, 2) * ( z * (-3.0e0*z + 8.0e0) - 6.0e0) + 1.0e0);
    return psi1r;
  }

  AMREX_GPU_HOST_DEVICE AMREX_INLINE      
  double dpsi1(z)
  {
    double dpsi1r = z*z * ( z * (-15.0e0*z + 32.0e0) - 18.0e0) + 1.0e0;
    return dpsi1r;
  }

  AMREX_GPU_HOST_DEVICE AMREX_INLINE  
  double ddpsi1(z)
  {
    double ddpsi1r = z * (z * (-60.0e0*z + 96.0e0) -36.0e0);
    return ddpsi1r;
  }

  // psi2  and its derivatives
  AMREX_GPU_HOST_DEVICE AMREX_INLINE    
  double psi2(z)
  {
    double psi2r = 0.5e0*z*z*( z* ( z * (-z + 3.0e0) - 3.0e0) + 1.0e0);
    return psi2r;
  }

  AMREX_GPU_HOST_DEVICE AMREX_INLINE  
  double dpsi2(z)
  {
    double dpsi2r = 0.5e0*z*( z*(z*(-5.0e0*z + 12.0e0) - 9.0e0) + 2.0e0);
    return dpsi2r;
  }

  AMREX_GPU_HOST_DEVICE AMREX_INLINE    
  double ddpsi2(z)
  {
    double ddpsi2r = 0.5e0*(z*( z * (-20.0e0*z + 36.0e0) - 18.0e0) + 2.0e0);
    return ddpsi2r;
  }

  // biquintic hermite polynomial function
  AMREX_GPU_HOST_DEVICE AMREX_INLINE      
  double h5(const Array1D<double, 36>& fi,
            const double w0t, const double w1t, const double w2t, const double w0mt,
            const double w1mt, const double w2mt, const double w0d, const double w1d,
            const double w2d, const double w0md, const double w1md, const double w2md)
  {
    double h5r;

    h5r =  (fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t 
            + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt 
            + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t 
            + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt 
            + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t 
            + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt 
            + fi(13) *w1d*w0t   + fi(14) *w1md*w0t 
            + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt 
            + fi(17) *w2d*w0t   + fi(18) *w2md*w0t 
            + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt 
            + fi(21) *w1d*w1t   + fi(22) *w1md*w1t 
            + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt 
            + fi(25) *w2d*w1t   + fi(26) *w2md*w1t 
            + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt 
            + fi(29) *w1d*w2t   + fi(30) *w1md*w2t 
            + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt 
            + fi(33) *w2d*w2t   + fi(34) *w2md*w2t 
            + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt);
    return h5r;
  }


  // cubic hermite polynomial functions
  // psi0 & derivatives
  AMREX_GPU_HOST_DEVICE AMREX_INLINE        
  double xpsi0(z)
  {
    double xpsi0r = z * z * (2.0e0*z - 3.0e0) + 1.0;
    return xpsi0r;
  }

  AMREX_GPU_HOST_DEVICE AMREX_INLINE          
  doubl xdpsi0(z)
  {
    double xdpsi0r = z * (6.0e0*z - 6.0e0);
    return xdpsi0r;
  }


  // psi1 & derivatives
  AMREX_GPU_HOST_DEVICE AMREX_INLINE            
  double xpsi1(z)
  {
    double xpsi1r = z * ( z * (z - 2.0e0) + 1.0e0);
    return xpsi1r;
  }

  AMREX_GPU_HOST_DEVICE AMREX_INLINE            
  double xdpsi1(z)
  {
    double xdpsi1r = z * (3.0e0*z - 4.0e0) + 1.0e0;
    return xdpsi1r;
  }

  // bicubic hermite polynomial function
  AMREX_GPU_HOST_DEVICE AMREX_INLINE              
  double h3(const Array1D<double, 36>& fi,
            const double w0t, const double w1t, const double w0mt,
            const double w1mt, const double w0d, const double w1d,
            const double w0md, const double w1md)
  {
    double h3r = (fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t 
                  + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt 
                  + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t 
                  + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt 
                  + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t 
                  + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt 
                  + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t 
                  + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt);
    return h3r;
  }

};
#endif
