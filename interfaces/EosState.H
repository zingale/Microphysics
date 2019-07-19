#ifndef EOS_STATE_H
#define EOS_STATE_H
#include "AMReX_REAL.H"

using namespace amrex;

enum class EosInput {rt, rh, tp, rp, re, ps, ph, th};
enum class EosIndex {temp=1, dens, ener, enth, entr, pres};
enum class EosError {general=1, input, iter_conv,
                     neg_e, neg_p, neg_h, neg_s,
                     iter_var, init, init_xn, out_of_bounds, not_implemented};

class EosState {
public:
  Real rho;
  Real T;
  Real p;
  Real e;
  Real h;
  Real s;
  Array1D<Real, nspec> xn;
  Array1D<Real, naux> aux;

  Real dpdT;
  Real dpdr;
  Real dedT;
  Real dedr;
  Real dhdT;
  Real dhdr;
  Real dsdT;
  Real dsdr;
  Real dpde;
  Real dpdr_e;

  Real cv;
  Real cp;
  Real xne;
  Real xnp;
  Real eta;
  Real pele;
  Real ppos;
  Real mu;
  Real mu_e;
  Real y_e;

#ifdef EXTRA_THERMO
  Array1D<Real, nspec> dedX;
  Array1D<Real, nspec> dpdX;
  Array1D<Real, nspec> dhdX;
#endif

  Real gam1;
  Real cs;

  Real abar;
  Real zbar;

#ifdef EXTRA_THERMO
  Real dpdA;
  Real dpdZ;
  Real dedA;
  Real dedZ;
#endif

  Real conductivity;
};
#endif
