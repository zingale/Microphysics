#ifndef NSE_UPDATE_H
#define NSE_UPDATE_H

#include <iostream>
#include <fstream>
#include <actual_network.H>

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <extern_parameters.H>

#include <burn_type.H>
#include <eos.H>

#ifdef NSE_TABLE
#include <nse_table.H>
#endif
#ifdef NSE_NET
#include <nse_solver.H>
#endif

using namespace amrex;

#ifdef SIMPLIFIED_SDC
///
/// update the state due to NSE changes for the simplified-SDC case
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void sdc_nse_burn(burn_t& state, const Real dt) {

    using namespace AuxZero;

    state.success = true;
    state.n_rhs = 0;
    state.n_jac = 0;

    // call the NSE table to get (dYe/dt)^n
    Real abar_out;
    Real dq_out;
    Real dyedt;
    Real X[NumSpec];

    Real ye_in = state.y[SFX+iye] / state.rho;

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

#if defined (NSE_TABLE)
    // get the current NSE state from the table

    nse_interp(T_in, state.rho, ye_in,
               abar_out, dq_out, dyedt, X);

#elif defined (NSE_NET)
    // solve for the NSE state directly

    // compute abar, dq_out (B/A), and dyedt

#endif

    // density and momentum have no reactive sources
    Real rho_old = state.y[SRHO];
    Real rhoe_old = state.y[SEINT];

    Real rho_bea_old = state.y[SFX+ibea];
    Real rho_abar_old = state.y[SFX+iabar];
    Real rho_ye_old = state.y[SFX+iye];

    // update the non-reactive state
    state.y[SRHO] += dt * state.ydot_a[SRHO];
    state.y[SMX] += dt * state.ydot_a[SMX];
    state.y[SMY] += dt * state.ydot_a[SMY];
    state.y[SMZ] += dt * state.ydot_a[SMZ];

    Real rho_new = state.y[SRHO];

    // predict the state at the midpoint in time
    Real rho_half = rho_old + 0.5_rt * dt * state.ydot_a[SRHO];

    // we are assuming that since we start out at NSE, that the energy
    // generation is ~ 0
    Real rhoe_half = rhoe_old + 0.5_rt * dt * state.ydot_a[SEINT];

    Real rho_ye_half = rho_ye_old + 0.5_rt * dt * state.ydot_a[SFX+iye] + 0.5 * dt * rho_old * dyedt;
    Real rho_abar_half = rho_abar_old + 0.5_rt * dt * state.ydot_a[SFX+iabar];

    // get the updated temperature here

    eos_re_t eos_state;
    eos_state.rho = rho_half;
    eos_state.e = rhoe_half / rho_half;
    eos_state.aux[iye] = rho_ye_half / rho_half;
    eos_state.aux[iabar] = rho_abar_half / rho_half;
    eos_state.T = 1.e6_rt;

    eos(eos_input_re, eos_state);

    // find the NSE state at the midpoint in time

    nse_interp(eos_state.T, eos_state.rho, eos_state.aux[iye],
               abar_out, dq_out, dyedt, X);

    // compute the sources at the midpoint in time

    Real dyedt_half = dyedt;


    // compute the energy release -- we need to remove the advective part

    Real rho_bea_star = rho_bea_old + 0.5_rt * dt * state.ydot_a[SFX+ibea];
    Real bea_star = rho_bea_star / rho_half;

    Real dbea = dq_out - bea_star;

    // convert the energy to erg / cm**3
    Real rho_enucdot_half  = 0.5_rt * (rho_old + rho_half) * dbea * C::MeV2eV * C::ev2erg * C::n_A / dt;

    Real rho_abar_star = rho_abar_old + 0.5_rt * dt * state.ydot_a[SFX+iabar];
    Real abar_star = rho_abar_star / rho_half;

    Real dabardt_half = (abar_out - abar_star) / dt;

    // finally update everything to the final time

    Real rhoe_new = rhoe_old + dt * state.ydot_a[SEINT] + dt * rho_enucdot_half;

    Real rho_ye_new = rho_ye_old + dt * state.ydot_a[SFX+iye] + dt * rho_half * dyedt_half;
    Real rho_abar_new = rho_abar_old + dt * state.ydot_a[SFX+iabar] + rho_half * dabardt_half;

    // get the final temperature

    eos_state.rho = rho_new;
    eos_state.e = rhoe_new / rho_new;
    eos_state.aux[iye] = rho_ye_new / rho_new;
    eos_state.aux[iabar] = rho_abar_new / rho_new;

    eos(eos_input_re, eos_state);

    // find the NSE state at the final time

    nse_interp(eos_state.T, eos_state.rho, eos_state.aux[iye],
               abar_out, dq_out, dyedt, X);


    // find the final energy release

    rho_bea_star = rho_bea_old + dt * state.ydot_a[SFX+ibea];
    bea_star = rho_bea_star / rho_new;

    dbea = dq_out - bea_star;

    // convert the energy to erg / cm**3
    Real rho_enucdot  = rho_half * dbea * C::MeV2eV * C::ev2erg * C::n_A / dt;


    // the new mass fractions are just those that come from the table
    // make sure they are normalized
    Real sum_X{0.0_rt};
    for (int n = 0; n < NumSpec; ++n) {
        X[n] = amrex::max(small_x, amrex::min(1.0_rt, X[n]));
        sum_X += X[n];
    }

    for (int n = 0; n < NumSpec; ++n) {
        X[n] /= sum_X;
    }

    for (int n = 0; n < NumSpec; ++n) {
        state.y[SFS+n] = eos_state.rho * X[n];
    }

    // aux data comes from the table (except Ye, which we predicted)

    state.y[SFX+iye] = eos_state.rho * eos_state.aux[iye];
    state.y[SFX+iabar] = eos_state.rho * abar_out;
    state.y[SFX+ibea] = eos_state.rho * dq_out;


    // density and momenta have already been updated

    // update the total and internal energy now

    state.y[SEINT] += dt * state.ydot_a[SEINT] + dt * rho_enucdot;
    state.y[SEDEN] += dt * state.ydot_a[SEDEN] + dt * rho_enucdot;

}

#else
///
/// update the state due to NSE changes for the Strang-split case
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_burn(burn_t& state, const Real dt) {

#if defined(NSE_TABLE)

  using namespace AuxZero;
  
  // use the NSE table
  Real abar_out;
  Real dq_out;
  Real dyedt;
  Real X[NumSpec];

  // get the energy -- we are assuming that rho, T are valid on input

  eos(eos_input_rt, state);

  // if we are doing drive_initial_convection, we want to use
  // the temperature that comes in through T_fixed

  Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

  // call the NSE table using the * state to get the t^{n+1}
  // source estimates.  The thermodynamnics here is specified
  // in terms of the auxillary composition, Ye, abar, and B/A

  nse_interp(T_in, state.rho, state.aux[iye],
             abar_out, dq_out, dyedt, X);

  // update Ye

  state.aux[iye] += dt * dyedt;

  // now get the composition from the table using the upated Ye

  nse_interp(T_in, state.rho, state.aux[iye],
             abar_out, dq_out, dyedt, X);


  // this is MeV / nucleon -- here aux has not yet been updated, so we
  // access the old binding energy
  Real deltaq = dq_out - state.aux[ibea];

  // under-relaxation / inertia (see Ma et el. 2013)
  deltaq = eta * deltaq;

  state.aux[ibea] += deltaq;
  state.aux[iabar] = abar_out;

#elif defined(NSE_NET)

  // solve for the NSE state -- first compute Ye
  state.y_e = 0.0_rt;
  for (int n = 0; n < NumSpec; ++n) {
      state.y_e += zion[n] * state.xn[n] * aion_inv[n];
  }

  auto nse_state = get_actual_nse_state(state);

  // for now, we'll take dyedt = 0.  Later we can evaluate this by
  // calling the RHS with the NSE state.
  Real dyedt = 0.0_rt;

  // compute the change in binding energy -- this is the energy release
  // we want to do sum {B dY}, where Y is the molar fraction.
  // this will be in MeV / nucleon
  Real deltaq = 0.0_rt;
  for (int n = 0; n < NumSpec; ++n) {
      deltaq += network::bion(n+1) * aion_inv[n] * (nse_state.xn[n] - state.xn[n]);
  }

  // under-relaxation / inertia (see Ma et el. 2013)
  deltaq = eta * deltaq;

#endif

  state.success = true;
  state.n_rhs = 0;
  state.n_jac = 0;

  // convert the energy to erg / g
  Real enuc = deltaq * C::MeV2eV * C::ev2erg * C::n_A;

  state.e = enuc + state.e;


  // store the new composition

#if defined(NSE_TABLE)
  for (int n = 0; n < NumSpec; ++n) {
    state.xn[n] = X[n];
  }
#elif defined(NSE_NET)
  for (int n = 0; n < NumSpec; ++n) {
    state.xn[n] = nse_state.xn[n];
  }
#endif


}
#endif

#endif
