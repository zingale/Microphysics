#ifndef _vode_rhs_H_
#define _vode_rhs_H_

#include <network.H>
#include <actual_network.H>
#include <actual_rhs.H>
#include <burn_type.H>
#include <extern_parameters.H>
#include <vode_type.H>

// The rhs routine provides the right-hand-side for the DVODE solver.
// This is a generic interface that calls the specific RHS routine in the
// network you're actually using.

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs (const Real time, burn_t& state, dvode_t& vode_state, RArray1D& ydot)
{

    // We are integrating a system of
    //
    // y(1:NumSpec) = dX/dt
    // y(net_itemp) = dT/dt
    // y(net_ienuc) = denuc/dt

    // Only do the burn if the incoming temperature is within the temperature
    // bounds. Otherwise set the RHS to zero and return.

    if (state.T <= EOSData::mintemp || state.T >= MAX_TEMP) {

        for (int n = 1; n <= VODE_NEQS; ++n) {
            ydot(n) = 0.0_rt;
        }

        return;

    }

    // Fix the state as necessary.    

    clean_state(vode_state);

    // Update the thermodynamics as necessary.

    update_thermodynamics(state, vode_state);

    // Call the specific network routine to get the RHS.

    vode_to_burn(vode_state, state);

    actual_rhs(state, ydot);

    // We integrate X, not Y
    for (int n = 1; n <= NumSpec; ++n) {
        ydot(n) *= aion[n-1];
    }

    // Allow temperature and energy integration to be disabled.
    if (!integrate_temperature) {
        ydot(net_itemp) = 0.0_rt;
    }

    if (!integrate_energy) {
        ydot(net_ienuc) = 0.0_rt;
    }

    // apply fudge factor:
    if (react_boost > 0.0_rt) {
        for (int n = 1; n <= VODE_NEQS; ++n) {
            ydot(n) *= react_boost;
        }
    }

    burn_to_vode(state, vode_state);

}



// Analytical Jacobian
template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac (burn_t& state, dvode_t& vode_state, MatrixType& pd)
{
    // NOTE: the time at which to evaluate the Jacobian is not
    // explicitly passed. VODE always evaluates the analytic
    // Jacobian at vode_state.tn, and we pass vode_state.

    // Only do the burn if the incoming temperature is within the temperature
    // bounds. Otherwise set the Jacobian to zero and return.

    if (state.T <= EOSData::mintemp || state.T >= MAX_TEMP) {

        for (int j = 1; j <= VODE_NEQS; ++j) {
            for (int i = 1; i <= VODE_NEQS; ++i) {
                pd(i,j) = 0.0_rt;
            }
        }

        return;

    }

    // Call the specific network routine to get the Jacobian.

    vode_to_burn(vode_state, state);

    actual_jac(state, pd);

    // We integrate X, not Y
    for (int j = 1; j <= NumSpec; ++j) {
        for (int i = 1; i <= VODE_NEQS; ++i) {
            pd.mul(j, i, aion[j-1]);
            pd.mul(i, j, aion_inv[j-1]);
        }
    }

    // apply fudge factor:
    if (react_boost > 0.0_rt) {
        for (int j = 1; j <= VODE_NEQS; ++j) {
            for (int i = 1; i <= VODE_NEQS; ++i) {
                pd.mul(i, j, react_boost);
            }
        }
    }

    // Allow temperature and energy integration to be disabled.
    if (!integrate_temperature) {
        for (int j = 1; j <= VODE_NEQS; ++j) {
            pd(net_itemp,j) = 0.0_rt;
        }
    }

    if (!integrate_energy) {
        for (int j = 1; j <= VODE_NEQS; ++j) {
            pd(net_ienuc,j) = 0.0_rt;
        }
    }

    burn_to_vode(state, vode_state);

}

#endif
