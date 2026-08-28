#ifndef MODIFIED_RATES_H
#define MODIFIED_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>
#include <temperature_table_rates.H>
#include <reaclib_rates.H>

namespace modified_rates {

    using namespace reaclib_rates;
    using namespace temp_tabular;

    template <int do_T_derivatives>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne19_3p_to_Mg22_beta_pos_modified(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

        // Ne19 + 3 p --> Mg22 (calls the underlying rate)
        // represents the sequence: Ne19(p,γ)Na20(p,γ)Mg21(β+)Na21(p,γ)Mg22

        rate_Ne19_p_to_Na20_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
    }

    template <int do_T_derivatives>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne20_p_p_to_Mg22_modified(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

        // Ne20 + p + p --> Mg22 (calls the underlying rate)
        // represents the sequence: Ne20(p,γ)Na21(p,γ)Mg22

        rate_Ne20_p_to_Na21_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
    }


    template <int do_T_derivatives, typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void
    fill_rates([[maybe_unused]] const tf_t& tfactors,
               [[maybe_unused]] T& rate_eval)
    {

        [[maybe_unused]] amrex::Real rate;
        [[maybe_unused]] amrex::Real drate_dT;

        {
            // Ne19_3p_to_Mg22_beta_pos_modified

            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ne19);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ne19);
            }
#endif
            rate_Ne19_3p_to_Mg22_beta_pos_modified<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
            rate_eval.screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified) = rate;
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                rate_eval.dscreened_rates_dT(k_Ne19_3p_to_Mg22_beta_pos_modified) = drate_dT;
            }
        }

        {
            // Ne20_p_p_to_Mg22_modified

            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ne20);
            }
#endif
            rate_Ne20_p_p_to_Mg22_modified<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
            rate_eval.screened_rates(k_Ne20_p_p_to_Mg22_modified) = rate;
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                rate_eval.dscreened_rates_dT(k_Ne20_p_p_to_Mg22_modified) = drate_dT;
            }
        }


    }
}
#endif
