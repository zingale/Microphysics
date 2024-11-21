#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

void actual_network_init();

const std::string network_name = "pynucastro-cxx";

namespace network
{

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real bion () {
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the binding energy of the element

        if constexpr (spec == N) {
            return 0.0_rt;
        }
        else if constexpr (spec == H1) {
            return 0.0_rt;
        }
        else if constexpr (spec == He4) {
            return 28.295662457999697_rt;
        }
        else if constexpr (spec == C12) {
            return 92.16173498399803_rt;
        }
        else if constexpr (spec == N13) {
            return 94.10522604799917_rt;
        }
        else if constexpr (spec == N14) {
            return 104.65860734799753_rt;
        }
        else if constexpr (spec == O16) {
            return 127.6193154119992_rt;
        }
        else if constexpr (spec == F18) {
            return 137.36950247599816_rt;
        }
        else if constexpr (spec == Ne20) {
            return 160.64482384000075_rt;
        }
        else if constexpr (spec == Ne21) {
            return 167.40598973999658_rt;
        }
        else if constexpr (spec == Na22) {
            return 174.14457080400098_rt;
        }
        else if constexpr (spec == Na23) {
            return 186.56435240400242_rt;
        }
        else if constexpr (spec == Mg24) {
            return 198.2570479679962_rt;
        }
        else if constexpr (spec == Al27) {
            return 224.95193723199915_rt;
        }
        else if constexpr (spec == Si28) {
            return 236.53684539599638_rt;
        }
        else if constexpr (spec == P31) {
            return 262.9161999600037_rt;
        }
        else if constexpr (spec == S32) {
            return 271.78016372399725_rt;
        }
        else if constexpr (spec == Ar36) {
            return 306.7167469519991_rt;
        }
        else if constexpr (spec == Ca40) {
            return 342.05218528000114_rt;
        }
        else if constexpr (spec == Ti44) {
            return 375.47496160800074_rt;
        }
        else if constexpr (spec == Cr48) {
            return 411.4679399359957_rt;
        }
        else if constexpr (spec == Mn51) {
            return 440.321747199996_rt;
        }
        else if constexpr (spec == Fe52) {
            return 447.6996182639923_rt;
        }
        else if constexpr (spec == Fe54) {
            return 471.76475446399854_rt;
        }
        else if constexpr (spec == Fe56) {
            return 492.2599506639962_rt;
        }
        else if constexpr (spec == Co55) {
            return 476.82912552799826_rt;
        }
        else if constexpr (spec == Co56) {
            return 486.91094362799777_rt;
        }
        else if constexpr (spec == Co57) {
            return 498.28746172798856_rt;
        }
        else if constexpr (spec == Ni56) {
            return 483.9956965919919_rt;
        }
        else if constexpr (spec == Ni58) {
            return 506.4596327920008_rt;
        }


        // Return zero if we don't recognize the species.
        return 0.0_rt;
    }

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real mion () {
        static_assert(spec >= 1 && spec <= NumSpec);

        constexpr amrex::Real A = NetworkProperties::aion(spec);
        constexpr amrex::Real Z = NetworkProperties::zion(spec);

        return (A - Z) * C::Legacy::m_n + Z * (C::Legacy::m_p + C::Legacy::m_e) - bion<spec>() * C::Legacy::MeV2gr;
    }

    // Legacy (non-templated) interfaces

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real bion (int spec) {
        using namespace Species;

        amrex::Real b = 0.0_rt;

        // Set the binding energy of the element
        constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                b = bion<n>();
            }
        });

        return b;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real mion (int spec) {
        using namespace Species;

        amrex::Real m = 0.0_rt;

        constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                m = mion<n>();
            }
        });

        return m;
    }
}

namespace Rates
{

    enum NetworkRates
    {
        k_N13_to_p_C12 = 1,
        k_O16_to_He4_C12 = 2,
        k_F18_to_He4_N14 = 3,
        k_Ne20_to_He4_O16 = 4,
        k_Na22_to_p_Ne21 = 5,
        k_Na22_to_He4_F18 = 6,
        k_Mg24_to_p_Na23 = 7,
        k_Mg24_to_He4_Ne20 = 8,
        k_Si28_to_p_Al27 = 9,
        k_Si28_to_He4_Mg24 = 10,
        k_S32_to_p_P31 = 11,
        k_S32_to_He4_Si28 = 12,
        k_Fe52_to_p_Mn51 = 13,
        k_Fe52_to_He4_Cr48 = 14,
        k_Co55_to_He4_Mn51 = 15,
        k_Ni56_to_p_Co55 = 16,
        k_Ni56_to_He4_Fe52 = 17,
        k_C12_to_He4_He4_He4 = 18,
        k_p_C12_to_N13 = 19,
        k_He4_C12_to_O16 = 20,
        k_He4_N14_to_F18 = 21,
        k_He4_O16_to_Ne20 = 22,
        k_He4_F18_to_Na22 = 23,
        k_He4_Ne20_to_Mg24 = 24,
        k_p_Ne21_to_Na22 = 25,
        k_p_Na23_to_Mg24 = 26,
        k_He4_Mg24_to_Si28 = 27,
        k_p_Al27_to_Si28 = 28,
        k_He4_Si28_to_S32 = 29,
        k_p_P31_to_S32 = 30,
        k_He4_Cr48_to_Fe52 = 31,
        k_p_Mn51_to_Fe52 = 32,
        k_He4_Mn51_to_Co55 = 33,
        k_He4_Fe52_to_Ni56 = 34,
        k_p_Co55_to_Ni56 = 35,
        k_C12_C12_to_p_Na23 = 36,
        k_C12_C12_to_He4_Ne20 = 37,
        k_He4_N13_to_p_O16 = 38,
        k_p_O16_to_He4_N13 = 39,
        k_C12_O16_to_p_Al27 = 40,
        k_C12_O16_to_He4_Mg24 = 41,
        k_O16_O16_to_p_P31 = 42,
        k_O16_O16_to_He4_Si28 = 43,
        k_He4_F18_to_p_Ne21 = 44,
        k_He4_Ne20_to_p_Na23 = 45,
        k_p_Ne21_to_He4_F18 = 46,
        k_p_Na23_to_He4_Ne20 = 47,
        k_He4_Mg24_to_p_Al27 = 48,
        k_p_Al27_to_He4_Mg24 = 49,
        k_He4_Si28_to_p_P31 = 50,
        k_p_P31_to_He4_Si28 = 51,
        k_He4_Cr48_to_p_Mn51 = 52,
        k_p_Mn51_to_He4_Cr48 = 53,
        k_He4_Fe52_to_p_Co55 = 54,
        k_p_Co55_to_He4_Fe52 = 55,
        k_He4_He4_He4_to_C12 = 56,
        k_C12_C12_to_Mg24_modified = 57,
        k_O16_O16_to_S32_modified = 58,
        k_C12_O16_to_Si28_modified = 59,
        k_Co55_to_p_Fe54 = 60,
        k_Co56_to_n_Co55 = 61,
        k_Co57_to_n_Co56 = 62,
        k_Co57_to_p_Fe56 = 63,
        k_Ni58_to_p_Co57 = 64,
        k_Ni58_to_He4_Fe54 = 65,
        k_p_Fe54_to_Co55 = 66,
        k_He4_Fe54_to_Ni58 = 67,
        k_p_Fe56_to_Co57 = 68,
        k_n_Co55_to_Co56 = 69,
        k_n_Co56_to_Co57 = 70,
        k_p_Co57_to_Ni58 = 71,
        k_He4_Mn51_to_p_Fe54 = 72,
        k_p_Fe54_to_He4_Mn51 = 73,
        k_He4_Fe54_to_p_Co57 = 74,
        k_p_Fe56_to_n_Co56 = 75,
        k_He4_Co55_to_p_Ni58 = 76,
        k_n_Co56_to_p_Fe56 = 77,
        k_p_Co56_to_n_Ni56 = 78,
        k_p_Co57_to_He4_Fe54 = 79,
        k_n_Ni56_to_p_Co56 = 80,
        k_p_Ni58_to_He4_Co55 = 81,
        k_He4_S32_to_Ar36_removed = 82,
        k_He4_S32_to_p_Cl35_removed = 83,
        k_p_Cl35_to_Ar36_removed = 84,
        k_Ar36_to_He4_S32_removed = 85,
        k_Ar36_to_p_Cl35_removed = 86,
        k_p_Cl35_to_He4_S32_removed = 87,
        k_He4_Ar36_to_Ca40_removed = 88,
        k_He4_Ar36_to_p_K39_removed = 89,
        k_p_K39_to_Ca40_removed = 90,
        k_Ca40_to_He4_Ar36_removed = 91,
        k_Ca40_to_p_K39_removed = 92,
        k_p_K39_to_He4_Ar36_removed = 93,
        k_He4_Ca40_to_Ti44_removed = 94,
        k_He4_Ca40_to_p_Sc43_removed = 95,
        k_p_Sc43_to_Ti44_removed = 96,
        k_Ti44_to_He4_Ca40_removed = 97,
        k_Ti44_to_p_Sc43_removed = 98,
        k_p_Sc43_to_He4_Ca40_removed = 99,
        k_He4_Ti44_to_Cr48_removed = 100,
        k_He4_Ti44_to_p_V47_removed = 101,
        k_p_V47_to_Cr48_removed = 102,
        k_Cr48_to_He4_Ti44_removed = 103,
        k_Cr48_to_p_V47_removed = 104,
        k_p_V47_to_He4_Ti44_removed = 105,
        k_n_Fe52_to_Fe53_removed = 106,
        k_n_Fe53_to_Fe54_removed = 107,
        k_Fe54_to_n_Fe53_removed = 108,
        k_Fe53_to_n_Fe52_removed = 109,
        k_n_Fe54_to_Fe55_removed = 110,
        k_n_Fe55_to_Fe56_removed = 111,
        k_Fe56_to_n_Fe55_removed = 112,
        k_Fe55_to_n_Fe54_removed = 113,
        k_n_Ni56_to_Ni57_removed = 114,
        k_n_Ni57_to_Ni58_removed = 115,
        k_Ni58_to_n_Ni57_removed = 116,
        k_Ni57_to_n_Ni56_removed = 117,
        k_Co56_to_Fe56 = 118,
        k_Co56_to_Ni56 = 119,
        k_Fe56_to_Co56 = 120,
        k_n_to_p = 121,
        k_Ni56_to_Co56 = 122,
        k_p_to_n = 123,
        k_S32_He4_to_Ar36_approx = 124,
        k_Ar36_to_S32_He4_approx = 125,
        k_Ar36_He4_to_Ca40_approx = 126,
        k_Ca40_to_Ar36_He4_approx = 127,
        k_Ca40_He4_to_Ti44_approx = 128,
        k_Ti44_to_Ca40_He4_approx = 129,
        k_Ti44_He4_to_Cr48_approx = 130,
        k_Cr48_to_Ti44_He4_approx = 131,
        k_Fe52_n_n_to_Fe54_approx = 132,
        k_Fe54_to_Fe52_n_n_approx = 133,
        k_Fe54_n_n_to_Fe56_approx = 134,
        k_Fe56_to_Fe54_n_n_approx = 135,
        k_Ni56_n_n_to_Ni58_approx = 136,
        k_Ni58_to_Ni56_n_n_approx = 137,
        NumRates = k_Ni58_to_Ni56_n_n_approx
    };

    // number of reaclib rates

    const int NrateReaclib = 117;

    // number of tabular rates

    const int NrateTabular = 6;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "N13_to_p_C12",  // 1,
        "O16_to_He4_C12",  // 2,
        "F18_to_He4_N14",  // 3,
        "Ne20_to_He4_O16",  // 4,
        "Na22_to_p_Ne21",  // 5,
        "Na22_to_He4_F18",  // 6,
        "Mg24_to_p_Na23",  // 7,
        "Mg24_to_He4_Ne20",  // 8,
        "Si28_to_p_Al27",  // 9,
        "Si28_to_He4_Mg24",  // 10,
        "S32_to_p_P31",  // 11,
        "S32_to_He4_Si28",  // 12,
        "Fe52_to_p_Mn51",  // 13,
        "Fe52_to_He4_Cr48",  // 14,
        "Co55_to_He4_Mn51",  // 15,
        "Ni56_to_p_Co55",  // 16,
        "Ni56_to_He4_Fe52",  // 17,
        "C12_to_He4_He4_He4",  // 18,
        "p_C12_to_N13",  // 19,
        "He4_C12_to_O16",  // 20,
        "He4_N14_to_F18",  // 21,
        "He4_O16_to_Ne20",  // 22,
        "He4_F18_to_Na22",  // 23,
        "He4_Ne20_to_Mg24",  // 24,
        "p_Ne21_to_Na22",  // 25,
        "p_Na23_to_Mg24",  // 26,
        "He4_Mg24_to_Si28",  // 27,
        "p_Al27_to_Si28",  // 28,
        "He4_Si28_to_S32",  // 29,
        "p_P31_to_S32",  // 30,
        "He4_Cr48_to_Fe52",  // 31,
        "p_Mn51_to_Fe52",  // 32,
        "He4_Mn51_to_Co55",  // 33,
        "He4_Fe52_to_Ni56",  // 34,
        "p_Co55_to_Ni56",  // 35,
        "C12_C12_to_p_Na23",  // 36,
        "C12_C12_to_He4_Ne20",  // 37,
        "He4_N13_to_p_O16",  // 38,
        "p_O16_to_He4_N13",  // 39,
        "C12_O16_to_p_Al27",  // 40,
        "C12_O16_to_He4_Mg24",  // 41,
        "O16_O16_to_p_P31",  // 42,
        "O16_O16_to_He4_Si28",  // 43,
        "He4_F18_to_p_Ne21",  // 44,
        "He4_Ne20_to_p_Na23",  // 45,
        "p_Ne21_to_He4_F18",  // 46,
        "p_Na23_to_He4_Ne20",  // 47,
        "He4_Mg24_to_p_Al27",  // 48,
        "p_Al27_to_He4_Mg24",  // 49,
        "He4_Si28_to_p_P31",  // 50,
        "p_P31_to_He4_Si28",  // 51,
        "He4_Cr48_to_p_Mn51",  // 52,
        "p_Mn51_to_He4_Cr48",  // 53,
        "He4_Fe52_to_p_Co55",  // 54,
        "p_Co55_to_He4_Fe52",  // 55,
        "He4_He4_He4_to_C12",  // 56,
        "C12_C12_to_Mg24_modified",  // 57,
        "O16_O16_to_S32_modified",  // 58,
        "C12_O16_to_Si28_modified",  // 59,
        "Co55_to_p_Fe54",  // 60,
        "Co56_to_n_Co55",  // 61,
        "Co57_to_n_Co56",  // 62,
        "Co57_to_p_Fe56",  // 63,
        "Ni58_to_p_Co57",  // 64,
        "Ni58_to_He4_Fe54",  // 65,
        "p_Fe54_to_Co55",  // 66,
        "He4_Fe54_to_Ni58",  // 67,
        "p_Fe56_to_Co57",  // 68,
        "n_Co55_to_Co56",  // 69,
        "n_Co56_to_Co57",  // 70,
        "p_Co57_to_Ni58",  // 71,
        "He4_Mn51_to_p_Fe54",  // 72,
        "p_Fe54_to_He4_Mn51",  // 73,
        "He4_Fe54_to_p_Co57",  // 74,
        "p_Fe56_to_n_Co56",  // 75,
        "He4_Co55_to_p_Ni58",  // 76,
        "n_Co56_to_p_Fe56",  // 77,
        "p_Co56_to_n_Ni56",  // 78,
        "p_Co57_to_He4_Fe54",  // 79,
        "n_Ni56_to_p_Co56",  // 80,
        "p_Ni58_to_He4_Co55",  // 81,
        "He4_S32_to_Ar36_removed",  // 82,
        "He4_S32_to_p_Cl35_removed",  // 83,
        "p_Cl35_to_Ar36_removed",  // 84,
        "Ar36_to_He4_S32_removed",  // 85,
        "Ar36_to_p_Cl35_removed",  // 86,
        "p_Cl35_to_He4_S32_removed",  // 87,
        "He4_Ar36_to_Ca40_removed",  // 88,
        "He4_Ar36_to_p_K39_removed",  // 89,
        "p_K39_to_Ca40_removed",  // 90,
        "Ca40_to_He4_Ar36_removed",  // 91,
        "Ca40_to_p_K39_removed",  // 92,
        "p_K39_to_He4_Ar36_removed",  // 93,
        "He4_Ca40_to_Ti44_removed",  // 94,
        "He4_Ca40_to_p_Sc43_removed",  // 95,
        "p_Sc43_to_Ti44_removed",  // 96,
        "Ti44_to_He4_Ca40_removed",  // 97,
        "Ti44_to_p_Sc43_removed",  // 98,
        "p_Sc43_to_He4_Ca40_removed",  // 99,
        "He4_Ti44_to_Cr48_removed",  // 100,
        "He4_Ti44_to_p_V47_removed",  // 101,
        "p_V47_to_Cr48_removed",  // 102,
        "Cr48_to_He4_Ti44_removed",  // 103,
        "Cr48_to_p_V47_removed",  // 104,
        "p_V47_to_He4_Ti44_removed",  // 105,
        "n_Fe52_to_Fe53_removed",  // 106,
        "n_Fe53_to_Fe54_removed",  // 107,
        "Fe54_to_n_Fe53_removed",  // 108,
        "Fe53_to_n_Fe52_removed",  // 109,
        "n_Fe54_to_Fe55_removed",  // 110,
        "n_Fe55_to_Fe56_removed",  // 111,
        "Fe56_to_n_Fe55_removed",  // 112,
        "Fe55_to_n_Fe54_removed",  // 113,
        "n_Ni56_to_Ni57_removed",  // 114,
        "n_Ni57_to_Ni58_removed",  // 115,
        "Ni58_to_n_Ni57_removed",  // 116,
        "Ni57_to_n_Ni56_removed",  // 117,
        "Co56_to_Fe56",  // 118,
        "Co56_to_Ni56",  // 119,
        "Fe56_to_Co56",  // 120,
        "n_to_p",  // 121,
        "Ni56_to_Co56",  // 122,
        "p_to_n",  // 123,
        "S32_He4_to_Ar36_approx",  // 124,
        "Ar36_to_S32_He4_approx",  // 125,
        "Ar36_He4_to_Ca40_approx",  // 126,
        "Ca40_to_Ar36_He4_approx",  // 127,
        "Ca40_He4_to_Ti44_approx",  // 128,
        "Ti44_to_Ca40_He4_approx",  // 129,
        "Ti44_He4_to_Cr48_approx",  // 130,
        "Cr48_to_Ti44_He4_approx",  // 131,
        "Fe52_n_n_to_Fe54_approx",  // 132,
        "Fe54_to_Fe52_n_n_approx",  // 133,
        "Fe54_n_n_to_Fe56_approx",  // 134,
        "Fe56_to_Fe54_n_n_approx",  // 135,
        "Ni56_n_n_to_Ni58_approx",  // 136,
        "Ni58_to_Ni56_n_n_approx"  // 137,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int H1_index = 1;
    constexpr int N_index = 0;
    constexpr int He4_index = 2;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
