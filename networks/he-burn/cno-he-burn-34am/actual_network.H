#ifndef actual_network_H
#define actual_network_H

#include <random>
#include <string_view>

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <extern_parameters.H>
#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex::literals;


inline
void actual_network_init() {


}

constexpr std::string_view network_name = "pynucastro-cxx";

namespace network
{

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real bion () {
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the binding energy of the element

        if constexpr (spec == H1) {
            return 0.0_rt;
        }
        else if constexpr (spec == He4) {
            return 28.295662457999697_rt;
        }
        else if constexpr (spec == C12) {
            return 92.16173498399803_rt;
        }
        else if constexpr (spec == C13) {
            return 97.10804378399916_rt;
        }
        else if constexpr (spec == N13) {
            return 94.10522604799917_rt;
        }
        else if constexpr (spec == N14) {
            return 104.65860734799753_rt;
        }
        else if constexpr (spec == N15) {
            return 115.49190414799887_rt;
        }
        else if constexpr (spec == O14) {
            return 98.73189611199996_rt;
        }
        else if constexpr (spec == O15) {
            return 111.95539521199862_rt;
        }
        else if constexpr (spec == O16) {
            return 127.6193154119992_rt;
        }
        else if constexpr (spec == O17) {
            return 131.76239561199873_rt;
        }
        else if constexpr (spec == O18) {
            return 139.8077658120019_rt;
        }
        else if constexpr (spec == F17) {
            return 128.21958437599824_rt;
        }
        else if constexpr (spec == F18) {
            return 137.36950247599816_rt;
        }
        else if constexpr (spec == F19) {
            return 147.80136567599766_rt;
        }
        else if constexpr (spec == Ne18) {
            return 132.14265544000227_rt;
        }
        else if constexpr (spec == Ne19) {
            return 143.7795235400008_rt;
        }
        else if constexpr (spec == Ne20) {
            return 160.64482384000075_rt;
        }
        else if constexpr (spec == Ne21) {
            return 167.40598973999658_rt;
        }
        else if constexpr (spec == Ne22) {
            return 177.77024384000106_rt;
        }
        else if constexpr (spec == Na21) {
            return 163.07672270400144_rt;
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
        else if constexpr (spec == Fe52) {
            return 447.6996182639923_rt;
        }
        else if constexpr (spec == Ni56) {
            return 483.9956965919919_rt;
        }


        // Return zero if we don't recognize the species.
        return 0.0_rt;
    }

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real mion () {
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the mass of the element

        if constexpr (spec == H1) {
            return 1.6735328377636005e-24_rt;
        }
        else if constexpr (spec == He4) {
            return 6.646479071584587e-24_rt;
        }
        else if constexpr (spec == C12) {
            return 1.99264687992e-23_rt;
        }
        else if constexpr (spec == C13) {
            return 2.1592578700883826e-23_rt;
        }
        else if constexpr (spec == N13) {
            return 2.1596537049448796e-23_rt;
        }
        else if constexpr (spec == N14) {
            return 2.3252651436495096e-23_rt;
        }
        else if constexpr (spec == N15) {
            return 2.490826682883827e-23_rt;
        }
        else if constexpr (spec == O14) {
            return 2.3261822098665393e-23_rt;
        }
        else if constexpr (spec == O15) {
            return 2.4913176568383383e-23_rt;
        }
        else if constexpr (spec == O16) {
            return 2.6560180592333686e-23_rt;
        }
        else if constexpr (spec == O17) {
            return 2.822772237905709e-23_rt;
        }
        else if constexpr (spec == O18) {
            return 2.9888307701990364e-23_rt;
        }
        else if constexpr (spec == F17) {
            return 2.8232643353472436e-23_rt;
        }
        else if constexpr (spec == F18) {
            return 2.989125964092377e-23_rt;
        }
        else if constexpr (spec == F19) {
            return 3.1547590653659727e-23_rt;
        }
        else if constexpr (spec == Ne18) {
            return 2.9899182681834417e-23_rt;
        }
        else if constexpr (spec == Ne19) {
            return 3.155336557821978e-23_rt;
        }
        else if constexpr (spec == Ne20) {
            return 3.3198227947612416e-23_rt;
        }
        else if constexpr (spec == Ne21) {
            return 3.4861102572650884e-23_rt;
        }
        else if constexpr (spec == Ne22) {
            return 3.651755410955497e-23_rt;
        }
        else if constexpr (spec == Na21) {
            return 3.486742553187392e-23_rt;
        }
        else if constexpr (spec == Na22) {
            return 3.652262279854593e-23_rt;
        }
        else if constexpr (spec == Na23) {
            return 3.817541002484691e-23_rt;
        }
        else if constexpr (spec == Mg24) {
            return 3.9828098739467446e-23_rt;
        }
        else if constexpr (spec == Al27) {
            return 4.480389861070653e-23_rt;
        }
        else if constexpr (spec == Si28) {
            return 4.6456779473820677e-23_rt;
        }
        else if constexpr (spec == P31) {
            return 5.14331418367544e-23_rt;
        }
        else if constexpr (spec == S32) {
            return 5.309087322384128e-23_rt;
        }
        else if constexpr (spec == Ar36) {
            return 5.972551377884467e-23_rt;
        }
        else if constexpr (spec == Ca40) {
            return 6.635944331004904e-23_rt;
        }
        else if constexpr (spec == Ti44) {
            return 7.299678247096977e-23_rt;
        }
        else if constexpr (spec == Cr48) {
            return 7.962953983065421e-23_rt;
        }
        else if constexpr (spec == Fe52) {
            return 8.626187166893794e-23_rt;
        }
        else if constexpr (spec == Ni56) {
            return 9.289408870379396e-23_rt;
        }


        // Return zero if we don't recognize the species.
        return 0.0_rt;
    }

    // Legacy (non-templated) interfaces

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real bion (int spec) {
        using namespace Species;

        amrex::Real b = 0.0_rt;

        // Set the binding energy of the element
        amrex::constexpr_for<1, NumSpec+1>([&] (auto n) {
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

        amrex::constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                m = mion<n>();
            }
        });

        return m;
    }
}

namespace Rates
{

    enum NetworkRates :
    std::uint8_t
    {
        k_N13_to_C13_reaclib = 1,
        k_O14_to_N14_reaclib = 2,
        k_O15_to_N15_reaclib = 3,
        k_N13_to_p_C12_reaclib = 4,
        k_N14_to_p_C13_reaclib = 5,
        k_O14_to_p_N13_reaclib = 6,
        k_O15_to_p_N14_reaclib = 7,
        k_O16_to_p_N15_reaclib = 8,
        k_O16_to_He4_C12_reaclib = 9,
        k_F17_to_p_O16_reaclib = 10,
        k_F18_to_p_O17_reaclib = 11,
        k_F18_to_He4_N14_reaclib = 12,
        k_F19_to_p_O18_reaclib = 13,
        k_F19_to_He4_N15_reaclib = 14,
        k_Ne18_to_p_F17_reaclib = 15,
        k_Ne18_to_He4_O14_reaclib = 16,
        k_Ne19_to_p_F18_reaclib = 17,
        k_Ne19_to_He4_O15_reaclib = 18,
        k_Ne20_to_p_F19_reaclib = 19,
        k_Ne20_to_He4_O16_reaclib = 20,
        k_Ne21_to_He4_O17_reaclib = 21,
        k_Ne22_to_He4_O18_reaclib = 22,
        k_Na21_to_p_Ne20_reaclib = 23,
        k_Na21_to_He4_F17_reaclib = 24,
        k_Na22_to_p_Ne21_reaclib = 25,
        k_Na22_to_He4_F18_reaclib = 26,
        k_Na23_to_p_Ne22_reaclib = 27,
        k_Na23_to_He4_F19_reaclib = 28,
        k_Mg24_to_p_Na23_reaclib = 29,
        k_Mg24_to_He4_Ne20_reaclib = 30,
        k_Al27_to_He4_Na23_reaclib = 31,
        k_Si28_to_p_Al27_reaclib = 32,
        k_Si28_to_He4_Mg24_reaclib = 33,
        k_P31_to_He4_Al27_reaclib = 34,
        k_S32_to_p_P31_reaclib = 35,
        k_S32_to_He4_Si28_reaclib = 36,
        k_C12_to_He4_He4_He4_reaclib = 37,
        k_p_C12_to_N13_reaclib = 38,
        k_He4_C12_to_O16_reaclib = 39,
        k_p_C13_to_N14_reaclib = 40,
        k_p_N13_to_O14_reaclib = 41,
        k_p_N14_to_O15_reaclib = 42,
        k_He4_N14_to_F18_reaclib = 43,
        k_p_N15_to_O16_reaclib = 44,
        k_He4_N15_to_F19_reaclib = 45,
        k_He4_O14_to_Ne18_reaclib = 46,
        k_He4_O15_to_Ne19_reaclib = 47,
        k_p_O16_to_F17_reaclib = 48,
        k_He4_O16_to_Ne20_reaclib = 49,
        k_p_O17_to_F18_reaclib = 50,
        k_He4_O17_to_Ne21_reaclib = 51,
        k_p_O18_to_F19_reaclib = 52,
        k_He4_O18_to_Ne22_reaclib = 53,
        k_p_F17_to_Ne18_reaclib = 54,
        k_He4_F17_to_Na21_reaclib = 55,
        k_p_F18_to_Ne19_reaclib = 56,
        k_He4_F18_to_Na22_reaclib = 57,
        k_p_F19_to_Ne20_reaclib = 58,
        k_He4_F19_to_Na23_reaclib = 59,
        k_p_Ne20_to_Na21_reaclib = 60,
        k_He4_Ne20_to_Mg24_reaclib = 61,
        k_p_Ne21_to_Na22_reaclib = 62,
        k_p_Ne22_to_Na23_reaclib = 63,
        k_p_Na23_to_Mg24_reaclib = 64,
        k_He4_Na23_to_Al27_reaclib = 65,
        k_He4_Mg24_to_Si28_reaclib = 66,
        k_p_Al27_to_Si28_reaclib = 67,
        k_He4_Al27_to_P31_reaclib = 68,
        k_He4_Si28_to_S32_reaclib = 69,
        k_p_P31_to_S32_reaclib = 70,
        k_He4_C12_to_p_N15_reaclib = 71,
        k_C12_C12_to_p_Na23_reaclib = 72,
        k_C12_C12_to_He4_Ne20_reaclib = 73,
        k_He4_N13_to_p_O16_reaclib = 74,
        k_He4_N14_to_p_O17_reaclib = 75,
        k_p_N15_to_He4_C12_reaclib = 76,
        k_He4_N15_to_p_O18_reaclib = 77,
        k_He4_O14_to_p_F17_reaclib = 78,
        k_He4_O15_to_p_F18_reaclib = 79,
        k_p_O16_to_He4_N13_reaclib = 80,
        k_He4_O16_to_p_F19_reaclib = 81,
        k_C12_O16_to_p_Al27_reaclib = 82,
        k_C12_O16_to_He4_Mg24_reaclib = 83,
        k_O16_O16_to_p_P31_reaclib = 84,
        k_O16_O16_to_He4_Si28_reaclib = 85,
        k_p_O17_to_He4_N14_reaclib = 86,
        k_p_O18_to_He4_N15_reaclib = 87,
        k_p_F17_to_He4_O14_reaclib = 88,
        k_He4_F17_to_p_Ne20_reaclib = 89,
        k_p_F18_to_He4_O15_reaclib = 90,
        k_He4_F18_to_p_Ne21_reaclib = 91,
        k_p_F19_to_He4_O16_reaclib = 92,
        k_He4_F19_to_p_Ne22_reaclib = 93,
        k_He4_Ne18_to_p_Na21_reaclib = 94,
        k_He4_Ne19_to_p_Na22_reaclib = 95,
        k_p_Ne20_to_He4_F17_reaclib = 96,
        k_He4_Ne20_to_p_Na23_reaclib = 97,
        k_p_Ne21_to_He4_F18_reaclib = 98,
        k_p_Ne22_to_He4_F19_reaclib = 99,
        k_p_Na21_to_He4_Ne18_reaclib = 100,
        k_He4_Na21_to_p_Mg24_reaclib = 101,
        k_p_Na22_to_He4_Ne19_reaclib = 102,
        k_p_Na23_to_He4_Ne20_reaclib = 103,
        k_p_Mg24_to_He4_Na21_reaclib = 104,
        k_He4_Mg24_to_p_Al27_reaclib = 105,
        k_p_Al27_to_He4_Mg24_reaclib = 106,
        k_He4_Si28_to_p_P31_reaclib = 107,
        k_p_P31_to_He4_Si28_reaclib = 108,
        k_He4_He4_He4_to_C12_reaclib = 109,
        k_C12_C12_to_n_Mg23_reaclib = 110,
        k_O16_O16_to_n_S31_reaclib = 111,
        k_C12_O16_to_n_Si27_reaclib = 112,
        k_He4_Na22_to_Al26_reaclib = 113,
        k_p_Mg24_to_Al25_reaclib = 114,
        k_He4_S32_to_Ar36_reaclib = 115,
        k_He4_S32_to_p_Cl35_reaclib = 116,
        k_p_Cl35_to_Ar36_reaclib = 117,
        k_Ar36_to_He4_S32_reaclib = 118,
        k_Ar36_to_p_Cl35_reaclib = 119,
        k_p_Cl35_to_He4_S32_reaclib = 120,
        k_He4_Ar36_to_Ca40_reaclib = 121,
        k_He4_Ar36_to_p_K39_reaclib = 122,
        k_p_K39_to_Ca40_reaclib = 123,
        k_Ca40_to_He4_Ar36_reaclib = 124,
        k_Ca40_to_p_K39_reaclib = 125,
        k_p_K39_to_He4_Ar36_reaclib = 126,
        k_He4_Ca40_to_Ti44_reaclib = 127,
        k_He4_Ca40_to_p_Sc43_reaclib = 128,
        k_p_Sc43_to_Ti44_reaclib = 129,
        k_Ti44_to_He4_Ca40_reaclib = 130,
        k_Ti44_to_p_Sc43_reaclib = 131,
        k_p_Sc43_to_He4_Ca40_reaclib = 132,
        k_He4_Ti44_to_Cr48_reaclib = 133,
        k_He4_Ti44_to_p_V47_reaclib = 134,
        k_p_V47_to_Cr48_reaclib = 135,
        k_Cr48_to_He4_Ti44_reaclib = 136,
        k_Cr48_to_p_V47_reaclib = 137,
        k_p_V47_to_He4_Ti44_reaclib = 138,
        k_He4_Cr48_to_Fe52_reaclib = 139,
        k_He4_Cr48_to_p_Mn51_reaclib = 140,
        k_p_Mn51_to_Fe52_reaclib = 141,
        k_Fe52_to_He4_Cr48_reaclib = 142,
        k_Fe52_to_p_Mn51_reaclib = 143,
        k_p_Mn51_to_He4_Cr48_reaclib = 144,
        k_He4_Fe52_to_Ni56_reaclib = 145,
        k_He4_Fe52_to_p_Co55_reaclib = 146,
        k_p_Co55_to_Ni56_reaclib = 147,
        k_Ni56_to_He4_Fe52_reaclib = 148,
        k_Ni56_to_p_Co55_reaclib = 149,
        k_p_Co55_to_He4_Fe52_reaclib = 150,
        k_F17_to_O17_weaktab = 151,
        k_O17_to_F17_weaktab = 152,
        k_F18_to_Ne18_weaktab = 153,
        k_F18_to_O18_weaktab = 154,
        k_Ne18_to_F18_weaktab = 155,
        k_O18_to_F18_weaktab = 156,
        k_F19_to_Ne19_weaktab = 157,
        k_Ne19_to_F19_weaktab = 158,
        k_Na21_to_Ne21_weaktab = 159,
        k_Ne21_to_Na21_weaktab = 160,
        k_Na22_to_Ne22_weaktab = 161,
        k_Ne22_to_Na22_weaktab = 162,
        k_S32_He4_to_Ar36_approx = 163,
        k_Ar36_to_S32_He4_approx = 164,
        k_Ar36_He4_to_Ca40_approx = 165,
        k_Ca40_to_Ar36_He4_approx = 166,
        k_Ca40_He4_to_Ti44_approx = 167,
        k_Ti44_to_Ca40_He4_approx = 168,
        k_Ti44_He4_to_Cr48_approx = 169,
        k_Cr48_to_Ti44_He4_approx = 170,
        k_Cr48_He4_to_Fe52_approx = 171,
        k_Fe52_to_Cr48_He4_approx = 172,
        k_Fe52_He4_to_Ni56_approx = 173,
        k_Ni56_to_Fe52_He4_approx = 174,
        k_C12_C12_to_Mg24_modified = 175,
        k_O16_O16_to_S32_modified = 176,
        k_C12_O16_to_Si28_modified = 177,
        k_He4_Na22_to_Si28_modified = 178,
        k_p_Mg24_to_Na22_He4_modified = 179,
        NumRates = k_p_Mg24_to_Na22_He4_modified
    };

    enum ScreenPairs :
    std::uint8_t
    {
        k_He4_Ne18 = 1,
        k_p_F18 = 2,
        k_p_C12 = 3,
        k_p_Na21 = 4,
        k_He4_O16 = 5,
        k_He4_F17 = 6,
        k_He4_N14 = 7,
        k_p_K39 = 8,
        k_He4_Ne19 = 9,
        k_p_F19 = 10,
        k_p_C13 = 11,
        k_He4_Ca40 = 12,
        k_p_Mg24 = 13,
        k_He4_Ar36 = 14,
        k_p_Sc43 = 15,
        k_p_Na22 = 16,
        k_He4_F18 = 17,
        k_He4_C12 = 18,
        k_p_O17 = 19,
        k_He4_Na21 = 20,
        k_He4_Ti44 = 21,
        k_p_N15 = 22,
        k_p_Ne20 = 23,
        k_He4_F19 = 24,
        k_He4_Mg24 = 25,
        k_C12_O16 = 26,
        k_p_O18 = 27,
        k_He4_Na22 = 28,
        k_p_V47 = 29,
        k_p_Ne21 = 30,
        k_He4_Cr48 = 31,
        k_He4_O17 = 32,
        k_C12_C12 = 33,
        k_p_Co55 = 34,
        k_He4_N15 = 35,
        k_He4_Ne20 = 36,
        k_p_Al27 = 37,
        k_p_Mn51 = 38,
        k_p_Ne22 = 39,
        k_p_Na23 = 40,
        k_He4_O18 = 41,
        k_He4_Fe52 = 42,
        k_He4_S32 = 43,
        k_p_N13 = 44,
        k_He4_O14 = 45,
        k_p_P31 = 46,
        k_He4_Na23 = 47,
        k_He4_Al27 = 48,
        k_He4_He4 = 49,
        k_p_O16 = 50,
        k_He4_Be8 = 51,
        k_p_N14 = 52,
        k_He4_Si28 = 53,
        k_p_Cl35 = 54,
        k_p_F17 = 55,
        k_O16_O16 = 56,
        k_He4_O15 = 57,
        k_He4_N13 = 58,
        NumScreenPairs = k_He4_N13
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "N13_to_C13_reaclib",  // 1,
        "O14_to_N14_reaclib",  // 2,
        "O15_to_N15_reaclib",  // 3,
        "N13_to_p_C12_reaclib",  // 4,
        "N14_to_p_C13_reaclib",  // 5,
        "O14_to_p_N13_reaclib",  // 6,
        "O15_to_p_N14_reaclib",  // 7,
        "O16_to_p_N15_reaclib",  // 8,
        "O16_to_He4_C12_reaclib",  // 9,
        "F17_to_p_O16_reaclib",  // 10,
        "F18_to_p_O17_reaclib",  // 11,
        "F18_to_He4_N14_reaclib",  // 12,
        "F19_to_p_O18_reaclib",  // 13,
        "F19_to_He4_N15_reaclib",  // 14,
        "Ne18_to_p_F17_reaclib",  // 15,
        "Ne18_to_He4_O14_reaclib",  // 16,
        "Ne19_to_p_F18_reaclib",  // 17,
        "Ne19_to_He4_O15_reaclib",  // 18,
        "Ne20_to_p_F19_reaclib",  // 19,
        "Ne20_to_He4_O16_reaclib",  // 20,
        "Ne21_to_He4_O17_reaclib",  // 21,
        "Ne22_to_He4_O18_reaclib",  // 22,
        "Na21_to_p_Ne20_reaclib",  // 23,
        "Na21_to_He4_F17_reaclib",  // 24,
        "Na22_to_p_Ne21_reaclib",  // 25,
        "Na22_to_He4_F18_reaclib",  // 26,
        "Na23_to_p_Ne22_reaclib",  // 27,
        "Na23_to_He4_F19_reaclib",  // 28,
        "Mg24_to_p_Na23_reaclib",  // 29,
        "Mg24_to_He4_Ne20_reaclib",  // 30,
        "Al27_to_He4_Na23_reaclib",  // 31,
        "Si28_to_p_Al27_reaclib",  // 32,
        "Si28_to_He4_Mg24_reaclib",  // 33,
        "P31_to_He4_Al27_reaclib",  // 34,
        "S32_to_p_P31_reaclib",  // 35,
        "S32_to_He4_Si28_reaclib",  // 36,
        "C12_to_He4_He4_He4_reaclib",  // 37,
        "p_C12_to_N13_reaclib",  // 38,
        "He4_C12_to_O16_reaclib",  // 39,
        "p_C13_to_N14_reaclib",  // 40,
        "p_N13_to_O14_reaclib",  // 41,
        "p_N14_to_O15_reaclib",  // 42,
        "He4_N14_to_F18_reaclib",  // 43,
        "p_N15_to_O16_reaclib",  // 44,
        "He4_N15_to_F19_reaclib",  // 45,
        "He4_O14_to_Ne18_reaclib",  // 46,
        "He4_O15_to_Ne19_reaclib",  // 47,
        "p_O16_to_F17_reaclib",  // 48,
        "He4_O16_to_Ne20_reaclib",  // 49,
        "p_O17_to_F18_reaclib",  // 50,
        "He4_O17_to_Ne21_reaclib",  // 51,
        "p_O18_to_F19_reaclib",  // 52,
        "He4_O18_to_Ne22_reaclib",  // 53,
        "p_F17_to_Ne18_reaclib",  // 54,
        "He4_F17_to_Na21_reaclib",  // 55,
        "p_F18_to_Ne19_reaclib",  // 56,
        "He4_F18_to_Na22_reaclib",  // 57,
        "p_F19_to_Ne20_reaclib",  // 58,
        "He4_F19_to_Na23_reaclib",  // 59,
        "p_Ne20_to_Na21_reaclib",  // 60,
        "He4_Ne20_to_Mg24_reaclib",  // 61,
        "p_Ne21_to_Na22_reaclib",  // 62,
        "p_Ne22_to_Na23_reaclib",  // 63,
        "p_Na23_to_Mg24_reaclib",  // 64,
        "He4_Na23_to_Al27_reaclib",  // 65,
        "He4_Mg24_to_Si28_reaclib",  // 66,
        "p_Al27_to_Si28_reaclib",  // 67,
        "He4_Al27_to_P31_reaclib",  // 68,
        "He4_Si28_to_S32_reaclib",  // 69,
        "p_P31_to_S32_reaclib",  // 70,
        "He4_C12_to_p_N15_reaclib",  // 71,
        "C12_C12_to_p_Na23_reaclib",  // 72,
        "C12_C12_to_He4_Ne20_reaclib",  // 73,
        "He4_N13_to_p_O16_reaclib",  // 74,
        "He4_N14_to_p_O17_reaclib",  // 75,
        "p_N15_to_He4_C12_reaclib",  // 76,
        "He4_N15_to_p_O18_reaclib",  // 77,
        "He4_O14_to_p_F17_reaclib",  // 78,
        "He4_O15_to_p_F18_reaclib",  // 79,
        "p_O16_to_He4_N13_reaclib",  // 80,
        "He4_O16_to_p_F19_reaclib",  // 81,
        "C12_O16_to_p_Al27_reaclib",  // 82,
        "C12_O16_to_He4_Mg24_reaclib",  // 83,
        "O16_O16_to_p_P31_reaclib",  // 84,
        "O16_O16_to_He4_Si28_reaclib",  // 85,
        "p_O17_to_He4_N14_reaclib",  // 86,
        "p_O18_to_He4_N15_reaclib",  // 87,
        "p_F17_to_He4_O14_reaclib",  // 88,
        "He4_F17_to_p_Ne20_reaclib",  // 89,
        "p_F18_to_He4_O15_reaclib",  // 90,
        "He4_F18_to_p_Ne21_reaclib",  // 91,
        "p_F19_to_He4_O16_reaclib",  // 92,
        "He4_F19_to_p_Ne22_reaclib",  // 93,
        "He4_Ne18_to_p_Na21_reaclib",  // 94,
        "He4_Ne19_to_p_Na22_reaclib",  // 95,
        "p_Ne20_to_He4_F17_reaclib",  // 96,
        "He4_Ne20_to_p_Na23_reaclib",  // 97,
        "p_Ne21_to_He4_F18_reaclib",  // 98,
        "p_Ne22_to_He4_F19_reaclib",  // 99,
        "p_Na21_to_He4_Ne18_reaclib",  // 100,
        "He4_Na21_to_p_Mg24_reaclib",  // 101,
        "p_Na22_to_He4_Ne19_reaclib",  // 102,
        "p_Na23_to_He4_Ne20_reaclib",  // 103,
        "p_Mg24_to_He4_Na21_reaclib",  // 104,
        "He4_Mg24_to_p_Al27_reaclib",  // 105,
        "p_Al27_to_He4_Mg24_reaclib",  // 106,
        "He4_Si28_to_p_P31_reaclib",  // 107,
        "p_P31_to_He4_Si28_reaclib",  // 108,
        "He4_He4_He4_to_C12_reaclib",  // 109,
        "C12_C12_to_n_Mg23_reaclib",  // 110,
        "O16_O16_to_n_S31_reaclib",  // 111,
        "C12_O16_to_n_Si27_reaclib",  // 112,
        "He4_Na22_to_Al26_reaclib",  // 113,
        "p_Mg24_to_Al25_reaclib",  // 114,
        "He4_S32_to_Ar36_reaclib",  // 115,
        "He4_S32_to_p_Cl35_reaclib",  // 116,
        "p_Cl35_to_Ar36_reaclib",  // 117,
        "Ar36_to_He4_S32_reaclib",  // 118,
        "Ar36_to_p_Cl35_reaclib",  // 119,
        "p_Cl35_to_He4_S32_reaclib",  // 120,
        "He4_Ar36_to_Ca40_reaclib",  // 121,
        "He4_Ar36_to_p_K39_reaclib",  // 122,
        "p_K39_to_Ca40_reaclib",  // 123,
        "Ca40_to_He4_Ar36_reaclib",  // 124,
        "Ca40_to_p_K39_reaclib",  // 125,
        "p_K39_to_He4_Ar36_reaclib",  // 126,
        "He4_Ca40_to_Ti44_reaclib",  // 127,
        "He4_Ca40_to_p_Sc43_reaclib",  // 128,
        "p_Sc43_to_Ti44_reaclib",  // 129,
        "Ti44_to_He4_Ca40_reaclib",  // 130,
        "Ti44_to_p_Sc43_reaclib",  // 131,
        "p_Sc43_to_He4_Ca40_reaclib",  // 132,
        "He4_Ti44_to_Cr48_reaclib",  // 133,
        "He4_Ti44_to_p_V47_reaclib",  // 134,
        "p_V47_to_Cr48_reaclib",  // 135,
        "Cr48_to_He4_Ti44_reaclib",  // 136,
        "Cr48_to_p_V47_reaclib",  // 137,
        "p_V47_to_He4_Ti44_reaclib",  // 138,
        "He4_Cr48_to_Fe52_reaclib",  // 139,
        "He4_Cr48_to_p_Mn51_reaclib",  // 140,
        "p_Mn51_to_Fe52_reaclib",  // 141,
        "Fe52_to_He4_Cr48_reaclib",  // 142,
        "Fe52_to_p_Mn51_reaclib",  // 143,
        "p_Mn51_to_He4_Cr48_reaclib",  // 144,
        "He4_Fe52_to_Ni56_reaclib",  // 145,
        "He4_Fe52_to_p_Co55_reaclib",  // 146,
        "p_Co55_to_Ni56_reaclib",  // 147,
        "Ni56_to_He4_Fe52_reaclib",  // 148,
        "Ni56_to_p_Co55_reaclib",  // 149,
        "p_Co55_to_He4_Fe52_reaclib",  // 150,
        "F17_to_O17_weaktab",  // 151,
        "O17_to_F17_weaktab",  // 152,
        "F18_to_Ne18_weaktab",  // 153,
        "F18_to_O18_weaktab",  // 154,
        "Ne18_to_F18_weaktab",  // 155,
        "O18_to_F18_weaktab",  // 156,
        "F19_to_Ne19_weaktab",  // 157,
        "Ne19_to_F19_weaktab",  // 158,
        "Na21_to_Ne21_weaktab",  // 159,
        "Ne21_to_Na21_weaktab",  // 160,
        "Na22_to_Ne22_weaktab",  // 161,
        "Ne22_to_Na22_weaktab",  // 162,
        "S32_He4_to_Ar36_approx",  // 163,
        "Ar36_to_S32_He4_approx",  // 164,
        "Ar36_He4_to_Ca40_approx",  // 165,
        "Ca40_to_Ar36_He4_approx",  // 166,
        "Ca40_He4_to_Ti44_approx",  // 167,
        "Ti44_to_Ca40_He4_approx",  // 168,
        "Ti44_He4_to_Cr48_approx",  // 169,
        "Cr48_to_Ti44_He4_approx",  // 170,
        "Cr48_He4_to_Fe52_approx",  // 171,
        "Fe52_to_Cr48_He4_approx",  // 172,
        "Fe52_He4_to_Ni56_approx",  // 173,
        "Ni56_to_Fe52_He4_approx",  // 174,
        "C12_C12_to_Mg24_modified",  // 175,
        "O16_O16_to_S32_modified",  // 176,
        "C12_O16_to_Si28_modified",  // 177,
        "He4_Na22_to_Si28_modified",  // 178,
        "p_Mg24_to_Na22_He4_modified"  // 179,
    };

    template<std::uint8_t rate>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real get_p_random() {
        return 0.0_rt;
    }

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    // A RatePair data table, 2D array of shape (NumNSERatePairs, 10).
    // Each row corresponds to a RatePair returned by _get_nse_rate_pairs().
    // Use enum names following NetworkSpecies and NetworkRates.

    //  - Columns 1–3: 1-based index of reactant nuclei (forward rate)
    //  - Columns 4–6: 1-based index of product nuclei (forward rate)
    //  - Columns 7-8: 1-based index of non (n, H1, He4) nuclei.
    //  - Column  9:   1-based index of the forward rate
    //  - Column  10:  1-based index of the reverse rate

    using namespace Species;
    using namespace Rates;

    constexpr int H1_index = 1;
    constexpr int N_index = -1;
    constexpr int He4_index = 2;
    constexpr int NumNSERatePairs = 56;

    inline AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        H1, C12, -1, N13, -1, -1, C12, N13, k_p_C12_to_N13_reaclib, k_N13_to_p_C12_reaclib,
        He4, C12, -1, O16, -1, -1, C12, O16, k_He4_C12_to_O16_reaclib, k_O16_to_He4_C12_reaclib,
        H1, C13, -1, N14, -1, -1, C13, N14, k_p_C13_to_N14_reaclib, k_N14_to_p_C13_reaclib,
        H1, N13, -1, O14, -1, -1, N13, O14, k_p_N13_to_O14_reaclib, k_O14_to_p_N13_reaclib,
        H1, N14, -1, O15, -1, -1, N14, O15, k_p_N14_to_O15_reaclib, k_O15_to_p_N14_reaclib,
        He4, N14, -1, F18, -1, -1, N14, F18, k_He4_N14_to_F18_reaclib, k_F18_to_He4_N14_reaclib,
        H1, N15, -1, O16, -1, -1, N15, O16, k_p_N15_to_O16_reaclib, k_O16_to_p_N15_reaclib,
        He4, N15, -1, F19, -1, -1, N15, F19, k_He4_N15_to_F19_reaclib, k_F19_to_He4_N15_reaclib,
        He4, O14, -1, Ne18, -1, -1, O14, Ne18, k_He4_O14_to_Ne18_reaclib, k_Ne18_to_He4_O14_reaclib,
        He4, O15, -1, Ne19, -1, -1, O15, Ne19, k_He4_O15_to_Ne19_reaclib, k_Ne19_to_He4_O15_reaclib,
        H1, O16, -1, F17, -1, -1, O16, F17, k_p_O16_to_F17_reaclib, k_F17_to_p_O16_reaclib,
        He4, O16, -1, Ne20, -1, -1, O16, Ne20, k_He4_O16_to_Ne20_reaclib, k_Ne20_to_He4_O16_reaclib,
        H1, O17, -1, F18, -1, -1, O17, F18, k_p_O17_to_F18_reaclib, k_F18_to_p_O17_reaclib,
        He4, O17, -1, Ne21, -1, -1, O17, Ne21, k_He4_O17_to_Ne21_reaclib, k_Ne21_to_He4_O17_reaclib,
        H1, O18, -1, F19, -1, -1, O18, F19, k_p_O18_to_F19_reaclib, k_F19_to_p_O18_reaclib,
        He4, O18, -1, Ne22, -1, -1, O18, Ne22, k_He4_O18_to_Ne22_reaclib, k_Ne22_to_He4_O18_reaclib,
        H1, F17, -1, Ne18, -1, -1, F17, Ne18, k_p_F17_to_Ne18_reaclib, k_Ne18_to_p_F17_reaclib,
        He4, F17, -1, Na21, -1, -1, F17, Na21, k_He4_F17_to_Na21_reaclib, k_Na21_to_He4_F17_reaclib,
        H1, F18, -1, Ne19, -1, -1, F18, Ne19, k_p_F18_to_Ne19_reaclib, k_Ne19_to_p_F18_reaclib,
        He4, F18, -1, Na22, -1, -1, F18, Na22, k_He4_F18_to_Na22_reaclib, k_Na22_to_He4_F18_reaclib,
        H1, F19, -1, Ne20, -1, -1, F19, Ne20, k_p_F19_to_Ne20_reaclib, k_Ne20_to_p_F19_reaclib,
        He4, F19, -1, Na23, -1, -1, F19, Na23, k_He4_F19_to_Na23_reaclib, k_Na23_to_He4_F19_reaclib,
        H1, Ne20, -1, Na21, -1, -1, Ne20, Na21, k_p_Ne20_to_Na21_reaclib, k_Na21_to_p_Ne20_reaclib,
        He4, Ne20, -1, Mg24, -1, -1, Ne20, Mg24, k_He4_Ne20_to_Mg24_reaclib, k_Mg24_to_He4_Ne20_reaclib,
        H1, Ne21, -1, Na22, -1, -1, Ne21, Na22, k_p_Ne21_to_Na22_reaclib, k_Na22_to_p_Ne21_reaclib,
        H1, Ne22, -1, Na23, -1, -1, Ne22, Na23, k_p_Ne22_to_Na23_reaclib, k_Na23_to_p_Ne22_reaclib,
        H1, Na23, -1, Mg24, -1, -1, Na23, Mg24, k_p_Na23_to_Mg24_reaclib, k_Mg24_to_p_Na23_reaclib,
        He4, Na23, -1, Al27, -1, -1, Na23, Al27, k_He4_Na23_to_Al27_reaclib, k_Al27_to_He4_Na23_reaclib,
        He4, Mg24, -1, Si28, -1, -1, Mg24, Si28, k_He4_Mg24_to_Si28_reaclib, k_Si28_to_He4_Mg24_reaclib,
        H1, Al27, -1, Si28, -1, -1, Al27, Si28, k_p_Al27_to_Si28_reaclib, k_Si28_to_p_Al27_reaclib,
        He4, Al27, -1, P31, -1, -1, Al27, P31, k_He4_Al27_to_P31_reaclib, k_P31_to_He4_Al27_reaclib,
        He4, Si28, -1, S32, -1, -1, Si28, S32, k_He4_Si28_to_S32_reaclib, k_S32_to_He4_Si28_reaclib,
        H1, P31, -1, S32, -1, -1, P31, S32, k_p_P31_to_S32_reaclib, k_S32_to_p_P31_reaclib,
        He4, N13, -1, H1, O16, -1, N13, O16, k_He4_N13_to_p_O16_reaclib, k_p_O16_to_He4_N13_reaclib,
        H1, N15, -1, He4, C12, -1, N15, C12, k_p_N15_to_He4_C12_reaclib, k_He4_C12_to_p_N15_reaclib,
        He4, O14, -1, H1, F17, -1, O14, F17, k_He4_O14_to_p_F17_reaclib, k_p_F17_to_He4_O14_reaclib,
        H1, O17, -1, He4, N14, -1, O17, N14, k_p_O17_to_He4_N14_reaclib, k_He4_N14_to_p_O17_reaclib,
        H1, O18, -1, He4, N15, -1, O18, N15, k_p_O18_to_He4_N15_reaclib, k_He4_N15_to_p_O18_reaclib,
        He4, F17, -1, H1, Ne20, -1, F17, Ne20, k_He4_F17_to_p_Ne20_reaclib, k_p_Ne20_to_He4_F17_reaclib,
        H1, F18, -1, He4, O15, -1, F18, O15, k_p_F18_to_He4_O15_reaclib, k_He4_O15_to_p_F18_reaclib,
        He4, F18, -1, H1, Ne21, -1, F18, Ne21, k_He4_F18_to_p_Ne21_reaclib, k_p_Ne21_to_He4_F18_reaclib,
        H1, F19, -1, He4, O16, -1, F19, O16, k_p_F19_to_He4_O16_reaclib, k_He4_O16_to_p_F19_reaclib,
        He4, F19, -1, H1, Ne22, -1, F19, Ne22, k_He4_F19_to_p_Ne22_reaclib, k_p_Ne22_to_He4_F19_reaclib,
        He4, Ne18, -1, H1, Na21, -1, Ne18, Na21, k_He4_Ne18_to_p_Na21_reaclib, k_p_Na21_to_He4_Ne18_reaclib,
        He4, Ne19, -1, H1, Na22, -1, Ne19, Na22, k_He4_Ne19_to_p_Na22_reaclib, k_p_Na22_to_He4_Ne19_reaclib,
        He4, Na21, -1, H1, Mg24, -1, Na21, Mg24, k_He4_Na21_to_p_Mg24_reaclib, k_p_Mg24_to_He4_Na21_reaclib,
        H1, Na23, -1, He4, Ne20, -1, Na23, Ne20, k_p_Na23_to_He4_Ne20_reaclib, k_He4_Ne20_to_p_Na23_reaclib,
        H1, Al27, -1, He4, Mg24, -1, Al27, Mg24, k_p_Al27_to_He4_Mg24_reaclib, k_He4_Mg24_to_p_Al27_reaclib,
        H1, P31, -1, He4, Si28, -1, P31, Si28, k_p_P31_to_He4_Si28_reaclib, k_He4_Si28_to_p_P31_reaclib,
        He4, He4, He4, C12, -1, -1, C12, -1, k_He4_He4_He4_to_C12_reaclib, k_C12_to_He4_He4_He4_reaclib,
        S32, He4, -1, Ar36, -1, -1, S32, Ar36, k_S32_He4_to_Ar36_approx, k_Ar36_to_S32_He4_approx,
        Ar36, He4, -1, Ca40, -1, -1, Ar36, Ca40, k_Ar36_He4_to_Ca40_approx, k_Ca40_to_Ar36_He4_approx,
        Ca40, He4, -1, Ti44, -1, -1, Ca40, Ti44, k_Ca40_He4_to_Ti44_approx, k_Ti44_to_Ca40_He4_approx,
        Ti44, He4, -1, Cr48, -1, -1, Ti44, Cr48, k_Ti44_He4_to_Cr48_approx, k_Cr48_to_Ti44_He4_approx,
        Cr48, He4, -1, Fe52, -1, -1, Cr48, Fe52, k_Cr48_He4_to_Fe52_approx, k_Fe52_to_Cr48_He4_approx,
        Fe52, He4, -1, Ni56, -1, -1, Fe52, Ni56, k_Fe52_He4_to_Ni56_approx, k_Ni56_to_Fe52_He4_approx
    };
}
#endif

#endif
