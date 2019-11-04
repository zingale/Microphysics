module eos_composition_module

  use eos_type_module, only : eos_t
  use network, only: nspec, aion, zion
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  type :: eos_comp_t
     real(rt) :: mu_e
     real(rt) :: abar
     real(rt) :: zbar
     real(rt) :: y_e
  end type eos_comp_t

  type :: eos_xderivs_t
    real(rt) :: dedX(nspec)
    real(rt) :: dpdX(nspec)
    real(rt) :: dhdX(nspec)
 end type eos_xderivs_t

contains

  ! Given a set of mass fractions, calculate quantities that depend
  ! on the composition like abar and zbar.

  subroutine composition(xn, state_comp)

    use amrex_constants_module, only: ONE
    use network, only: aion_inv, zion

    implicit none

    real(rt), intent(in) :: xn(nspec)
    type (eos_comp_t), intent(out) :: state_comp

    !$gpu

    ! Calculate abar, the mean nucleon number,
    ! zbar, the mean proton number,
    ! mu, the mean molecular weight,
    ! mu_e, the mean number of nucleons per electron, and
    ! y_e, the electron fraction.

    state_comp % mu_e = ONE / (sum(xn(:) * zion(:) * aion_inv(:)))
    state_comp % y_e = ONE / state_comp % mu_e

    state_comp % abar = ONE / (sum(xn(:) * aion_inv(:)))
    state_comp % zbar = state_comp % abar / state_comp % mu_e

  end subroutine composition


  ! Compute thermodynamic derivatives with respect to xn(:)

  subroutine composition_derivatives(state, state_comp, state_xderivs)

    use amrex_constants_module, only: ZERO
    use network, only: aion, aion_inv, zion

    implicit none

    type (eos_t), intent(in) :: state
    type (eos_comp_t), intent(in) :: state_comp
    type (eos_xderivs_t), intent(out) :: state_xderivs

    !$gpu

#ifdef EXTRA_THERMO
    state_xderivs % dpdX(:) = state % dpdA * (state_comp % abar * aion_inv(:))   &
                                        * (aion(:) - state_comp % abar) &
                            + state % dpdZ * (state_comp % abar * aion_inv(:))   &
                                        * (zion(:) - state_comp % zbar)

    state_xderivs % dEdX(:) = state % dedA * (state_comp % abar * aion_inv(:))   &
                                        * (aion(:) - state_comp % abar) &
                            + state % dedZ * (state_comp % abar * aion_inv(:))   &
                                        * (zion(:) - state_comp % zbar)

    if (state % dPdr .ne. ZERO) then

       state_xderivs % dhdX(:) = state_xderivs % dedX(:) &
            + (state % p / state % rho**2 - state % dedr) &
            *  state_xderivs % dPdX(:) / state % dPdr

    endif
#endif

  end subroutine composition_derivatives

end module eos_composition_module
