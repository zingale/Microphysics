! the burner_aux module contains thermodynamic state variables that are
! needed in the RHS and Jacobian routines

module burner_aux_module

  use amrex_fort_module, only : rt => amrex_real
  use network, only : nspec

  implicit none

  real(kind=rt), save :: sdc_rhoX_pass(nspec)
  real(kind=rt), save :: sdc_rhoh_pass
  real(kind=rt), save :: p0_pass  

end module burner_aux_module
