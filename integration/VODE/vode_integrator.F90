! Common variables and routines for burners
! that use VODE for their integration.

module vode_integrator_module

  use eos_type_module, only: eos_input_rt
  use network
  use vode_rpar_indices
  use vode_type_module
  use burn_type_module
  use cuvode_parameters_module
  use amrex_fort_module, only: rt => amrex_real

  implicit none
  
contains

  subroutine vode_integrator_init()

    implicit none

  end subroutine vode_integrator_init


  ! Main interface
  subroutine vode_integrator(state_in, state_out, dt, time, status)

    !$acc routine seq

    use vode_rpar_indices
    use extern_probin_module, only: jacobian, use_jacobian_caching, &
         burner_verbose, &
         burning_mode, burning_mode_factor, &
         call_eos_in_rhs, dt_crit, ode_max_steps
    use cuvode_module, only: dvode
    use eos_module, only: eos
    use eos_type_module, only: eos_t, copy_eos_t
    use cuvode_types_module, only: dvode_t
    use amrex_constants_module, only: ZERO, ONE
    use integrator_scaling_module, only: temp_scale, ener_scale, inv_ener_scale
    use integration_data, only: integration_status_t
    use temperature_integration_module, only: self_heat
#ifndef CUDA
    use amrex_error_module, only: amrex_error
#endif

    implicit none

    ! Input arguments

    type (burn_t), intent(in   ) :: state_in
    type (burn_t), intent(inout) :: state_out
    real(rt),    intent(in   ) :: dt, time
    type (integration_status_t), intent(inout) :: status

    ! Local variables

    type (eos_t) :: eos_state_in, eos_state_temp
    type (dvode_t) :: dvode_state

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.

    real(rt) :: ener_offset
    real(rt) :: t_enuc, t_sound, limit_factor

    logical :: integration_failed
    real(rt), parameter :: failure_tolerance = 1.e-2_rt

    !$gpu

    dvode_state % jacobian = jacobian

    if (use_jacobian_caching) then
       dvode_state % JSV = 1
    else
       dvode_state % JSV = -1
    endif

    integration_failed = .false.

    ! Set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of
    ! steps allowed.

    dvode_state % atol(1:nspec)   = status % atol_spec ! mass fractions
    dvode_state % atol(net_itemp) = status % atol_temp ! temperature
    dvode_state % atol(net_ienuc) = status % atol_enuc ! energy generated

    dvode_state % rtol(1:nspec)   = status % rtol_spec ! mass fractions
    dvode_state % rtol(net_itemp) = status % rtol_temp ! temperature
    dvode_state % rtol(net_ienuc) = status % rtol_enuc ! energy generated

    ! We want VODE to re-initialize each time we call it.

    dvode_state % istate = 1

    ! Set the maximum number of steps allowed.

    dvode_state % MXSTEP = ode_max_steps

    ! Start off by assuming a successful burn.

    state_out % success = .true.

    ! Initialize the integration time.
    dvode_state % T = ZERO
    dvode_state % TOUT = dt

    ! Convert our input burn state into an EOS type.

    call burn_to_eos(state_in, eos_state_in)

    ! We assume that (rho, T) coming in are valid, do an EOS call
    ! to fill the rest of the thermodynamic variables.

    call eos(eos_input_rt, eos_state_in)

    ! Convert the EOS state data into the form VODE expects.

    call eos_to_vode(eos_state_in, dvode_state)

    ener_offset = eos_state_in % e * inv_ener_scale

    dvode_state % y(net_ienuc) = ener_offset

    ! Pass through whether we are doing self-heating.

    if (self_heat) then
       dvode_state % rpar(irp_self_heat) = ONE
    else
       dvode_state % rpar(irp_self_heat) = -ONE
    endif

    ! Copy in the zone size.

    dvode_state % rpar(irp_dx) = state_in % dx

    ! Set the sound crossing time.

    dvode_state % rpar(irp_t_sound) = state_in % dx / eos_state_in % cs

    ! Set the time offset -- this converts between the local integration 
    ! time and the simulation time

    dvode_state % rpar(irp_t0) = time

    ! If we are using the dT_crit functionality and therefore doing a linear
    ! interpolation of the specific heat in between EOS calls, do a second
    ! EOS call here to establish an initial slope.

    dvode_state % rpar(irp_Told) = eos_state_in % T

    if (dT_crit < 1.0e19_rt) then

       call copy_eos_t(eos_state_temp, eos_state_in)
       eos_state_temp % T = eos_state_in % T * (ONE + sqrt(epsilon(ONE)))

       call eos(eos_input_rt, eos_state_temp)

       dvode_state % rpar(irp_dcvdt) = (eos_state_temp % cv - eos_state_in % cv) / &
                                       (eos_state_temp % T - eos_state_in % T)
       dvode_state % rpar(irp_dcpdt) = (eos_state_temp % cp - eos_state_in % cp) / &
                                       (eos_state_temp % T - eos_state_in % T)

    endif

    ! Call the integration routine.
    call dvode(dvode_state)

    ! If we are using hybrid burning and the energy release was negative (or we failed),
    ! re-run this in self-heating mode.

    if ( burning_mode == 2 .and. &
         (dvode_state % y(net_ienuc) - ener_offset < ZERO .or. &
          dvode_state % istate < 0) ) then

       dvode_state % rpar(irp_self_heat) = ONE

       dvode_state % istate = 1

       dvode_state % T = ZERO
       dvode_state % TOUT = dt

       call eos_to_vode(eos_state_in, dvode_state)

       dvode_state % rpar(irp_Told) = eos_state_in % T

       if (dT_crit < 1.0e19_rt) then

          dvode_state % rpar(irp_dcvdt) = (eos_state_temp % cv - eos_state_in % cv) / &
                                          (eos_state_temp % T - eos_state_in % T)
          dvode_state % rpar(irp_dcpdt) = (eos_state_temp % cp - eos_state_in % cp) / &
                                          (eos_state_temp % T - eos_state_in % T)

       endif

       dvode_state % y(net_ienuc) = ener_offset

       ! Call the integration routine.
       call dvode(dvode_state)

    endif

    ! VODE does not always fail even though it can lead to unphysical states,
    ! so add some sanity checks that trigger a retry even if VODE thinks
    ! the integration was successful.

    if (dvode_state % istate < 0) then
       integration_failed = .true.
    end if

    if (dvode_state % y(net_itemp) < ZERO) then
       integration_failed = .true.
    end if

    if (any(dvode_state % y(1:nspec) < -failure_tolerance)) then
       integration_failed = .true.
    end if

    if (any(dvode_state % y(1:nspec) > 1.e0_rt + failure_tolerance)) then
       integration_failed = .true.
    end if

    ! If we failed, print out the current state of the integration.

    if (integration_failed) then
#ifndef CUDA
       print *, 'ERROR: integration failed in net'
       print *, 'istate = ', dvode_state % istate
       print *, 'time = ', dvode_state % T
       print *, 'dens = ', state_in % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'temp current = ', dvode_state % y(net_itemp) * temp_scale
       print *, 'xn current = ', dvode_state % y(1:nspec) * aion(1:nspec)
       print *, 'energy generated = ', (dvode_state % y(net_ienuc) - ener_offset) * ener_scale
#endif

       state_out % success = .false.
       return
    endif

    ! Subtract the energy offset
    dvode_state % y(net_ienuc) = dvode_state % y(net_ienuc) - ener_offset

    ! Store the final data, and then normalize abundances.
    call vode_to_burn(dvode_state, state_out)

    ! get the number of RHS calls and jac evaluations from the VODE
    ! work arrays
    state_out % n_rhs = dvode_state % NFE
    state_out % n_jac = dvode_state % NJE

    state_out % time = dvode_state % t

    ! For burning_mode == 3, limit the burning.

    if (burning_mode == 3) then

       t_enuc = eos_state_in % e / max(abs(state_out % e - state_in % e) / max(dt, 1.e-50_rt), 1.e-50_rt)
       t_sound = state_in % dx / eos_state_in % cs

       limit_factor = min(1.0e0_rt, burning_mode_factor * t_enuc / t_sound)

       state_out % e = state_in % e + limit_factor * (state_out % e - state_in % e)
       state_out % xn(:) = state_in % xn(:) + limit_factor * (state_out % xn(:) - state_in % xn(:))

    endif

    call normalize_abundances_burn(state_out)

    ! set the integration time for any diagnostics
    state_out % time = time + dt

#ifndef CUDA
    if (burner_verbose) then

       ! Print out some integration statistics, if desired.
       print *, 'integration summary: '
       print *, 'dens: ', state_out % rho, ' temp: ', state_out % T, &
                ' energy released: ', state_out % e - state_in % e
       print *, 'number of steps taken: ', dvode_state % NST
       print *, 'number of f evaluations: ', dvode_state % NFE
    endif
#endif
    
  end subroutine vode_integrator

end module vode_integrator_module
