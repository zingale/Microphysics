! Setup a 3D grid of smoothly varying rho, T, and user-defined X.  Then
! call react_state() on the grid and output the results.

program test_react

  use BoxLib
  use bl_constants_module
  use bl_types
  use bl_space
  use f2kcli
  use box_util_module
  use ml_layout_module
  use multifab_module
  use variables
  use probin_module, only: dens_min, dens_max, &
                           temp_min, temp_max, test_set, tmax, run_prefix, &
                           small_temp, small_dens
  use runtime_init_module
  use burn_type_module
  use actual_burner_module
  use microphysics_module
  use eos_type_module, only : eos_get_small_temp, eos_get_small_dens
  use network
  use util_module
  use variables
  use fabio_module
  use build_info_module

  ! For temperature sorting
  use linked_list_module

  implicit none

  ! Conventional fluid state multifabs
  type(multifab) , allocatable :: s(:)

  real(kind=dp_t) :: dx(1, MAX_SPACEDIM)

  logical :: pmask(MAX_SPACEDIM)

  type(ml_layout) :: mla
  type(ml_boxarray) :: mba

  integer :: i, j, n
  integer :: ii, jj, kk
  integer :: nrho, nT, nX

  integer :: dm, nlevs

  integer :: n_rhs_min, n_rhs_max, n_rhs_avg

  type(plot_t) :: pf

  real(kind=dp_t), pointer :: sp(:,:,:,:)

  real(kind=dp_t), allocatable :: state(:,:,:,:)

  ! Begin: Declarations for Temperature Sorting
  real(kind=dp_t), allocatable :: state_flat_sort(:,:)
  integer               :: flatlen
  real(kind=dp_t)       :: tscratch
  type(cell_t)          :: trunk
  type(cell_t), pointer :: aptr, bptr
  integer               :: pf_irho, pf_itemp, pf_ispec, pf_ispec_old, pf_irodot, pf_irho_hnuc
  ! End: Declarations for Temperature Sorting
  
  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: domlo(MAX_SPACEDIM), domhi(MAX_SPACEDIM)

  type (burn_t) :: burn_state_in, burn_state_out

  real (kind=dp_t) :: dens_zone, temp_zone
  real (kind=dp_t) :: dlogrho, dlogT
  real (kind=dp_t), allocatable :: xn_zone(:, :)

  real (kind=dp_t) :: sum_X

  character (len=256) :: out_name
  
  call boxlib_initialize()
  call bl_prof_initialize(on = .true.)

  call runtime_init(.true.)

  ! initialize a grid -- since we are not doing anything hydroy, set the
  ! periodic mask to true
  pmask(:) = .true.
  call read_a_hgproj_grid(mba, test_set)
  call ml_layout_build(mla, mba, pmask)

  nlevs = mla % nlevel
  if (nlevs /= 1) then
     call bl_error("ERROR: only 1 level of refinement currently supported")
  endif

  dm = mla % dim
  if (dm /= 3) then
     call bl_error("ERROR: we require dm = 3")
  endif

  ! we don't care about dx -- we have no physical size
  dx(1,:) = ONE

  ! microphysics
  call microphysics_init(small_temp=small_temp, small_dens=small_dens)

  call eos_get_small_temp(small_temp)
  print *, "small_temp = ", small_temp

  call eos_get_small_dens(small_dens)
  print *, "small_dens = ", small_dens

  ! we'll store everything in a multifab -- inputs and outputs
  call init_variables(pf)

  allocate(s(nlevs))

  do n = 1,nlevs
    call multifab_build(s(n), mla%la(n), pf % n_plot_comps, 0)
  end do

  nrho = extent(mla%mba%pd(1),1)
  nT = extent(mla%mba%pd(1),2)
  nX = extent(mla%mba%pd(1),3)

  allocate(state(0:nrho-1, 0:nT-1, 0:nX-1, pf % n_plot_comps))

  ! Allocate the flattened/sorted state
  flatlen = nrho * nT * nX
  allocate(state_flat_sort(flatlen, pf % n_plot_comps))
  
  dlogrho = (log10(dens_max) - log10(dens_min))/(nrho - 1)
  dlogT   = (log10(temp_max) - log10(temp_min))/(nT - 1)

  ! read from the input file to get all the species data
  domlo = lwb(get_pd(get_layout(s(1))))
  domhi = upb(get_pd(get_layout(s(1))))

  allocate(xn_zone(nspec, 0:nX-1))   ! this assumes that lo(3) = 0

  call get_xn(xn_zone, domlo(3), domhi(3))

  ! normalize -- just in case
  do kk = domlo(3), domhi(3)
     sum_X = sum(xn_zone(:, kk))
     xn_zone(:, kk) = xn_zone(:, kk)/sum_X
  enddo

  n = 1  ! single level assumption

  n_rhs_avg = 0
  n_rhs_max = -100000000
  n_rhs_min = 100000000

  do i = 1, nfabs(s(n))
     sp => dataptr(s(n), i)

     lo = lwb(get_box(s(n), i))
     hi = upb(get_box(s(n), i))

     ! Create incoming state
     do kk = lo(3), hi(3)
        do jj = lo(2), hi(2)
           do ii = lo(1), hi(1)
              
              temp_zone = 10.0_dp_t**(log10(temp_min) + dble(jj)*dlogT)
              dens_zone = 10.0_dp_t**(log10(dens_min) + dble(ii)*dlogrho)
              
              state(ii, jj, kk, pf % irho) = dens_zone
              state(ii, jj, kk, pf % itemp) = temp_zone

              do j = 1, nspec
                 state(ii, jj, kk, pf % ispec_old + j - 1) = max(xn_zone(j, kk), 1.e-10_dp_t)
              enddo

           enddo
        enddo
     enddo

     write(*,*) 'Incoming state created'

     ! Create a linked list containing the sorted incoming state
     do kk = lo(3), hi(3)
        do jj = lo(2), hi(2)
           do ii = lo(1), hi(1)
              ! write(*,*) '(i,j,k): ', ii, jj, kk
              tscratch = state(ii, jj, kk, pf % itemp)
              if ( associated(trunk % first_child) ) then
                 ! There exists a child cell
                 aptr => trunk % first_child
                 do while ( associated(aptr) )
                    if ( tscratch .ge. aptr % data % T ) then
                       ! Add before current cell
                       call aptr % new_before(bptr)
                       if ( .not. associated(bptr % prev) ) then
                          trunk % first_child => bptr
                       end if
                       bptr % data % T = tscratch
                       bptr % data % ni = ii
                       bptr % data % nj = jj
                       bptr % data % nk = kk
                       ! write(*,*) 'Making new cell BEFORE'
                       exit
                    else if ( .not. associated(aptr % next) ) then
                       ! Add after current cell
                       call aptr % new_after(bptr)
                       if ( .not. associated(bptr % next) ) then
                          trunk % last_child => bptr
                       end if
                       bptr % data % T = tscratch
                       bptr % data % ni = ii
                       bptr % data % nj = jj
                       bptr % data % nk = kk
                       ! write(*,*) 'Making new cell AFTER'
                       exit
                    else
                       aptr => aptr % next
                    end if
                 end do
              else
                 ! Create first child cell
                 call trunk % create_only_child(bptr)
                 bptr % data % T = tscratch
                 bptr % data % ni = ii
                 bptr % data % nj = jj
                 bptr % data % nk = kk
                 ! write(*,*) 'Created first child cell'
              end if
           end do
        end do
     end do

     write(*,*) 'Sorted linked list created'
     
     ! Create flattened and sorted state
     ii = 1
     aptr => trunk % first_child
     do while ( associated(aptr) )
        ! write(*,*) 'i=', ii
        state_flat_sort(ii, :) = state(aptr % data % ni, aptr % data % nj, aptr % data % nk, :)
        aptr => aptr % next
        ii = ii + 1
     end do
     if ( ii-1 /= flatlen ) then
        write(*,*) 'ERROR: state_flat_sort not filled with number of zones in state!'
     else
        write(*,*) 'Filled state_flat_sort'
     end if

     ! Carry out calculations on incoming state

     ! Check incoming state
     ! write(*,*) '(rho, T): X in flattened/sorted state:'
     ! do kk = 1, flatlen
     !    write(*,*) 'kk=', kk
     !    write(*,*) '(', state_flat_sort(kk, pf % irho), ', ', state_flat_sort(kk, pf % itemp), '):'
     !    do j = 1, nspec
     !       write(*,*) state_flat_sort(kk, pf % ispec_old + j - 1)
     !    enddo
     ! end do
     
     ! Some indices to copyin 
     pf_irho      = pf % irho
     pf_itemp     = pf % itemp
     pf_ispec     = pf % ispec
     pf_ispec_old = pf % ispec_old
     pf_irodot    = pf % irodot
     pf_irho_hnuc = pf % irho_hnuc

     write(*,*) 'Entering burn loop'

     !$OMP PARALLEL DO PRIVATE(kk,j) &
     !$OMP PRIVATE(burn_state_in, burn_state_out) &
     !$OMP REDUCTION(+:n_rhs_avg) REDUCTION(MAX:n_rhs_max) REDUCTION(MIN:n_rhs_min) &
     !$OMP SCHEDULE(DYNAMIC,1)

     !$acc data copyin(flatlen, tmax, pf_irho, pf_itemp, pf_ispec, pf_ispec_old) &
     !$acc      copyin(pf_irodot, pf_irho_hnuc) &
     !$acc      copy(state_flat_sort) 

     !$acc parallel reduction(+:n_rhs_avg) reduction(max:n_rhs_max) reduction(min:n_rhs_min)
     !$acc loop gang vector &
     !$acc private(burn_state_in, burn_state_out, kk, j)

     do kk = 1, flatlen
        
        burn_state_in % rho = state_flat_sort(kk, pf_irho)
        burn_state_in % T = state_flat_sort(kk, pf_itemp)

        do j = 1, nspec
           burn_state_in % xn(j) = state_flat_sort(kk, pf_ispec_old + j - 1) 
        enddo

        call normalize_abundances_burn(burn_state_in)

        ! the integrator doesn't actually care about the initial internal
        ! energy.
        burn_state_in % e = ZERO

        call actual_burner(burn_state_in, burn_state_out, tmax, ZERO)

        do j = 1, nspec
           state_flat_sort(kk, pf_ispec + j - 1) = burn_state_out % xn(j)
        enddo
                            
        do j=1, nspec
           ! an explicit loop is needed here to keep the GPU happy
           state_flat_sort(kk, pf_irodot + j - 1) = &
                (burn_state_out % xn(j) - burn_state_in % xn(j)) / tmax
        enddo

        state_flat_sort(kk, pf_irho_hnuc) = &
             burn_state_in % rho * (burn_state_out % e - burn_state_in % e) / tmax
              
        n_rhs_avg = n_rhs_avg + burn_state_out % n_rhs
        n_rhs_min = min(n_rhs_min, burn_state_out % n_rhs)
        n_rhs_max = max(n_rhs_max, burn_state_out % n_rhs)

     enddo
     !$acc end parallel
     !$acc end data

     !$OMP END PARALLEL DO

     ! Unflatten/Unsort outgoing state
     write(*,*) 'Finished burning, unflattening/unsorting state'
     ii = 1
     aptr => trunk % first_child
     do while ( associated(aptr) )
        ! write(*,*) 'i=', ii
        state(aptr % data % ni, aptr % data % nj, aptr % data % nk, :) = state_flat_sort(ii, :)
        aptr => aptr % next
        ii = ii + 1
     end do
     if ( ii-1 /= flatlen ) then
        write(*,*) 'ERROR: state not filled with number of zones in flattened/sorted state!'
     else
        write(*,*) 'Unflattened/unsorted state'
     end if

     ! Store outgoing state
     sp(:,:,:,:) = state(:,:,:,:)
  enddo

  ! note: integer division
  n_rhs_avg = n_rhs_avg/(nT*nrho*nX)

  print *, "RHS stats:"
  print *, "  min: ", n_rhs_min
  print *, "  avg: ", n_rhs_avg
  print *, "  max: ", n_rhs_max

  ! output
  out_name = trim(run_prefix) // "test_react." // trim(integrator_dir)

  call fabio_ml_multifab_write_d(s, mla%mba%rr(:,1), trim(out_name), names=pf%names)
  
  call write_job_info(out_name, mla%mba)


  ! if you (or a subroutine) built it, destroy it!
  
  ! Deallocate flattened/sorted state
  deallocate(state_flat_sort)
  ! Deallocate linked list
  call trunk % delete_children
  
  do n = 1,nlevs
    call destroy(s(n))
  end do
  
  call destroy(mla)

  deallocate(s)
  deallocate(xn_zone)

  call runtime_close()

    
  call microphysics_finalize()

  ! end boxlib
  call bl_prof_glean("bl_prof_res")
  call bl_prof_finalize()
  call boxlib_finalize()

end program test_react
