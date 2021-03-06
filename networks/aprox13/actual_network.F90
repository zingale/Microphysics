module actual_network

  use amrex_fort_module, only : rt => amrex_real
  use network_properties

  implicit none

  integer, parameter :: ihe4  = 1
  integer, parameter :: ic12  = 2
  integer, parameter :: io16  = 3
  integer, parameter :: ine20 = 4
  integer, parameter :: img24 = 5
  integer, parameter :: isi28 = 6
  integer, parameter :: is32  = 7
  integer, parameter :: iar36 = 8
  integer, parameter :: ica40 = 9
  integer, parameter :: iti44 = 10
  integer, parameter :: icr48 = 11
  integer, parameter :: ife52 = 12
  integer, parameter :: ini56 = 13

  real(rt)        , allocatable :: bion(:), mion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: bion, mion
#endif

  character (len=32), parameter :: network_name = "aprox13"

  ! Some fundamental physical constants

  real(rt)        , parameter :: avo = 6.0221417930e23_rt
  real(rt)        , parameter :: c_light = 2.99792458e10_rt

  real(rt)        , parameter :: ev2erg  = 1.60217648740e-12_rt
  real(rt)        , parameter :: mev2erg = ev2erg*1.0e6_rt
  real(rt)        , parameter :: mev2gr  = mev2erg/c_light**2

  real(rt)        , parameter :: mn = 1.67492721184e-24_rt
  real(rt)        , parameter :: mp = 1.67262163783e-24_rt
  real(rt)        , parameter :: me = 9.1093821545e-28_rt

  ! Conversion factor for the nuclear energy generation rate.

  real(rt)        , parameter :: enuc_conv2 = -avo*c_light*c_light

  !$acc declare create(aion, zion, nion, bion, mion)

  ! Rates data

  integer, parameter :: nrates = 67
  integer, parameter :: num_rate_groups = 2

  integer, parameter :: ir3a   = 1
  integer, parameter :: irg3a  = 2
  integer, parameter :: ircag  = 3
  integer, parameter :: iroga  = 4
  integer, parameter :: ir1212 = 5
  integer, parameter :: ir1216 = 6
  integer, parameter :: ir1616 = 7
  integer, parameter :: iroag  = 8
  integer, parameter :: irnega = 9
  integer, parameter :: irneag = 10
  integer, parameter :: irmgga = 11
  integer, parameter :: irmgag = 12
  integer, parameter :: irsiga = 13
  integer, parameter :: irmgap = 14
  integer, parameter :: iralpa = 15
  integer, parameter :: iralpg = 16
  integer, parameter :: irsigp = 17
  integer, parameter :: irsiag = 18
  integer, parameter :: irsga  = 19
  integer, parameter :: irsiap = 20
  integer, parameter :: irppa  = 21
  integer, parameter :: irppg  = 22
  integer, parameter :: irsgp  = 23
  integer, parameter :: irsag  = 24
  integer, parameter :: irarga = 25
  integer, parameter :: irsap  = 26
  integer, parameter :: irclpa = 27
  integer, parameter :: irclpg = 28
  integer, parameter :: irargp = 29
  integer, parameter :: irarag = 30
  integer, parameter :: ircaga = 31
  integer, parameter :: irarap = 32
  integer, parameter :: irkpa  = 33
  integer, parameter :: irkpg  = 34
  integer, parameter :: ircagp = 35
  integer, parameter :: ircaag = 36
  integer, parameter :: irtiga = 37
  integer, parameter :: ircaap = 38
  integer, parameter :: irscpa = 39
  integer, parameter :: irscpg = 40
  integer, parameter :: irtigp = 41
  integer, parameter :: irtiag = 42
  integer, parameter :: ircrga = 43
  integer, parameter :: irtiap = 44
  integer, parameter :: irvpa  = 45
  integer, parameter :: irvpg  = 46
  integer, parameter :: ircrgp = 47
  integer, parameter :: ircrag = 48
  integer, parameter :: irfega = 49
  integer, parameter :: ircrap = 50
  integer, parameter :: irmnpa = 51
  integer, parameter :: irmnpg = 52
  integer, parameter :: irfegp = 53
  integer, parameter :: irfeag = 54
  integer, parameter :: irniga = 55
  integer, parameter :: irfeap = 56
  integer, parameter :: ircopa = 57
  integer, parameter :: ircopg = 58
  integer, parameter :: irnigp = 59
  integer, parameter :: irr1   = 60
  integer, parameter :: irs1   = 61
  integer, parameter :: irt1   = 62
  integer, parameter :: iru1   = 63
  integer, parameter :: irv1   = 64
  integer, parameter :: irw1   = 65
  integer, parameter :: irx1   = 66
  integer, parameter :: iry1   = 67

  character (len=16), save :: ratenames(nrates)

#ifdef REACT_SPARSE_JACOBIAN
  ! Shape of Jacobian in Compressed Sparse Row format
  integer, parameter   :: NETWORK_SPARSE_JAC_NNZ = 107
  integer, allocatable :: csr_jac_col_index(:), csr_jac_row_count(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: csr_jac_col_index, csr_jac_row_count
#endif
#endif

contains

  subroutine actual_network_init

    implicit none

    call network_properties_init()
    allocate(bion(nspec))
    allocate(mion(nspec))

    ! Set the binding energy of the element (MeV)
    bion(ihe4)  =  28.29603e0_rt
    bion(ic12)  =  92.16294e0_rt
    bion(io16)  = 127.62093e0_rt
    bion(ine20) = 160.64788e0_rt
    bion(img24) = 198.25790e0_rt
    bion(isi28) = 236.53790e0_rt
    bion(is32)  = 271.78250e0_rt
    bion(iar36) = 306.72020e0_rt
    bion(ica40) = 342.05680e0_rt
    bion(iti44) = 375.47720e0_rt
    bion(icr48) = 411.46900e0_rt
    bion(ife52) = 447.70800e0_rt
    bion(ini56) = 484.00300e0_rt

    ! Set the mass
    mion(:) = nion(:) * mn + zion(:) * (mp + me) - bion(:) * mev2gr

    !$acc update device(aion, zion, nion, bion, mion)

    ratenames(ir3a)   = 'r3a  '
    ratenames(irg3a)  = 'rg3a '
    ratenames(ircag)  = 'rcag '
    ratenames(ir1212) = 'r1212'
    ratenames(ir1216) = 'r1216'
    ratenames(ir1616) = 'r1616'
    ratenames(iroga)  = 'roga '
    ratenames(iroag)  = 'roag '
    ratenames(irnega) = 'rnega'
    ratenames(irneag) = 'rneag'
    ratenames(irmgga) = 'rmgga'
    ratenames(irmgag) = 'rmgag'
    ratenames(irsiga) = 'rsiga'
    ratenames(irmgap) = 'rmgap'
    ratenames(iralpa) = 'ralpa'
    ratenames(iralpg) = 'ralpg'
    ratenames(irsigp) = 'rsigp'
    ratenames(irsiag) = 'rsiag'
    ratenames(irsga)  = 'rsga '
    ratenames(irsiap) = 'rsiap'
    ratenames(irppa)  = 'rppa '
    ratenames(irppg)  = 'rppg '
    ratenames(irsgp)  = 'rsgp '
    ratenames(irsag)  = 'rsag '
    ratenames(irarga) = 'rarga'
    ratenames(irsap)  = 'rsap '
    ratenames(irclpa) = 'rclpa'
    ratenames(irclpg) = 'rclpg'
    ratenames(irargp) = 'rargp'
    ratenames(irarag) = 'rarag'
    ratenames(ircaga) = 'rcaga'
    ratenames(irarap) = 'rarap'
    ratenames(irkpa)  = 'rkpa '
    ratenames(irkpg)  = 'rkpg '
    ratenames(ircagp) = 'rcagp'
    ratenames(ircaag) = 'rcaag'
    ratenames(irtiga) = 'rtiga'
    ratenames(ircaap) = 'rcaap'
    ratenames(irscpa) = 'rscpa'
    ratenames(irscpg) = 'rscpg'
    ratenames(irtigp) = 'rtigp'
    ratenames(irtiag) = 'rtiag'
    ratenames(ircrga) = 'rcrga'
    ratenames(irtiap) = 'rtiap'
    ratenames(irvpa)  = 'rvpa '
    ratenames(irvpg)  = 'rvpg '
    ratenames(ircrgp) = 'rcrgp'
    ratenames(ircrag) = 'rcrag'
    ratenames(irfega) = 'rfega'
    ratenames(ircrap) = 'rcrap'
    ratenames(irmnpa) = 'rmnpa'
    ratenames(irmnpg) = 'rmnpg'
    ratenames(irfegp) = 'rfegp'
    ratenames(irfeag) = 'rfeag'
    ratenames(irniga) = 'rniga'
    ratenames(irfeap) = 'rfeap'
    ratenames(ircopa) = 'rcopa'
    ratenames(ircopg) = 'rcopg'
    ratenames(irnigp) = 'rnigp'
    ratenames(irr1)   = 'r1   '
    ratenames(irs1)   = 's1   '
    ratenames(irt1)   = 't1   '
    ratenames(iru1)   = 'u1   '
    ratenames(irv1)   = 'v1   '
    ratenames(irw1)   = 'w1   '
    ratenames(irx1)   = 'x1   '
    ratenames(iry1)   = 'y1   '

#ifdef REACT_SPARSE_JACOBIAN
    ! Set CSR format metadata for Jacobian
    allocate(csr_jac_col_index(NETWORK_SPARSE_JAC_NNZ))
    allocate(csr_jac_row_count(nspec + 3)) ! neq + 1

    csr_jac_col_index = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, &
                         1, 2, 3, 14, &
                         1, 2, 3, 4, 14, &
                         1, 2, 3, 4, 5, 14, &
                         1, 2, 3, 4, 5, 6, 14, &
                         1, 2, 3, 5, 6, 7, 14, &
                         1, 3, 6, 7, 8, 14, &
                         1, 7, 8, 9, 14, &
                         1, 8, 9, 10, 14, &
                         1, 9, 10, 11, 14, &
                         1, 10, 11, 12, 14, &
                         1, 11, 12, 13, 14, &
                         1, 12, 13, 14, &
                         1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, &
                         1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

    csr_jac_row_count = [1, 15, 19, 24, 30, &
                         37, 44, 50, 55, 60, &
                         65, 70, 75, 79, 93, 108]
#endif

  end subroutine actual_network_init



  subroutine actual_network_finalize

    implicit none

    call network_properties_finalize

    if (allocated(bion)) then
       deallocate(bion)
    endif
    if (allocated(mion)) then
       deallocate(mion)
    endif

  end subroutine actual_network_finalize

end module actual_network
