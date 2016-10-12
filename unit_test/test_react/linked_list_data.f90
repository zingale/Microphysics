module linked_list_data
  
  use bl_types

  implicit none

  type :: cell_data_t
     integer :: ni, nj, nk
     real(kind=dp_t) :: T
   contains
     procedure :: init => init
     procedure :: terminate => terminate
  end type cell_data_t

contains

  subroutine init(cell_data)
    class(cell_data_t) :: cell_data
    cell_data % T = 0.0d0
    cell_data % ni = 0
    cell_data % nj = 0
    cell_data % nk = 0
  end subroutine init
  
  subroutine terminate(cell_data)
    class(cell_data_t) :: cell_data
    ! STUB
  end subroutine terminate
  
end module linked_list_data
