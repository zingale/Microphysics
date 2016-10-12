module linked_list_module

  use linked_list_data, only: cell_data_t
  
  implicit none

  type :: cell_t
     type(cell_t), pointer :: next=>NULL(), prev=>NULL()
     type(cell_t), pointer :: first_child=>NULL(), last_child=>NULL()
     type(cell_data_t)     :: data
   contains
     procedure :: new_after  => new_after
     procedure :: new_before => new_before
     procedure :: create_only_child => create_only_child
     procedure :: append_new_child  => append_new_child
     procedure :: prepend_new_child => prepend_new_child
     procedure :: delete_cell       => delete_cell
     procedure :: delete_children   => delete_children
  end type cell_t

contains

  subroutine new_after(cell, cell_new)
    class(cell_t), target :: cell
    type(cell_t), pointer :: cell_new
    type(cell_t), pointer :: cell_ptr_next
    cell_ptr_next => cell % next
    nullify(cell % next)
    allocate(cell % next)
    call cell % next % data % init
    cell_new => cell % next
    cell_new % next => cell_ptr_next
    cell_new % prev => cell
    if ( associated(cell_ptr_next) ) then
       cell_new % next % prev => cell_new
    end if
  end subroutine new_after

  subroutine new_before(cell, cell_new)
    class(cell_t), target :: cell
    type(cell_t), pointer :: cell_new
    type(cell_t), pointer :: cell_ptr_prev
    cell_ptr_prev => cell % prev
    nullify(cell % prev)
    allocate(cell % prev)
    call cell % prev % data % init
    cell_new => cell % prev
    cell_new % prev => cell_ptr_prev
    cell_new % next => cell
    if ( associated(cell_ptr_prev) ) then
       cell_new % prev % next => cell_new
    end if
  end subroutine new_before

  ! NOT TESTED:
  ! subroutine insert_after(cell, cell_insert)
  !   class(cell_t), target :: cell
  !   type(cell_t), target  :: cell_insert
  !   cell_insert % next => cell % next
  !   cell_insert % prev => cell
  !   cell % next => cell_insert
  ! end subroutine insert_after

  ! NOT TESTED:
  ! subroutine insert_before(cell, cell_insert)
  !   class(cell_t), target :: cell
  !   type(cell_t), target  :: cell_insert
  !   cell_insert % next => cell
  !   cell_insert % prev => cell % prev
  !   cell % prev => cell_insert
  ! end subroutine insert_before

  subroutine create_only_child(cell, cell_new)
    ! Creates the first child of an otherwise childless cell
    class(cell_t) :: cell
    type(cell_t), pointer :: cell_new
    allocate(cell % first_child)
    call cell % first_child % data % init
    cell_new => cell % first_child
    cell % last_child => cell_new
  end subroutine create_only_child

  subroutine append_new_child(cell, cell_new)
    class(cell_t) :: cell
    type(cell_t), pointer :: cell_new
    if ( associated(cell % last_child) ) then
       call cell % last_child % new_after(cell_new)
       cell % last_child => cell_new
    else
       call cell % create_only_child(cell_new)
    end if
  end subroutine append_new_child

  subroutine prepend_new_child(cell, cell_new)
    class(cell_t) :: cell
    type(cell_t), pointer :: cell_new
    if ( associated(cell % first_child) ) then
       call cell % first_child % new_before(cell_new)
       cell % first_child => cell_new
    else
       call cell % create_only_child(cell_new)
    end if
  end subroutine prepend_new_child

  subroutine delete_cell(cell, cell_next)
    ! Delete this cell's data and its child cells
    class(cell_t) :: cell
    type(cell_t), pointer :: cell_next
!    write(*,*) 'Entered delete_cell'
    call cell % data % terminate
    if ( associated(cell % prev) ) then
       cell % prev % next => cell % next
    end if
    if ( associated(cell % next) ) then
       cell % next % prev => cell % prev
    end if
    if ( associated(cell % first_child) ) then
       call cell % delete_children
    end if
    cell_next => cell % next
!    write(*,*) 'Leaving delete_cell'
  end subroutine delete_cell

  subroutine delete_children(cell)
    ! Delete and deallocate all child cells
    class(cell_t) :: cell
    type(cell_t), pointer    :: cell_ptr
    type(cell_t), pointer    :: cell_next
!    write(*,*) 'Entered delete_children'
!    write(*,*) ''
    cell_ptr => cell % first_child
    do while ( associated(cell_ptr) )
!       write(*,*) 'Deleting cell cell_ptr'
       call cell_ptr % delete_cell(cell_next)
!       write(*,*) 'Deallocating cell_ptr'
       deallocate(cell_ptr)
!       write(*,*) 'Deallocated cell_ptr'
!       write(*,*) ''
       cell_ptr => cell_next
    end do
  end subroutine delete_children

end module linked_list_module
