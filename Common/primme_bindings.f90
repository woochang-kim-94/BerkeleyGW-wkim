module primme_bindings_m

#ifdef USEPRIMME

  use, intrinsic :: iso_c_binding
  use global_m
  implicit none

  !include 'primme_f90.inc'

  private
    integer, parameter :: primme_int=c_int64_t
    integer, parameter :: primme_real=c_double

    type primme_t
      type(c_ptr) :: ptr = c_null_ptr
      type(c_ptr) :: fp = c_null_ptr
    contains
      procedure :: init => primme_init
      procedure :: free => primme_free
      !
      procedure :: solve_eigs_real => primme_solve_eigs_real
      procedure :: solve_eigs_cplx => primme_solve_eigs_cplx
      generic :: solve_eigs => solve_eigs_real, solve_eigs_cplx
      !
      procedure :: set_int4 => primme_set_int4
      procedure :: set_int8 => primme_set_int8
      procedure :: set_real => primme_set_real
      procedure :: set_ptr => primme_set_ptr
      procedure :: set_sub => primme_set_sub
      generic :: set => set_int4, set_int8, set_real, set_ptr, set_sub
      !
      procedure :: set_method_int4 => primme_set_method_int4
      procedure :: set_method_int8 => primme_set_method_int8
      generic :: set_method => set_method_int4, set_method_int8
    end type primme_t

  public :: &
    primme_t, primme_int, primme_real

  interface
    integer(c_int) function primme_set_member_c_ptr(primme_ptr, label, val) bind(C, name='primme_set_member')
      import c_ptr, c_int
      type(c_ptr), value :: primme_ptr
      integer(c_int), value :: label
      type(c_ptr), value :: val
    end function primme_set_member_c_ptr
  end interface

  interface
    integer(c_int) function primme_set_member_c_funptr(primme_ptr, label, val) bind(C, name='primme_set_member')
      import c_ptr, c_int, c_funptr
      type(c_ptr), value :: primme_ptr
      integer(c_int), value :: label
      type(c_funptr), value :: val
    end function primme_set_member_c_funptr
  end interface

  interface
    integer(c_int) function primme_set_method_c(val, primme_ptr) bind(C, name='primme_set_method')
      import c_ptr, primme_int, c_int
      integer(primme_int), value :: val
      type(c_ptr), value :: primme_ptr
    end function primme_set_method_c
  end interface

contains

  subroutine primme_die(ierr, sub, label)
    integer(c_int), intent(in):: ierr
    character(len=*), intent(in) :: sub
    integer(c_int), intent(in), optional :: label

    character(len=256) :: errstr

    write(errstr,'(a,i0,a,a)') &
      'Got ierr=', ierr, ' after a call to PRIMME`s ', trim(sub)
    if (present(label)) &
      write(errstr,'(2a,i0)') trim(errstr), ' while setting variable ', label
    call die(errstr)

  end subroutine primme_die

  subroutine primme_init(this)
    class(primme_t), intent(inout) :: this

    call primme_initialize_f77(this%ptr)

  end subroutine primme_init


  subroutine primme_free(this)
    class(primme_t), intent(inout) :: this

    call primme_free_f77(this%ptr)

  end subroutine primme_free


  subroutine primme_solve_eigs_real(this, evals, evecs, rnorms)
    class(primme_t), intent(inout) :: this
    real(kind=kind(1d0)), intent(inout) :: evals(:)
    real(kind=kind(1d0)), intent(inout) :: evecs(:,:)
    real(kind=kind(1d0)), intent(inout) :: rnorms(:)

    integer(c_int) :: ierr

    call dprimme_f77(evals, evecs, rnorms, this%ptr, ierr)
    if (ierr/=0) call primme_die(ierr, 'dprimme_f77')

  end subroutine primme_solve_eigs_real


  subroutine primme_solve_eigs_cplx(this, evals, evecs, rnorms)
    class(primme_t), intent(inout) :: this
    real(kind=kind(1d0)), intent(inout) :: evals(:)
    complex(kind=kind((1d0,1d0))), intent(inout) :: evecs(:,:)
    real(kind=kind(1d0)), intent(inout) :: rnorms(:)

    integer(c_int) :: ierr

    call zprimme_f77(evals, evecs, rnorms, this%ptr, ierr)
    if (ierr/=0) call primme_die(ierr, 'zprimme_f77')

  end subroutine primme_solve_eigs_cplx


  subroutine primme_set_int4(this, label, val)
    class(primme_t), intent(inout) :: this
    integer(c_int), intent(in) :: label
    integer(kind=4), intent(in) :: val

    integer(c_int) :: ierr

    call primme_set_member_f77(this%ptr, label, int(val, kind=primme_int), ierr)
    if (ierr/=0) call primme_die(ierr, 'primme_set_member_f77', label)

  end subroutine primme_set_int4


  subroutine primme_set_int8(this, label, val)
    class(primme_t), intent(inout) :: this
    integer(c_int), intent(in) :: label
    integer(kind=8), intent(in) :: val

    integer(c_int) :: ierr

    call primme_set_member_f77(this%ptr, label, int(val, kind=primme_int), ierr)
    if (ierr/=0) call primme_die(ierr, 'primme_set_member_f77', label)

  end subroutine primme_set_int8


  subroutine primme_set_real(this, label, val)
    class(primme_t), intent(inout) :: this
    integer(c_int), intent(in) :: label
    real(kind=kind(1d0)), intent(in) :: val

    integer(c_int) :: ierr

    call primme_set_member_f77(this%ptr, label, real(val, kind=primme_real), ierr)
    if (ierr/=0) call primme_die(ierr, 'primme_set_member_f77', label)

  end subroutine primme_set_real


  subroutine primme_set_ptr(this, label, val)
    class(primme_t), intent(inout) :: this
    integer(c_int), intent(in) :: label
    type(c_ptr), intent(in) :: val

    integer(c_int) :: ierr

    ierr = primme_set_member_c_ptr(this%ptr, label, val)
    if (ierr/=0) call primme_die(ierr, 'primme_set_member_c_ptr', label)

  end subroutine primme_set_ptr


  subroutine primme_set_sub(this, label, val)
    class(primme_t), intent(inout) :: this
    integer(c_int), intent(in) :: label
    external :: val

    integer(c_int) :: ierr

    ierr = primme_set_member_c_funptr(this%ptr, label, C_FUNLOC(val))
    if (ierr/=0) call primme_die(ierr, 'primme_set_member_c_funptr', label)

  end subroutine primme_set_sub

  subroutine primme_set_method_int4(this, val)
    class(primme_t), intent(inout) :: this
    integer(kind=4), intent(in) :: val

    integer(c_int) :: ierr

    ierr = primme_set_method_c(int(val, kind=primme_int), this%ptr)
    if (ierr/=0) call primme_die(ierr, 'primme_set_method_c')

  end subroutine primme_set_method_int4

  subroutine primme_set_method_int8(this, val)
    class(primme_t), intent(inout) :: this
    integer(kind=8), intent(in) :: val

    integer(c_int) :: ierr

    ierr = primme_set_method_c(int(val, kind=primme_int), this%ptr)
    if (ierr/=0) call primme_die(ierr, 'primme_set_method_c')

  end subroutine primme_set_method_int8

#endif

end module primme_bindings_m
