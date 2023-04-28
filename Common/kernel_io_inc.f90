!=========================================================================
!
! Included from file kernel_io.F90.
! You might be expected to understand this. --FHJ
!
!=========================================================================

#ifdef READ
  #define READWRITE(x) read ## x
  #define INTENT out
#else
  #define READWRITE(x) write ## x
  #define INTENT in
#endif

#ifdef FORMATTED
  #define READWRITE_FORM(x) READWRITE(_format ## x)
  #define FORMAT , *
#elif defined BINARY
  #define READWRITE_FORM(x) READWRITE(_binary ## x)
  #define FORMAT
#else
  #define READWRITE_FORM(x) READWRITE(x)
#endif

#define READWRITE_FORM_KERNELHEADER READWRITE_FORM(_kernel_header)
#define READWRITE_KERNELHEADER READWRITE(_kernel_header)

#if defined FORMATTED || defined BINARY 

!> Defines a subroutine with the template {read,write}_{formatted,binary}_kernel_header
subroutine READWRITE_FORM_KERNELHEADER(iunit, kernel)
  integer, intent(in) :: iunit !< unit number
  type(kernel_header_t), intent(INTENT) :: kernel !< kernel_header_t type

  integer :: ik

  PUSH_SUB(READWRITE_FORM_KERNELHEADER)

  call READWRITE(_mf_header)(iunit, kernel%mf)
#ifdef READ
  if (kernel%mf%sheader/='KER') then
    write(0,*) 'ERROR: header mismatch (got "'//kernel%mf%sheader//'", expected "KER")'
    call die('Input file is not from a kernel calculation (header="'//kernel%mf%sheader//'")', &
      only_root_writes=.true.)
  endif
#endif

  if (peinf%inode==0) then
    ! General information 
    READWRITE()(iunit FORMAT) kernel%iscreen, kernel%icutv, kernel%ecuts

    ! Variables specific to kernel files: kpts
    READWRITE()(iunit FORMAT) kernel%nk
#ifdef READ
    SAFE_ALLOCATE(kernel%kpts, (3,kernel%nk))
#endif
    do ik = 1, kernel%nk
      READWRITE()(iunit FORMAT) kernel%kpts(1:3, ik)
    enddo

    ! Variables specific to kernel files: everything else
    READWRITE()(iunit FORMAT) kernel%ns, kernel%nspinor, kernel%nvb, kernel%ncb, kernel%n1b, kernel%n2b
    READWRITE()(iunit FORMAT) kernel%theory, kernel%nmat, kernel%storage, kernel%nblocks
    ! Empty records: you can use these slots in the future to extend the file format
    READWRITE()(iunit FORMAT)
    READWRITE()(iunit FORMAT)
    READWRITE()(iunit FORMAT)
    READWRITE()(iunit FORMAT)
    READWRITE()(iunit FORMAT)
  endif
#if defined READ && defined MPI
  if (peinf%npes > 1) then
    ! General information 
    call MPI_BCAST(kernel%iscreen, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kernel%icutv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kernel%ecuts, 1, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)

    ! Variables specific to kernel files: kpts
    call MPI_BCAST(kernel%nk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    if (peinf%inode>0) then
      SAFE_ALLOCATE(kernel%kpts, (3,kernel%nk))
    endif
    call MPI_BCAST(kernel%kpts, 3*kernel%nk, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)

    ! Variables specific to kernel files: everything else
    call MPI_BCAST(kernel%ns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kernel%nspinor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kernel%nvb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kernel%ncb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kernel%n1b, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kernel%n2b, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

    call MPI_BCAST(kernel%theory, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kernel%nmat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kernel%storage, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kernel%nblocks, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

  endif
#endif

  POP_SUB(READWRITE_FORM_KERNELHEADER)

end subroutine READWRITE_FORM_KERNELHEADER

#else

!> Automatically call {read,write}_{formatted,binary}_kernel_header
subroutine READWRITE_KERNELHEADER(iunit, kernel)
  integer, intent(in) :: iunit !< unit number
  type(kernel_header_t), intent(INTENT) :: kernel !< kernel_header_t type

  character(len=16) :: fmt_str
  logical :: is_fmt = .false.

  PUSH_SUB(READWRITE_KERNELHEADER)

  if (peinf%inode==0) then
    inquire(unit=iunit, form=fmt_str)
    if (TRUNC(fmt_str)=='FORMATTED') then
      is_fmt = .true.
    else if (TRUNC(fmt_str)/='UNFORMATTED') then
      call die('Unknown value for formatted string: '//TRUNC(fmt_str), &
        only_root_writes=.true.)
    endif
  endif
  ! FHJ: No need to bcast is_fmt because only inode==0 reads the file.

  if (is_fmt) then
    call READWRITE(_format_kernel_header)(iunit, kernel)
  else
    call READWRITE(_binary_kernel_header)(iunit, kernel)
  endif

  POP_SUB(READWRITE_KERNELHEADER)
 
end subroutine READWRITE_KERNELHEADER

#endif

! these undefs prevent lots of cpp warnings
#undef READWRITE_KERNELHEADER
#undef READWRITE_FORM_KERNELHEADER
#undef READWRITE_FORM
#undef READWRITE
#undef FORMAT
#undef INTENT
