!>===================================================================
!!
!! Module misc_m
!!
!! Routines:
!!
!! 1. checknorm()        Originally By (SIB)     Last Modified 5/02/2012 (BAB)
!!
!!    Checks normalization of a wavefunction, one spin at a time.
!!    By normalization we mean that sum_l { |z(l,m)|^2 } = 1 for each m.
!!    It aborts if norm is off by more than TOL_Small from 1.
!!
!! 2. compute_norm()      Originally by BAB   5/02/12
!!
!!    Computes the norm of an input wavefunction.
!!
!! 3. get_volume()       Originally By (SIB)     Last Modified 6/12/2008 (JRD)
!!
!!    This assumes that b is a symmetric matrix.  It sets
!!    vol = (2pi)^3 / square_root(|det(b)|)
!!    This makes sense if b is the matrix of dot products of the recip
!!    lattice vectors, so vol is the real space volume.
!!
!! 4. findvector()       Originally By (SIB)     Last Modified 6/12/2008 (JRD)
!!
!!    Looks for the vector in the list of vectors
!!    gvec%components(1:3,1:gvec%ng).  If found, iout is its index.  Otherwise
!!    iout is zero.
!!
!! 5. invert_matrix()    Originally By (SIB)     Last Modified 6/12/2008 (JRD)
!!
!!    Inverts 3x3 matrix.
!!
!! 6. invert_matrix_int()   Originally by DAS  12/28/11  Last Modified 8/27/2012 (BAB)
!!
!!    Like invert_matrix, but for integer input. Dies if output is not integers.
!!
!! 7. compute_det()      Originally By (BAB)     Last Modified 8/27/2012 (BAB)
!!
!!    Computes determinant of 3x3 matrix.
!!
!! 8. compute_det_int()      Originally By (BAB)     Last Modified 8/27/2012 (BAB)
!!
!!    Like compute_det, but for integer input. Does not die if output not integer
!!
!! 9. compute_cofac()        Originally By (BAB)      Last Modified 9/4/2012 (BAB)
!!
!!    Computes matrix of cofactors for input matrix.
!!
!! 10. compute_cofac_int()        Originally By (BAB)      Last Modified 9/4/2012 (BAB)
!!
!!    Computes matrix of cofactors for input matrix, but for integer input.
!!
!! 11. procmem()          Originally By (gsm)     Last Modified 4/14/2009 (gsm)
!!
!!     Determines the amount of free memory per processor
!!     from the proc file system
!!
!! 12. sizeof_scalar()    Originally By (DAS)     Last Modified 1/25/2011 (DAS)
!!
!!     Return the size of the SCALAR type, for memory estimation.
!!
!! 13. voigt()            Originally By (gsm)     Last Modified 1/31/2011 (gsm)
!!
!!     Returns Voigt function (convolution of Gaussian and Lorentzian).
!!     Based on the rational approximation to the complex error function
!!     from A. K. Hui, B. H. Armstrong and A. A. Wray,
!!     "Rapid computation of the Voigt and complex error functions,"
!!     Journal of Quantitative Spectroscopy and Radiative Transfer,
!!     Volume 19, Issue 5, Pages 509 - 516, Year 1978.
!!
!! 14. k_range()        Originally By gsm       Last Modified 8/18/2010 (gsm)
!!
!!    Translates k-point kpt(1:3) to [0,1) interval. Returns G-vector gpt(1:3)
!!    that brings kpt(1:3) to [0,1) interval. The interval is satisfied within
!!    a given tolerance, the actual interval is [-tol,1-tol).
!!
!! 15. bse_index()         Originally by DAS   4/19/12
!!
!!    Defines the mapping of k-point, conduction band, valence band, and spin
!!    into a single index (traditionally called 'ikcvs'), as used in BSE codes.
!!
!! 16. c_to_f_string()
!!
!!    Converts a C string to a Fortran string.
!!
!! 17. get_host_name()
!!
!!    Get the hostname using a call to POSIX`s gethostname.
!!
!!===================================================================

#include "f_defs.h"

module misc_m

  use, intrinsic  :: iso_c_binding
  use global_m
  use blas_m
  implicit none

  private

  public ::        &
    checknorm,     &
    compute_norm,  &
    get_volume,    &
    findvector,    &
    invert_matrix, &
    invert_matrix_int, &
    compute_det,   &
    compute_det_int,   &
    procmem,       &
    sizeof_scalar, &
    voigt,         &
    k_range,       &
    bse_index

contains

!> Checking normalization of only one spin component at a time
  subroutine checknorm(filename,iband,ik,ng,ispin,nspinor,wfn)
    character (len=*), intent(in) :: filename
    integer, intent(in) :: iband,ik,ng,ispin,nspinor
    SCALAR, intent(in) :: wfn(:,:) !< (ng,nspin*nspinor)

    real(DP) :: xnorm

    PUSH_SUB(checknorm)

    if (.not.peinf%check_norms) return
    call compute_norm(xnorm,ispin,ng,nspinor,wfn)
    if(abs(xnorm - 1.0d0) > TOL_Small) then
      write(0,555) TRUNC(filename),abs(xnorm-1.0d0),iband,ispin,ik
555   format(1x,'Wavefunction is not normalized in file',1x,a,/,&
        3x,'abs(norm - 1) =',f10.7,/,&
        3x,'iband =',i6,1x,'ispin =',i2,1x,'ik =',i6,/)
      call die("Incorrect normalization.")
    endif

    POP_SUB(checknorm)
    
    return
  end subroutine checknorm

!=====================================================================

  subroutine compute_norm(xnorm,ispin,ng,nspinor,wfn)
    real(DP), intent(out) :: xnorm
    integer, intent(in) :: ispin,ng,nspinor
    SCALAR, intent(in) :: wfn(:,:) !< (ng,nspin*nspinor)

    integer :: ispinor
    real(DP) :: vnorm(nspinor)

    PUSH_SUB(compute_norm)

    do ispinor=1,nspinor
      vnorm(ispinor) = blas_nrm2(ng, wfn(:,ispin*ispinor), 1)
    enddo
    xnorm = sqrt(sum(vnorm(:)**2))

    POP_SUB(compute_norm)
    
    return
  end subroutine compute_norm

!=====================================================================

  subroutine get_volume(vol,b)
    real(DP), intent(out) :: vol
    real(DP), intent(in)  :: b(3,3)
    
    PUSH_SUB(get_volume)
    
    vol = b(1,1)*(b(2,2)*b(3,3) - b(2,3)**2) &
      + 2*b(1,2)*b(2,3)*b(3,1) &
      - b(2,2)*b(1,3)**2 - b(3,3)*b(1,2)**2
    vol = sqrt(abs(vol))
    vol = ((2.0d0*PI_D)**3)/vol
    
    POP_SUB(get_volume)

    return
  end subroutine get_volume

!=====================================================================

  subroutine findvector(iout,kk,gvec)
    integer, intent(out) :: iout
    integer, intent(in)  :: kk(3)
    type (gspace), intent(in) :: gvec
    
! no push/pop since called too frequently

    iout=((kk(1)+gvec%FFTgrid(1)/2)*gvec%FFTgrid(2)+kk(2)+gvec%FFTgrid(2)/2)* &
      gvec%FFTgrid(3)+kk(3)+gvec%FFTgrid(3)/2+1
    if (iout .ge. 1 .and. iout .le. gvec%nFFTgridpts) then
      iout=gvec%index_vec(iout)
      if (iout .ge. 1 .and. iout .le. gvec%ng) then
        if (any(kk(1:3) /= gvec%components(1:3, iout))) iout = 0
      else
        iout = 0
      endif
    else
      iout = 0
    endif

    return
  end subroutine findvector

!=====================================================================

  subroutine invert_matrix(mat, inv)
    real(DP), intent(in) :: mat(:,:)  !< (3,3)
    real(DP), intent(out) :: inv(:,:) !< (3,3)
    
    real(DP) :: aa(3,3), det
    
    PUSH_SUB(invert_matrix)
    
    call compute_cofac(mat, aa)
    call compute_det(mat, det)
    if (abs(det) .lt. TOL_ZERO) call die('Cannot invert singular matrix.')
    
    inv(1:3, 1:3) = aa(1:3, 1:3) / det
    
    POP_SUB(invert_matrix)
    
    return
  end subroutine invert_matrix

!================================================================================

  subroutine invert_matrix_int(mat, inv)
    integer, intent(in)  :: mat(:,:) !< (3,3)
    integer, intent(out) :: inv(:,:) !< (3,3)

    integer :: aa(3,3), det
    
    PUSH_SUB(invert_matrix_int)
    
    call compute_cofac_int(mat, aa)
    call compute_det_int(mat, det)
    if (det == 0) call die('Cannot invert singular matrix.')

    inv(1:3, 1:3) = aa(1:3, 1:3) / det
    
    if (any(inv(1:3, 1:3) * det /= aa(1:3, 1:3))) then
      write(0,*) 'determinant = ', det
      call die('Inverse of this integer matrix is not an integer matrix.')
    endif
    
    POP_SUB(invert_matrix_int)
    
    return
  end subroutine invert_matrix_int

!=====================================================================

  subroutine compute_det(mat, det)
    real(DP), intent(in) :: mat(:,:)  !< (3,3)
    real(DP), intent(out) :: det
    
    real(DP) :: aa(3,3)

    PUSH_SUB(compute_det)
    
!> Compute matrix of cofactors

    call compute_cofac(mat, aa)

!> Compute determinant

    det = sum(mat(1, 1:3) * aa(1:3, 1))
    
    POP_SUB(compute_det)
    
    return
  end subroutine compute_det

!================================================================================

  subroutine compute_det_int(mat, det)
    integer, intent(in) :: mat(:,:)  !< (3,3)
    integer, intent(out) :: det
    
    integer :: aa(3,3)

    PUSH_SUB(compute_det_int)
    
!> Compute matrix of cofactors

    call compute_cofac_int(mat, aa)

!> Compute determinant

    det = sum(mat(1, 1:3) * aa(1:3, 1))
    
    POP_SUB(compute_det_int)
    
    return
  end subroutine compute_det_int

!================================================================================

  subroutine compute_cofac(mat, aa)
    real(DP), intent(in) :: mat(:,:)  !< (3,3)
    real(DP), intent(out) :: aa(3,3)
    
    PUSH_SUB(compute_cofac)
    
!> Compute matrix of cofactors

    aa(1,1) =  mat(2,2) * mat(3,3) - mat(2,3) * mat(3,2)
    aa(2,1) = -mat(2,1) * mat(3,3) + mat(2,3) * mat(3,1)
    aa(3,1) =  mat(2,1) * mat(3,2) - mat(2,2) * mat(3,1)
    aa(1,2) = -mat(1,2) * mat(3,3) + mat(1,3) * mat(3,2)
    aa(2,2) =  mat(1,1) * mat(3,3) - mat(1,3) * mat(3,1)
    aa(3,2) = -mat(1,1) * mat(3,2) + mat(1,2) * mat(3,1)
    aa(1,3) =  mat(1,2) * mat(2,3) - mat(1,3) * mat(2,2)
    aa(2,3) = -mat(1,1) * mat(2,3) + mat(1,3) * mat(2,1)
    aa(3,3) =  mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)
    
    POP_SUB(compute_cofac)
    
    return
  end subroutine compute_cofac

!================================================================================

  subroutine compute_cofac_int(mat, aa)
    integer, intent(in) :: mat(:,:)  !< (3,3)
    integer, intent(out) :: aa(3,3)
    
    PUSH_SUB(compute_cofac_int)
    
!> Compute matrix of cofactors

    aa(1,1) =  mat(2,2) * mat(3,3) - mat(2,3) * mat(3,2)
    aa(2,1) = -mat(2,1) * mat(3,3) + mat(2,3) * mat(3,1)
    aa(3,1) =  mat(2,1) * mat(3,2) - mat(2,2) * mat(3,1)
    aa(1,2) = -mat(1,2) * mat(3,3) + mat(1,3) * mat(3,2)
    aa(2,2) =  mat(1,1) * mat(3,3) - mat(1,3) * mat(3,1)
    aa(3,2) = -mat(1,1) * mat(3,2) + mat(1,2) * mat(3,1)
    aa(1,3) =  mat(1,2) * mat(2,3) - mat(1,3) * mat(2,2)
    aa(2,3) = -mat(1,1) * mat(2,3) + mat(1,3) * mat(2,1)
    aa(3,3) =  mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)
    
    POP_SUB(compute_cofac_int)
    
    return
  end subroutine compute_cofac_int

!================================================================================

  subroutine procmem(mem, ranks_per_node, nfreq_group)
    real(DP), intent(out) :: mem !< Memory per MPI rank (lower bound)
    integer, intent(out) :: ranks_per_node !< # of MPI ranks per node (upper bound)
    integer, intent(in), optional :: nfreq_group

    integer, parameter :: host_len=64
    integer :: ierr,eof,info,iunit,m,n,p,i,j,pagesize
    real(DP) :: x,y,mac_m,mac_n
!    integer :: ntot
!    real(DP) :: xtot,ytot ! we do not use the total memory actually
    character(len=256) :: s, filename, my_hostname
    character(len=host_len), allocatable :: hostnames(:)

    PUSH_SUB(procmem)

!-----------------------------------------------------
!> determine the amount of free memory per node in kB

    m=0
    iunit=14
    x=0

    call open_file(unit=iunit,file='/proc/meminfo',form='formatted',iostat=ierr,status='old')
    if (ierr.eq.0) then
      eof=0
      do while(eof.eq.0)
        read(iunit,'(a)',iostat=eof)s
        if (s(1:7).eq."MemFree") then
          read(s(9:),*)n
          m=m+n
        endif
!        if (s(1:8).eq."MemTotal") then
!          read(s(10:),*)ntot
!        endif
        if (s(1:6).eq."Cached") then
          read(s(8:),*)n
          m=m+n
        endif
      enddo
      x=dble(m)/dble(peinf%npes)
      call close_file(iunit)
    endif

    if(m == 0) then
      !> this is for Mac OS
      !! total memory is accessible instead from sysctl -n hw.usermem
      write(filename,'(a,i9.9)') 'vm_stat_', peinf%inode
      call system("vm_stat > " + TRUNC(filename) + " 2> /dev/null")
      !> Fortran 2008 would use execute_command_line instead
      !! even if the command failed, still open file in order to delete it
      call open_file(unit=iunit,file=TRUNC(filename),form='formatted',iostat=ierr,status='old')
      if (ierr.eq.0) then
        eof=0
        do while(eof.eq.0)
          read(iunit,'(a)',iostat=eof)s
          if (s(1:45).eq."Mach Virtual Memory Statistics: (page size of") then
            read(s(46:),*)pagesize ! in bytes
          endif
          if (s(1:11).eq."Pages free:") then
            read(s(12:),*) mac_n
            mac_m = mac_m + mac_n
          endif
          if (s(1:18).eq."Pages speculative:") then
            read(s(19:),*) mac_n
            mac_m = mac_m + mac_n
          endif
        enddo
        call close_file(iunit, delete = .true.)
        x = mac_m * dble(pagesize) / dble(peinf%npes * 1024) ! to kB
      endif
    endif

!> === Example output from vm_stat ===
!! Mach Virtual Memory Statistics: (page size of 4096 bytes)
!! Pages free:                           2886.
!! Pages active:                       139635.
!! Pages inactive:                      66906.
!! Pages speculative:                    2376.
!! Pages wired down:                    50096.
!! "Translation faults":            123564742.
!! Pages copy-on-write:              10525831.
!! Pages zero filled:                53274329.
!! Pages reactivated:                  739514.
!! Pageins:                           2282166.
!! Pageouts:                           306727.
!! Object cache: 25 hits of 522230 lookups (0% hit rate)

    if(m == 0 .and. mac_m == 0) then ! BSD
      !> http://mario79t.wordpress.com/2008/08/29/memory-usage-on-freebsd/
      !! -bash-2.05b$ sysctl vm.stats.vm.v_free_count
      !! vm.stats.vm.v_free_count: 29835
      !! -bash-2.05b$ sysctl vm.stats.vm.v_page_count
      !! vm.stats.vm.v_page_count: 124419
      !! -bash-2.05b$ sysctl hw.pagesize
      !! hw.pagesize: 4096
      write(filename,'(a,i9.9)') 'sysctl_', peinf%inode
      call system("sysctl -a > " + TRUNC(filename) + " 2> /dev/null")
      !> Fortran 2008 would use execute_command_line instead
      !! even if the command failed, still open file in order to delete it
      call open_file(unit=iunit,file=TRUNC(filename),form='formatted',iostat=ierr,status='old')
      if (ierr.eq.0) then
        eof=0
        do while(eof.eq.0)
          read(iunit,'(a)',iostat=eof)s
          if (s(1:12).eq."hw.pagesize:") then
            read(s(13:),*)pagesize ! in bytes
          endif
          if (s(1:25).eq."vm.stats.vm.v_free_count:") then
            read(s(26:),*) mac_n
            mac_m = mac_m + mac_n
          endif
          if (s(1:26).eq."vm.stats.vm.v_cache_count:") then
            read(s(27:),*) mac_n
            mac_m = mac_m + mac_n
          endif
        enddo
        call close_file(iunit, delete = .true.)
        x = mac_m * dble(pagesize) / dble(peinf%npes * 1024) ! to kB
      endif
    endif

!    xtot=dble(ntot)/dble(peinf%npes)
    
#ifdef MPI
    call MPI_Allreduce(x,y,1,MPI_REAL_DP,MPI_SUM,MPI_COMM_WORLD,mpierr)
!    call MPI_Allreduce(xtot,ytot,1,MPI_REAL_DP,MPI_SUM,MPI_COMM_WORLD,mpierr)
#else
    y=x
!    ytot=xtot
#endif

!----------------------------------------------
! Determine the number of processors per node
! FHJ: we do this in parallel. Each MPI rank compares with all other hostnames,
! and figures out the number of ranks running on the same node. At the end, we
! pick the upper bound for ranks_per_node.

    if(present(nfreq_group)) then
      SAFE_ALLOCATE(hostnames, (peinf%npes_orig))
    else
      SAFE_ALLOCATE(hostnames, (peinf%npes))
    endif
    my_hostname = get_host_name()
    hostnames(peinf%inode+1) = my_hostname(1:host_len)
    
#ifdef MPI
    call MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, hostnames(1), &
      host_len, MPI_BYTE, MPI_COMM_WORLD, mpierr)
    ranks_per_node = 0
    do i = 1, peinf%npes
      ! FHJ: no need to trim, all strings have the same length
      if (hostnames(peinf%inode+1)==hostnames(i)) ranks_per_node = ranks_per_node + 1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE, ranks_per_node, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
#else
    ranks_per_node = 1
#endif

!-----------------------------------
!> report the available memory in B

    if (ranks_per_node>1) y = y / dble(ranks_per_node)
    mem=y*1024.0d0

!-----------------------------------
!> warn if zero memory

    if (mem .lt. TOL_Small .and. peinf%inode .eq. 0) then
      write(0,666)
666   format(1x,'WARNING: estimation of memory available is zero, probably failed.',/)
    endif

    POP_SUB(procmem)
  
    return
  end subroutine procmem

!================================================================================

!> for memory estimation, tell what size of SCALAR type is
  integer function sizeof_scalar()

#ifdef NOSIZEOF
!> if no sizeof, take a guess

#ifdef CPLX
    sizeof_scalar = 16
#else
    sizeof_scalar = 8
#endif

#else

    SCALAR :: dummy
    sizeof_scalar = sizeof(dummy)

#endif

  end function sizeof_scalar

!================================================================================

  real(DP) function voigt(x, sigma, gamma)
    real(DP), intent(in) :: x, sigma, gamma

    real(DP), parameter :: a0 = 122.607931777104326d0
    real(DP), parameter :: a1 = 214.382388694706425d0
    real(DP), parameter :: a2 = 181.928533092181549d0
    real(DP), parameter :: a3 =  93.155580458138441d0
    real(DP), parameter :: a4 =  30.180142196210589d0
    real(DP), parameter :: a5 =   5.912626209773153d0
    real(DP), parameter :: a6 =   0.564189583562615d0
    real(DP), parameter :: b0 = 122.607931773875350d0
    real(DP), parameter :: b1 = 352.730625110963558d0
    real(DP), parameter :: b2 = 457.334478783897737d0
    real(DP), parameter :: b3 = 348.703917719495792d0
    real(DP), parameter :: b4 = 170.354001821091472d0
    real(DP), parameter :: b5 =  53.992906912940207d0
    real(DP), parameter :: b6 =  10.479857114260399d0

    complex(DPC) :: z, zh, f

    PUSH_SUB(voigt)

    if (sigma .lt. TOL_Zero .or. gamma.lt.-TOL_Zero) &
     call die('Voigt function invalid broadening')

    z = CMPLX(abs(x), gamma) / (sqrt(2.0d0) * sigma)
    zh = CMPLX(IMAG(z), -dble(z))
    f = ((((((a6*zh + a5)*zh + a4)*zh + a3)*zh + a2)*zh + a1)*zh + a0) / &
     (((((((zh + b6)*zh + b5)*zh + b4)*zh + b3)*zh + b2)*zh + b1)*zh + b0)
    if (x .lt. 0.0d0) f = conjg(f)
    voigt = dble(f) / (sqrt(2.0d0 * PI_D) * sigma)

    POP_SUB(voigt)

    return
  end function voigt

!================================================================================

  subroutine k_range(kpt, gpt, tol)
    real(DP), intent(inout) :: kpt(3)
    integer, intent(out) :: gpt(3)
    real(DP), intent(in) :: tol

    integer :: ii

    ! no push_sub, called too frequently
    
    do ii = 1, 3
      gpt(ii) = 0
      do while (kpt(ii) .lt. -tol)
        gpt(ii) = gpt(ii) + 1
        kpt(ii) = kpt(ii) + 1.0d0
      enddo
      do while (kpt(ii) .ge. 1.0d0 - tol)
        gpt(ii) = gpt(ii) - 1
        kpt(ii) = kpt(ii) - 1.0d0
      enddo
    enddo

    return
  end subroutine k_range

!=====================================================================

  integer function bse_index(ik, ic, iv, is, xct, ncband, nvband)
    integer, intent(in) :: ik, ic, iv, is
    type(xctinfo), intent(in) :: xct
    integer, optional, intent(in) :: ncband !< default is xct%ncb_fi
    integer, optional, intent(in) :: nvband !< default is xct%nvb_fi

    integer :: ncband_, nvband_

    ! The optionals are needed for the parallelization scheme sometimes, to be set to 1.
    if(present(ncband)) then
      ncband_ = ncband
    else
      ncband_ = xct%ncb_fi
    endif

    if(present(nvband)) then
      nvband_ = nvband
    else
      nvband_ = xct%nvb_fi
    endif

    ! no push_sub, called too frequently

    bse_index = is + (iv - 1 + (ic - 1 + (ik - 1)*ncband_)*nvband_)*xct%nspin
    return
  end function bse_index

end module misc_m
