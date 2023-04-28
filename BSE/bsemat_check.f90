!===============================================================================
!
! bsemat_check       Originally By FHJ (02/10/2013)
!
!     This utility can be used to:
!       1) check the hermiticity of a single BSE kernel file, or
!       2) compare the hermiticity and the equivalence of a restricted (TDA)
!          BSE kernel with an extended one.
!
!     NOTE: all matrices are read and *stored in memory*, so avoid large files!
!
!===============================================================================

#include "f_defs.h"

program bsemat_check

  use global_m
  use misc_m
  implicit none

  character*256 :: fname_tda, fname_ext
  character :: ftype
  integer, parameter :: F_TDA=20, F_EXT=21
  integer :: nargs
  integer :: nk, nc, nv, ns, nb, bse_sz
  integer :: nk_ext, nc_ext, nv_ext, ns_ext
  integer :: ik, ik_, kk1(3), kk2(3), nmatrix
  SCALAR, allocatable :: bsemat_tda(:,:,:), bsemat_ext(:,:,:)
  type(xctinfo) :: xct_tda, xct_ext
  character(len=*), parameter :: fmt1 = &
    "('  -  (matrix ',i1,')  MAX | ',a,' | = ',es8.2e2)"

  nargs = command_argument_count()
  if(nargs == 2 .or. nargs == 3) then
    call get_command_argument(1, ftype)
    call get_command_argument(2, fname_tda)
    if (nargs==3) call get_command_argument(3, fname_ext)
  else
    call die('Usage: bsemat_check.x "X"|"D" bsemat_tda bsemat_extended OR "X"|"D" bsemat')
  endif
  write(6,'(a)') 'TDA kernel file: '//TRUNC(fname_tda)
  if (nargs==3) then
    write(6,'(a)') 'Extended kernel file: '//TRUNC(fname_ext)
  endif
  if (ftype=="X") then
    write(6,'(a)') 'Matrix type: Exchange'
    nmatrix = 1
  elseif (ftype=="D") then
    write(6,'(a)') 'Matrix type: Direct'
    nmatrix = 3
  else
    call die('Unknown matrix type. Options are "X" or "D".')
  endif
  write(6,'(a)')

  call open_file(F_TDA, file=TRUNC(fname_tda), form='unformatted', status='old')
  read(F_TDA) nk, nc, nv, ns
  nb = nv + nc
  ! FHJ: technically, we should use ncb_co for kernel stuff, but the bse_index
  ! function thinks we are dealing with the absorption program.
  xct_tda%nvb_fi = nv
  xct_tda%ncb_fi = nc
  xct_tda%nspin = ns
  bse_sz = (ns*nv*nc*nk)
  SAFE_ALLOCATE(bsemat_tda, (bse_sz,bse_sz,nmatrix))

  if (nargs==2) then
    do ik = 1,nk
      read(F_TDA)
    enddo
  else
    call open_file(F_EXT, file=TRUNC(fname_ext), form='unformatted', status='old')
    read(F_EXT) nk_ext, nc_ext, nv_ext, ns_ext
    if (nk/=nk_ext) call die("Number of k-points don't match.")
    if (ns/=ns_ext) call die("Number of spins don't match.")
    if (nv_ext/=nc_ext) call die("File does not contain an extended BSE kernel.")
    if (nb/=nc_ext) call die("Number of bands are not compatible between files.")
    do ik = 1,nk
      read(F_TDA) ik_, kk1(1:3)
      read(F_EXT) ik_, kk2(1:3)
      if (any(kk1/=kk2)) call die("K-points do not match.")
    enddo
    xct_ext%nvb_fi = nb
    xct_ext%ncb_fi = nb
    xct_ext%nspin = ns
    bse_sz = (ns*nb*nb*nk)
    SAFE_ALLOCATE(bsemat_ext, (bse_sz,bse_sz,nmatrix))
  endif

  write(6,'(a,a)') 'Parameters for TDA kernel file:'
  write(6,'(a,i8)') 'nk=',nk
  write(6,'(a,i8)') 'nc=',nc
  write(6,'(a,i8)') 'nv=',nv
  write(6,'(a,i8)') 'ns=',ns
  write(6,'(a)')
  write(6,'(a)') 'Reading TDA BSE matrices.'
  call read_bsemat(F_TDA, xct_tda, bsemat_tda)
  if (nargs==3) then
    write(6,'(a)') 'Reading extended BSE matrices.'
    call read_bsemat(F_EXT, xct_ext, bsemat_ext)
  endif
  write(6,'(a)')

  write(6,'(a)') 'Checking hermiticity of TDA kernel file.'
  call check_hermiticity(xct_tda, bsemat_tda)
  write(6,'(a)')
  if (nargs==3) then
    write(6,'(a)') 'Checking hermiticity of extended kernel file.'
    call check_hermiticity(xct_ext, bsemat_ext)
    write(6,'(a)')
    write(6,'(a)') 'Checking equivalence of TDA and extended kernels:'
    call check_equivalence()
    write(6,'(a)')
    SAFE_DEALLOCATE(bsemat_ext)
    call close_file(F_EXT)
  endif
  SAFE_DEALLOCATE(bsemat_tda)
  call close_file(F_TDA)
  write(6,'(a)') 'All done.'

contains

  !> Reads bsemat file F_BSE and stores kernel matrices in bsemat (it, itp, imatrix)
  subroutine read_bsemat(F_BSE, xct, bsemat)
    integer, intent(in) :: F_BSE
    type(xctinfo), intent(in) :: xct
    SCALAR, intent(out) :: bsemat(:,:,:)

    integer :: nv_, nc_, ik_, ic_, iv_, imatrix
    integer :: ik, ic, iv, is, it
    integer :: ikp, icp, ivp, isp, itp
    SCALAR, allocatable :: buf(:,:,:,:,:) !<isp, is, ivp, icp, ikp

    PUSH_SUB(bsemat_check.read_bsemat)

    nv_ = xct%nvb_fi
    nc_ = xct%ncb_fi
    SAFE_ALLOCATE(buf, (ns,ns,nv_,nc_,nk))
    do ik = 1, nk
      do imatrix = 1, nmatrix
        do ic = 1, nc_
          do iv = 1, nv_
            read(F_BSE) ik_, ic_, iv_, buf(:,:,:,:,:)
            do is = 1, ns
              it = bse_index(ik, ic, iv, is, xct)
              do ikp = 1, nk
                do icp = 1, nc_
                  do ivp = 1, nv_
                    do isp = 1, ns
                      itp = bse_index(ikp, icp, ivp, isp, xct)
                      bsemat(it, itp, imatrix) = buf(isp, is, ivp, icp, ikp)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    SAFE_DEALLOCATE(buf)

    POP_SUB(bsemat_check.read_bsemat)

  end subroutine read_bsemat

  !> Calculates the deviation from the hermiticity via K - K^H.
  subroutine check_hermiticity(xct, bsemat)
    type(xctinfo), intent(in) :: xct
    SCALAR, intent(in) :: bsemat(:,:,:)

    real(DP) :: diff_max
    integer :: imatrix, it, itp

    PUSH_SUB(bsemat_check.check_hermiticity)

    bse_sz = ns*xct%nvb_fi*xct%ncb_fi*nk
    do imatrix = 1, nmatrix
      diff_max = 0d0
      do itp = 1, bse_sz
        do it = itp, bse_sz
          diff_max = max(diff_max, abs(bsemat(it,itp,imatrix) &
            - MYCONJG(bsemat(itp,it,imatrix))))
        enddo
      enddo
      write(6,fmt1) imatrix, 'K - K^H', diff_max
    enddo

    POP_SUB(bsemat_check.check_hermiticity)

  end subroutine check_hermiticity

  !> Checks the equivalence of the TDA kernel matrices with the extended ones.
  !! For the extended matrices, both vc,v`c` and cv,c`v` transitions are
  !! considered.
  subroutine check_equivalence()
    real(DP) :: diff_max, diff_max2
    integer :: imatrix, it, itp
    integer :: it_ext, it_ext2, itp_ext, itp_ext2
    integer :: ik, ic, iv, is
    integer :: ikp, icp, ivp, isp

    PUSH_SUB(bsemat_check.check_equivalence)

    do imatrix = 1, nmatrix
      diff_max = 0d0
      diff_max2 = 0d0
      do ikp = 1, nk
        do icp = 1, nc
          do ivp = 1, nv
            do isp = 1, ns
              itp = bse_index(ikp, icp, ivp, isp, xct_tda)
              itp_ext = bse_index(ikp, nv+icp, ivp, isp, xct_ext)
              itp_ext2 = bse_index(ikp, ivp, nv+icp, isp, xct_ext)
              do ik = 1, nk
                do ic = 1, nc
                  do iv = 1, nv
                    do is = 1, ns
                      it = bse_index(ik, ic, iv, is, xct_tda)
                      it_ext = bse_index(ik, nv+ic, iv, is, xct_ext)
                      it_ext2 = bse_index(ik, iv, nv+ic, is, xct_ext)
                      diff_max = max(diff_max, abs(bsemat_tda(it,itp,imatrix) &
                        - bsemat_ext(it_ext,itp_ext,imatrix)))
                      diff_max2 = max(diff_max2, abs(bsemat_tda(it,itp,imatrix) &
                        - MYCONJG(bsemat_ext(it_ext2,itp_ext2,imatrix))))
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      write(6,fmt1) imatrix, "K_TDA(vc,v`c`) - K_EXT(vc,v`c`)  ", diff_max
      write(6,fmt1) imatrix, "K_TDA(vc,v`c`) - K_EXT(cv,c`v`)^*", diff_max2
    enddo

    POP_SUB(bsemat_check.check_equivalence)

  end subroutine check_equivalence

end program bsemat_check
