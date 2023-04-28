!================================================================================
!
! Routines:
!
! 1. trunc_cell_wire()        Originally By MLT       Last Modified 6/12/2008 (JRD)
!
!    Calculate Coulomb interaction with wire boundary conditions. V_coul
!    is calculated using truncation of periodic images following this method:
!    S. Ismail-Beigi, PRB 73, 233103 (2006).
!
!    The Coulomb potential is calculated through 2-D FFT on the extended cell.
!
!    n_in_wire is a (positive integer) parameter for real-space resolution.
!    If the energy cutoff is too small, then the FFT grid may not be very
!    dense and the Coulomb singularity at r=0 may be poorly handled.
!    Larger values of n_in_wire will add extra points in between the original
!    FFT grid. If the original grid has N points covering a space of length
!    L, then the grid used to do 2-dimensional FFTs will have N * n_in_wire
!    points covering a space of length L.
!
!    When V_coul(x,y,G_3+qz) is calculated, we add a shift on the grid in
!    order to avoid the singularity at x = y = 0.
!
!    Output is truncated Vcoul(G_1,G_2,G_3+qz) in reciprocal space such that,
!    without truncation, its value would be:
!
!             Vcoul(ig) = 8 Pi/(q + G)*2
!
!    for each G-vector  G = gvec%components(:,isrtq(ig)) in the list.
!
!    Calculation is distributed over processors but array vcoul returns
!    with global value.
!
!================================================================================

#include "f_defs.h"

module trunc_cell_wire_m

  use bessel_m
  use fft_parallel_m
  use fftw_m
  use global_m
  use misc_m
  implicit none

  private

  public :: &
    trunc_cell_wire

contains

subroutine trunc_cell_wire(gvec,verbose,inode,npes,bdot,qz,ncoul,isrtq,vcoul)
  type (gspace), intent(in)  :: gvec         !< Reciprocal space structure
  logical,       intent(in)  :: verbose      !< Flag for extra output on unit 6
  integer,       intent(in)  :: inode, npes  !< MPI parameters: processor rank, number of processors
  real(DP),      intent(in)  :: bdot(3,3)    !< Metric matrix in reciprocal space, assumed to be symmetric 
                                             !! and such that bdot(1,3) = bdot(2,3) = 0
  real(DP),      intent(in)  :: qz           !< Length of q-vector along z (z component only, not the entire vector!)
  integer,       intent(in)  :: ncoul        !< Number of G-vectors
  integer,       intent(in)  :: isrtq(ncoul) !< Indices of G-vectors
  real(DP),      intent(out) :: vcoul(ncoul) !< Coulomb potential

  integer :: i1, i2, j1, j2, j3, l1, l2, l3, ig, &
    Nfft(3), dNfft(3), dkmax(3), i, j, Gzmin, Gzmax, NSize(3)
  real(DP) :: gpq_z, r_len, rr(3), t_len, tt(3), adot(3,3), &
    scale, dscale, phase, vimag, vdummy
  complex(DPC) :: vtemp
  real(DP), allocatable :: xdummy(:)
  complex(DPC), allocatable :: fftbox_2D(:,:,:)
  integer, allocatable :: inv_indx(:,:,:)

  PUSH_SUB(trunc_cell_wire)
  
  ! Check metric matrix in reciprocal space
  
  if (abs(bdot(1,3)) .gt. TOL_Zero .or. abs(bdot(2,3)) .gt. TOL_Zero) &
    call die('For wire truncation, the 1st and 2nd lattice vectors must be perpendicular to the 3rd.', only_root_writes = .true.)

  ! Initialize FFT grids.
  
  call setup_FFT_sizes(gvec%FFTgrid,Nfft,scale)
  dkmax(1) = gvec%FFTgrid(1) * n_in_wire
  dkmax(2) = gvec%FFTgrid(2) * n_in_wire
  dkmax(3) = 1
  call setup_FFT_sizes(dkmax,dNfft,dscale)

! Initialize real-space metric. If everything goes well, value of
! rr(3) should remain always zero.

  rr = 0.0d0
  tt = 0.0d0
  call invert_matrix(bdot, adot)
  adot = adot * 4.d0 * PI_D * PI_D
  
  do i=1,2
    do j=1,2
      adot(i,j)=adot(i,j)/(1.d0 * dNfft(i) * dNfft(j))
    enddo
  enddo

! Initialize boxes.

  SAFE_ALLOCATE(fftbox_2D, (dNfft(1),dNfft(2),1))
  
  if (verbose) then
    if (inode .eq. 0) then
      write(6,555) dNfft(1:2), n_in_wire
    endif
  endif
555 format(/,1x,"Cell truncation.",/,2x,"Building truncated Coulomb", &
      1x,"potential using 2-D FFTs.",/,2x,"2-D FFT grid =",2i5,/,2x, &
      "Extra points between two points in a 2-D FFT grid =",i3,/)

! Invert G-index.

  SAFE_ALLOCATE(inv_indx, (Nfft(1),Nfft(2),Nfft(3)))
  call fft_map_s(ncoul,Nfft,isrtq,gvec%components,inv_indx)

! Find Gzmin & Gzmax.

  Gzmin = 0
  Gzmax = 0
  do ig = 1, ncoul
    if (gvec%components(3,isrtq(ig)) .lt. Gzmin) &
      Gzmin = gvec%components(3,isrtq(ig))
    if (gvec%components(3,isrtq(ig)) .gt. Gzmax) &
      Gzmax = gvec%components(3,isrtq(ig))
  enddo

  vcoul = 0.0d0
  vimag = 0.0d0
  do j3 = Gzmin, Gzmax
    if (mod(j3 - Gzmin, npes) /= inode) cycle

    l3 = j3 + 1
    if (j3 < 0) l3 = Nfft(3) + l3
    gpq_z = abs(j3 + qz) * sqrt( bdot(3,3) )

! For each G_3 plane, calculate the potential V_trunc(r1,r2,G_3+qz).
! The potential is not zero only inside the Wigner-Seitz cell.
! Because of the logarithmic singularity, we must separate the case
! G_3 + qz = 0.

    do i2 = 1, dNfft(2)
      rr(2) = dble(i2 - 1) + trunc_shift(2)
      do i1 = 1, dNfft(1)
        rr(1) = dble(i1 - 1) + trunc_shift(1)
        
        r_len = INF
        
        do l2 = -ncell+1, ncell
          tt(2) = rr(2) - dble(l2 * dNfft(2))
          do l1 = -ncell+1, ncell
            tt(1) = rr(1) - dble(l1 * dNfft(1))
            t_len = dot_product(tt,matmul(adot,tt))
            if (t_len < r_len) r_len = t_len
          enddo
        enddo
        
        r_len = sqrt(r_len)

! JRD/MJ: You can split the log into two terms and the divergent gpq_z term
! is zero if you use the framework of Ismail-Beigi

        if (abs(gpq_z) < TOL_Small) then
          fftbox_2D(i1,i2,1) = -log( r_len )
        else
          fftbox_2D(i1,i2,1) = dbesk0( gpq_z * r_len )
        endif
        
      enddo
    enddo

! Do a two-dimensional Fourier transform:
! V_trunc(r1,r2,G_3+qz) -> V_trunc(G_1,G_2,G_3+qz)

! Do FFT

    NSize(:)=dNfft(:)
    NSize(3)=1
    call do_fft(fftbox_2D,NSize,-1)

! Collect components of V_trunc for G-vectors in Coulomb list.

    do j2 = - Nfft(2)/2, Nfft(2) - Nfft(2)/2 - 1
      l2 = j2 + 1
      if (j2 < 0) l2 = Nfft(2) + l2
      i2 = j2 + 1
      if (j2 < 0) i2 = dNfft(2) + i2
      do j1 = - Nfft(1)/2, Nfft(1) - Nfft(1)/2 - 1
        l1 = j1 + 1
        if (j1 < 0) l1 = Nfft(1) + l1
        i1 = j1 + 1
        if (j1 < 0) i1 = dNfft(1) + i1
        ig = inv_indx(l1,l2,l3)
        if (ig == 0) cycle

! (gsm) [2010-06-17] there was a bug here
! the singularity of the Coulomb potential was shifted from
! the origin of the coordinate system by half a grid step

!            vcoul(ig) = dble(fftbox_2D(i1,i2))

        phase = dble(j1) * trunc_shift(1) / dble(dNfft(1)) + dble(j2) * trunc_shift(2) / dble(dNfft(2))
        phase = 2.0d0 * PI_D * phase
        vtemp = fftbox_2D(i1,i2,1)
        vtemp = vtemp * CMPLX(cos(phase), -sin(phase))
        vcoul(ig) = dble(vtemp)
        vdummy = abs(IMAG(vtemp))
        if (vdummy.gt.vimag) vimag=vdummy
      enddo
    enddo

  enddo ! j3

#ifdef MPI
  if (npes .gt. 1) then
    vdummy = vimag
    call MPI_Reduce(vdummy, vimag, 1, MPI_REAL_DP, MPI_MAX, 0, MPI_COMM_WORLD, mpierr)
  endif
#endif
  if (vimag.gt.TOL_Small) &
    call die("The Coulomb interaction was incorrectly computed as complex: most likely a problem with your FFT library.", &
    only_root_writes = .true.)

  call destroy_fftw_plans()
  SAFE_DEALLOCATE(fftbox_2D)
  SAFE_DEALLOCATE(inv_indx)
  
! Rescale with the area of the unit cell in the xy plane.

  scale = adot(1,1)*adot(2,2) - adot(1,2)*adot(2,1)
  scale = 4.d0 * sqrt(scale)
  vcoul = vcoul * scale

! Global reduction if there is more than one processor.

#ifdef MPI
  if (npes .gt. 1) then
    SAFE_ALLOCATE(xdummy, (ncoul))
    xdummy = vcoul
    call MPI_Allreduce(xdummy, vcoul, ncoul, MPI_REAL_DP, MPI_SUM, MPI_COMM_WORLD, mpierr)
    SAFE_DEALLOCATE(xdummy)
  endif
#endif

  POP_SUB(trunc_cell_wire)

  return
end subroutine trunc_cell_wire

end module trunc_cell_wire_m
