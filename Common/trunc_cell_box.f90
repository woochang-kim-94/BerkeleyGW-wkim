!================================================================================
!
! Routines:
!
! 1. trunc_cell_box()        Originally By JRD       Last Modified 3/01/2009 (JRD)
!
!    Calculate Coulomb interaction with Cell Box Truncation.
!
!    The Coulomb potential is calculated through FFT on the extended cell.
!
!    n_in_box is a (positive integer) parameter for real-space resolution.
!    If the energy cutoff is too small, then the FFT grid may not be very
!    dense and the Coulomb singularity at r=0 may be poorly handled.
!    Larger values of n_in_box will add extra points in between the original
!    FFT grid. If the original grid has N points covering a space of length
!    L, then the grid used to do 3-dimensional FFTs will have N * n_in_box
!    points covering a space of length L.
!
!    When V_coul(x,y,z) is calculated, we add a shift on the grid in
!    order to avoid the singularity at x = y = z = 0.
!
!    Output is truncated Vcoul(G_1,G_2,G_3) in reciprocal space such that,
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

module trunc_cell_box_m

  use global_m
  use fft_parallel_m
  use fftw_m
  use misc_m

  implicit none

  private

  public ::          &
    trunc_cell_box

contains

subroutine trunc_cell_box(gvec,verbose,bdot,ncoul,isrtq,vcoul,coulomb_mod)

  type (gspace),      intent(in)    :: gvec         !< Reciprocal space structure
  logical,            intent(in)    :: verbose      !< Flag for extra output on unit 6
  real(DP),           intent(in)    :: bdot(3,3)    !< Metric matrix in reciprocal space
  integer,            intent(in)    :: ncoul        !< Number of G-vectors
  integer,            intent(in)    :: isrtq(ncoul) !< Indices of G-vectors
  real(DP),           intent(out)   :: vcoul(ncoul) !< Coulomb potential
  !> Optional coulomb modifier parameters. Used in TDDFT and hybrid functional calculations
  !! they allow one to change the bare coulomb interaction. 
  type(coulomb_modifier_t), optional, intent(in) :: coulomb_mod

  integer :: i1, i2, i3, j1, j2, j3, l1, l2, l3, ig, &
    Nfft(3), dNfft(3), dkmax(3), i, j
  real(DP) :: r_len, rr(3), t_len, tt(3), adot(3,3), &
    scale, dscale, phase, vimag, vdummy, fac1, fac2, &
    screening_length
  real(DP) :: b(3,3)
  complex(DPC) :: vtemp
  complex(DPC), allocatable :: fftbox_3D(:,:,:)
  complex(DPC), allocatable :: dummy(:,:,:)
  integer, allocatable :: inv_indx(:,:,:)

  PUSH_SUB(trunc_cell_box)

! Initialize FFT grids.

  call setup_FFT_sizes(gvec%FFTgrid,Nfft,scale)
  dkmax(1) = gvec%FFTgrid(1) * n_in_box
  dkmax(2) = gvec%FFTgrid(2) * n_in_box
  dkmax(3) = gvec%FFTgrid(3) * n_in_box
  call setup_FFT_sizes(dkmax,dNfft,dscale)

! Initialize real-space metric.

  rr = 0.0d0
  call invert_matrix(bdot, adot)
  adot = adot * 4.0d0 * PI_D * PI_D

  do i =1,3
    do j=1,3
      adot(i,j) = adot(i,j) / dble(dNfft(i) * dNfft(j))
    enddo
  enddo

! Initialize boxes.

  SAFE_ALLOCATE(fftbox_3D, (dNfft(1),dNfft(2),dNfft(3))) 
  SAFE_ALLOCATE(dummy, (dNfft(1),dNfft(2),dNfft(3)))

  if (verbose) then
    if (peinf%inode .eq. 0) then
      write(6,555) dNfft, n_in_box
    endif
  endif
555 format(/,1x,"Cell truncation.",/,2x,"Building truncated Coulomb", &
      1x,"potential using 3-D FFTs.",/,2x,"3-D FFT grid =",3i5,/,2x, &
      "Extra points between two points in a 3-D FFT grid =",i3,/)

! Rescale with the volume of the unit cell.

  b = adot
  scale = b(1,1)*(b(2,2)*b(3,3) - b(2,3)**2) &
    + 2*b(1,2)*b(2,3)*b(3,1) &
    - b(2,2)*b(1,3)**2 - b(3,3)*b(1,2)**2
  scale = 2.0d0 * dsqrt(scale)
  if (present(coulomb_mod)) then
    screening_length = coulomb_mod%screening_length*BOHR
    fac1 = coulomb_mod%short_range_frac_fock * scale
    fac2 = coulomb_mod%long_range_frac_fock * scale
  else
    fac1 = scale
    fac2 = 0.0d0
    screening_length = 0.0d0
  endif


! Invert G-index

  SAFE_ALLOCATE(inv_indx, (Nfft(1),Nfft(2),Nfft(3)))
  call fft_map_s(ncoul,Nfft,isrtq,gvec%components,inv_indx)

! For each G_3 plane, calculate the potential V_trunc(r1,r2,r3).
! The potential is not zero only inside the Wigner-Seitz cell.

  fftbox_3D = (0.0d0,0.0d0)
  do i3 = 1, dNfft(3)
    if (mod(i3-1,peinf%npes) /= peinf%inode) cycle
    rr(3) = dble(i3 - 1) + trunc_shift(3)
    do i2 = 1, dNfft(2)
      rr(2) = dble(i2 - 1) + trunc_shift(2)
      do i1 = 1, dNfft(1)
        rr(1) = dble(i1 - 1) + trunc_shift(1)
        
        r_len = INF
        
        do l3 = -ncell+1, ncell
          tt(3) = rr(3) - dble(l3 * dNfft(3))
          do l2 = -ncell+1, ncell
            tt(2) = rr(2) - dble(l2 * dNfft(2))
            do l1 = -ncell+1, ncell
              tt(1) = rr(1) - dble(l1 * dNfft(1))
              t_len = dot_product(tt,matmul(adot,tt))
              if (t_len < r_len) r_len = t_len
            enddo
          enddo
        enddo

        r_len = sqrt(r_len)
        
        !fftbox_3D(i1,i2,i3) = 1.0d0 / r_len
!#BEGIN_INTERNAL_ONLY
        if (fac1 > 0.0d0 .and. fac2 > 0.0d0) then
          fftbox_3D(i1,i2,i3) = (fac1*erfc(r_len*screening_length) + fac2*erf(r_len*screening_length))/ r_len
        elseif (screening_length > 0.0d0) then
          if (fac1 > 0.0d0) then
            fftbox_3D(i1,i2,i3) = fac1*erfc(r_len*screening_length)/ r_len
          elseif (fac2 > 0.0d0) then
            fftbox_3D(i1,i2,i3) = fac2*erf(r_len*screening_length)/ r_len
          endif
        else
!#END_INTERNAL_ONLY
          fftbox_3D(i1,i2,i3) = fac1/ r_len
!#BEGIN_INTERNAL_ONLY
        endif
!#END_INTERNAL_ONLY

        
      enddo
    enddo
  enddo

#ifdef MPI
  call MPI_Allreduce(fftbox_3D(1,1,1),dummy(1,1,1),dNfft(1)*dNfft(2)*dNfft(3), &
    MPI_COMPLEX_DPC,MPI_SUM,MPI_COMM_WORLD,mpierr)
  fftbox_3D=dummy
#endif

!---------------------------------------
! gsm: Below is equivalent to fft_r2g_s except we multiply by a phase,
!      collect the real part and check that the imaginary part is smaller
!      than the tolerance. Calling fft_r2g_s directly and then altering
!      the result would cost extra memory and code execution.

! Do a three-dimensional Fourier transform:
! V_trunc(r1,r2,r3) -> V_trunc(G_1,G_2,G_3)

  call do_FFT(fftbox_3D,dNfft,-1)

! Collect components of V_trunc for G-vectors in Coulomb list.

  vcoul = 0.0d0
  vimag = 0.0d0
  do j3 = - Nfft(3)/2, Nfft(3) - Nfft(3)/2 - 1
    l3 = j3 + 1
    if (j3 < 0) l3 = Nfft(3) + l3
    i3 = j3 + 1
    if (j3 < 0) i3 = dNfft(3) + i3
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

!            vcoul(ig) = dble(fftbox_3D(i1,i2,i3))

        phase = dble(j1) * trunc_shift(1) / dble(dNfft(1)) &
          + dble(j2) * trunc_shift(2) / dble(dNfft(2)) &
          + dble(j3) * trunc_shift(3) / dble(dNfft(3))
        phase = 2.0d0 * PI_D * phase
        vtemp = fftbox_3D(i1,i2,i3)
        vtemp = vtemp * CMPLX(cos(phase), -sin(phase))
        vcoul(ig) = dble(vtemp)
        vdummy = abs(IMAG(vtemp))
        if (vdummy.gt.vimag) vimag=vdummy
      enddo
    enddo
  enddo

  if (vimag.gt.TOL_Small) &
    call die("The Coulomb interaction was incorrectly computed as complex: most likely a problem with your FFT library.", &
    only_root_writes = .true.)

  call destroy_fftw_plans()
  SAFE_DEALLOCATE(fftbox_3D)
  SAFE_DEALLOCATE(dummy)
  SAFE_DEALLOCATE(inv_indx)

  POP_SUB(trunc_cell_box)

  return
end subroutine trunc_cell_box

end module trunc_cell_box_m
