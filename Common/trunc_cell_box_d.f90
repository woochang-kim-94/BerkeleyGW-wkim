!================================================================================
!
! Routines:
!
! 1. trunc_cell_box_d()        Originally By JRD       Last Modified 9/09/2009 (gsm)
!
!    Modification of trunc_cell_box with parallel FFT.
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

module trunc_cell_box_d_m

  use global_m
  use fft_parallel_m
  use fftw_m
  use misc_m

  implicit none

  private

  public ::          &
    trunc_cell_box_d

contains

subroutine trunc_cell_box_d(gvec,verbose,inode,npes,bdot,ncoul,isrtq,vcoul,coulomb_mod)

  type (gspace),      intent(in)    :: gvec         !< Reciprocal space structure
  logical,            intent(in)    :: verbose      !< Flag for extra output on unit 6
  integer,            intent(in)    :: inode, npes  !< MPI parameters: processor rank, number of processors
  real(DP),           intent(in)    :: bdot(3,3)    !< Metric matrix in reciprocal space
  integer,            intent(in)    :: ncoul        !< Number of G-vectors
  integer,            intent(in)    :: isrtq(ncoul) !< Indices of G-vectors
  real(DP),           intent(out)   :: vcoul(ncoul) !< Coulomb potential
  !> Optional coulomb modifier parameters. Used in TDDFT and hybrid functional calculations
  !! they allow one to change the bare coulomb interaction. 
  type(coulomb_modifier_t), optional, intent(in) :: coulomb_mod

  integer :: i1, i2, i3, j1, j2, j3, l1, l2, l3, ig, &
    Nfft(3), dNfft(3), dkmax(3), Nrod, Nplane, i, j, k, NSize(3)
  real(DP) :: r_len, rr(3), t_len, tt(3), adot(3,3), &
    scale, dscale, phase, vimag, vdummy, screening_length, &
    fac1, fac2
  real(DP) :: b(3,3)
  complex(DPC) :: vtemp
  real(DP), allocatable :: xdummy(:)
  complex(DPC), allocatable :: fftbox_temp(:,:,:)
  complex(DPC), allocatable :: fftbox_2D(:,:,:)
  complex(DPC), allocatable :: fftbox_1D(:,:)
  complex(DPC), allocatable :: buffer_2D(:,:,:)
  complex(DPC), allocatable :: buffer_1D(:,:,:)
  integer, allocatable :: inv_indx(:,:,:)
  real(DP) :: erf, erfc


  PUSH_SUB(trunc_cell_box_d)

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

  call fft_set_p(dNfft,Nplane,Nrod)
  SAFE_ALLOCATE(fftbox_2D, (dNfft(1),dNfft(2),Nplane))
  SAFE_ALLOCATE(fftbox_1D, (dNfft(3),Nrod))
  SAFE_ALLOCATE(buffer_2D, (Nrod,Nplane,npes))
  SAFE_ALLOCATE(buffer_1D, (Nrod,Nplane,npes))

  if (verbose) then
    if (inode .eq. 0) then
      write(6,555) dNfft, n_in_box
    endif
  endif
555 format(/,1x,"Cell truncation.",/,2x,"Building truncated Coulomb", &
      1x,"potential using 2-D/1-D FFTs.",/,2x,"2-D/1-D FFT grid =", &
      3i5,/,2x,"Extra points between two points in a 2-D/1-D FFT", &
      1x,"grid =",i3,/)

! Rescale with the volume of the unit cell.

  b = adot
  scale = b(1,1)*(b(2,2)*b(3,3) - b(2,3)**2) &
    + 2*b(1,2)*b(2,3)*b(3,1) &
    - b(2,2)*b(1,3)**2 - b(3,3)*b(1,2)**2
  scale = 2.0d0 * dsqrt(scale)
  !vcoul = vcoul * scale
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

  fftbox_2D(:,:,:)=(0.0d0,0.0d0)
  do i3 = Nplane*inode+1, min(Nplane*(inode+1),dNfft(3))
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
        
        !fftbox_2D(i1,i2,i3-Nplane*inode) = 1.0d0 / r_len
!#BEGIN_INTERNAL_ONLY
        if (fac1 > 0.0d0 .and. fac2 > 0.0d0) then
          fftbox_2D(i1,i2,i3-Nplane*inode) = (fac1*erfc(r_len*screening_length) + fac2*erf(r_len*screening_length))/ r_len
        elseif (screening_length > 0.0d0) then
          if (fac1 > 0.0d0) then
            fftbox_2D(i1,i2,i3-Nplane*inode) = fac1*erfc(r_len*screening_length)/ r_len
          elseif (fac2 > 0.0d0) then
            fftbox_2D(i1,i2,i3-Nplane*inode) = fac2*erf(r_len*screening_length)/ r_len
          endif
        else
!#END_INTERNAL_ONLY
          fftbox_2D(i1,i2,i3-Nplane*inode) = fac1/ r_len
!#BEGIN_INTERNAL_ONLY
        endif
!#END_INTERNAL_ONLY
        
      enddo
    enddo
  enddo

!---------------------------------------
! gsm: Below is equivalent to fft_r2g_p/_s except we multiply by a phase,
!      collect the real part and check that the imaginary part is smaller
!      than the tolerance. Calling fft_r2g_p directly and then altering
!      the result would cost extra memory and code execution.

! Do two-dimensional Fourier transforms:
! V_trunc(r1,r2,r3) -> V_trunc(G_1,G_2,r3)

  NSize(:)=dNfft(:)
  NSize(3)=1
  SAFE_ALLOCATE(fftbox_temp,(NSize(1),NSize(2),NSize(3)))
  do i = 1, Nplane
    fftbox_temp(:,:,1)=fftbox_2D(:,:,i)
    call do_FFT(fftbox_temp,NSize,-1)
    fftbox_2D(:,:,i)=fftbox_temp(:,:,1)
  enddo
  SAFE_DEALLOCATE(fftbox_temp)

! Transfer data from fftbox_2D to fftbox_1D

  buffer_2D(:,:,:)=(0.0d0,0.0d0)
  do k = 1, Nplane
    do i2 = 1, dNfft(2)
      do i1 = 1, dNfft(1)
        i = (i2-1)*dNfft(1)+i1-1
        j = i/Nrod
        buffer_2D(i-j*Nrod+1,k,j+1) = fftbox_2D(i1,i2,k)
      enddo
    enddo
  enddo
#ifdef MPI
  call MPI_Alltoall(buffer_2D(1,1,1),Nrod*Nplane,MPI_COMPLEX_DPC, &
    buffer_1D(1,1,1),Nrod*Nplane,MPI_COMPLEX_DPC,MPI_COMM_WORLD,mpierr)
#else
  buffer_1D(:,:,:)=buffer_2D(:,:,:)
#endif
  do i3 = 1, dNfft(3)
    j = (i3-1)/Nplane
    do i = 1, Nrod
      fftbox_1D(i3,i) = buffer_1D(i,i3-j*Nplane,j+1)
    enddo
  enddo
      
! Do one-dimensional Fourier transforms:
! V_trunc(G_1,G_2,r3) -> V_trunc(G_1,G_2,G_3)

  NSize(:)=1
  NSize(1)=dNfft(3)
  SAFE_ALLOCATE(fftbox_temp,(NSize(1),NSize(2),NSize(3)))
  do i = 1, Nrod
    fftbox_temp(:,1,1)=fftbox_1D(:,i)
    call do_FFT(fftbox_temp,NSize,-1)
    fftbox_1D(:,i)=fftbox_temp(:,1,1)
  enddo
  SAFE_DEALLOCATE(fftbox_temp)

! Collect components of V_trunc for G-vectors in Coulomb list.

  vcoul = 0.0d0
  vimag = 0.0d0
! XXX THREAD?
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
        i = (i2-1)*dNfft(1)+i1-1
        j = i/Nrod

! (gsm) [2010-06-17] there was a bug here
! the singularity of the Coulomb potential was shiftd from
! the origin of the coordinate system by half a grid step

!            if (j == inode) vcoul(ig) = &
!            dble(fftbox_1D(i3,i-j*Nrod+1))

        if (j == inode) then
          phase = dble(j1) * trunc_shift(1) / dble(dNfft(1)) &
            + dble(j2) * trunc_shift(2) / dble(dNfft(2)) &
            + dble(j3) * trunc_shift(3) / dble(dNfft(3))
          phase = 2.0d0 * PI_D * phase
          vtemp = fftbox_1D(i3,i-j*Nrod+1)
          vtemp = vtemp * CMPLX(cos(phase), -sin(phase))
          vcoul(ig) = dble(vtemp)
          vdummy = abs(IMAG(vtemp))
          if (vdummy.gt.vimag) vimag=vdummy
        endif
      enddo
    enddo
  enddo

#ifdef MPI
  vdummy = vimag
  call MPI_Reduce(vdummy, vimag, 1, MPI_REAL_DP, MPI_MAX, 0, MPI_COMM_WORLD, mpierr)
#endif
  if (vimag.gt.TOL_Small) &
    call die("The Coulomb interaction was incorrectly computed as complex: most likely a problem with your FFT library.", &
    only_root_writes = .true.)

  call destroy_fftw_plans()
  SAFE_DEALLOCATE(fftbox_2D)
  SAFE_DEALLOCATE(fftbox_1D)
  SAFE_DEALLOCATE(buffer_2D)
  SAFE_DEALLOCATE(buffer_1D)
  SAFE_DEALLOCATE(inv_indx)

! Global reduction if there is more than one processor.

#ifdef MPI
  SAFE_ALLOCATE(xdummy, (ncoul))
  xdummy = vcoul
  call MPI_Allreduce(xdummy, vcoul, ncoul, MPI_REAL_DP, MPI_SUM, MPI_COMM_WORLD, mpierr)
  SAFE_DEALLOCATE(xdummy)
#endif

  POP_SUB(trunc_cell_box_d)

  return
end subroutine trunc_cell_box_d

end module trunc_cell_box_d_m
