!==============================================================================
!
! Module:
!
! (1) read_wannier_m        Originally By DAS  April 2012
!
!     Reads binary checkpoint file "seedname.chk" produced by Wannier90 v1.2
!     to obtain the coefficients relating the Wannier to Bloch functions.
!     Reference: wannier90-1.2/src/parameters.F90 routine param_write_chkpt
!
!==============================================================================

#include "f_defs.h"

module read_wannier_m

  use global_m
  implicit none

  private

  public :: read_wannier

contains

  subroutine read_wannier(seedname, crys, kg, iwann, nvband, coefficients)
    character(len=*), intent(in) :: seedname !< read file seedname.chk
    type(crystal), intent(in) :: crys !< for checking lattice vectors
    type(grid), intent(in) :: kg !< for assigning k-points
    integer, intent(in) :: iwann !< index of Wannier function to read
    integer, intent(in) :: nvband !< number of valence bands in plotxct calculation
    complex(DPC), intent(out) :: coefficients(:,:) !< u_matrix relating Bloch to Wannier. (band, k-point)
    
    integer, parameter :: chk_unit = 73
    integer :: ib,iw,ii,jj,ik,ik2,num_bands,num_kpts,num_wann,ikmatch
    character :: filename*50, header*33, chkpt1*20
    logical :: have_disentangled
    real(DP) :: real_lattice(3,3), recip_lattice(3,3)
    real(DP), allocatable :: kpt_latt(:,:)
    complex(DPC), allocatable :: u_matrix(:,:,:)
    
    PUSH_SUB(read_wannier)
    
    filename = trim(seedname)//'.chk'
    write(6,'(a)') "Reading file " // trim(filename)
    call open_file(unit=chk_unit,file=filename,form='unformatted',status='old')
    
    if(iwann <= 0) then
      call die("iwann cannot be negative", only_root_writes = .true.)
    endif

    read(chk_unit) header                                       ! Date and time
    write(6,'(a)') header
    
    read(chk_unit) num_bands                                    ! Number of bands
    if(num_bands < nvband) then
      if(peinf%inode == 0) then
        write(0,'(a,9f11.6)') 'read ', num_bands, ' bands but need at least ', nvband, ' = number of valence bands'
      endif
      call die("not enough bands in " // trim(filename), only_root_writes = .true.)
    endif

    read(chk_unit) ! num_exclude_bands (don`t care)             ! Number of excluded bands  
    read(chk_unit) ! (exclude_bands(ib),ib=1,num_exclude_bands) ! Excluded bands 

    read(chk_unit) ((real_lattice(ii,jj),ii=1,3),jj=1,3)        ! Real lattice (a.u.)
    if(any(abs(real_lattice(:,:) - crys%avec(:,:) * crys%alat) > TOL_Small)) then
      if(peinf%inode == 0) then
        write(0,'(a,9f11.6)') 'read lattice vectors ', real_lattice(:,:)
        write(0,'(a,9f11.6)') 'should be ', crys%avec(:,:) * crys%alat
      endif
      call die("incorrect real-space lattice vectors", only_root_writes = .true.)
    endif

    read(chk_unit) ((recip_lattice(ii,jj),ii=1,3),jj=1,3)       ! Reciprocal lattice (a.u.)
    if(any(abs(real_lattice(:,:) - crys%bvec(:,:) * crys%alat) > TOL_Small)) then
      if(peinf%inode == 0) then
        write(0,'(a,9f11.6)') 'read lattice vectors ', recip_lattice(:,:)
        write(0,'(a,9f11.6)') 'should be ', crys%bvec(:,:) * crys%blat
      endif
      call die("incorrect reciprocal lattice vectors", only_root_writes = .true.)
    endif

    read(chk_unit) num_kpts                                     ! Number of k-points
    if(num_kpts < kg%nf) then
      if(peinf%inode == 0) write(0,'(a,i6,a,i6)') 'Need ', kg%nf, 'kpoints but file has only ', num_kpts
      call die("Not enough kpoints in " // trim(filename), only_root_writes = .true.)
    endif

    read(chk_unit) !(mp_grid(ii),ii=1,3) (don`t care)           ! M-P grid
    SAFE_ALLOCATE(kpt_latt, (3, num_kpts))
    read(chk_unit) ((kpt_latt(ii,ik),ii=1,3),ik=1,num_kpts)     ! K-points in lattice vectors
    read(chk_unit) !nntot (don`t care)                          ! Number of nearest k-point neighbours

    read(chk_unit) num_wann                                     ! Number of wannier functions
    if(iwann > num_wann) then
      if(peinf%inode == 0) write(0,'(2a)') 'Requested Wannier function ', iwann, ' > number of Wannier functions ', num_wann
      call die("iwann out of bounds", only_root_writes = .true.)
    endif
    read(chk_unit) chkpt1                                       ! Position of checkpoint
    if(trim(chkpt1) /= 'post_wann') then
      if(peinf%inode == 0) write(0,'(2a)') 'chkpt1 = ', trim(chkpt1)
      call die("Wannier90 checkpoint position must be post_wann.", only_root_writes = .true.)
    endif
    read(chk_unit) have_disentangled                            ! Whether a disentanglement has been performed
    if (have_disentangled) then
      call die("disentanglement not implemented")
!      read(chk_unit) omega_invariant     ! Omega invariant
! lwindow, ndimwin and U_matrix_opt 
!      read(chk_unit) ((lwindow(i,nkp),i=1,num_bands),nkp=1,num_kpts)
!      read(chk_unit) (ndimwin(nkp),nkp=1,num_kpts)
!      read(chk_unit) (((u_matrix_opt(i,j,nkp),i=1,num_bands),j=1,num_wann),nkp=1,num_kpts)
    endif
    SAFE_ALLOCATE(u_matrix, (num_wann, num_wann, num_kpts))
    read(chk_unit) (((u_matrix(ib,iw,ik),ib=1,num_wann),iw=1,num_wann),ik=1,num_kpts) ! U_matrix
!    read(chk_unit) ((((m_matrix(i,j,k,l),i=1,num_wann),j=1,num_wann),k=1,nntot),l=1,num_kpts) ! M_matrix
!    read(chk_unit) ((wannier_centres(i,j),i=1,3),j=1,num_wann)
!    read(chk_unit) (wannier_spreads(ib),i=1,num_wann)
    call close_file(chk_unit)

    ! relate wannier90 kpts to those in plotxct
    do ik = 1, kg%nf
      ikmatch = 0
      do ik2 = 1, num_kpts
        if(all(abs(kg%f(1:3, ik) - kpt_latt(1:3, ik2)) < TOL_Small)) then
          ikmatch = ik2
          exit
        endif
      enddo
      if(ikmatch > 0) then
        do ib = 1, nvband
          coefficients(ib, ik) = u_matrix(nvband - ib + 1, iwann, ikmatch)
        enddo
      else
        if(peinf%inode == 0) write(0,*) 'no match for kpoint ', kg%f(1:3, ik)
        call die("need kpoint not present in " // trim(filename), only_root_writes = .true.)
      endif
    enddo

    SAFE_DEALLOCATE(u_matrix)
    SAFE_DEALLOCATE(kpt_latt)
    
    POP_SUB(read_wannier)
    return
  end subroutine read_wannier
  
end module read_wannier_m
