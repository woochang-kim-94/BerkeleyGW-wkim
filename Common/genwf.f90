#include "f_defs.h"

module genwf_m

  use global_m
  use blas_m
  use gmap_m
  use susymmetries_m
  use random_m
  use input_utils_m
  use sort_m
  use misc_m

  implicit none

  private

  public :: genwf

contains

  !> Generate the wavefunctions at a particular k-point of the full Brillouin
  !! zone by rotating the wavefunctions at a symmetry-related k-point of the
  !! irreducible Brillouin zone.
  subroutine genwf(crys,gvec,kg,syms,wfn,ik,ik2,nspin,nband,work,intwfn,iwriteint,is_cond)
    type (crystal), intent(in) :: crys
    type (gspace), intent(in) :: gvec
    type (grid), intent(in) :: kg
    type (symmetry), intent(in) :: syms
    type (wavefunction), intent(out) :: wfn
    integer, intent(in) :: ik  !< FIXME very hard to understand what is ik vs. ik2
    integer, intent(in) :: ik2  !< Index of the k-point in the full BZ.
    integer, intent(in) :: nspin  !< Number of spins
    integer, intent(in) :: nband  !< Number of bands
    type (work_genwf), intent(inout) :: work
    type (int_wavefunction), intent(in) :: intwfn
    integer, intent(in) :: iwriteint  !< 0 = write/read wavefunction from disc
                                      !! 1 = use MPI communication
    logical, intent(in) :: is_cond
    
    integer :: iwrite
    character :: filename*20, wfnname*10
    integer :: irk, iunit
    integer :: hh, ii, jj, kk, eof, ig
    
    integer, allocatable :: isorti(:)
    real(DP), allocatable :: ekin(:)
#ifdef CPLX
    SCALAR, allocatable :: zvtemp(:,:)
    complex(DP) :: umtrx(2,2)
#endif
    
    PUSH_SUB(genwf)
  
    if (iwriteint .eq. 0) then
      if(peinf%inode.lt.10000) then
        if(is_cond) then
          write(filename,'(a,i4.4)') 'INT_CWFN_', peinf%inode
          iunit=128+(2*peinf%inode)+1
          wfnname = "conduction"
        else
          write(filename,'(a,i4.4)') 'INT_VWFNQ_', peinf%inode
          iunit=128+(2*peinf%inode)+2
          wfnname = "valence"
        endif
      else
        call die('genwf: cannot use more than 9999 nodes')
      endif
    endif
  
  !-----------------------------------------------------------------------
  !     Start looking for the right k-point in unit iunit
  
    if(ik.ne.work%ikold) then
      if (iwriteint .eq. 0) then
        call open_file(iunit,file=filename,form='unformatted',status='old')
        
        eof=0
        read(iunit) irk,work%ng,work%nb,work%ns,work%nspinor
        do while((irk.ne.kg%indr(ik2)).and.(eof.eq.0))
          read(iunit)
          read(iunit,iostat=eof) irk,work%ng,work%nb,work%ns,work%nspinor
        enddo
        if(eof.ne.0) then
          write(0,*) 'BUG: PE ', peinf%inode, ' could not find the ', trim(wfnname), &
            'wavefunctions for k-point ', ik2
          write(0,*) '  (equivalent to k-point ', kg%indr(ik2), 'in the IBZ) in file ', filename
          call die("genwf wfns missing")
        endif
      else
        iwrite=0
        do ii=1, peinf%ikt(peinf%inode+1)
          if(peinf%ik(peinf%inode+1,ii).eq.ik) then
            iwrite=ii
            work%ng=intwfn%ng(ii)
            work%nb=nband
            work%ns=nspin
            work%nspinor=intwfn%nspinor
          endif
        enddo
      endif
      
      if(work%ikold.ne.0) then
        SAFE_DEALLOCATE_P(work%cg)
        SAFE_DEALLOCATE_P(work%ph)
        SAFE_DEALLOCATE_P(work%ind)
        SAFE_DEALLOCATE_P(work%isort)
      endif
      SAFE_ALLOCATE(work%cg, (work%ng,work%nb,work%ns*work%nspinor))
      SAFE_ALLOCATE(work%ind, (work%ng))
      SAFE_ALLOCATE(work%ph, (work%ng))
      SAFE_ALLOCATE(work%isort, (gvec%ng))
    endif
  
    wfn%ng=work%ng
    wfn%nband=work%nb
    wfn%nspin=work%ns
    wfn%nspinor=work%nspinor
    if (work%ns.ne.nspin) then
      write(0,*) 'spin number mismatch in file ', filename, nspin, work%ns
      call die("genwf spin number mismatch")
    endif
    
    SAFE_ALLOCATE(wfn%cg, (wfn%ng,wfn%nband,wfn%nspin*wfn%nspinor))
    SAFE_ALLOCATE(wfn%isort, (gvec%ng))
    if(ik.ne.work%ikold) then
  !       Read the wavefunctions for the rk-kpoint
      if (iwriteint .eq. 0) then
        read(iunit) (work%isort(ii),ii=1,gvec%ng), &
          (((work%cg(ii,jj,kk),ii=1,wfn%ng),jj=1,wfn%nband), kk=1,wfn%nspin*wfn%nspinor)
      else
        work%isort(:)=intwfn%isort(:,iwrite)
        work%cg(1:wfn%ng,:,:)=intwfn%cgk(1:wfn%ng,:,:,iwrite)
      endif
  
  !       Compute inverse index array of Fourier components around rk-kpoint
      SAFE_ALLOCATE(isorti, (gvec%ng))
      isorti(:)=0
      do ii=1,wfn%ng
        isorti(work%isort(ii))=ii
      enddo
  
  !       Compute index array of Fourier components around fk-kpoint
      SAFE_ALLOCATE(ekin, (gvec%ng))
      call kinetic_energies(gvec, crys%bdot, ekin, qvec = kg%f(:, ik2))
      call sortrx(gvec%ng, ekin, work%isort, gvec = gvec%components)
      SAFE_DEALLOCATE(ekin)
  
  !       Find ind and ph relating wavefunctions in fk to rk-kpoint
      work%ind=0
      work%ph=ZERO
      call gmap(gvec,syms,wfn%ng,kg%itran(ik2), &
        kg%kg0(:,ik2),work%isort,isorti,work%ind,work%ph,.true.)
      SAFE_DEALLOCATE(isorti)

  !       Compute and renormalize wavefunctions
      do kk=1,wfn%nspin
        do jj=1,wfn%nband
          do hh=1,wfn%nspinor
            do ii=1,wfn%ng
              if (work%ind(ii) .gt. 0) then
                wfn%cg(ii,jj,kk*hh)=work%ph(ii)*work%cg(work%ind(ii),jj,kk*hh)
              else
                wfn%cg(ii,jj,kk*hh)=ZERO
              endif
            enddo
          enddo
          call checknorm('wfn%cg',jj,ik,wfn%ng,kk,wfn%nspinor,wfn%cg(1:wfn%ng,jj,:))
        enddo
      enddo
      work%cg=wfn%cg

  !   In spinor case, we must rotate spinors according to spinor rotation matrix umtrx
#ifdef CPLX
      if (work%nspinor.eq.2) then
        do ig=1,work%ng
          SAFE_ALLOCATE(zvtemp, (work%nb,work%ns*work%nspinor))
          call susymmetries(crys%bvec(:,:),syms%mtrx(:,:,kg%itran(ik2)),umtrx,kg%itran(ik2))
          zvtemp = matmul(work%cg(ig,:,:),umtrx)
          work%cg(ig,:,:) = zvtemp
          SAFE_DEALLOCATE(zvtemp)
        enddo
      endif
#endif

      if(iwriteint == 0) call close_file(iunit)
    endif
    
    wfn%cg=work%cg
    wfn%isort=work%isort
    
    work%ikold=ik
    
    POP_SUB(genwf)
    
    return
  end subroutine genwf
  
end module genwf_m
