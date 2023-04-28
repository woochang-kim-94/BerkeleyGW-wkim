! a header is needed here

#include "f_defs.h"

program kptsummedspectralfcn

  use global_m
  implicit none

  integer :: nkpts,nval,nfreq_ext,nfreq_F,i_en, &
    ii,jj,ll,ik,iband,iw,nfreq_Fadd,nlines,&
    iispin,iiband,iiw,nfreq,ind,&
    i_method,ntime,it,nfreq_zero,diverge_int,i_shift
  real(DP) :: ecor,bx,Resx,Rech,Resig,mefact1,intfact,diff_small,freq_interv,&
    Imsx,Imch,Imsig,Imeqp0,kx,ky,kz,Ew,mefact2,diff, &
    alpha,gamma,gamma_big,Eqp,E_f,slope,beta,betap,&
    w0,ff,timerange,timeint,eta,num,freqsm,band,kpt
  complex(DPC) :: imag,Z_qp
  character :: Aname*15
  real(DP), allocatable :: lda(:,:),wts(:),freq(:),vxc(:,:),time_array(:),&
    Resigma(:,:,:),Imsigma(:,:,:),Imsigma2(:,:,:),&
    freq_ext(:),convolt(:),convol1(:),convol2(:),&
    convol3(:),convol4(:),convolint(:),Awb(:),&
    Awbwme(:),Gauss(:),Reeqp0(:,:),dResigma(:,:,:),&
    dImsigma(:,:,:),Reeqp1(:,:),daImsigma(:,:,:),&
    d3Imsigma(:,:,:),func_lda(:,:,:),func_qp(:,:,:),&
    d2Imsigma(:,:,:),matr_el(:,:)
                     
  complex(DPC), allocatable :: Cst(:),Cs_bandt(:,:),Cs_band_kt(:,:,:),&
    Cqpt(:),Cqp_bandt(:,:),Cqp_band_kt(:,:,:),&
    Aqp(:),Atot(:),Aqp_band(:,:),Atot_band(:,:),&
    Aqp_band_k(:,:,:),Atot_band_k(:,:,:),&
    Aconv_band_k(:,:,:),Cs(:),Cs_band(:,:),Cs_band_k(:,:,:)
  integer, allocatable :: flda_ind(:,:),fqp_ind(:,:)
  
  nlines=35
  ntime = 25001
  timerange = 1.25d2
  timeint = 1d-2
  i_method = 1  ! =1 for 1st order, =2 for full exponentiation
  w0 = 2d0
  Z_qp = CMPLX(0.0d0,0.0d0)
  diverge_int = 2
  E_f = 0.113d0 !Fermi energy for Si from paratec
  freq_interv = 0.1d0
  intfact=(20.0d0*1.325d0)
  mefact1=2.5d0
  mefact2=2.5d0
  nkpts=1
  nfreq_Fadd=5
  nfreq_F=34   !number of frequencies to get to Fermi energy
  nfreq=51     !number of frequencies we have in spectrum.dat
  nfreq_ext=120 !more frequencies needed for convolution smoothing
  nval=1
  imag=CMPLX(0.0d0,1.0d0)
  
  allocate(wts(nkpts))
  allocate(lda(nkpts,nval))
  allocate(fqp_ind(nkpts,nval))
  allocate(flda_ind(nkpts,nval))
  allocate(vxc(nkpts,nval))
  allocate(Cst(ntime))
  allocate(Cs_bandt(ntime,nval))
  allocate(Cs_band_kt(ntime,nval,nkpts))
  allocate(Cqpt(ntime))
  allocate(Cqp_bandt(ntime,nval))
  allocate(Cqp_band_kt(ntime,nval,nkpts))
  allocate(Cs(nfreq))
  allocate(Cs_band(nfreq,nval))
  allocate(Cs_band_k(nfreq,nval,nkpts))
  allocate(Aqp(nfreq))
  allocate(Aqp_band(nfreq,nval))
  allocate(Aqp_band_k(nfreq,nval,nkpts))
  allocate(Atot(nfreq))
  allocate(Atot_band(nfreq,nval))
  allocate(Atot_band_k(nfreq,nval,nkpts))
  allocate(Aconv_band_k(nfreq,nval,nkpts))
  allocate(convolint(nfreq_F+nfreq_Fadd))
  allocate(Awb(nfreq_F+nfreq_Fadd))
  allocate(Awbwme(nfreq_F+nfreq_Fadd))
  allocate(convolt(nfreq))
  allocate(convol1(nfreq))
  allocate(convol2(nfreq))
  allocate(convol3(nfreq))
  allocate(convol4(nfreq))
  allocate(freq(nfreq))
  allocate(freq_ext(nfreq_ext))
  allocate(Gauss(nfreq_ext))
  allocate(Resigma(nfreq,nval,nkpts))
  allocate(Reeqp0(nkpts,nval))
  allocate(Reeqp1(nkpts,nval))
  allocate(Imsigma(nfreq,nval,nkpts))
  allocate(Imsigma2(nfreq,nval,nkpts))
  allocate(dResigma(nfreq,nval,nkpts))
  allocate(dImsigma(nfreq,nval,nkpts))
  allocate(daImsigma(nfreq,nval,nkpts))
  allocate(d2Imsigma(nfreq,nval,nkpts))
  allocate(d3Imsigma(nfreq,nval,nkpts))
  allocate(func_lda(nfreq,nval,nkpts))
  allocate(func_qp(nfreq,nval,nkpts))
  allocate(time_array(ntime))
  allocate(matr_el(nkpts,nval))
  
  Reeqp1=0.0d0
  Reeqp0=0.0d0
  Gauss = 0.0d0
  convolint=0.0d0 
  convolt=0.0d0 
  convol1=0.0d0 
  convol2=0.0d0 
  convol3=0.0d0 
  convol4=0.0d0
  Awb = 0.0d0
  Awbwme = 0.0d0 
  wts = 0.0d0
  freq = 0.0d0
  lda = 0.0d0
  vxc = 0.0d0
  Cs = CMPLX(0.0d0,0.0d0)
  Cs_band = CMPLX(0.0d0,0.0d0)
  Cs_band_k = CMPLX(0.0d0,0.0d0)
  Cst = CMPLX(0.0d0,0.0d0)
  Cs_bandt = CMPLX(0.0d0,0.0d0)
  Cs_band_kt = CMPLX(0.0d0,0.0d0)
  Cqpt = CMPLX(0.0d0,0.0d0)
  Cqp_bandt = CMPLX(0.0d0,0.0d0)
  Cqp_band_kt = CMPLX(0.0d0,0.0d0)
  Aqp = CMPLX(0.0d0,0.0d0)
  Aqp_band = CMPLX(0.0d0,0.0d0)
  Aqp_band_k = CMPLX(0.0d0,0.0d0)
  Atot = CMPLX(0.0d0,0.0d0)
  Atot_band = CMPLX(0.0d0,0.0d0)
  Atot_band_k = CMPLX(0.0d0,0.0d0)
  Aconv_band_k = CMPLX(0.0d0,0.0d0)
  Resigma = 0.0d0
  Imsigma = 0.0d0
  Imsigma2 = 0.0d0
  dResigma = 0.0d0
  dImsigma = 0.0d0
  daImsigma = 0.0d0
  d2Imsigma = 0.0d0
  d3Imsigma = 0.0d0
  func_lda = 0.0d0
  func_qp = 0.0d0
  time_array = 0.0d0
  matr_el = 0.0d0
  
  call open_file(unit=15,file='matrixelements',status='old',form='formatted')
  
  do ii=1,nkpts
    do jj=1,nval
      read(15,*) kpt,band,matr_el(ii,jj)
    enddo
  enddo
  
  call close_file(15)
  
  call open_file(unit=10,file='weights',status='old',form='formatted')
  
  do ii=1,nkpts
    read(10,*) wts(ii)
  enddo
  
  call close_file(10)
  
  call open_file(unit=20,file='sig_hp.log',status='old',form='formatted')
  call open_file(unit=21,file='energies',status='replace',form='formatted')
  
  do ii=1,nlines
    read(20,*)
  enddo
  
  do ik=1,nkpts
    do ii=1,4
      read(20,*)
    enddo
    do jj=1,nval
      read(20,200) ll,lda(ik,jj),ecor,bx,Resx,Rech,Resig,vxc(ik,jj),Reeqp0(ik,jj)
      read(20,201) Imsx,Imch,Imsig,Imeqp0
    enddo
    read(20,*)
  enddo
  
  do jj=1,nval
    do ik=1,nkpts
      write(21,*) Reeqp0(ik,jj)
    enddo
  enddo
  
  call close_file(20)
  call close_file(21)
  
  call open_file(unit=30,file='spectrum.dat',status='old',form='formatted')
  
  freqsm =1d6
  
  do ik=1,nkpts
    do iband=1,nval
      read(30,*)
      do iw=1,nfreq
        if(ik .eq. 1 .and. iband .eq. 1) then
          read(30,202) kx,ky,kz,iispin,iiband,iiw,freq(iw),Resigma(iw,iband,ik),&
            Imsigma(iw,iband,ik),Imsigma2(iw,iband,ik) 
          if(abs(freq(iw)) .lt. freqsm) then
            freqsm = abs(freq(iw)) 
            nfreq_zero = iw         
          endif
        else  
          read(30,202) kx,ky,kz,iispin,iiband,iiw,Ew,Resigma(iw,iband,ik),&
            Imsigma(iw,iband,ik),Imsigma2(iw,iband,ik)
        endif
      enddo
      read(30,*)
    enddo
  enddo
  
  write(*,*) nfreq_zero,freq(nfreq_zero)
  
  call close_file(30)
  
  do it=1,ntime
    time_array(it) = -timerange + (it-1)*timeint
  enddo
  
!Get the value on our frequency grid closest to the lda energy for each kpoint
!and band. We will need this in evaluating sigma and its derivatives evaluated
!at the lda energy (Aryasetiawan) 

  do ik=1,nkpts
    do iband=1,nval
      ind=0
      diff = 0.0d0
      diff_small = 1d6
      do iw=1,nfreq
        diff=lda(ik,iband)-freq(iw)
        diff = abs(diff)
        if(diff .le. diff_small) then
          ind=iw
          diff_small = diff
        endif
      enddo
      flda_ind(ik,iband) = ind
    enddo
  enddo
  
!calculate numerical derivative of Sigma for use in cumulant expansion
!expressions

!Derivatives for QP part
  do ik = 1,nkpts
    do iband = 1,nval
      do iw=1,nfreq-2
        dResigma(iw+1,iband,ik) = (Resigma(iw+2,iband,ik)-&
          Resigma(iw,iband,ik))/(2*freq_interv)
        dImsigma(iw+1,iband,ik) = (Imsigma2(iw+2,iband,ik)-&
          Imsigma2(iw,iband,ik))/(2*freq_interv)
      enddo
    enddo
  enddo
  
!Derivatives for satellite part
  do ik = 1,nkpts
    do iband = 1,nval
      do iw=1,nfreq-2
        daImsigma(iw+1,iband,ik) = (abs(Imsigma2(iw+2,iband,ik))- &
          abs(Imsigma2(iw,iband,ik)))/(2*freq_interv)
!      write(*,*) (iw+1),1d0/PI_D*daImsigma(iw+1,iband,ik)
      enddo
    enddo
  enddo

  do ik = 1,nkpts
    do iband = 1,nval
      do iw=1,nfreq-4
        d2Imsigma(iw+2,iband,ik) = (daImsigma(iw+3,iband,ik)- &
          daImsigma(iw+1,iband,ik))/(2*freq_interv)
      enddo
    enddo
  enddo
  
  do ik = 1,nkpts
    do iband = 1,nval
      do iw=1,nfreq-6
        d3Imsigma(iw+3,iband,ik) = (d2Imsigma(iw+4,iband,ik)- &
          d2Imsigma(iw+3,iband,ik))/(2*freq_interv)
      enddo
    enddo
  enddo
  
!find the linear approximation to the self-energy correction evaluated at
!the QP energy

  do ik =1,nkpts
    do iband=1,nval
      Z_qp = CMPLX(dResigma(flda_ind(ik,iband),iband,ik),dImsigma(flda_ind(ik,iband),iband,ik))/ &
        (1-CMPLX(dResigma(flda_ind(ik,iband),iband,ik),dImsigma(flda_ind(ik,iband),iband,ik)))
      Reeqp1(ik,iband) = Reeqp0(ik,iband) + real(Z_qp)*(Reeqp0(ik,iband)-lda(ik,iband))
      write(*,*) ik,iband,Reeqp1(ik,iband),Reeqp0(ik,iband)
    enddo
  enddo
  
!Get the value on our frequency grid closest to the QP energy for each kpoint
!and band. We will need this in evaluating sigma and its derivatives evaluated
!at the QP energy (Hedin Ambladh). 

  do ik=1,nkpts
    do iband=1,nval
      ind=0
      diff = 0.0d0
      diff_small = 1d6
      do iw=1,nfreq
        diff=Reeqp1(ik,iband)-freq(iw)
        diff = abs(diff)
        if(diff .le. diff_small) then
          ind=iw
          diff_small = diff
        endif
      enddo
      fqp_ind(ik,iband) = ind
    enddo
  enddo

!Construct the function for the limit when w->0 of Cs(w)

  do ik = 1,nkpts
    do iband = 1,nval
      do iw=1,nfreq
        if((iw .gt. 3) .and. (iw .lt. nfreq-3)) then
          func_lda(iw,iband,ik) = d2Imsigma(flda_ind(ik,iband),iband,ik)/2 !+ freq(iw)*&
!                                d3Imsigma(flda_ind(ik,iband),iband,ik)/6
          func_qp(iw,iband,ik) =  d2Imsigma(fqp_ind(ik,iband),iband,ik)/2 !+ freq(iw)*&
!                                d3Imsigma(fqp_ind(ik,iband),iband,ik)/6
        else
          func_lda(iw,iband,ik) = 0.0d0 !We have the derivatives everywhere but the endpoints, so we set
          func_qp(iw,iband,ik) = 0.0d0  !the function to zero here. This works because the function will
          !be multiplied by a gaussian far from it`s center and the product 
          !will hence be zero anyhow
        endif
      enddo
    enddo
  enddo


!Calculate the QP and satellite spectral functions, evaluated at LDA (i_en=1)
! and QP (i_en = 2) energies

  do i_en=1,2 
!Here`s Cs(w), or As(w)to first order
!We shift the frequency grid here so that we can use all the data we have 
!for sigma in the calculations. Corresponding corrections are made in the 
!convolution in the first order treatment, and in the fourier transform
!to time in the full treatment
    do iband=1,nval
      do ik=1,nkpts
        if(i_en .eq. 1) then
          beta = (1d0/PI_D)*abs(Imsigma2(flda_ind(ik,iband),iband,ik))
          betap = (1d0/PI_D)*daImsigma(flda_ind(ik,iband),iband,ik)
!        write(*,*) "beta",beta,betap
        else
          beta = (1d0/PI_D)*abs(Imsigma2(fqp_ind(ik,iband),iband,ik))
          betap = (1d0/PI_D)*daImsigma(fqp_ind(ik,iband),iband,ik)
        endif
        do iw=1,nfreq
          
          if(i_en .eq. 1) then   
            ff=freq(iw)-freq(flda_ind(ik,iband))
            if(ff+2*lda(ik,iband) .lt. E_f) then
              num=((1d0/PI_D)*abs(Imsigma2(iw,iband,ik))-beta-ff*betap)
              Cs_band_k(iw,iband,ik) = num/ff**2                                                             
            else
              num=(-beta-ff*betap)
              Cs_band_k(iw,iband,ik) = num/ff**2
            endif
            
          else
            ff=freq(iw)-freq(fqp_ind(ik,iband))
            if(ff+2*Reeqp1(ik,iband) .lt. E_f) then
              num=((1d0/PI_D)*abs(Imsigma2(iw,iband,ik))-beta-ff*betap)
              Cs_band_k(iw,iband,ik) = num/ff**2                          
            else
              num=(-beta-ff*betap)
              Cs_band_k(iw,iband,ik) = num/ff**2
            endif
            !         write(*,*) ff,real(Cs_band_k(iw,iband,ik))          
            
          endif
        enddo
!! The next 5 lines are a bit of a hack to deal with the w = 0 divergence
!! It just connects two points separated by 2*diverge_int with a straight
!! line
        i_shift=floor(freq(flda_ind(ik,iband))/freq_interv)
        write(*,*) i_shift,freq(flda_ind(ik,iband))
        slope=(real(Cs_band_k(nfreq_zero+i_shift+diverge_int,iband,ik)) &
          -real(Cs_band_k(nfreq_zero+i_shift-diverge_int,iband,ik)))/real(2*diverge_int)
        do ii=1,1+2*diverge_int
          Cs_band_k(nfreq_zero+i_shift-(diverge_int+1)+ii,iband,ik) = &
            Cs_band_k(nfreq_zero+i_shift-diverge_int,iband,ik)+slope*real((ii-1))
        enddo
        Cs_band(:,iband) = Cs_band(:,iband) + Cs_band_k(:,iband,ik)*wts(ik)
      enddo
      Cs(:) = Cs(:) + Cs_band(:,iband)
    enddo
    
!Here`s the QP spectral function
    alpha = 0.0d0
    gamma = 0.0d0
    gamma_big = 0.0d0
    Eqp = 0.0d0
    do iband=1,nval
      do ik=1,nkpts
        if(i_en .eq. 1) then
          alpha = dImsigma(flda_ind(ik,iband),iband,ik)
          gamma = -dResigma(flda_ind(ik,iband),iband,ik)
          gamma_big = Imsigma2(flda_ind(ik,iband),iband,ik) 
          Eqp = lda(ik,iband)
        else
          alpha = dImsigma(fqp_ind(ik,iband),iband,ik)
          gamma = -dResigma(fqp_ind(ik,iband),iband,ik)
          gamma_big = Imsigma2(fqp_ind(ik,iband),iband,ik) 
          Eqp = Reeqp1(ik,iband)
        endif
        do iw=1,nfreq
          Aqp_band_k(iw,iband,ik) =(1d0/PI_D)*exp(-gamma)*(abs(gamma_big)*cos(alpha)-&
            (freq(iw)-Eqp)*sin(alpha))/((freq(iw)-Eqp)**2+abs(gamma_big)**2)
        enddo
        Aqp_band(:,iband) = Aqp_band(:,iband) + Aqp_band_k(:,iband,ik)*wts(ik)
      enddo
      Aqp(:) = Aqp(:) + Aqp_band(:,iband)
    enddo
    
    if(i_method .eq. 1) then  ! 1st order approximation

!Convolute the QP and satellite spectral functions to get the 
!total spectrum
      do iband=1,nval
        do ik=1,nkpts
          if(i_en .eq. 1) then
            alpha = dImsigma(flda_ind(ik,iband),iband,ik)
            gamma = -dResigma(flda_ind(ik,iband),iband,ik)
            gamma_big = Imsigma2(flda_ind(ik,iband),iband,ik) 
            Eqp = lda(ik,iband)
          else
            alpha = dImsigma(fqp_ind(ik,iband),iband,ik)
            gamma = -dResigma(fqp_ind(ik,iband),iband,ik)
            gamma_big = Imsigma2(fqp_ind(ik,iband),iband,ik) 
            Eqp = Reeqp1(ik,iband)
          endif
          do iw=1,nfreq_F+nfreq_Fadd
            do ii=1,nfreq_F+nfreq_Fadd
              Aconv_band_k(iw,iband,ik)=Aconv_band_k(iw,iband,ik) + freq_interv* &
                (1d0/PI_D)*exp(-gamma)*(abs(gamma_big)*cos(alpha)-(freq_interv*(iw-ii))*sin(alpha))/ &
                ((freq_interv*(iw-ii))**2+abs(gamma_big)**2)*Cs_band_k(ii,iband,ik)
            enddo
            Atot_band_k(iw,iband,ik)=Aqp_band_k(iw,iband,ik)+&
              (8.0d0/3.0d0)*Aconv_band_k(iw,iband,ik)
          enddo
          Atot_band(:,iband)= Atot_band(:,iband) + Atot_band_k(:,iband,ik)*wts(ik)*matr_el(ik,iband)
        enddo
        Atot(:) = Atot(:) + Atot_band(:,iband)
      enddo
      
    else     !Full exponentiation
!Fourier transform Cs(w) to time domain
      do iband=1,nval
        do ik=1,nkpts
          do it=1,ntime
            do iw=1,nfreq
              Cs_band_kt(it,iband,ik) = Cs_band_kt(it,iband,ik) + exp(-imag*freq(iw)*time_array(it))* &
                Cs_band_k(iw,iband,ik)*freq_interv
            enddo
            if(i_en .eq. 1) then
              Cs_band_kt(it,iband,ik) = Cs_band_kt(it,iband,ik)*exp(imag*lda(ik,iband)*time_array(it))
            else
              Cs_band_kt(it,iband,ik) = Cs_band_kt(it,iband,ik)*exp(imag*Reeqp1(ik,iband)*time_array(it))
            endif
 !         write(*,*) time_array(it),real(Cs_band_kt(it,iband,ik))
          enddo
          Cs_bandt(:,iband) = Cs_bandt(:,iband) + Cs_band_kt(:,iband,ik)*wts(ik) 
        enddo
        Cst(:) = Cst(:) + Cs_bandt(:,iband)
      enddo
      
      do it=1,ntime
!  write(*,*) time_array(it),real(exp(Cs_band_kt(it,1,1)))
      enddo

!Construct Cqp(t)
      do iband=1,nval
        do ik=1,nkpts
          if(i_en .eq. 1) then
            alpha = dImsigma(flda_ind(ik,iband),iband,ik)
            gamma = -dResigma(flda_ind(ik,iband),iband,ik)
            eta = abs(Imsigma2(flda_ind(ik,iband),iband,ik))
            Eqp = Reeqp0(ik,iband)
          else
            alpha = dImsigma(fqp_ind(ik,iband),iband,ik)
            gamma = -dResigma(fqp_ind(ik,iband),iband,ik)
            eta = abs(Imsigma2(fqp_ind(ik,iband),iband,ik))
            Eqp = Reeqp1(ik,iband)
          endif
!        write(*,*) "eta",eta,gamma,Eqp,alpha,flda_ind(ik,iband)
!       write(*,*) "alpha",alpha,Eqp
          do it=1,ntime
            Cqp_band_kt(it,iband,ik) = -gamma+imag*alpha*sign(1.0d0,time_array(it)) &
              -imag*Eqp*time_array(it)-abs(eta)*abs(time_array(it))
 !         write(*,*) time_array(it),real(Cqp_band_kt(it,iband,ik)),aimag(Cqp_band_kt(it,iband,ik))
          enddo
        enddo
      enddo

!do it=1,ntime
!  write(*,*) time_array(it),real(exp(Cs_band_kt(it,1,1)+Cqp_band_kt(it,1,1)))
!enddo
       
!Fourier transform exp(Cqp(t)+Cs(t)) back to frequency space to get spectral
!function. Note that we get Atot_band_k on the extended grid to see second
!plasmon (we`ll see if we actually see it)
      do iband=1,nval
        do ik=1,nkpts
          do iw=1,nfreq
            do it=1,ntime
              Atot_band_k(iw,iband,ik) = Atot_band_k(iw,iband,ik)+&
                exp(imag*freq(iw)*time_array(it))*exp(Cs_band_kt(it,iband,ik)+Cqp_band_kt(it,iband,ik))*timeint
              if(iw .eq. 1 .and. iband .eq. 2) then
!              write(*,*) time_array(it),real(exp(Cs_band_kt(it,iband,ik)+Cqp_band_kt(it,iband,ik)))
              endif
            enddo
            Atot_band_k(iw,iband,ik) = Atot_band_k(iw,iband,ik)/(2d0*PI_D)
 !         write(*,*) freq(iw),real(Atot_band_k(iw,iband,ik)),aimag(Atot_band_k(iw,iband,ik))
          enddo
 !       write(*,*)
          Atot_band(:,iband) = Atot_band(:,iband) + Atot_band_k(:,iband,ik)*wts(ik)*matr_el(ik,iband)
        enddo
        Atot(:) = Atot(:) + Atot_band(:,iband)
      enddo
    endif
  

    iband = 0
    freq = freq - E_f


!Write out the QP, satellite, and total spectral functions

    if(i_en .eq. 1) then
      write(Aname,102) "A_e_qp_bnd",iband,"_k00"
    else
      write(Aname,102) "A_E_qp_bnd",iband,"_k00"
    endif
    call open_file(unit=30,file=Aname,status='replace',form='formatted')  
    do iw=1,nfreq
      write(30,100) freq(iw),real(Aqp(iw))
    enddo
    call close_file(30)
    
    do iband=1,nval
      if(i_en .eq. 1) then
        write(Aname,102) "A_e_qp_bnd",iband,"_k00"
      else
        write(Aname,102) "A_E_qp_bnd",iband,"_k00"
      endif
      call open_file(unit=30,file=Aname,status='replace',form='formatted')  
      do iw=1,nfreq
        write(30,100) freq(iw),real(Aqp_band(iw,iband))
      enddo
      call close_file(30)
      do ik=1,nkpts
        if(i_en .eq. 1) then
          write(Aname,103) "A_e_qp_bnd",iband,"_k",ik
        else
          write(Aname,103) "A_E_qp_bnd",iband,"_k",ik
        endif
        call open_file(unit=30,file=Aname,status='replace',form='formatted')  
        do iw=1,nfreq
          write(30,100) freq(iw),real(Aqp_band_k(iw,iband,ik))
        enddo
        call close_file(30)
      enddo
    enddo
    
    iband=0
    
    if(i_en .eq. 1) then
      write(Aname,102) "A_e_sa_bnd",iband,"_k00"
    else
      write(Aname,102) "A_E_sa_bnd",iband,"_k00"
    endif
    call open_file(unit=30,file=Aname,status='replace',form='formatted')  
    do iw=1,nfreq
      write(30,100) freq(iw),real(Cs(iw))
    enddo
    call close_file(30)
    
    do iband=1,nval
      if(i_en .eq. 1) then
        write(Aname,102) "A_e_sa_bnd",iband,"_k00"
      else
        write(Aname,102) "A_E_sa_bnd",iband,"_k00"
      endif
      call open_file(unit=30,file=Aname,status='replace',form='formatted')  
      do iw=1,nfreq
        write(30,100) freq(iw),real(Cs_band(iw,iband))
      enddo
      call close_file(30)
      do ik=1,nkpts
        if(i_en .eq. 1) then
          write(Aname,103) "A_e_sa_bnd",iband,"_k",ik
        else
          write(Aname,103) "A_E_sa_bnd",iband,"_k",ik
        endif
        call open_file(unit=30,file=Aname,status='replace',form='formatted')  
        do iw=1,nfreq
          write(30,100) freq(iw),real(Cs_band_k(iw,iband,ik))
        enddo
        call close_file(30)
      enddo
    enddo
    
    iband=0
    
    if(i_en .eq. 1) then
      write(Aname,102) "A_e_to_bnd",iband,"_k00"
    else
      write(Aname,102) "A_E_to_bnd",iband,"_k00"
    endif
    call open_file(unit=30,file=Aname,status='replace',form='formatted')  
    do iw=1,nfreq_F+nfreq_Fadd
      write(30,100) freq(iw),real(Atot(iw))
    enddo
    call close_file(30)
    
    do iband=1,nval
      if(i_en .eq. 1) then
        write(Aname,102) "A_e_to_bnd",iband,"_k00"
      else
        write(Aname,102) "A_E_to_bnd",iband,"_k00"
      endif
      call open_file(unit=30,file=Aname,status='replace',form='formatted')  
      do iw=1,nfreq_F+nfreq_Fadd
        write(30,100) freq(iw),real(Atot_band(iw,iband))
      enddo
      call close_file(30)
      do ik=1,nkpts
        if(i_en .eq. 1) then
          write(Aname,103) "A_e_to_bnd",iband,"_k",ik
        else
          write(Aname,103) "A_E_to_bnd",iband,"_k",ik
        endif
        call open_file(unit=30,file=Aname,status='replace',form='formatted')  
        do iw=1,nfreq_F+nfreq_Fadd
          write(30,100) freq(iw),real(Atot_band_k(iw,iband,ik))
        enddo
        call close_file(30)
      enddo
    enddo

!Smoothing the data here

    do ii=1,nfreq_ext
      freq_ext(ii) = -90.0d0 + freq_interv*(ii-1)
      Gauss(ii) = 1/sqrt(2 * PI_D * 0.25)*exp(-freq_ext(ii)**2 /(2*0.25))
    enddo
    
    do ii=1,nfreq_F+nfreq_Fadd
      do jj=1,nfreq_F+nfreq_Fadd
        convol1(ii) = convol1(ii) + freq_interv*Atot_band(jj,1)*Gauss(451+ii-jj)
        convol2(ii) = convol2(ii) + freq_interv*Atot_band(jj,2)*Gauss(451+ii-jj)
        convol3(ii) = convol3(ii) + freq_interv*Atot_band(jj,3)*Gauss(451+ii-jj)
        convol4(ii) = convol4(ii) + freq_interv*Atot_band(jj,4)*Gauss(451+ii-jj)
      enddo
    enddo
    
    convolt = convol1 + convol2 + convol3 + convol4
    
    do ii=1,nfreq_F+nfreq_Fadd
      do jj=ii,nfreq_F+nfreq_Fadd
        convolint(ii) = convolint(ii) + freq_interv*convolt(jj)
      enddo
    enddo
    
    do ii=1,nfreq_F+nfreq_Fadd
      Awb(ii) = convolt(ii) + convolint(ii)/intfact
    enddo
    
    convolt=0.0d0
    convolint=0.0d0
    
    convolt = mefact1*convol1 + mefact2*convol2 + convol3 + convol4
    
    do ii=1,nfreq_F+nfreq_Fadd
      do jj=ii,nfreq_F+nfreq_Fadd
        convolint(ii) = convolint(ii) + freq_interv*convolt(jj)
      enddo
    enddo
    
    do ii=1,nfreq_F+nfreq_Fadd
      Awbwme(ii) = convolt(ii) + convolint(ii)/intfact
    enddo
    
    convolt = convol1 + convol2 + convol3 + convol4
    if(i_en .eq. 1) then
      call open_file(unit=30,file='I_e_curr',status='replace',form='formatted')  
      call open_file(unit=40,file='I_e_wb',status='replace',form='formatted')
      call open_file(unit=50,file='I_e_wbwme',status='replace',form='formatted')
    else
      call open_file(unit=30,file='I_E_curr',status='replace',form='formatted')  
      call open_file(unit=40,file='I_E_wb',status='replace',form='formatted')
      call open_file(unit=50,file='I_E_wbwme',status='replace',form='formatted')
    endif
    do ii=1,nfreq_F+nfreq_Fadd
      write(30,100) freq(ii),convolt(ii)
    enddo
    do ii=1,nfreq_F+nfreq_Fadd
      write(40,100) freq(ii),real(Awb(ii))
    enddo
    do ii=1,nfreq_F+nfreq_Fadd
      write(50,100) freq(ii),real(Awbwme(ii))
    enddo
    
    call close_file(30)
    call close_file(40)
    call close_file(50)
    
    if(i_en .eq. 1) then
      freq = freq + E_f
    endif
  enddo
  
  deallocate(Reeqp0)
  deallocate(Reeqp1)
  deallocate(Gauss)
  deallocate(wts)
  deallocate(lda)
  deallocate(fqp_ind)
  deallocate(flda_ind)
  deallocate(vxc)
  deallocate(Cs)
  deallocate(Cs_band)
  deallocate(Cs_band_k)
  deallocate(time_array)
  deallocate(Cst)
  deallocate(Cs_bandt)
  deallocate(Cs_band_kt)
  deallocate(Cqpt)
  deallocate(Cqp_bandt)
  deallocate(Cqp_band_kt)
  deallocate(Aqp)
  deallocate(Aqp_band)
  deallocate(Aqp_band_k)
  deallocate(Atot)
  deallocate(Atot_band)
  deallocate(Atot_band_k)
  deallocate(Aconv_band_k)
  deallocate(convolint)
  deallocate(Awb)
  deallocate(Awbwme)
  deallocate(convolt)
  deallocate(convol1)
  deallocate(convol2)
  deallocate(convol3)
  deallocate(convol4)
  deallocate(freq_ext)
  deallocate(freq)
  deallocate(Resigma)
  deallocate(Imsigma)
  deallocate(Imsigma2)
  deallocate(dResigma)
  deallocate(dImsigma)
  deallocate(daImsigma)
  deallocate(d2Imsigma)
  deallocate(d3Imsigma)
  deallocate(func_qp)
  deallocate(func_lda)
  
  stop
  
100 format(f6.2,7x,f12.8)
102 format(a10,i1,a4)
103 format(a10,i1,a2,i2.2)
  
200 format(i4,8f12.6)
201 format(40x,3f12.6,12x,f12.6)
202 format(3f12.5,2x,3i4,2x,4f12.5)

end program kptsummedspectralfcn
