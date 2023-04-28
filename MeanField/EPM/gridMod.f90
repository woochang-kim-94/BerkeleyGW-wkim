!*************************************************************************
Module gridMod
!*************************************************************************

  use SysParams, only : double

  Type GridT
     complex(double), pointer :: V(:,:,:)
     complex(double), pointer :: Psi(:,:,:)
     complex(double), pointer :: VPsi(:,:,:)
     complex(double), pointer :: VPsiG(:,:,:)
     complex(double), pointer :: Work(:,:,:)
     integer                  :: NumGVec,k
     integer,         pointer :: GVecIndex(:,:)
     integer,         pointer :: GInvIndex(:,:,:)
     integer,         pointer :: BasisIndex(:,:)
     integer,    dimension(3) :: Size,Length
     integer(double)          :: PlanF,PlanB
     real(double)             :: fac
  end Type GridT

end Module gridMod

