
subroutine A51(j1,j2,j3,j4,j5,za,zb)
  use, intrinsic :: iso_c_binding
  use constants_mod
  use epinv_mod
  use epinv2_mod
  use mxpart_mod
  use scale_mod
  use sprods_com_mod

  implicit none

  integer :: j1,j2,j3,j4,j5
  complex(dp) :: za(mxpart,mxpart),zb(mxpart,mxpart)
  complex(dp) :: Vcc,Fcc,Vsc,Fsc,l12,l23,L0,L1,Lsm1,A5lom
  complex(dp) :: lnrat,zab2
  real(dp) :: s123

  zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

  s123=s(j1,j2)+s(j2,j3)+s(j3,j1)
  A5lom=-za(j3,j4)*zab2(j3,j1,j2,j5)/(za(j1,j2)*za(j2,j3)*s123)

  l12=lnrat(musq,-s(j1,j2))
  l23=lnrat(musq,-s(j2,j3))
  Vcc=
     & -(epinv*epinv2+epinv*l12+half*pow(l12,2))
     & -(epinv*epinv2+epinv*l23+half*pow(l23,2))
     & -two*(epinv+l23)-four

  Fcc=zab2(j3,j1,j2,j5)/(za(j1,j2)*za(j2,j3)*s123)
     & *(za(j3,j4)*Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)
     & +two*za(j3,j1)*zab2(j4,j2,j3,j1)
     &   *L0(-s(j2,j3),-s123)/s123)

  Vsc =half*(epinv+l23)+one
  Fsc =-za(j3,j4)*za(j3,j1)*zb(j1,j5)
     & /(za(j1,j2)*za(j2,j3))*L0(-s(j2,j3),-s123)/s123
     & +half*za(j3,j1)**2*zb(j1,j5)*zab2(j4,j2,j3,j1)
     & /(za(j1,j2)*za(j2,j3))*L1(-s(j2,j3),-s123)/pow(s123,2)

  A51=(Vcc+Vsc)*A5Lom+Fcc+Fsc

  return
end subroutine A51
