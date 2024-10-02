function A51(j1,j2,j3,j4,j5,za,zb) result(res)
  use, intrinsic :: iso_c_binding
  use types
  use constants_mod
  use mxpart_mod
  use sprods_com_mod
  use scale_mod
  use epinv_mod
  use epinv2_mod

  implicit none
  complex(c_double) :: res
  complex(c_double), dimension(mxpart,mxpart), intent(in) :: za, zb
  integer, intent(in) :: j1, j2, j3, j4, j5

  interface
    function A51_wrapper(j1,j2,j3,j4,j5,za,zb) bind(C, name="A51_wrapper")
      use, intrinsic :: iso_c_binding
      import :: mxpart
      complex(c_double) :: A51_wrapper
      integer(c_int), value :: j1, j2, j3, j4, j5
      complex(c_double), dimension(mxpart,mxpart), intent(in) :: za, zb
    end function A51_wrapper
  end interface

  res = A51_wrapper(j1,j2,j3,j4,j5,za,zb)
end function A51
