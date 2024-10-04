
function A51(j1,j2,j3,j4,j5,za,zb)
  use, intrinsic :: iso_c_binding
  use constants_mod
  use mxpart_mod
  use sprods_com_mod
  use scale_mod
  use epinv_mod
  use epinv2_mod
  
  implicit none

  integer(c_int), value :: j1, j2, j3, j4, j5
  complex(c_double_complex) :: za(mxpart,mxpart), zb(mxpart,mxpart)
  complex(c_double_complex) :: A51

  interface
    function A51_wrapper(j1,j2,j3,j4,j5,za,zb) bind(C, name="A51_wrapper")
      import :: c_int, c_double_complex
      integer(c_int), value :: j1, j2, j3, j4, j5
      complex(c_double_complex), dimension(mxpart,mxpart) :: za, zb
      complex(c_double_complex) :: A51_wrapper
    end function A51_wrapper
  end interface

  A51 = A51_wrapper(j1, j2, j3, j4, j5, za, zb)
end function A51
