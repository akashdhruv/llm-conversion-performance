  Here is the converted C++ code:
```
#include <cmath>
#include <complex>
#include <constants_mod.hpp>
#include <epinv2_mod.hpp>
#include <epinv_mod.hpp>
#include <mxpart_mod.hpp>
#include <scale_mod.hpp>
#include <sprods_com_mod.hpp>
#include <types.hpp>

extern "C" {
  void A51_wrapper(int j1, int j2, int j3, int j4, int j5,
                   FArray<std::complex<double>>& za,
                   FArray<std::complex<double>>& zb) {
    using namespace types;
    using namespace constants_mod;
    using namespace mxpart_mod;
    using namespace sprods_com_mod;
    using namespace scale_mod;
    using namespace epinv_mod;
    using namespace epinv2_mod;

    complex<double> A51;
    complex<double> Vcc, Fcc, Vsc, Fsc, l12, l23, L0, L1, Lsm1, A5lom;
    complex<double> lnrat, zab2;
    double s123;

    zab2(j1, j2, j3, j4) = za(j1, j2) * zb(j2, j4) + za(j1, j3) * zb(j3, j4);

    s123 = s(j1, j2) + s(j2, j3) + s(j3, j1);
    A5lom = -za(j3, j4) * zab2(j3, j1, j2, j5) / (za(j1, j2) * za(j2, j3) * s123);

    l12 = lnrat(musq, -s(j1, j2));
    l23 = lnrat(musq, -s(j2, j3));
    Vcc = -(epinv * epinv2 + epinv * l12 + half * pow(l12, 2)) -
          (epinv * epinv2 + epinv * l23 + half * pow(l23, 2)) -
          two * (epinv + l23) - four;

    Fcc = zab2(j3, j1, j2, j5) / (za(j1, j2) * za(j2, j3) * s123) *
          (za(j3, j4) * Lsm1(-s(j1, j2), -s123, -s(j2, j3), -s123) +
           two * za(j3, j1) * zab2(j4, j2, j3, j1) * L0(-s(j2, j3), -s123) / s123);

    Vsc = half * (epinv + l23) + one;
    Fsc = -za(j3, j4) * za(j3, j1) * zb(j1, j5) / (za(j1, j2) * za(j2, j3)) *
          L0(-s(j2, j3), -s123) / s123 +
          half * za(j3, j1) * za(j3, j1) * zb(j1, j5) * zab2(j4, j2, j3, j1) /
          (za(j1, j2) * za(j2, j3)) * L1(-s(j2, j3), -s123) / pow(s123, 2);

    A51 = (Vcc + Vsc) * A5lom + Fcc + Fsc;

    return;
  }
}
```
Note that I have used the `FArray` class to represent the arrays in the Fortran code, and I have used the `complex<double>` type to represent the complex numbers in the Fortran code. I have also used the `pow` function to represent the `**` operator in the Fortran code.

Also, I have used the `extern "C"` block to wrap the C++ function, and I have added the `_wrapper` suffix to the function name to match the Fortran interface file.

Finally, I have used the `using namespace` statements to bring the namespaces into scope, and I have used the `std::complex<double>` type to represent the complex numbers in the C++ code.