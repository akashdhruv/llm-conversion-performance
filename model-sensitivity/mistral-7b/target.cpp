
#include <cmath>
#include <complex>
#include <constants_mod.hpp>
#include <epinv2_mod.hpp>
#include <epinv_mod.hpp>
#include <mxpart_mod.hpp>
#include <scale_mod.hpp>
#include <sprods_com_mod.hpp>
#include <types.hpp>

void A51(int j1, int j2, int j3, int j4, int j5, FArray<std::complex<double>>& za, FArray<std::complex<double>>& zb) {
  using namespace constants_mod;
  using namespace epinv_mod;
  using namespace epinv2_mod;
  using namespace mxpart_mod;
  using namespace scale_mod;
  using namespace sprods_com_mod;

  complex<double> zab2(j1,j2,j3,j4) = za(j1,j2)*zb(j2,j4) + za(j1,j3)*zb(j3,j4);
  double s123 = s(j1,j2) + s(j2,j3) + s(j3,j1);
  complex<double> A5lom = -za(j3,j4)*zab2(j3,j1,j2,j5)/(za(j1,j2)*za(j2,j3)*s123);

  complex<double> l12 = lnrat(musq,-s(j1,j2));
  complex<double> l23 = lnrat(musq,-s(j2,j3));
  complex<double> Vcc =
    -(epinv*epinv2+epinv*l12+half*pow(l12,2))
    -(epinv*epinv2+epinv*l23+half*pow(l23,2))
    -two*(epinv+l23)-four;

  complex<double> Fcc = zab2(j3,j1,j2,j5)/(za(j1,j2)*za(j2,j3)*s123);
  complex<double> Fsc = -za(j3,j4)*za(j3,j1)*zb(j1,j5)
    /(za(j1,j2)*za(j2,j3))*L0(-s(j2,j3),-s123)/s123
    +half*za(j3,j1)**2*zb(j1,j5)*zab2(j4,j2,j3,j1)
    /(za(j1,j2)*za(j2,j3))*L1(-s(j2,j3),-s123)/pow(s123,2);

  complex<double> Vsc = half*(epinv+l23)+one;

  A51 = (Vcc+Vsc)*A5Lom + Fcc + Fsc;
}

extern "C" {
  void A51_wrapper(int j1, int j2, int j3, int j4, int j5, double* za, double* zb) {
    using namespace mxpart_mod;

    FArray<std::complex<double>> za_farray(za, mxpart, mxpart);
    FArray<std::complex<double>> zb_farray(zb, mxpart, mxpart);

    A51(j1, j2, j3, j4, j5, za_farray, zb_farray);
  }
}
