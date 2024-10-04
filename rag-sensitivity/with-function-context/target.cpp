
#include <complex>
#include <constants_mod.hpp>
#include <epinv2_mod.hpp>
#include <epinv_mod.hpp>
#include <mxpart_mod.hpp>
#include <scale_mod.hpp>
#include <sprods_com_mod.hpp>
#include <types.hpp>
#include <FArray.hpp>

using namespace types;
using namespace constants_mod;
using namespace mxpart_mod;
using namespace sprods_com_mod;
using namespace scale_mod;
using namespace epinv_mod;
using namespace epinv2_mod;

std::complex<double> A51(int j1, int j2, int j3, int j4, int j5, FArray2D<std::complex<double>>& za, FArray2D<std::complex<double>>& zb) {
    std::complex<double> Vcc, Fcc, Vsc, Fsc, l12, l23, A5lom;
    std::complex<double> zab2(mxpart, mxpart);
    double s123;

    zab2(j1,j2) = za(j1,j2)*zb(j2,j4) + za(j1,j3)*zb(j3,j4);

    s123 = s(j1,j2) + s(j2,j3) + s(j3,j1);
    A5lom = -za(j3,j4) * zab2(j3,j1) / (za(j1,j2) * za(j2,j3) * s123);

    l12 = lnrat(musq,-s(j1,j2));
    l23 = lnrat(musq,-s(j2,j3));

    Vcc = -(epinv * epinv2 + epinv * l12 + half * std::pow(l12, 2))
         -(epinv * epinv2 + epinv * l23 + half * std::pow(l23, 2))
         -two * (epinv + l23) - four;

    Fcc = zab2(j3,j1) / (za(j1,j2) * za(j2,j3) * s123)
         * (za(j3,j4) * Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)
         + two * za(j3,j1) * zab2(j4,j2) * L0(-s(j2,j3),-s123) / s123);

    Vsc = half * (epinv + l23) + one;
    Fsc = -za(j3,j4) * za(j3,j1) * zb(j1,j5)
         / (za(j1,j2) * za(j2,j3)) * L0(-s(j2,j3),-s123) / s123
         + half * std::pow(za(j3,j1), 2) * zb(j1,j5) * zab2(j4,j2)
         / (za(j1,j2) * za(j2,j3)) * L1(-s(j2,j3),-s123) / std::pow(s123, 2);

    return (Vcc + Vsc) * A5lom + Fcc + Fsc;
}

extern "C" {
    std::complex<double> A51_wrapper(int j1, int j2, int j3, int j4, int j5, std::complex<double>* zaf, std::complex<double>* zbf) {
        FArray2D<std::complex<double>> za(zaf, mxpart, mxpart);
        FArray2D<std::complex<double>> zb(zbf, mxpart, mxpart);
        return A51(j1, j2, j3, j4, j5, za, zb);
    }
}
