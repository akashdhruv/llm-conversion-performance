
#include <cmath>
#include <complex>
#include <constants_mod.hpp>
#include <epinv2_mod.hpp>
#include <epinv_mod.hpp>
#include <mxpart_mod.hpp>
#include <scale_mod.hpp>
#include <sprods_com_mod.hpp>
#include <types.hpp>
#include <FArray.hpp>

std::complex<double> A51(int j1, int j2, int j3, int j4, int j5, FArray2D<std::complex<double>>& za, FArray2D<std::complex<double>>& zb) {
    using namespace constants_mod;
    using namespace mxpart_mod;
    using namespace sprods_com_mod;
    using namespace scale_mod;
    using namespace epinv_mod;
    using namespace epinv2_mod;

    std::complex<double> A51_return;
    
    FArray2D<std::complex<double>> zab2(j1, j2, j3, j4, za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4));
    
    double s123 = s(j1,j2) + s(j2,j3) + s(j3,j1);
    std::complex<double> A5lom = -za(j3,j4) * zab2(j3,j1,j2,j5) / (za(j1,j2) * za(j2,j3) * s123);
    
    std::complex<double> l12 = lnrat(musq, -s(j1,j2));
    std::complex<double> l23 = lnrat(musq, -s(j2,j3));
    std::complex<double> Vcc = -(epinv*epinv2 + epinv*l12 + 0.5*l12*l12) 
                             - (epinv*epinv2 + epinv*l23 + 0.5*l23*l23)
                             - 2.0*(epinv + l23) - 4.0;
    
    std::complex<double> Fcc = zab2(j3,j1,j2,j5) / (za(j1,j2) * za(j2,j3) * s123) *
                               (za(j3,j4) * Lsm1(-s(j1,j2), -s123, -s(j2,j3), -s123)
                               + 2.0*za(j3,j1)*zab2(j4,j2,j3,j1) * L0(-s(j2,j3), -s123) / s123);
    
    std::complex<double> Vsc = 0.5*(epinv + l23) + 1.0;
    std::complex<double> Fsc = -za(j3,j4) * za(j3,j1) * zb(j1,j5) / (za(j1,j2) * za(j2,j3)) * L0(-s(j2,j3), -s123) / s123
                             + 0.5*za(j3,j1)*za(j3,j1) * zb(j1,j5) * zab2(j4,j2,j3,j1) / (za(j1,j2) * za(j2,j3)) * L1(-s(j2,j3), -s123) / (s123 * s123);
    
    A51_return = (Vcc + Vsc) * A5lom + Fcc + Fsc;
    
    return A51_return;
}

extern "C" {
    void A51_wrapper(int j1, int j2, int j3, int j4, int j5, std::complex<double>* fza, std::complex<double>* fzb, std::complex<double>* result) {
        using namespace types;
        
        FArray2D<std::complex<double>> za(fza, mxpart, mxpart);
        FArray2D<std::complex<double>> zb(fzb, mxpart, mxpart);
        *result = A51(j1, j2, j3, j4, j5, za, zb);
    }
}
