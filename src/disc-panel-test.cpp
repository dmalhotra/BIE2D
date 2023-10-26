#include <bie2d.hpp>
using namespace sctl;

int main()
{
    using Real = double;
    using Long = long long;
    constexpr Integer Order = 16;
    
    Long Ndisc = 3;
    Real r = 0.5;
    Vector<Real> xy(2*Ndisc);
    
    xy[0] = -0.51;
    xy[1] = 0;

    xy[2] = 0.51;
    xy[3] = 0;

    xy[4] = 0;
    xy[5] = 1;

    DiscPanelLst<Real, Order> discs;
    discs.Init(xy, r);

    for (int i = 0; i < Ndisc; i++) {
        auto X = discs.SurfCoord(i);
        for (int j = 0; j < X.Dim(); j += discs.CoordDim()) {
            std::cout << X[j] << " " << X[j+1] << std::endl;
        }
    }

    auto nearlist = discs.GetNearList();
    for (auto& pair : nearlist) {
        std::cout << "(j,k) = " << pair.disc_idx0 << " " << pair.disc_idx1 << std::endl;
        std::cout << "j range = " << pair.panel_idx_range0[0] << " " << pair.panel_idx_range0[1] << std::endl;
        std::cout << "k range = " << pair.panel_idx_range1[0] << " " << pair.panel_idx_range1[1] << std::endl;
    }
}
