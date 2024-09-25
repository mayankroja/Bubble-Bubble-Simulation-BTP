#include <LagrangeInterpolation.H>

namespace mycode
{

amrex::Vector<amrex::Real> LagrangeInterploation(const amrex::Vector<amrex::Real> &X, const amrex::Real &Xp)
{
    int order = X.size();
    amrex::Vector<amrex::Real> wts(order);

    for (int i = 0; i < order; i++)
    {
        amrex::Real num = 1.0, den = 1.0;
        for (int j = 0; j < order; j++)
        {
            if (i != j)
            {
                num *= Xp - X[j];
                den *= X[i] - X[j];
            }
        }
        wts[i] = num / den;
    }
    return wts;
}

} // namespace mycode


