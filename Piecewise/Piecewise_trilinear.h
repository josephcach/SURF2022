#ifndef PIECEWISE_TRILINEAR_H
#define PIECEWISE_TRILINEAR_H

#include <array>
#include <vector>
#include <functional>

class Piecewise_trilinear
{
private:

public:
    Piecewise_trilinear() {}
    ~Piecewise_trilinear() {}

    static std::vector<std::vector<std::vector<double> > > GetCorners(const int N,
        const std::function<double(double, double, double)>& fun);

    static double Interpolate(const int N, const std::vector<std::vector<std::vector<double> > >& corners, 
        const double x, const double y, const double z );

};

#endif
