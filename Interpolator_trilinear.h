#ifndef INTERPOLATOR_TRILINEAR_H
#define INTERPOLATOR_TRILINEAR_H

#include <array>
#include <vector>

class Interpolator_trilinear
{
private:
    /* data */
public:
    Interpolator_trilinear();
    ~Interpolator_trilinear();

    static std::array<double,8> GetCoefficients(const std::array<double, 8>& corner_values);
    
    static std::vector<double> InterpolateWCoeffs(const std::array<double, 8>& coeffs, 
    const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z);
    
    static std::vector<double> Interpolate(const std::array<double, 8>& cube_corner_values, 
        const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z);
};


#endif

