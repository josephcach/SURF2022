#include "Interpolator_trilinear.h"
#include <exception>


Interpolator_trilinear::Interpolator_trilinear(/* args */)
{
}

Interpolator_trilinear::~Interpolator_trilinear()
{
}


std::vector<double> Interpolator_trilinear::Interpolate(const std::array<double, 8>& cube_corner_values, 
        const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z)
{
    if (x.size() != y.size() || x.size() != z.size())
        throw std::invalid_argument("different sized x, y, z.");

    std::vector<double> result(x.size(), 0.0);
    double xd, yd, zd;
    double xinterp[4];
    for (long long iter = 0; iter < x.size(); iter++)
    {
        //for each point do trilinear interpolation using the given 8 cube_corner_values
         xd = (x[iter]+1)/2;
         yd = (y[iter]+1)/2;
         zd = (z[iter] +1)/2;

        for(int i=0; i<4;i++){
            xinterp[i] = cube_corner_values[i]*(1-xd)+cube_corner_values[i+4]*xd;
        }
        result[iter] =  (xinterp[0]*(1-yd)+xinterp[2]*yd)*(1-zd)+(xinterp[1]*(1-yd)+xinterp[3]*yd)*zd;;//write the result of the interpolation into this vector.
    }
    
    return result;
}

std::array<double,8> Interpolator_trilinear::GetCoefficients(const std::array<double, 8>& corner_values){
    std::array<double,8> coefs;
    std::array<std::array<int,8>,8> scalars = {{
        {-1,-1,-1,-1,-1,-1,-1,-1},
        {1,1,1,1,-1,-1,-1,-1},
        {1,1,-1,-1,1,1,-1,-1},
        {1,-1,1,-1,1,-1,1,-1},
        {-1,-1,1,1,1,1,-1,-1},
        {-1,1,-1,1,1,-1,1,-1},
        {-1,1,1,-1,-1,1,1,-1},
        {1,-1,-1,1,-1,1,1,-1}
        }};

        for(int i=0;i<8;i++){
            for(int j=0;j<8;j++){
                coefs[i]+=corner_values[j]*scalars[i][j];
            }
            coefs[i] /= -8;
        }
    return coefs;
}

std::vector<double> Interpolator_trilinear::InterpolateWCoeffs(const std::array<double, 8>& coeffs, 
    const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z){
        const double N = x.size();
        if (N!=y.size() || N!= z.size()){
            throw std::invalid_argument("different sized x, y, z.");
        }

        std::vector<double> result(N);
        for(int i=0;i<N;i++){
            result[i] = coeffs[0]+coeffs[1]*x[i]+coeffs[2]*y[i]+coeffs[3]*z[i]+coeffs[4]*x[i]*y[i]
                       +coeffs[5]*x[i]*z[i]+coeffs[6]*y[i]*z[i]+coeffs[7]*x[i]*y[i]*z[i];
        }
        return result;
    }
