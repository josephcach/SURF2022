#include "Interpolator_chebyshev.h"
#include <exception>
#include <cmath>
Interpolator_chebyshev::Interpolator_chebyshev()
{
}

Interpolator_chebyshev::~Interpolator_chebyshev()
{
}

/*std::vector<std::vector<double>>  Interpolator_chebyshev::GetPoints(int Nx, int Ny, int Nz){
    
}
*/

double Interpolator_chebyshev::ChebPol(int deg, double x){
    if(deg == 0){
        return 1.;
    }
    else if(deg == 1){
        return x;
    }

    double T0 = 1.;
    double T1 = x;
    double Ti;
    for(int i=2; i<=deg;i++){
        Ti = 2*x*T1-T0;
        T0 = T1;
        T1 = Ti;
    }
    return Ti;
}

std::vector<double>  Interpolator_chebyshev::GetPoints(int N){
    std::vector<double> points(N);
    for(int i=0;i<N;i++){
        points[i] = std::cos(M_PI*double(2*i+1)/(double(2*N)));
    }
    return points;
}

std::vector<std::vector<std::vector<double> > > Interpolator_chebyshev::GetCoefficients(const int N){
    std::vector<std::vector<std::vector<double> > > coefs(N, std::vector<std::vector<double> >(N, std::vector<double>(N, 0))); 
    std::vector<std::vector<std::vector<double> > > fvals(N, std::vector<std::vector<double> >(N, std::vector<double>(N, 0)));
    std::vector<double> points = GetPoints(N);
    
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(int k=0;k<N;k++){
                fvals[i][j][k] = function(points[i],points[j],points[k]);
            }
        }
    }
    
    //for each coeff
    for(int i1=0;i1<N;i1++){
        for(int j1=0;j1<N;j1++){
            for(int k1=0;k1<N;k1++){
                //for each summand per coeff
                for(int i2=0;i2<N;i2++){
                    for(int j2=0;j2<N;j2++){
                        for(int k2=0;k2<N;k2++){
                            coefs[i1][j1][k1] += fvals[i2][j2][k2]*ChebPol(i1,points[i2])*ChebPol(j1,points[j2])*ChebPol(k1,points[k2]);
                        }
                    }
                }
                coefs[i1][j1][k1] /= (std::pow(N,3));
                if(i1!=0){
                    coefs[i1][j1][k1] *= 2;
                }
                if(j1!=0){
                    coefs[i1][j1][k1] *= 2;
                }
                if(k1!=0){
                    coefs[i1][j1][k1] *= 2;
                }
            }
        }
    }
    return coefs;
}

double Interpolator_chebyshev::InterpolateWCoeffs(std::vector<std::vector<std::vector<double> > > coefs, int N, double x, double y, double z){
    double result = 0;
     for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(int k=0;k<N;k++){
                result += coefs[i][j][k]*ChebPol(i,x)*ChebPol(j,y)*ChebPol(k,z);
            }
        }
     }
    return result; 
}

double Interpolator_chebyshev::function(double x,double y, double z){
    int k=1;
    double mag = std::sqrt(std::pow(x,2)+std::pow(y,2)+std::pow(z,2));
    return std::cos(k*mag);
}

