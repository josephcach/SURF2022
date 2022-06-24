#include <iostream>
#include "Interpolator_trilinear.h"
#include <cstdlib>
#include <cmath>
#include <random>


const long long N = 10000;
const  int Np = 10;
    
double function(double x, double y, double z, std::array<double,3> x0)
{
    double k = 1;
    double mag = std::sqrt(std::pow(x-x0[0],2)+std::pow(y-x0[1],2)+std::pow(z-x0[2],2));
    return cos(k*mag);
}


int main()
{
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-1.0,1.0);
    
    std::vector<double> x(N,0.0), y(N, 0.0), z(N,0.0);
    std::array<double,8> corner_values;
    std::vector<double> error(N,0.0);
    std::array<double,3> x0;
    std::array<std::array<double,3>,8> corners = 
    {{{-1,-1,-1},{-1,-1,1},{-1,1,-1},{-1,1,1},{1,-1,-1},{1,-1,1},{1,1,-1},{1,1,1}}};
    
    for(int j=0; j<Np;j++){

        for(int i=0; i<3;i++){
            x0[i]=dis(gen);
        }

    //fill corner values using the 'function'
        for(int i=0;i<8;i++){
            corner_values[i] = function(corners[i][0],corners[i][1],corners[i][2],x0);
        }
        //std::array<double,8> coefs = Interpolator_trilinear::GetCoefficients(corner_values);

        for(int i=0;i<N;i++){
            x[i] = dis(gen);
            y[i] = dis(gen);
            z[i] = dis(gen);
        }
    //And call the Interpolator on them (or multiple)
    //You also instantiate your interpolation scheme

        std::vector<double> result = Interpolator_trilinear::Interpolate(corner_values, x, y, z);
        //std::vector<double> result = Interpolator_trilinear::InterpolateWCoeffs(coefs, x,y,z);
        ///Compute the error using 'function' again and the positions x y z in a simple for loop.
        double total_err = 0;
        double avg_err;
        for(int i=0; i<N;i++){
            error[i] = std::abs(result[i]-function(x[i],y[i],z[i],x0));
            total_err += error[i];
        }
        avg_err = total_err / N; 
    //Print x0 and mean error.
        printf("[%10f, %10f, %10f] -->  %10f\n",x0[0],x0[1],x0[2],avg_err);
    }

   
    return 0;
}
