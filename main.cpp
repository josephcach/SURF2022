#include <stdio.h>
#include <iostream>
#include "Interpolator_trilinear.cpp"



double function(double x, double y, double z)
{
    return x*x;
}


int main()
{
    long long N = 10000;
    //this is your driver code
    //here you generate N random points in [-1, 1]^3 as three vectors of size N: x, y, z
    std::vector<double> x(N,0.0), y(N, 0.0), z(N,0.0);
    std::array<double,8> corner_values;
    std::vector<double> error(N,0.0);
    std::array<std::array<double,3>,8> corners = 
    {{{-1,-1,-1},{-1,-1,1},{-1,1,-1},{-1,1,1},{1,-1,-1},{1,-1,1},{1,1,-1},{1,1,1}}};

    //fill corner values using the 'function'
    for(int i=0;i<8;i++){
        corner_values[i] = function(corners[i][0],corners[i][1],corners[i][2]);
    }

    for(int i=0;i<N;i++){
        x[i] = (rand() % 2000000)/1000000 -1;
        y[i] = (rand() % 2000000)/1000000 -1;
        z[i] = (rand() % 2000000)/1000000 -1;
    }
    //And call the Interpolator on them (or multiple)
    //You also instantiate your interpolation scheme
    
    std::vector<double> result = Interpolator_trilinear::Interpolate(corner_values, x, y, z);
   
    ///Compute the error using 'function' again and the positions x y z in a simple for loop.

    for(int i=0; i<N;i++){
        error[i] = abs(result[i]-function(x[i],y[i],z[i]));
        printf("%f\n",error[i]);
    }
    return 0;
}
