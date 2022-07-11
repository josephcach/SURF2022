#include "Piecewise_trilinear.h"
#include <cstdlib>
#include <cmath>
#include <stdio.h>


double function(double x, double y, double z){
    return x*x;
}

const int N = 100;
int main()
{
   std::vector<std::vector<std::vector<double> > > corners = Piecewise_trilinear::GetCorners(N,&function);
   double result = Piecewise_trilinear::Interpolate(N, corners, .05,.05,.05);
   printf("%f", result); 
   return 0; 
}