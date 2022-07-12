#include "Piecewise_trilinear.h"
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <random>

double function(double x, double y, double z){
    int k = 1;
    double mag = std::sqrt(std::pow(x,2)+std::pow(y,2)+std::pow(z,2));
    return std::cos(k*mag);
}

const int N = 101;
const int num_targets = 10000000;
int main()
{
    std::vector<std::vector<std::vector<double> > > corners = Piecewise_trilinear::GetCorners(N,&function);
    double result;
    double total_error;
    std::random_device rd;
    std::mt19937 gen(1);
    std::uniform_real_distribution<double> dis(-1.0,1.0);
    double x,y,z;
    for(int i=0; i<num_targets;i++){
        x = dis(gen);
        y = dis(gen);
        z = dis(gen);
        result = Piecewise_trilinear::Interpolate(N, corners, x, y, z);
        total_error += std::abs(result - function(x,y,z));
    }
    printf("Test Complete.\nGrid Size: %d\nNum Targets: %d\nAverage Error: %10f", N-1, num_targets, total_error/num_targets);
    return 0; 
}