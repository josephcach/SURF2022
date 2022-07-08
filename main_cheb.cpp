#include "Interpolator_chebyshev.h"
#include <random>


const long long N = 10000;
const int Nmesh = 5;

double function(double x,double y, double z){
    int k=1;
    double mag = std::sqrt(std::pow(x,2)+std::pow(y,2)+std::pow(z,2));
    return std::cos(k*mag);
}

int main(){ 

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-1.0,1.0);
    std::vector<std::vector<std::vector<double> > > coefs;
    coefs = Interpolator_chebyshev::GetCoefficients(Nmesh,&function);
    double x,y,z;
    double error;
    double result;
    double total_err = 0;
    

    for(int i=0;i<N;i++){
            x = dis(gen);
            y = dis(gen);
            z = dis(gen);
            result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,Nmesh,x,y,z);
            error = std::abs(result - function(x,y,z));
            total_err += error;
            printf("[%10f, %10f, %10f] -->  %10f\n",x,y,z,error);
    }
    printf("Average Error over %llu points: %f",N,total_err/N);
return 0;
}

