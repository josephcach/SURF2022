#include <gtest/gtest.h>
#include "Piecewise_trilinear.h"
#include <functional>
#include <random>


std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(-1.0,1.0);

TEST(PWTrilinearTest,Degree0Error_Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return 1.0;
    };
    int N = 101;
    std::vector<std::vector<std::vector<double> > > corners = Piecewise_trilinear::GetCorners(N,func);
    double x = dis(gen);
    double y = dis(gen);
    double z = dis(gen);
    double result = Piecewise_trilinear::Interpolate(N, corners, x,y,z);
    EXPECT_NEAR(result, func(x,y,z),std::pow(10,-14));
}

TEST(PWTrilinearTest, Degree1x_Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return x;
    };
    int N = 101;
    std::vector<std::vector<std::vector<double> > > corners = Piecewise_trilinear::GetCorners(N,func);;
    double x = dis(gen); double y = dis(gen); double z = dis(gen);
    double result  = Piecewise_trilinear::Interpolate(N, corners, x,y,z);
    EXPECT_NEAR(result, func(x,y,z),std::pow(10,-14));
}

TEST(PWTrilinearTest, Degree1y_Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return y;
    };
    int N = 101;
    std::vector<std::vector<std::vector<double> > > corners = Piecewise_trilinear::GetCorners(N,func);;
    double x = dis(gen); double y = dis(gen); double z = dis(gen);
    double result  = Piecewise_trilinear::Interpolate(N, corners, x,y,z);
    EXPECT_NEAR(result, func(x,y,z),std::pow(10,-14));
}

TEST(PWTrilinearTest, Degree1z_Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return z;
    };
    int N = 101;
    std::vector<std::vector<std::vector<double> > > corners = Piecewise_trilinear::GetCorners(N,func);;
    double x = dis(gen); double y = dis(gen); double z = dis(gen);
    double result  = Piecewise_trilinear::Interpolate(N, corners, x,y,z);
    EXPECT_NEAR(result, func(x,y,z),std::pow(10,-14));
}

TEST(PWTrilinearTest, Error_Convergence1_Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
    double mag = std::sqrt(std::pow(x,2)+std::pow(y,2)+std::pow(z,2));
    return std::cos(mag);
    };
    int N1 = 51; int N2 = 101; int N3 = 201;
    int num_targets = 1000000;
    double total_error1, total_error2, total_error3;
    std::vector<std::vector<std::vector<double> > > corners1 = Piecewise_trilinear::GetCorners(N1,func);
    std::vector<std::vector<std::vector<double> > > corners2 = Piecewise_trilinear::GetCorners(N2,func);
    std::vector<std::vector<std::vector<double> > > corners3 = Piecewise_trilinear::GetCorners(N3,func);
    for(int i = 0; i<num_targets; i++){
        double x = dis(gen); double y = dis(gen); double z = dis(gen);
        total_error1 += std::abs(Piecewise_trilinear::Interpolate(N1, corners1,x,y,z)-func(x,y,z));
        total_error2 += std::abs(Piecewise_trilinear::Interpolate(N2, corners2,x,y,z)-func(x,y,z));
        total_error3 += std::abs(Piecewise_trilinear::Interpolate(N3, corners3,x,y,z)-func(x,y,z));
    }

    EXPECT_NEAR(total_error1/(4*num_targets),total_error2/num_targets,std::pow(10,-6));
    EXPECT_NEAR(total_error2/(4*num_targets),total_error3/num_targets,std::pow(10,-6));
}