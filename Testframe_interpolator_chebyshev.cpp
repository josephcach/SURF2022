#include <gtest/gtest.h>
#include "Interpolator_chebyshev.h"
#include <functional>
#include <random>
    constexpr int N = 5;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-1.0,1.0);

TEST(ChebyshevTest, Degree0Coeff__Test){

    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return 1.0;
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
   
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(int k=0;k<N;k++){
                if(i+j+k!=0){
                EXPECT_DOUBLE_EQ(0.0,coefs[i][j][k]);
                }
            }
        }
    }
}

TEST(ChebyshevTest, Degree0Error__Test){

    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return 1.0;
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_DOUBLE_EQ(result,1.0);
}

TEST(ChebyshevTest, Degree1x__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return x;
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_DOUBLE_EQ(result,x);
}

TEST(ChebyshevTest, Degree1y__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return y;
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_DOUBLE_EQ(result,y);
}

TEST(ChebyshevTest, Degree1z__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return z;
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_DOUBLE_EQ(result,z);
}

TEST(ChebyshevTest, Degree4x__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::pow(x,4);
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_DOUBLE_EQ(result,std::pow(x,4));
}

TEST(ChebyshevTest, Degree4y__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::pow(y,4);
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_DOUBLE_EQ(result,std::pow(y,4));
}

TEST(ChebyshevTest, Degree4z__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::pow(z,4);
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_DOUBLE_EQ(result,std::pow(z,4));
}

TEST(ChebyshevTest, ChebPoints__Test){
    // Tests angular equidistance of chebyshev points
    std::vector<double> points = Interpolator_chebyshev::GetPoints(N);
    double theta1;
    double theta2;
    double theta3;
    for(int i=0; i<N-1; i++){
        theta1 = std::acos(points[i]);
        theta2 = std::acos(points[i+1]);
        theta3 = std::acos(points[i+2]);
        EXPECT_DOUBLE_EQ(std::abs(theta1-theta2),std::abs(theta2-theta3));
    }
}

int main(){
testing::InitGoogleTest();

return RUN_ALL_TESTS();
}