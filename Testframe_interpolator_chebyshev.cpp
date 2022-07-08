#include <gtest/gtest.h>
#include "Interpolator_chebyshev.h"
#include <functional>
#include <random>


constexpr int N = 5;
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(-1.0,1.0);

//change all to expect_near

TEST(ChebyshevTest, Degree0Coeff__Test){

    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return 1.0;
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
   
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(int k=0;k<N;k++){
                if(i+j+k!=0){
                EXPECT_NEAR(0.0,coefs[i][j][k],std::pow(10,-14));
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
    EXPECT_NEAR(result,1.0,std::pow(10,-14));
}

TEST(ChebyshevTest, Degree1x__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return x;
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_NEAR(result,x,std::pow(10,-14));
}

TEST(ChebyshevTest, Degree1y__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return y;
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_NEAR(result,y,std::pow(10,-14));
}

TEST(ChebyshevTest, Degree1z__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return z;
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_NEAR(result,z,std::pow(10,-14));
}


TEST(ChebyshevTest, Degree2x__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::pow(x,2);
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_NEAR(result,std::pow(x,2),std::pow(10,-14));
}

TEST(ChebyshevTest, Degree2y__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::pow(y,2);
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_NEAR(result,std::pow(y,2),std::pow(10,-14));
}

TEST(ChebyshevTest, Degree2z__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::pow(z,2);
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_NEAR(result, std::pow(z,2),std::pow(10,-14));
}

TEST(ChebyshevTest, Degree3x__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::pow(x,3);
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_NEAR(result,std::pow(x,3),std::pow(10,-14));
}

TEST(ChebyshevTest, Degree3y__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::pow(y,3);
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_NEAR(result,std::pow(y,3),std::pow(10,-14));
}

TEST(ChebyshevTest, Degree3z__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::pow(z,3);
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_NEAR(result,std::pow(z,3),std::pow(10,-14));
}

TEST(ChebyshevTest, Degree4x__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::pow(x,4);
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_NEAR(result,std::pow(x,4),std::pow(10,-14));
}

TEST(ChebyshevTest, Degree4y__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::pow(y,4);
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_NEAR(result,std::pow(y,4),std::pow(10,-14));
}

TEST(ChebyshevTest, Degree4z__Test){
    std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::pow(z,4);
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    double x=dis(gen); double y = dis(gen); double z = dis(gen);
    double result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,x,y,z);
    EXPECT_NEAR(result,std::pow(z,4),std::pow(10,-14));
}

TEST(ChebyshevTest, ChebPointsError__Test1){
    // Tests accuracy at cheb points for obscure function
     std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::exp(std::cos(x*y*z));
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    std::vector<double> points = Interpolator_chebyshev::GetPoints(N);
    double result;
    double error;
    for(int i=0; i<N; i++){
        for(int j=0; j<N;j++){
            for(int k=0; k<N;k++){
                result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,points[i],points[j],points[k]);
                error = result - func(points[i],points[j],points[k]);
                EXPECT_NEAR(error,0.0,std::pow(10,-14));
            }
        }
      
    }
}

TEST(ChebyshevTest, ChebPointsError__Test2){
    // Tests accuracy at cheb points for obscure function
     std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::atan(x*x +std::exp(20*y)*std::pow(z,3));
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    std::vector<double> points = Interpolator_chebyshev::GetPoints(N);
    double result;
    double error;
    for(int i=0; i<N; i++){
        for(int j=0; j<N;j++){
            for(int k=0; k<N;k++){
                result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,points[i],points[j],points[k]);
                error = result - func(points[i],points[j],points[k]);
                EXPECT_NEAR(error,0.0,std::pow(10,-14));
            }
        }
      
    }
}

TEST(ChebyshevTest, ChebPointsError__Test3){
    // Tests accuracy at cheb points for obscure function
     std::function<double(double, double, double)> func = [](double x, double y, double z){
        return std::sqrt(std::exp(std::sin(1000*x)+std::cos(1000*y)+std::tan(1000*z)));
    };
    std::vector<std::vector<std::vector<double> > > coefs = Interpolator_chebyshev::GetCoefficients(N,func);
    std::vector<double> points = Interpolator_chebyshev::GetPoints(N);
    double result;
    double error;
    for(int i=0; i<N; i++){
        for(int j=0; j<N;j++){
            for(int k=0; k<N;k++){
                result = Interpolator_chebyshev::InterpolateWCoeffs(coefs,N,points[i],points[j],points[k]);
                error = result - func(points[i],points[j],points[k]);
                EXPECT_NEAR(error,0.0,std::pow(10,-14));
            }
        }
      
    }
}
